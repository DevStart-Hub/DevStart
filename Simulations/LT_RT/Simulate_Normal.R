# =============================================================================
# SIMULATE INTERACTION TRIAL DATA - NORMAL DISTRIBUTION
# =============================================================================
# This script simulates normally distributed dependent variable data with:
# - Random intercepts and slopes for subjects
# - Fixed effects for categorical conditions (Spoon vs. Hammer)
# - Fixed effects for continuous variable (trial number)
# - Interaction effects between condition and trial number
# - Realistic measurement error
#
# The simulation includes a validation loop to ensure statistical power
# and model assumptions are met before saving the final dataset.
# =============================================================================

# Load Required Libraries -----------------------------------------------------
library(tidyverse) # Data manipulation and visualization
library(easystats) # Easy statistical modeling
library(lmerTest) # Mixed-effect models
library(GGally) # Plotting relationships between variables
library(faux) # Simulating correlated variables

# Set Global Options ----------------------------------------------------------
options("scipen" = 10, "digits" = 4) # Prevent scientific notation, display up to 4 digits
# set.seed(8675309) # Uncomment to set a seed for reproducibility

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Sample Size Parameters -----------------------------------------------------
n_subjects <- 20 # Number of participants in the study
trialn <- 20     # Number of trials per condition per participant

# Random Effects Parameters --------------------------------------------------
# These control individual differences between participants
random_intercept_sd <- 150           # Standard deviation of baseline differences between subjects
random_slope_sd <- 8                 # Standard deviation of learning rate differences between subjects
random_intercept_slope_cor <- -0.2   # Correlation between baseline and learning rate
                                     # Negative correlation: higher baseline = slower learning

# Fixed Effects Parameters ---------------------------------------------------
# These are the population-level effects we want to detect
grand_mean_dv <- 400                 # Overall average of the dependent variable across all conditions

# Effect sizes (adjusted to create realistic opposing effects)
fixed_categorical_effect <- -50      # Main effect of condition (Spoon vs. Hammer)
                                     # When effect-coded: Spoon = +25, Hammer = -25
fixed_continuous_effect <- -8        # Learning effect: DV decreases by 8 units per trial
interaction_effect <- -4             # Interaction: differential learning rates between conditions
                                     # Adds/subtracts 4 units per trial depending on condition

# Error Parameters ------------------------------------------------------------
residual_sd <- 30                    # Within-subject measurement error (trial-to-trial noise)

# =============================================================================
# GENERATE SUBJECT-LEVEL RANDOM EFFECTS
# =============================================================================

# Create Correlated Random Effects -------------------------------------------
# This section generates subject-specific deviations from the population means.
# Each subject gets their own intercept (baseline performance) and slope 
# (learning rate), which are correlated according to the specified correlation.

subjects <- faux::rnorm_multi(
  n = n_subjects,                              # Number of subjects to generate
  vars = 2,                                    # Two random effects: intercept and slope
  r = random_intercept_slope_cor,              # Correlation between the two effects
  mu = 0,                                      # Both effects centered at zero (deviations from grand mean)
  sd = c(random_intercept_sd, random_slope_sd), # Different variability for intercept vs. slope
  varnames = c("random_intercept", "random_slope_trial") # Descriptive names
) %>%
  mutate(subject_id = 1:n_subjects)            # Add unique identifier for each subject

# =============================================================================
# GENERATE TRIAL-LEVEL DATA STRUCTURE
# =============================================================================

# Create All Possible Trial Combinations -------------------------------------
# This creates the full factorial design: every subject experiences every 
# combination of condition and trial number.

trial_data <- crossing(
  subject_id = subjects$subject_id,                    # All subject IDs
  categorical_condition = c("Spoon", "Hammer"),        # Both experimental conditions
  trial_number = 1:trialn                             # All trial numbers (1 to 20)
) 
# Total observations = n_subjects × 2 conditions × trialn trials = 20 × 2 × 20 = 800

# Merge Subject Random Effects with Trial Data -------------------------------
# This joins the subject-specific random effects with the trial structure,
# so each trial inherits the appropriate subject's random intercept and slope.
trial_data <- left_join(trial_data, subjects, by = "subject_id")

# =============================================================================
# SIMULATION VALIDATION LOOP
# =============================================================================
# This loop ensures that the simulated data meets our statistical requirements:
# - Significant interaction effect (p < 0.05)
# - Normal residuals in the mixed-effects model
# If these conditions aren't met, new data is generated automatically.

# Initialize Validation Variables ---------------------------------------------
P1 <- -Inf  # p-value for interaction effect in linear model
P2 <- -Inf  # p-value for normality test in mixed-effects model
iter <- 1   # Iteration counter

# Validation Loop -------------------------------------------------------------
while(P1 < 0.05 | P2 < 0.05) {
  print(paste("Simulation attempt:", iter))
  iter <- iter + 1
  
  # ==========================================================================
  # CALCULATE DEPENDENT VARIABLE
  # ==========================================================================
  
  simulated_data <- trial_data %>%
    mutate(
      # Effect Coding for Categorical Variable -------------------------------
      # Convert categorical condition to numeric codes for analysis
      # Effect coding ensures that the intercept represents the grand mean
      categorical_coded = recode(
        categorical_condition,
        "Spoon" = -0.5,    # Spoon condition gets -0.5
        "Hammer" = 0.5     # Hammer condition gets +0.5
      ),
      
      # Calculate Component Effects ------------------------------------------
      
      # Fixed effect of categorical variable (condition main effect)
      fixed_categorical = fixed_categorical_effect * categorical_coded,
      
      # Interaction effect between condition and trial number
      # This allows the learning rate to differ between conditions
      interaction_term = categorical_coded * trial_number * interaction_effect,
      
      # Random measurement error for each trial
      random_error = rnorm(nrow(.), 0, residual_sd),
      
      # Combine All Effects into Final DV ------------------------------------
      # The dependent variable is the sum of all systematic and random effects
      dependent_variable = grand_mean_dv +                              # Grand mean
        (random_intercept) +                                           # Subject's baseline deviation
        (random_slope_trial * trial_number) +                         # Subject's learning rate × trial
        (fixed_continuous_effect * trial_number) +                    # Population learning effect
        fixed_categorical +                                            # Condition main effect
        interaction_term +                                             # Condition × trial interaction
        random_error                                                   # Random measurement noise
    )
  
  # Rescale to Realistic Range -----------------------------------------------
  # Ensure the dependent variable falls within a plausible range (400-1980)
  # This maintains the relative differences while constraining absolute values
  simulated_data$dependent_variable <- rescale(
    simulated_data$dependent_variable,
    to = c(400, 1980)
  )
  
  # ==========================================================================
  # MODEL FITTING AND VALIDATION
  # ==========================================================================
  
  # Linear Model (Ignores Subject-Level Dependencies) ----------------------
  # This model treats all observations as independent
  # Used primarily to check the interaction effect
  mod_lm <- lm(
    dependent_variable ~ categorical_condition * trial_number,
    data = simulated_data
  )
  
  # Mixed-Effects Model (Accounts for Subject Dependencies) ----------------
  # This model includes random intercepts and slopes for each subject
  # More appropriate for the hierarchical structure of the data
  mod_mixed <- lmer(
    dependent_variable ~ categorical_condition * trial_number + 
      (1 + trial_number | subject_id),
    data = simulated_data
  )
  
  # Extract Validation Metrics ----------------------------------------------
  P1 <- parameters(mod_lm)[4, 'p']        # p-value for interaction in linear model
  P2 <- check_normality(mod_mixed)$p       # p-value for normality of residuals in mixed model
  
  # Continue loop if validation criteria not met
}

# =============================================================================
# SAVE VALIDATED DATA
# =============================================================================

# Export Final Dataset -------------------------------------------------------
write.csv(simulated_data, "CONTENT\\Simulation\\simulatedNormal.csv", row.names = FALSE)

# =============================================================================
# MODEL DIAGNOSTICS
# =============================================================================

# Check Model Assumptions ----------------------------------------------------
# Generate diagnostic plots to verify model fit and assumptions
check_model(mod_lm)      # Diagnostics for linear model
check_model(mod_mixed)   # Diagnostics for mixed-effects model

# =============================================================================
# DATA VISUALIZATION
# =============================================================================

# Distribution Visualization -------------------------------------------------
# Show the overall distribution of the dependent variable by condition
ggplot(
  simulated_data,
  aes(x = dependent_variable, fill = categorical_condition)
) +
  geom_density(alpha = 0.3, color = "transparent") +  # Filled density curves by condition
  geom_density(color = "black", fill = "transparent") + # Overall density outline
  labs(
    title = "Distribution of Simulated Dependent Variable",
    x = "Dependent Variable",
    y = "Density",
    fill = "Condition"
  ) +
  theme_modern()

# Individual Subject Trajectories --------------------------------------------
# Plot each subject's performance across trials, separated by condition
# This visualization shows both individual differences and condition effects
ggplot(
  simulated_data,
  aes(x = trial_number, y = dependent_variable, color = categorical_condition)
) +
  geom_line() +    # Connect points within each condition
  geom_point() +   # Show individual data points
  labs(
    title = "Individual Subject Trajectories: Performance Across Trials",
    subtitle = "Random Intercepts and Slopes with Categorical and Interaction Effects",
    x = "Trial Number",
    y = "Dependent Variable (DV)",
    color = "Condition"
  ) +
  theme_minimal() +
  facet_wrap(vars(subject_id), scales = 'free') +  # Separate panel for each subject
  theme(
    strip.text = element_text(size = 8),           # Smaller subject ID labels
    legend.position = "bottom"                     # Move legend to bottom
  )

# =============================================================================
# MODEL PREDICTIONS AND EXPECTED VALUES
# =============================================================================

# Generate Model-Based Predictions -------------------------------------------
# Use the fitted mixed-effects model to estimate expected values across
# the range of trial numbers and categorical conditions, including confidence intervals

Est <- estimate_expectation(
  mod_mixed,
  get_datagrid(mod_mixed, by = c('categorical_condition', 'trial_number'))
)

# Plot Expected Trajectories with Uncertainty --------------------------------
ggplot(
  Est,
  aes(
    x = trial_number,
    y = Predicted,
    color = categorical_condition,
    fill = categorical_condition
  )
) +
  geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE), alpha = 0.4) + # Confidence bands
  geom_line(size = 1.2) +  # Predicted trajectories
  labs(
    title = "Model Predictions: Expected Performance Trajectories",
    subtitle = "Ribbons show ±1 Standard Error",
    x = "Trial Number",
    y = "Predicted Dependent Variable",
    color = "Condition",
    fill = "Condition"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")