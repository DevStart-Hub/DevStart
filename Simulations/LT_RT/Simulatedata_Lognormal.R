# =============================================================================
# SIMULATE LOGNORMAL REACTION TIME DATA WITH INTERACTIONS
# =============================================================================
# This script simulates lognormally distributed reaction time data with:
# - Random intercepts and slopes for subjects
# - Fixed effects for categorical conditions 
# - Fixed effects for continuous trial number (learning/practice effect)
# - Interaction effects between condition and trial number
# - Realistic measurement error and data constraints
#
# The simulation creates realistic reaction time patterns with individual
# differences and experimental effects, suitable for mixed-effects modeling.
# =============================================================================

# Load Libraries -----------------------------------------------------------
# Load necessary libraries for data manipulation, visualization, and statistical modeling.
library(tidyverse)    # For data manipulation (dplyr) and plotting (ggplot2)
library(easystats)    # Provides functions like estimate_expectation and check_model
library(lmerTest)     # For fitting linear mixed-effects models
library(GGally)       # For visualizing pairwise relationships
library(faux)         # For simulating correlated variables

# Set Options --------------------------------------------------------------
# Configure how numbers are displayed in console output
options("scipen" = 10, "digits" = 4) # Prevent scientific notation; show 4 decimal places


# Variables ----------------------------------------------------------------
# Define parameters for simulating reaction time data
n_subjects  <- 20                        # Number of subjects in the study
ntrials = 20                             # Trials per subject
random_intercept_sd <- 0.3               # SD of random intercepts (log scale)
random_slope_sd <- 0.02                  # SD of random slopes for trial number
random_intercept_slope_cor <- -0.2       # Correlation between intercept & slope

# Parameters on log scale (will be exponentiated to get ms RTs)
grand_mean_log_rt <- 6                   # Grand mean on log scale (~400ms when exponentiated)
fixed_categorical_effect <- -0.1         # Effect of categorical condition (e.g., Complex vs Simple)
fixed_continuous_effect <- -0.02         # Continuous effect of trial number (learning/practice)
interaction_effect <- -0.034             # Interaction between condition and trial number
residual_sd <- 0.1                     # Residual error on log scale


# =============================================================================
# CREATE SUBJECT-TRIAL DATA STRUCTURE
# =============================================================================

# Create subject - trials - Event ----------------------------------------------------------

# If the trial sequence file doesn't exist, generate it; otherwise load it
if (!file.exists("trial_sequences.csv")) {
  # Generate all combinations of subjects and trials
  SubTrials = crossing(
    subject_id = 1:n_subjects,
    trial_number = 1:ntrials)

  SubTrials$categorical_condition <- NA

  # Loop through each subject to assign a constrained-randomized condition sequence
  for (S in 1:n_subjects) {
    # Keep generating until no condition repeats more than 3 times in a row
    repeat {
      seq <- sample(c("Complex", "Simple"), ntrials, replace = TRUE)
      if (all(rle(seq)$lengths <= 3)) break # Exit loop once valid sequence is found
    }
    # Assign generated sequence to current subject
    SubTrials$categorical_condition[SubTrials$subject_id == S] <- seq
  }
  # Save to disk for future use
  write.csv(SubTrials, "Simulations/LT_RT/trial_sequences.csv", row.names = FALSE)
} else {
  # Load existing trial sequences
  SubTrials <- read.csv("Simulations/LT_RT/trial_sequences.csv")
}


# =============================================================================
# GENERATE SUBJECT-LEVEL RANDOM EFFECTS
# =============================================================================

# Random Effects for Subjects ----------------------------------------------
# Simulate subject-level random effects (intercept and slope) with specified correlation
subjects <- faux::rnorm_multi(
  n = n_subjects,
  vars = 2,
  r = random_intercept_slope_cor,
  mu = 0,
  sd = c(random_intercept_sd, random_slope_sd),
  varnames = c("random_intercept", "random_slope_trial")
) %>%
  mutate(subject_id = 1:n_subjects)


# =============================================================================
# MERGE TRIAL AND SUBJECT DATA
# =============================================================================

# Trials ------------------------------------------------------------------

# Merge trial-level data with subject-level random effects
trial_data <- left_join(SubTrials, subjects, by = "subject_id")
trial_data$subject_id <- as.factor(trial_data$subject_id)


# =============================================================================
# SIMULATE DEPENDENT VARIABLE (REACTION TIMES)
# =============================================================================

# Dependent Variable Calculation (Lognormal Reaction Times) ----------------

##### Simulate Data ####
simulated_data <- trial_data %>%
  mutate(
    # Code categorical condition as numeric for modeling (-0.5 or 0.5)
    categorical_coded = recode(categorical_condition, "Complex" = 0.5, "Simple" = -0.5),
    
    # Compute fixed effect components
    fixed_categorical = fixed_categorical_effect * categorical_coded,
    interaction_term = categorical_coded * trial_number * interaction_effect,
    
    # Add random error term
    random_error = rnorm(nrow(.), 0, residual_sd),
    
    # Calculate log reaction time using full model formula
    log_rt = grand_mean_log_rt + 
      (random_intercept) +
      (random_slope_trial * trial_number) +
      (fixed_continuous_effect * trial_number) +
      fixed_categorical +
      interaction_term +
      random_error,
    
    # Convert log RT to milliseconds (lognormal distribution)
    SaccadicRT = exp(log_rt),
    
    # Remove implausibly fast responses (< threshold)
    treshold = rnorm(nrow(.), 210, 10),
    SaccadicRT = ifelse(SaccadicRT < treshold, NA, SaccadicRT)
  )


# =============================================================================
# FIT STATISTICAL MODELS
# =============================================================================

# Models ------------------------------------------------------------------
simulated_data$stand_TrialN = datawizard::standardise(simulated_data$trial_number)

# Fit a linear mixed-effects model on log-transformed reaction times
mod_rt <- lmer(SaccadicRT ~ categorical_condition * stand_TrialN + 
                 (1 + stand_TrialN | subject_id), data = simulated_data)


# Fit a generalized linear mixed model using Gamma distribution
mod_gam = glmer(SaccadicRT ~ categorical_condition * stand_TrialN + 
                  (1 + stand_TrialN | subject_id),
                family = Gamma(link = 'log'), data = simulated_data)


# =============================================================================
# MODEL DIAGNOSTICS AND VALIDATION
# =============================================================================

#### Check models ####

# Output summaries and diagnostics
summary(mod_rt)
check_model(mod_rt)

summary(mod_gam)
check_model(mod_gam)



# =============================================================================
# SAVE SIMULATED DATA
# =============================================================================

#### Save data ####
# Export simulated data to CSV
write.csv(simulated_data, "Simulations/LT_RT/simulatedLognormal.csv", row.names = FALSE)



# =============================================================================
# VISUALIZE SIMULATED DATA
# =============================================================================

# Visualize the Simulated Data --------------------------------------------
# Plot individual subject lines colored by condition
ggplot(simulated_data, aes(x = trial_number, y = SaccadicRT, color = categorical_condition)) +
  geom_line() +
  geom_point() +
  labs(title = "Simulated Reaction Times (ms)",
       x = "Trial Number",
       y = "Reaction Time (ms)",
       color = "Condition") +
  theme_minimal() +
  facet_wrap(vars(subject_id), scales = 'free')


# Histogram showing density of reaction times by condition
ggplot(simulated_data, aes(x = SaccadicRT, fill = categorical_condition)) +
  geom_density(alpha = 0.3, color = "black") +
  labs(title = "Distribution of Reaction Times",
       x = "Reaction Time (ms)",
       y = "Count",
       fill = "Condition") +
  theme_minimal()




# =============================================================================
# MODEL PREDICTIONS AND CONTRASTS
# =============================================================================

# Estimated Expectation ---------------------------------------------------

# Estimate expected values from the Gamma model
Pred = estimate_expectation(mod_gam, by = c('stand_TrialN', 'categorical_condition'), transform = TRUE)

# Plot predicted means with confidence intervals
ggplot(Pred, aes(x = stand_TrialN, y = Predicted, color = categorical_condition, fill = categorical_condition)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .5, color = 'transparent') +
  geom_line(lwd = 1.3)+
  theme_modern()

# Test contrasts between conditions at different trial numbers
estimate_contrasts(mod_gam, contrast = 'stand_TrialN', by = 'categorical_condition', transform = TRUE)
