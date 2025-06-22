# =============================================================================
# COMBINE DATASETS SCRIPT
# =============================================================================
# This script combines simulated reaction time and looking time datasets into 
# a single comprehensive dataset for analysis. It also adds demographic variables
# and creates visualizations of the combined data.
#
# Input files:
# - simulatedLognormal.csv (reaction time data)
# - simulatedNormal.csv (looking time data)
#
# Output:
# - Dataset.csv (combined dataset)
# - CombinedDensity.png (visualization)
# =============================================================================

# Load Required Libraries -----------------------------------------------------
library(tidyverse)  # Data manipulation and visualization
library(patchwork)  # Combining multiple plots
library(lmerTest)   # Linear mixed-effects models  
library(easystats)  # Easy statistical modeling and reporting

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

# Load Reaction Time Data -----------------------------------------------------
# Load the lognormally distributed reaction time data from the simulation
ReactionTime <- vroom::vroom("CONTENT\\Simulation\\simulatedLognormal.csv") %>% 
  rename(
    ReactionTime = reaction_time,  # Rename for clarity
    Id = subject_id,              # Standardize ID column name
    Event = categorical_condition, # Standardize condition column name
    TrialN = trial_number         # Standardize trial number column name
  ) %>% 
  mutate(Event = fct_rev(Event)) %>%  # Reverse factor levels for consistent ordering
  select(Id, Event, TrialN, ReactionTime)  # Keep only relevant columns

# Load Looking Time Data ------------------------------------------------------
# Load the normally distributed looking time data from the simulation
LookingTime <- vroom::vroom("CONTENT\\Simulation\\simulatedNormal.csv") %>% 
  rename(
    LookingTime = dependent_variable,  # Rename for clarity
    Id = subject_id,                  # Standardize ID column name
    Event = categorical_condition,     # Standardize condition column name
    TrialN = trial_number             # Standardize trial number column name
  ) %>% 
  select(Id, Event, TrialN, LookingTime)  # Keep only relevant columns

# =============================================================================
# DATA COMBINATION AND CLEANING
# =============================================================================

# Combine Datasets ------------------------------------------------------------
# Merge reaction time and looking time data by participant, condition, and trial
Df <- left_join(ReactionTime, LookingTime, by = c("Id", "Event", "TrialN")) %>% 
  mutate(
    # Recode event labels to be more meaningful
    Event = fct_recode(Event, 
                       NoReward = "Hammer",  # Hammer condition = No Reward
                       Reward = "Spoon")     # Spoon condition = Reward
  )

# Handle Missing Data ---------------------------------------------------------
# If reaction time is missing (NA), also set looking time to missing
# This maintains data consistency across both measures
Df[is.na(Df$ReactionTime), ]$LookingTime <- NA

# =============================================================================
# ADD DEMOGRAPHIC VARIABLES
# =============================================================================

# Add Socioeconomic Status (SES) ---------------------------------------------
# Randomly assign SES levels to each participant (one per participant)
# This creates between-subject variation in socioeconomic status
Df <- Df %>%
  group_by(Id) %>%
  mutate(SES = sample(c("low", "medium", "high"), 1)) %>%  # Random SES assignment
  ungroup()

# Quick check of SES distribution (optional visualization)
# ggplot(Df, aes(x = SES)) + geom_bar()  # Uncomment to see SES distribution

# =============================================================================
# SAVE COMBINED DATASET
# =============================================================================

# Export Final Dataset -------------------------------------------------------
# Save the combined and cleaned dataset for further analysis
write_csv(Df, "resources\\Stats\\Dataset.csv")

# =============================================================================
# STATISTICAL MODELS
# =============================================================================

# Looking Time Model ----------------------------------------------------------
# Linear mixed-effects model for normally distributed looking time data
# Random intercepts and slopes for each participant
mod_l <- lmer(LookingTime ~ Event * TrialN + (1 + TrialN | Id), data = Df)

# Reaction Time Model ---------------------------------------------------------
# Generalized linear mixed-effects model for reaction time data
# Uses Gamma distribution with log link (appropriate for reaction times)
Df$Stand_TrialN <- datawizard::standardise(Df$TrialN)  # Standardize trial numbers
mod_gam <- glmer(ReactionTime ~ Event * Stand_TrialN + 
                   (1 + Stand_TrialN | Id),
                 family = Gamma(link = 'log'), 
                 data = Df)

# =============================================================================
# DATA VISUALIZATION
# =============================================================================

# Create Density Plots -------------------------------------------------------

# Looking Time Distribution
D1 <- ggplot(Df, aes(x = LookingTime, fill = Event)) +
  geom_density(color = 'transparent', alpha = 0.3) +      # Filled density curves
  geom_density(fill = 'transparent', color = 'black', lwd = 2) +  # Overall outline
  labs(
    title = "Distribution of Looking Times",
    x = "Looking Time (ms)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Reaction Time Distribution  
D2 <- ggplot(Df, aes(x = ReactionTime, fill = Event)) +
  geom_density(color = 'transparent', alpha = 0.3) +      # Filled density curves
  geom_density(fill = 'transparent', color = 'black', lwd = 2) +  # Overall outline
  labs(
    title = "Distribution of Reaction Times", 
    x = "Reaction Time (ms)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine and Save Plots ------------------------------------------------------
# Use patchwork to combine both density plots side by side
combined_plot <- D1 + D2
combined_plot

# Save the combined visualization
ggsave("CONTENT\\Simulation\\CombinedDensity.png", 
       plot = combined_plot,
       width = 20, 
       height = 10,  # Reduced height since plots are side by side
       dpi = 300)

# =============================================================================
# MODEL VISUALIZATION AND EXPLORATION
# =============================================================================

# Visualize Model Relationships ----------------------------------------------
# These plots show estimated relationships from the fitted models

# Looking time model visualizations
# plot(estimate_relation(mod_l), by = c('TrialN', 'Event'))      # Trial effects by Event
# plot(estimate_relation(mod_l), by = c('Event', 'TrialN'))      # Event effects by Trial
# plot(estimate_expectation(mod_l, by = c('TrialN', 'Event'), transform = T))  # Expected values

# Uncomment the lines above to generate model visualization plots