# --- Load required libraries ---
library(lme4)         # For linear mixed-effects models
library(lmerTest)     # For p-values in lmer models
library(easystats)    # For model checking and prediction utilities
library(tidyverse)    # For data manipulation and plotting
library(gt)           # For creating tables
library(patchwork)    # For combining ggplot objects

# --- Load and preprocess data ---
df = read.csv("resources/Stats/Dataset.csv")

# Ensure categorical variables are factors
df$Id    = factor(df$Id)
df$SES   = factor(df$SES)
df$Event = factor(df$Event)

# Standardize trial number for modeling
df$StandardTrialN = standardize(df$TrialN)

# --- Fit mixed-effects models ---
# Linear mixed model for LookingTime
mod <- lmer(
  LookingTime ~ StandardTrialN * Event + (1 + StandardTrialN | Id),
  data = df
)
check_model(mod) # Model diagnostics

# Generalized linear mixed model (Gamma) for ReactionTime
mod_gamma <- glmer(
  SaccadicRT ~ StandardTrialN * Event + (1 + StandardTrialN | Id),
  data   = df,
  family = Gamma(link = "log")
)
check_model(mod_gamma) # Model diagnostics

# Save model check plots
ggsave("Paper/Models/ModelChecks.png", width = 18, height = 12)

# --- Generate predictions for plotting ---
# Create prediction grid and get predicted values for linear model
Grid1 = get_datagrid(mod, by = c('Event', 'StandardTrialN'))
Pred1 = as.data.frame(get_predicted(mod, Grid1, ci = .95))
Pred1 = bind_cols(Grid1, Pred1)
Pred1$TrialN = unstandardize(Pred1$StandardTrialN, reference = df$TrialN)

# Create prediction grid and get predicted values for gamma model
Grid2 = get_datagrid(mod_gamma, by = c('Event', 'StandardTrialN'))
Pred2 = as.data.frame(get_predicted(mod_gamma, Grid2, ci = .95))
Pred2 = bind_cols(Grid2, Pred2)
Pred2$TrialN = unstandardize(Pred2$StandardTrialN, reference = df$TrialN)

# --- Plot model predictions ---
# Plot for linear model (LookingTime)
Linear = ggplot(Pred1, aes(x = TrialN, y = Predicted, color = Event, fill = Event)) +
  geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE), alpha = 0.4, color = 'transparent') +
  geom_line(lwd = 1.3, show.legend = FALSE) +
  labs(x = NULL, y = "Predicted looking time (ms)") +
  theme_classic(base_size = 35) +
  theme(legend.position = 'none') +
  scale_x_continuous(breaks = seq(0, 20, 2))

# Plot for gamma model (ReactionTime)
Gamma = ggplot(Pred2, aes(x = TrialN, y = Predicted, color = Event, fill = Event)) +
  geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE), alpha = 0.4, color = 'transparent') +
  geom_line(lwd = 1.3, show.legend = FALSE) +
  labs(x = "Trial Number", y = "Predicted saccadic reaction time (ms)") +
  theme_classic(base_size = 35) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = seq(0, 20, 2))

# --- Combine and save plots ---
(Linear | Gamma) +
  plot_layout(guides = 'collect') +
  plot_annotation(theme = theme(legend.position = 'bottom'))

ggsave("Paper/Models/ModelPredictionsWide.png", width = 28, height = 14, dpi = 300)
