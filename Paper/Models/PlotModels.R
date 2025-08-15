library(lme4)
library(lmerTest)
library(easystats)
library(tidyverse)

df = read.csv("resources/Stats/Dataset.csv")
df$Id = factor(df$Id) # make sure subject_id is a factor
df$SES = factor(df$SES) # make sure subject_id is a factor
df$Event = factor(df$Event) # make sure subject_id is a factor

df$StandardTrialN = standardize(df$TrialN) # create standardize trial column


mod_gamma <- glmer(
  ReactionTime ~ StandardTrialN * Event + (1 + StandardTrialN| Id),
  data   = df,
  family = Gamma(link = "log")
)
check_model(mod_gamma)
ggsave("Paper/Models/ModelChecks.png", width = 18, height = 12)


Grid = get_datagrid(mod_gamma, by = c('Event', 'StandardTrialN'))
Pred = as.data.frame(get_predicted(mod_gamma, Grid, ci  = .95))
Pred = bind_cols(Grid, Pred)
Pred$TrialN = unstandardize(Pred$StandardTrialN, reference  = df$TrialN) # create standardize trial column


ggplot(Pred, aes(x = TrialN, y = Predicted, color = Event, fill=Event))+
  geom_ribbon(aes(ymin = Predicted -SE, ymax = Predicted + SE), alpha = 0.4, color = 'transparent') + # Confidence bands
  geom_line(lwd = 1.3, show.legend = FALSE)+
  labs(x = "Trial Number", y = "Predicted Looking Time (ms)")+
  theme_classic(base_size = 25)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(breaks = seq(0,20,2))
ggsave("Paper/Models/ModelPredictions.png", width = 18, height = 10,dpi=300)

