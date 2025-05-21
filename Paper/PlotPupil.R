# Libraries and files --------------------------------------------------------------------

library(PupillometryR) # Library to process pupil signal
library(tidyverse) # Library to wrangle dataframes
library(patchwork)

csv_files = list.files(
  path = "resources\\Pupillometry\\RAW",
  pattern = "\\.csv$", # regex pattern to match .csv files
  full.names = TRUE # returns the full file paths
)


# Prepare data --------------------------------------------------------------------

## Event settings --------------------------------------------------------------------

# Settings to cut events
Fs = 300 # framerate
Step = 1000 / Fs

Pre_stim = 100 # pre stimulus time (100ms)
Post_stim = 2000 # post stimulus time (2000ms)
Pre_stim_samples = Pre_stim / Step # pre stimulus in samples
Post_stim_samples = Post_stim / Step # post stimulus in samples

# Time vector based on the event duration
Time = seq(
  from = -Pre_stim,
  by = Step,
  length.out = (Pre_stim + Post_stim) / Step
) # time vector

# Events to keep
Events_to_keep = c('Circle', 'Square')



# Rawplot all ------------------------------------------------------------

Raw_plot = vroom::vroom(csv_files) |> 
  filter(time<78000)


# 1) build a little data.frame of your event windows
event_windows <- Raw_plot %>%
  filter(Event %in% Events_to_keep) %>%
  mutate(
    xmin = time - 100,
    xmax = time + 2000
  )

P0 =  ggplot() +
  # shaded windows
  geom_rect( data = event_windows, inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Event),
    alpha = 0.5
  ) +
  # then your original lines and vlines
  geom_line(data = Raw_plot, aes(x = time, y = R_P,   color = 'Right'), size = 1.1 ) +
  geom_line(data = Raw_plot, aes(x = time, y = L_P,   color = 'Left'), size = 1.1 ) +
  geom_vline(data = event_windows, aes(xintercept = time), size = 1.3 , linetype = 'dashed') +

  theme_classic(base_size = 35) +
  ylim(2, 4.8) +
  labs(color = 'Eye', fill = 'Event',
       linetype = 'Events', x = 'Time (ms)', y = 'Pupil size') +
  scale_color_manual(values = c(Right = '#4A6274', Left = '#E2725A')) +
  theme(legend.position = 'top') +
  guides(
    color = guide_legend(override.aes = list(lwd = 20)),
    fill     = guide_legend(override.aes = list(size = 16, alpha = 1))
  ) +
  facet_grid(
    rows = vars(Subject),
    labeller = labeller(Subject = as_labeller(function(x) {
      paste0("S", sub("^Subject_", "", x))
    }))
  )



P0
ggsave("Paper/1_Raw.svg", P0,
        width = 28, height = 12, dpi = 300)





## Event fixes --------------------------------------------------------------------

List_of_subject_dataframes = list() # Empty list to be filled with dataframes

# Loop for each subject
for (sub in 1:length(csv_files)) {
  Raw_data = read.csv(csv_files[sub]) # Raw data
  Events = filter(Raw_data, Event %in% Events_to_keep) %>% # filter data
    mutate(TrialN = seq(n()))

  # Loop for each event
  for (trial in 1:nrow(Events)) {
    # Extract the information
    Event = Events[trial, ]$Event
    TrialN = Events[trial, ]$TrialN

    # Event onset information
    Onset = Events[trial, ]$time
    Onset_index = which.min(abs(Raw_data$time - Onset))

    # Find the rows to update based on pre post samples
    rows_to_update = seq(
      Onset_index - Pre_stim_samples,
      Onset_index + Post_stim_samples - 1
    )

    # Repeat the values of interest for all the rows
    Raw_data[rows_to_update, 'time'] = Time
    Raw_data[rows_to_update, 'Event'] = Event
    Raw_data[rows_to_update, 'TrialN'] = TrialN
  }

  # Filter only events of interest
  Trial_data = Raw_data %>%
    filter(Event %in% Events_to_keep)

  # Add daframe to list
  List_of_subject_dataframes[[sub]] = Trial_data
}

# Combine the list of dataframes into 1 dataframe
Trial_data = bind_rows(List_of_subject_dataframes)


# Pre-processing -----------------------------------------------------------------

## Filter Out Trials with all NA -----------------------------------------------------------------

Trial_data = Trial_data %>%
  group_by(Subject, TrialN) %>%
  filter(!all(is.na(L_P) & is.na(R_P))) %>%
  ungroup()


## Make PupillometryR Data -----------------------------------------------------------------
PupilR_data = make_pupillometryr_data(
  data = Trial_data,
  subject = Subject,
  trial = TrialN,
  time = time,
  condition = Event
)

### Plot ------------------------------------------------------------------
plot(PupilR_data, pupil = L_P, group = 'condition', geom = 'line') +
  theme_classic(base_size = 40) +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 20)))


## Regress Data -----------------------------------------------------------------

Regressed_data = regress_data(data = PupilR_data, pupil1 = L_P, pupil2 = R_P)


## Calculate Mean Pupil -----------------------------------------------------------------

Mean_data = calculate_mean_pupil_size(
  data = Regressed_data,
  pupil1 = L_P,
  pupil2 = R_P
)

### Plot --------------------------------------------------------------------

P1 = plot(Mean_data, pupil = mean_pupil, group = 'condition', geom = 'line') +
  theme_classic(base_size = 40) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(), # remove x-axis title
    axis.text.x = element_blank(), # remove x-axis labels
    axis.ticks.x = element_blank(), # remove x-axis ticks
    axis.title.y = element_blank())



## Lowpass Filter -----------------------------------------------------------------

filtered_data = filter_data(
  data = Mean_data,
  pupil = mean_pupil,
  filter = 'median',
  degree = 11
)

### Plot --------------------------------------------------------------------

plot(filtered_data, pupil = mean_pupil, group = 'condition', geom = 'line') +
  theme_bw(base_size = 40) +
  theme(
    legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(lwd = 20)))


## Downsample -----------------------------------------------------------------

NewHz = 20

timebinSize = 1 / NewHz

Downsampled_data = downsample_time_data(
  data = filtered_data,
  pupil = mean_pupil,
  timebin_size = timebinSize,
  option = 'median'
)

# Plot --------------------------------------------------------------------

P2 = plot(
  Downsampled_data,
  pupil = mean_pupil,
  group = 'condition',
  geom = 'line'
) +
  labs(y = 'Pupil size')+
  theme_classic(base_size = 40) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(), # remove x-axis title
    axis.text.x = element_blank(), # remove x-axis labels
    axis.ticks.x = element_blank()
  )

  ## Calculate Missing Data -----------------------------------------------------------------

  Missing_data = calculate_missing_data(Downsampled_data, mean_pupil)


## Clean Missing Data -----------------------------------------------------------------

Clean_data = clean_missing_data(
  Downsampled_data,
  pupil = mean_pupil,
  trial_threshold = .75,
  subject_trial_threshold = .75
)


## Detect Blinks -----------------------------------------------------------------

Blink_data = detect_blinks_by_velocity(
  Clean_data,
  mean_pupil,
  threshold = 0.1,
  extend_forward = 50,
  extend_back = 50
)


## Interpolate Data -----------------------------------------------------------------

Int_data = interpolate_data(
  data = Blink_data,
  pupil = mean_pupil,
  type = 'linear'
)


# Baseline correction -----------------------------------------------------

Base_data = baseline_data(
  data = Int_data,
  pupil = mean_pupil,
  start = -100,
  stop = 0
)

P3 = plot(Base_data, pupil = mean_pupil, group = 'condition', geom = 'line') +
  theme_classic(base_size = 40) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(y = "Baseline Corrected\nPupil Size")

# Remove the baseline
Final_data = subset_data(Base_data, start = 0)


### Final plot --------------------------------------------------------------

P4 = plot(Final_data, pupil = mean_pupil, group = 'condition',geom = 'line') +
  xlim(-100, 2000) +
  labs(x = 'Time (ms)') +
  theme_classic(base_size = 40) +
  theme(legend.position = 'bottom',
    legend.title = element_blank(),
    axis.title.y = element_blank())+
  guides(color = guide_legend(override.aes = list(lwd = 20)))

P4



# Plot toal ---------------------------------------------------------------




var_titles <- c("Mean\nPupil Size",
  "Smoothed\nPupil Size",
  "Baseline Corrected\nPupil Size")

(P1 / P2  / P4) +
  plot_annotation(
    tag_levels = list(var_titles),
    tag_prefix = "",
    tag_suffix = "",
    tag_sep    = " " ) &
  theme(
      plot.tag = element_text(
        family = "Times",                   # font family
        face   = "bold.italic",             # bold + italic
        size   = 30,                        # visible size
        hjust  = 0,                         # left-align tags
        vjust  = 1                          # nudge them down a bit
      )
    )

ggsave("Paper/PlotPupil.svg", 
       width = 28, height = 31, dpi = 300)






# MGCV Plot -------------------------------------------------------------------

library("mgcv")
library(easystats)


Final_data$Subject = as.factor(Final_data$Subject)
Final_data$Event = as.factor(Final_data$Event)

data_avg_time = Final_data %>%
  group_by(Event, time) %>%
  summarise(mean_pupil = mean(mean_pupil, na.rm=TRUE))




# Additive model
additive_model = bam(mean_pupil ~ Event
                     + s(time, by=Event, k=20)
                     + s(time, Subject, bs='fs', m=1),
                     data=Final_data)



# Data grid
datagrid = get_datagrid(additive_model, length = 100, include_random = T)

# Estimate expectation and uncertainty (Predicted and SE)
Est = estimate_expectation(additive_model, datagrid, exclude=c("s(time,Subject)"))


# Plot predictions with confidence intervals and the observed data
P5 = ggplot() +
  # Real data line
  geom_line(data = data_avg_time, aes(x=time, y=mean_pupil, color=Event), size=1.5) +
  
  # Predicted ribbons
  geom_ribbon(data = Est, aes(x=time, ymin = CI_low, ymax = CI_high,
                              fill = Event), alpha = 0.2) +
  
  # Predicted lines
  geom_line(data = Est, aes(x=time, y=Predicted, color=Event), lwd=4, linetype = "dashed") +
  
  labs(y = "Pupil Size", x = "Time (ms)") +
  theme_classic(base_size = 40) +
  theme(legend.position = 'bottom',
    legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(lwd = 20))) 

ggsave("Paper/MGCV.svg", P5,
       width = 28, height = 18, dpi = 300)




# TtoalPlot --------------------------------------------------------------


PP4 = P4+
  theme(legend.position = 'none')

PP5 =  P5 + 
  theme(plot.margin = unit(c(1, 1, 1, 8), "cm"),
  legend.position = 'none')



(P1 / P2 / PP4) | PP5

ggsave("Paper/PlotPupilFinal.svg", 
       width = 28, height = 12, dpi = 300)

