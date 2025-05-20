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



# Plot raw ---------------------------------------------------------------
Raw_plot = read.csv(csv_files[1])
P0 = ggplot(Raw_plot, aes(x = time, y = R_P)) +
      geom_line(aes(y = R_P, color = 'Pupil Right'), lwd = 1.2) +
      geom_line(aes(y = L_P, color = 'Pupil Left'), lwd = 1.2) +
      geom_vline(data = Raw_plot |> dplyr::filter(Event %in% Events_to_keep ), aes(xintercept = time, linetype = Event), lwd = 1.3) +
      
      theme_bw(base_size = 35) +
      ylim(1, 6) +
      labs(color= 'Signal', y='Pupil size')+
      scale_color_manual(
        values = c('Pupil Right' = '#4A6274', 'Pupil Left' = '#E2725A') )  +
      theme(
        legend.position = 'bottom'  ) +
      guides(
        color = guide_legend(override.aes = list(lwd = 20)),
        linetype = guide_legend(override.aes = list(lwd = 1.2))
      )

P0





# rawplot all ------------------------------------------------------------

Raw_plot = vroom::vroom(csv_files) |> 
  filter(time<78000)

P0 = ggplot(Raw_plot, aes(x = time, y = R_P)) +
      geom_line(aes(y = R_P, color = 'Right'), lwd = 1.2) +
      geom_line(aes(y = L_P, color = 'Left'), lwd = 1.2) +
      geom_vline(data = Raw_plot |> dplyr::filter(Event %in% Events_to_keep ), aes(xintercept = time, linetype = Event), lwd = 1.3) +
      
      theme_classic(base_size = 35) +
      ylim(2, 4.8) +
      labs(color= 'Eye', y ="", linetype = 'Events', x = 'Time (ms)')+
      scale_color_manual(
        values = c('Right' = '#4A6274', 'Left' = '#E2725A') )  +
      theme(
        legend.position = 'bottom',
        strip.text = element_text(size = 11)) +
      guides(
        color = guide_legend(override.aes = list(lwd = 10)),
        linetype = guide_legend(override.aes = list(lwd = 2))
      )+
      facet_grid(rows = vars(Subject)) 
      # facet_wrap(~Subject, ncol = 1)

P0




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
  theme_classic(base_size = 45) +
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
  theme_classic(base_size = 45) +
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
  theme_bw(base_size = 45) +
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
  theme_classic(base_size = 45) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(), # remove x-axis title
    axis.text.x = element_blank(), # remove x-axis labels
    axis.ticks.x = element_blank(), # remove x-axis ticks
    axis.title.y = element_blank())

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
  theme_classic(base_size = 45) +
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
  theme_classic(base_size = 45) +
  theme(legend.position = 'bottom',
    legend.title = element_blank(),
    axis.title.y = element_blank())+
  guides(color = guide_legend(override.aes = list(lwd = 10)))

P4



# Plot toal ---------------------------------------------------------------

# P1 /  P2 /  P3 /  P4 +
#   plot_annotation(
#     tag_levels = c("Mean\nPupil Size","Smoothed\nPupil Size", "Windowed\nPupil Size", "Baseline Corrected\nPupil Size"),
#     theme = theme(
#       plot.tag = element_text(size = 18, face = "italic"),
#       plot.tag.position = c(0, 1)
#     )
#   )


var_titles <- c("Raw\nPupil Size",
  '',
  "Mean\nPupil Size",
  "Smoothed\nPupil Size",
  "Baseline Corrected\nPupil Size")

(P0 / (P1 / P2  / P4)) +
  plot_layout(height = c(1, 1.5))+
  plot_annotation(
    tag_levels = list(var_titles),
    tag_prefix = "",
    tag_suffix = "",
    tag_sep    = " " ) &
  theme(
      plot.tag = element_text(
        family = "Times",                   # font family
        face   = "bold.italic",             # bold + italic
        size   = 21,                        # visible size
        hjust  = 0,                         # left-align tags
        vjust  = 1                          # nudge them down a bit
      )
    )


# make a separator plot
sep <- ggplot() +
  geom_hline(yintercept = 0, size = 1, linetype = 'dashed') +   # adjust size for thickness
  theme_void()                              # no axes, no background

# now stitch P0, sep, and the rest together
layout <- (P0 / sep / (P1 / P2 / P4)) +
  plot_layout(heights = c(1, .1, 1.5)) +  # sep is only 2% of total height
  plot_annotation(
    tag_levels = list(var_titles),
    tag_prefix = "",
    tag_suffix = "",
    tag_sep    = " "
  ) & theme(
    plot.tag = element_text(
      family = "Times",
      face   = "bold.italic",
      size   = 21,
      hjust  = 0,
      vjust  = 1
    )
  )

ggsave("Paper/PlotPupil.png", layout,
       width = 20, height = 28, dpi = 300)



