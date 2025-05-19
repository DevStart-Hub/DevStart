# In your main analysis script:
library(tidyverse)
library(signal)
library(truncnorm)
library(datawizard) # For change_scale function used in scale_signal
library(patchwork)

# Source my functions
source("Simulations/PupilSimulationFunctions.R")


# Read original data -------------------------------------------------------------------

pup = vroom::vroom("resources/FromFixationToData/DATA/RAW/Adult1.csv") |> 
  mutate(time = (time - min(time))+1) |> 
  rowwise() |>
  mutate(pupil = mean(c(R_P, L_P), rm.na = TRUE),
          GazeX = mean(c(R_X, L_X), rm.na = TRUE),
          GazeY = mean(c(R_Y, L_Y), rm.na = TRUE)) |> 
  ungroup()

Event = pup |>
  dplyr::filter(
    !is.na(Event) & Event != "Preferential" ) |>
  mutate(indexes = row_number()) |>
  select(time, indexes, Event)

# Plot real data
OriginalData = ggplot(pup, aes(x = time, y = R_P)) +
  geom_line(color = 'red') +
  geom_line(aes(y = L_P), color = 'blue') +
  geom_vline(
    data = Event,
    aes(xintercept = time, color = Event),
    linetype = 'dashed'
  ) +
  labs(title = "Original Pupil Data", y= "Pupil") +
  theme_classic(base_size = 20) +
  ylim(0, 6)+
  theme(legend.position = 'none' )




# Event specific settings ------------------------------------------------

# Define parameter mappings for each event type
param_map <- data.frame(
  EventType = c("Circle", "Square", "Reward", "NoReward", "Fixation"),
  npar_mean = c(6, 2, 2, 2, 2),
  tmax_mean = c(1300, 1300, 2300, 2300, 1000), 
  reduction_mean = c(0.8, 0.02, 0.2, 0.1, 00.1),
  stringsAsFactors = FALSE
)

# Compute standard deviations as 5% of the mean
param_map <- param_map %>%
  mutate(
    npar_sd = npar_mean * 0.05,
    tmax_sd = tmax_mean * 0.05,
    reduction_sd = reduction_mean * 0.05
  )


# Join Event with parameters
Event <- Event %>%
  left_join(param_map, by = c("Event" = "EventType"))


# Simulate data ----------------------------------------------------------

sim_data = list()
sub=1

for (sub in 1:6) {
  sim_data[[sub]] <- generate_pupil_data(
    fs = 300,
    segment_length = max(pup$time, na.rm = TRUE), # T should be replaced by your desired total duration in ms
    event_onsets = Event$time, # Using the 'time' column from Event
    event_names = Event$Event,

    prf_npar_means = Event$npar_mean, # Use the per-event npar means
    prf_npar_sds = Event$npar_sd, # Use the per-event npar sds

    prf_tmax_means = Event$tmax_mean, # Use the per-event tmax means
    prf_tmax_sds = Event$tmax_sd, # Use the per-event tmax sds

    # The 'scale_evoked_means' and 'scale_evoked_sds' parameters can be mapped to 'reduction_mean' and 'reduction_sd'
    # if you want your amplitude scaling to come from reduction. This sets the amplitude scaling factor per event.
    scale_evoked_means = Event$reduction_mean,
    scale_evoked_sds = Event$reduction_sd,

    baseline_lowpass = 0.6,
    num_spurious_events = sample(3:6, 1),
    noise_amp = 0.1,
    scale_signal = c(2.6, 4.4)
  )

  ### Add blinks
  sim_data[[sub]]$pupilBlinks <- insert_blinks(
    sim_data[[sub]]$pupil,
    fs = 300,
    num_blinks = 37,
    min_duration = 100,
    max_duration = 300,
    phase_length = 40,
    phase_mean = 0.92,
    phase_sd = 0.12,
    min_gap = 400
  )

  # Add subject info to the dataframe
  sim_data[[sub]]$Subject <- paste('Subject', sub, sep = '_')

  ### Difference between Left and Right eyes
  sim_data[[sub]]$L_P = sim_data[[sub]]$pupilBlinks
  sim_data[[sub]]$R_P = sim_data[[sub]]$pupilBlinks
  sim_data[[sub]]$L_P = sim_data[[sub]]$pupilBlinks *
    rnorm(1, mean = 1, sd = 0.05)
  sim_data[[sub]]$R_P = sim_data[[sub]]$pupilBlinks *
    rnorm(1, mean = 0.97, sd = 0.05)

  ### Add NAs segments independently for each eye
  sim_data[[sub]]$L_P <- insert_na_segments(
    sim_data[[sub]]$L_P,
    fs = 300,
    num_segments = sample(2:3, 1),
    min_duration = 1000,
    max_duration = 2500
  )
  sim_data[[sub]]$R_P <- insert_na_segments(
    sim_data[[sub]]$R_P,
    fs = 300,
    num_segments = sample(2:3, 1),
    min_duration = 1000,
    max_duration = 2500
  )
  
  ## Add events column in the data ------------------------------------------

  Event2 = sim_data[[sub]] %>% dplyr::filter(Event %in% c('Circle', 'Square'))

  if (sub == 1) {
    sim_data[[sub]][
      sim_data[[sub]]$time > (Event2[1, ]$time - 200) &
        sim_data[[sub]]$time < (Event2[1, ]$time + 2500),
      c('L_P', 'R_P')
    ] = NA
  } else if (sub == 6) {
    sim_data[[sub]][
      sim_data[[sub]]$time > (Event2[10, ]$time - 200) &
        sim_data[[sub]]$time < (Event2[10, ]$time + 2500),
      c('L_P', 'R_P')
    ] = NA
  } else if (sub == 3) {
    sim_data[[sub]][
      sim_data[[sub]]$time > (Event2[3, ]$time - 50) &
        sim_data[[sub]]$time < (Event2[3, ]$time + 1800),
      c('L_P', 'R_P')
    ] = NA
    sim_data[[sub]][
      sim_data[[sub]]$time > (Event2[6, ]$time - 50) &
        sim_data[[sub]]$time < (Event2[6, ]$time + 1800),
      c('L_P', 'R_P')
    ] = NA
    sim_data[[sub]][
      sim_data[[sub]]$time > (Event2[8, ]$time - 50) &
        sim_data[[sub]]$time < (Event2[8, ]$time + 1800),
      c('L_P', 'R_P')
    ] = NA
  }

  write.csv(sim_data[[sub]],
            paste("/resources/Pupillometry/Raw/",
            "Pupil_Subject", sub, ".csv", sep= ''))
}




# Combine the data -------------------------------------------------------

Pup = do.call(rbind, sim_data) |>
  select(Subject, Time_ms, L_P, R_P, Event)


# Final plot -------------------------------------------------------------

# Subject specific plots
SimulatedData1 =
  Pup |>
  dplyr::filter(Subject == 'Subject_1') |>
  ggplot(aes(x = Time_ms, y = L_P, color = Subject)) +
  geom_line(color = 'red') +
  geom_line(aes(y = R_P), color = 'blue') +
  geom_vline(
    data = Event,
    aes(xintercept = time, color = Event),
    linetype = 'dashed'
  ) +
  labs(title = "Simulated Pupil Data", y = "Pupil") +
  theme_classic(base_size = 20) +
  ylim(0, 6) +
  theme(legend.position = 'none')


SimulatedAll = ggplot(Pup, aes(x = Time_ms, y = L_P, color = Subject)) +
  geom_line(color = 'red') +
  geom_line(aes(y = R_P), color = 'blue') +
  geom_vline(
    data = Event,
    aes(xintercept = time, color = Event),
    linetype = 'dashed'
  ) +
  labs(y = "Pupil") +
  theme_minimal(base_size = 20) +
  facet_wrap(~Subject, ncol = 1) +
  ylim(0, 6) +
  theme(legend.position = 'top')

 
(OriginalData + SimulatedData1) / SimulatedAll
