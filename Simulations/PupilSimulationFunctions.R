#%% LbrLibrariesaries and functions
library(tidyverse)
library(signal)
library(truncnorm)


#%% SECTION 1: CORE PUPILLOMETRY FUNCTIONS

#' Pupil Response Function (PRF) Kernel
#'
#' Creates a gamma-shaped kernel that models the canonical pupillary response to events.
#' The shape follows a gamma function, which has been shown to fit the pupillary
#' response well in empirical studies.
#'
#' @param duration Duration of the kernel in milliseconds (default: 4000)
#' @param fs Sampling rate in Hz (default: 1000)
#' @param npar Shape parameter that controls the rise time (default: 10.1)
#' @param tmax Time-to-peak parameter in milliseconds (default: 930.0)
#' @return Vector containing the pupil response kernel
pupil_kernel <- function(
  duration = 4000,
  fs = 1000,
  npar = 10.1,
  tmax = 930.0
) {
  # Calculate number of samples based on duration and sampling rate
  n <- as.integer(duration / 1000 * fs)

  # Create time vector
  t <- seq(0, duration, length.out = n)

  # Calculate normalization factor (maximum amplitude of the unnormalized kernel)
  hmax <- exp(-npar) * tmax^npar

  # Create gamma-shaped kernel normalized to have a peak of 1
  h <- (t^npar * exp(-npar * t / tmax)) / hmax

  return(h)
}


#' Generate Simulated Pupil Data
#'
#' Creates realistic pupil time series data with:
#'   - Baseline low-frequency fluctuations
#'   - Event-related pupil responses with customizable parameters
#'   - Optional spurious (non-event-related) pupil dilations
#'   - Gaussian noise
#'
#' Special behavior:
#'   - If for an event both `scale_evoked_means[i] == 0` and `scale_evoked_sds[i] == 0`,
#'     that event is still marked in the `Event` column but elicits **no** evoked response
#'     (its weight is forced to zero, avoiding NA propagation).
#'
#' @param fs Sampling rate in Hz (default: 1000)
#' @param segment_length Total duration of simulated data in milliseconds (default: 30000)
#' @param event_onsets Vector of event onset times in milliseconds (default: c(5000,10000,15000))
#' @param event_names Vector of event names or single name to auto-enumerate (default: "Event")
#' @param scale_signal Optional global scaling of signal as c(min, max) (default: NULL)
#' @param noise_amp Amplitude of additive Gaussian noise (default: 0.1)
#' @param baseline_lowpass Cutoff frequency for baseline fluctuations in Hz (default: 0.2)
#' @param num_spurious_events Number of random non-event-related dilations (default: 0)
#' @param scale_evoked_means Amplitude scaling factor(s); single value or vector matching events
#'                          (default: 0.5).  A value of 0 with matching `scale_evoked_sds=0`
#'                          produces no evoked response for that event.
#' @param scale_evoked_sds Standard deviation(s) for amplitude; single value or vector matching events
#'                        (default: 0.05).  Must be 0 to explicitly suppress response when mean=0.
#' @param prf_npar_means Shape parameter(s) for event(s); single value or vector matching events
#'                     (default: 10.1)
#' @param prf_npar_sds Standard deviation(s) for shape parameter(s); single value or vector matching events
#'                   (default: 0.5)
#' @param prf_tmax_means Time-to-peak parameter(s) in ms; single value or vector matching events
#'                     (default: 930)
#' @param prf_tmax_sds Standard deviation(s) for time-to-peak; single value or vector matching events
#'                   (default: 50)
#' @return Data frame with columns: Time_ms, pupil, Event
generate_pupil_data <- function(
  fs = 300,
  segment_length = 30000,
  event_onsets = c(5000, 10000, 15000),
  event_names = "Event",
  scale_signal = NULL,
  noise_amp = 0.1,
  baseline_lowpass = 0.2,
  num_spurious_events = 0,
  scale_evoked_means = 0.5,
  scale_evoked_sds = 0.05,
  prf_npar_means = 10.1,
  prf_npar_sds = 0.5,
  prf_tmax_means = 930,
  prf_tmax_sds = 50
) {
  nevents <- length(event_onsets)
  # Validate lengths
  if (
    !(length(prf_npar_means) %in% c(1, nevents)) ||
    !(length(prf_npar_sds)   %in% c(1, nevents)) ||
    !(length(prf_tmax_means) %in% c(1, nevents)) ||
    !(length(prf_tmax_sds)   %in% c(1, nevents)) ||
    !(length(scale_evoked_means) %in% c(1, nevents)) ||
    !(length(scale_evoked_sds)   %in% c(1, nevents)) ||
    !(length(event_names) %in% c(1, nevents))
  ) stop("Lengths of parameter vectors must either be 1 or match number of events")

  # Expand singletons
  if (length(prf_npar_means)==1) prf_npar_means <- rep(prf_npar_means, nevents)
  if (length(prf_npar_sds)==1)   prf_npar_sds   <- rep(prf_npar_sds, nevents)
  if (length(prf_tmax_means)==1) prf_tmax_means <- rep(prf_tmax_means, nevents)
  if (length(prf_tmax_sds)==1)   prf_tmax_sds   <- rep(prf_tmax_sds, nevents)
  if (length(scale_evoked_means)==1) scale_evoked_means <- rep(scale_evoked_means, nevents)
  if (length(scale_evoked_sds)==1)   scale_evoked_sds   <- rep(scale_evoked_sds, nevents)
  if (length(event_names)==1) event_names <- paste0(event_names, "_", seq_len(nevents))

  # Time vector
  n <- as.integer(segment_length/1000 * fs)
  time <- seq(1, segment_length, 1000/fs)[1:n]

  # Baseline fluctuations
  slack <- as.integer(0.5 * n)
  base_noise <- runif(n + 2*slack)
  bfilt <- butter(2, baseline_lowpass/(fs/2), type="low")
  x0_full <- filtfilt(bfilt, base_noise)
  x0 <- x0_full[(slack+1):(slack+n)]
  x0 <- (x0 - mean(x0))/sd(x0)

  # Sample shape & timing
  npars <- tmaxs <- numeric(nevents)
  for (i in seq_len(nevents)) {
    npars[i] <- if (prf_npar_sds[i]==0) prf_npar_means[i] else rnorm(1, prf_npar_means[i], prf_npar_sds[i])
    tmaxs[i] <- if (prf_tmax_sds[i]==0) prf_tmax_means[i] else rnorm(1, prf_tmax_means[i], prf_tmax_sds[i])
  }
  maxdur <- max(tmaxs)*5

  # Sample amplitudes with zero‐response guard
  delta_weights <- numeric(nevents)
  for (i in seq_len(nevents)) {
    if (scale_evoked_means[i]==0 && scale_evoked_sds[i]==0) {
      delta_weights[i] <- 0
    } else {
      a <- (0 - scale_evoked_means[i]) / scale_evoked_sds[i]
      delta_weights[i] <- rtruncnorm(1, a=a, b=Inf,
                                     mean=scale_evoked_means[i],
                                     sd=scale_evoked_sds[i])
    }
  }

  sy <- numeric(n)
  event_column <- rep(NA_character_, n)

  # Add evoked responses
  for (i in seq_along(event_onsets)) {
    kernel <- pupil_kernel(duration=maxdur, fs=fs,
                           npar=npars[i], tmax=tmaxs[i])
    onset_ix <- as.integer(event_onsets[i]/1000 * fs) + 1
    if (onset_ix <= n) event_column[onset_ix] <- event_names[i]
    end_ix <- min(n, onset_ix + length(kernel)-1)
    sy[onset_ix:end_ix] <- sy[onset_ix:end_ix] +
      kernel[1:(end_ix-onset_ix+1)] * delta_weights[i]
  }

  # Spurious events if requested
  if (num_spurious_events>0) {
    sp_ix <- sample(seq_len(n-100), num_spurious_events)
    mean_npar <- mean(prf_npar_means)
    mean_tmax <- mean(prf_tmax_means)
    sp_npars <- rnorm(num_spurious_events, mean_npar, mean(prf_npar_sds))
    sp_tmaxs <- rnorm(num_spurious_events, mean_tmax, mean(prf_tmax_sds))
    sp_weights <- rtruncnorm(num_spurious_events,
                              a=(0-mean(scale_evoked_means))/mean(scale_evoked_sds),
                              b=Inf,
                              mean=mean(scale_evoked_means),
                              sd=mean(scale_evoked_sds))
    for (j in seq_along(sp_ix)) {
      k_sp <- pupil_kernel(duration=maxdur, fs=fs,
                           npar=sp_npars[j], tmax=sp_tmaxs[j])
      st <- sp_ix[j]
      en <- min(n, st+length(k_sp)-1)
      sy[st:en] <- sy[st:en] + k_sp[1:(en-st+1)] * sp_weights[j]
    }
  }

  # Combine, noise, optional scaling
  sy <- x0 + sy + rnorm(n, 0, noise_amp)
  if (!is.null(scale_signal) && length(scale_signal)==2) {
    sy <- datawizard::change_scale(sy, scale_signal)
  }

  data.frame(time=time, pupil=sy, Event=event_column)
}



#%% SECTION 2: ARTIFACT SIMULATION FUNCTIONS

#' Insert Realistic Blinks into Pupil Data
#'
#' Adds NA values to simulate eye blinks with realistic characteristics:
#'   - Closing phase (decreasing pupil size) - optional
#'   - Complete blink (NA values)
#'   - Reopening phase (increasing pupil size) - optional
#'
#' @param pupil Vector of pupil measurements
#' @param fs Sampling rate in Hz (default: 300)
#' @param num_blinks Number of blinks to insert (default: 30)
#' @param min_duration Minimum blink duration in ms (default: 100)
#' @param max_duration Maximum blink duration in ms (default: 300)
#' @param phase_length Duration of closing/reopening phases in ms (default: 40, set to 0 for immediate transitions)
#' @param phase_mean Mean proportion of baseline for closing/reopening (default: 0.95)
#' @param phase_sd SD of proportion for closing/reopening (default: 0.1)
#' @param min_gap Minimum spacing between blinks in ms (default: 500)
#' @return Pupil vector with blinks inserted as NA values
insert_blinks <- function(
  pupil,
  fs = 300,
  num_blinks = 30,
  min_duration = 100,
  max_duration = 300,
  phase_length = 40,
  phase_mean = 0.95,
  phase_sd = 0.1,
  min_gap = 500
) {
  # Determine if we should use gradual transitions or immediate jumps
  # When phase_length > 0, create realistic transitions between data and blinks
  # When phase_length = 0, insert blinks as immediate data gaps (no transition)
  has_transition_phases <- phase_length > 0

  # Convert time values from milliseconds to sample counts based on sampling rate
  # If phase_length is 0, phase_samples will be 0 (no transition samples)
  phase_samples <- floor(phase_length / 1000 * fs)

  # Calculate the total duration of the recording in milliseconds
  total_time_ms <- length(pupil) / fs * 1000

  # Distribute blinks evenly throughout the recording by dividing into segments
  # We exclude the first and last 1000ms to avoid placing blinks at the very edges
  segment_length <- floor((total_time_ms - 2000) / num_blinks)

  # Create a vector to store the onset times for each blink (in milliseconds)
  blink_onsets <- numeric(num_blinks)

  # Generate realistic onset times for each blink, one per segment
  # This ensures blinks are distributed throughout the recording instead of clustered
  for (i in 1:num_blinks) {
    # Calculate the start and end of the current segment in milliseconds
    segment_start <- 1000 + (i - 1) * segment_length
    segment_end <- segment_start + segment_length

    # Place the blink randomly within this segment, respecting constraints
    if (has_transition_phases) {
      # With transition phases, we need to leave room for:
      # - phase_length before the blink (for closing)
      # - min_gap before that (minimum time between blinks)
      # - phase_length after the blink (for reopening)
      # - max_duration for the blink itself
      blink_time <- runif(
        1,
        min = segment_start + phase_length + min_gap,
        max = segment_end - phase_length - max_duration
      )
    } else {
      # With immediate transitions, we only need space for:
      # - min_gap (minimum time between blinks)
      # - max_duration for the blink itself
      blink_time <- runif(
        1,
        min = segment_start + min_gap,
        max = segment_end - max_duration
      )
    }
    # Store the onset time for this blink
    blink_onsets[i] <- blink_time
  }

  # Process each blink and modify the pupil data accordingly
  for (b in blink_onsets) {
    # Convert the blink onset time from milliseconds to sample index
    start_samp <- floor(b / 1000 * fs) + 1

    # Select a random duration for this blink between min and max
    # (each blink has a different duration for realism)
    blink_dur_samp <- floor(runif(1, min_duration, max_duration) / 1000 * fs)
    end_samp <- start_samp + blink_dur_samp

    # Skip this blink if it would extend beyond the valid data range
    if (has_transition_phases) {
      # With transitions, we need room for the closing and reopening phases too
      if (
        (start_samp - phase_samples) < 1 ||
          (end_samp + phase_samples) > length(pupil)
      ) {
        next # Skip to the next blink
      }
      # Get the baseline pupil value just before the closing phase starts
      baseline_val <- pupil[start_samp - phase_samples]
      # Skip if the baseline value is not a valid number (e.g., NA, NaN, Inf)
      if (!is.finite(baseline_val)) next
    } else {
      # Without transitions, we only need to check if the blink itself fits
      if (start_samp < 1 || end_samp > length(pupil)) {
        next # Skip to the next blink
      }
    }

    # If using transition phases, create the closing phase (pupil gradually decreases)
    if (has_transition_phases) {
      # Generate random scaling factors that determine how much pupil size
      # decreases during closing and increases during reopening
      # These factors vary for each blink to create natural variation
      closing_factor <- rtruncnorm(
        1,
        a = 0, # Lower bound: never negative
        b = 1, # Upper bound: never larger than baseline
        mean = phase_mean, # Target is phase_mean (default: 95% of baseline)
        sd = phase_sd # With some random variation
      )
      reopening_factor <- rtruncnorm(
        1,
        a = 0,
        b = 1,
        mean = phase_mean,
        sd = phase_sd
      )

      # Get the baseline pupil value before the blink starts
      baseline_val <- pupil[start_samp - phase_samples]

      # Create the closing phase: pupil size gradually decreases from baseline to a smaller value
      # This simulates the eyelid beginning to close
      pupil[(start_samp - phase_samples):(start_samp - 1)] <-
        seq(
          from = baseline_val, # Start at the baseline value
          to = baseline_val * closing_factor, # End at a scaled-down value
          length.out = phase_samples # Create a smooth transition over phase_samples
        )
    }

    # Set the main blink phase to NA (representing missing data when eye is fully closed)
    pupil[start_samp:end_samp] <- NA

    # If using transition phases, create the reopening phase (pupil gradually increases)
    if (has_transition_phases) {
      # Calculate the sample indices for the reopening phase
      reopen_start <- end_samp + 1
      reopen_end <- reopen_start + phase_samples - 1

      # Only create reopening phase if it fits within the data range
      if (reopen_end <= length(pupil)) {
        # Create the reopening phase: pupil size gradually increases back to baseline
        # This simulates the eyelid opening again
        pupil[reopen_start:reopen_end] <-
          seq(
            from = baseline_val * reopening_factor, # Start at a scaled-down value
            to = baseline_val, # Return to the baseline value
            length.out = phase_samples # Create a smooth transition over phase_samples
          )
      }
    }
  }

  # Return the modified pupil vector with blinks inserted
  return(pupil)
}


#' Insert Missing Data Segments
#'
#' Adds blocks of NA values to simulate periods of missing data
#' (e.g., from tracker loss, participant movement, etc.)
#'
#' @param pupil Vector of pupil measurements
#' @param fs Sampling rate in Hz (default: 300)
#' @param num_segments Number of NA segments to insert (default: 10)
#' @param min_duration Minimum segment duration in ms (default: 1000)
#' @param max_duration Maximum segment duration in ms (default: 2000)
#' @return Pupil vector with NA segments inserted
insert_na_segments <- function(
  pupil,
  fs = 300,
  num_segments = 10,
  min_duration = 1000,
  max_duration = 2000
) {
  # Get total length of pupil data
  total_samples <- length(pupil)

  # Divide recording into equal segments, avoiding the edges
  segment_samples <- floor((total_samples - 2000) / num_segments)

  # Insert one NA block in each segment
  for (i in 1:num_segments) {
    # Maximum valid start position that ensures the gap fits in this segment
    max_start <- segment_samples - max_duration / 1000 * fs

    # Skip if segment is too small to fit a gap
    if (max_start < 1) next

    # Random position within this segment
    start_pos <- sample(1:max_start, 1) + (i - 1) * segment_samples

    # Random duration between min and max (ms to samples)
    na_duration <- floor(runif(1, min_duration, max_duration) / 1000 * fs)

    # Ensure we don't exceed data length
    end_pos <- min(start_pos + na_duration - 1, total_samples)

    # Insert the NA block
    pupil[start_pos:end_pos] <- NA
  }

  return(pupil)
}

