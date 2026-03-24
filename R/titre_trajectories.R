#' Simulate Antibody Titre Trajectories for an Individual
#'
#' Generates a simulated time course of true and measured antibody titres
#' for an individual, given sampling times, infection/vaccination events,
#' and baseline serostatus. The function applies a rise-and-decay model
#' with measurement error and limit of detection.
#'
#' @param sampling_times Numeric vector of sampling times (in days) relative to
#'  the start of the study.
#' @param treatment_group Integer (0 or 1) indicating treatment group.
#'  (0 = Placebo, 1 = Vaccine).
#' @param subject_id Identifier for the subject (numeric or character).
#' @param perm_rise Numeric scalar. Permanent rise in titre following an
#'  infection. It must have the same length as the number of infections.
#' @param temp_rise Numeric scalar. Initial temporary rise in titre
#'  following an infection. It must have the same length as the number of
#'  infections.
#' @param temp_decay Numeric scalar. Decay rate of the temporary titre rise.
#' @param meas_sd A non-negative real value that corresponds to the standard
#'   deviation of a normal measurement model of log titres.
#' @param infection_times Numeric vector of times of infection events relative
#'   to enrolment. Defaults to `NULL` if no infections.
#' @param vac_times Numeric vector of times of vaccination events relative to
#'   enrolment. Defaults to `NULL` if no vaccinations.

#'
#' @returns A `data.frame` with columns:
#' \describe{
#'   \item{subject_id}{Subject identifier.}
#'   \item{time}{Sampling time.}
#'   \item{true}{Simulated true log-titre at each time.}
#'   \item{meas}{Observed log-titre with measurement error and LOD applied.}
#' }
#' @export
#'
#' @examples
#' sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#' simulate_titre_trajectory(sampling_times   = sampling_times,
#'                           treatment_group  = 0, # Placebo
#'                           subject_id       = 1,
#'                           infection_times  = c(300, 900),
#'                           perm_rise        = 6,
#'                           temp_rise        = 2,
#'                           temp_decay       = 0.003,
#'                           meas_sd          = 0.3)
simulate_titre_trajectory <- function(sampling_times,
                                      treatment_group,
                                      subject_id,
                                      perm_rise  = NULL,
                                      temp_rise  = NULL,
                                      temp_decay = NULL,
                                      meas_sd,
                                      infection_times = NULL,
                                      vac_times       = NULL)
{
  if (length(perm_rise) != length(infection_times)) {
    stop("Vectors 'perm_rise' and 'infection_times' must have the same length.")
  }

  if (length(temp_rise) != length(infection_times)) {
    stop("Vectors 'temp_rise' and 'infection_times' must have the same length.")
  }

  if (length(temp_decay) != length(temp_decay)) {
    stop("Vectors 'temp_decay' and 'infection_times' must have the same length.")
  }


  df <- data.frame(subject_id = subject_id,
                   time       = sampling_times)

  log_titre_vals <- simulate_true_titre_DENV(
    sampling_times  = sampling_times,
    perm_rise       = perm_rise,
    temp_rise       = temp_rise,
    temp_decay      = temp_decay,
    treatment_group = treatment_group,
    infection_times = infection_times,
    vac_times       = vac_times)

  df$true <- log_titre_vals

  df$meas <- simulate_observed_titres(log_titre_vals,
                                      measurement_error = meas_sd,
                                      LOD = 1)

  df
}

#' Simulate Observed Titres with Gaussian Measurement Error
#'
#' This function simulates observed titre values by adding Gaussian (normal)
#' measurement error to true underlying titres. Any observed value falling below
#' a specified limit of detection (LOD) is reported as zero.
#'
#' @param true_titre Numeric vector. True underlying titre values.
#' @param measurement_error Numeric scalar. Standard deviation of the Gaussian measurement error.
#' @param LOD Numeric scalar. Limit of detection below which observed titres are set to zero.
#'
#' @return A numeric vector of observed titres after adding measurement error
#'   and applying the LOD rule.
#'
#' @details
#' Observations are generated as:
#' \deqn{obs = true\_titre + \epsilon,\quad \epsilon \sim \mathcal{N}(0, \sigma^2)}
#' where `measurement_error` is the standard deviation \eqn{\sigma} of the noise.
#' Any value less than `LOD` is returned as 0, mimicking assay detection limits.
#'
#' @examples
#' set.seed(123)
#' true <- c(5, 10, 2, 1)
#' simulate_titre_observation(true_titre = true,
#'                            measurement_error = 1.5,
#'                            LOD = 1)
#'
#' @export
simulate_observed_titres <- function(true_titre, measurement_error, LOD)
{
  obs_titre <- true_titre + stats::rnorm(n    = length(true_titre),
                                         mean = 0,
                                         sd   = measurement_error)

  obs_titre <- ifelse(obs_titre < LOD, 0, obs_titre)

  obs_titre
}

#' Simulate true titres
#'
#' This function generates simulated log-transformed antibody titre values
#' at specified sampling times, given infection (and potentially vaccination)
#' events and immune response parameters. This function assumes that the
#' individual is seronegative at time 0. Time 0 can represent the start of a
#' study or the time at which the individual became susceptible.
#'
#' @inheritParams simulate_titre_trajectory
#'
#' @return A numeric vector of length equal to \code{sampling_times},
#'  containing simulated log antibody titre values.
#'
#' @export
#'
#' @examples
#' sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#' simulate_true_titre(sampling_times   = sampling_times,
#'                     perm_rise        = 6,
#'                     temp_rise        = 2,
#'                     temp_decay       = 0.003,
#'                     treatment_group  = 0, # Placebo
#'                     infection_times  = c(300, 900))
simulate_true_titre_DENV <- function(sampling_times,
                                     perm_rise,
                                     temp_rise,
                                     temp_decay,
                                     treatment_group,
                                     infection_times = NULL,
                                     vac_times       = NULL)
{
  n_inf     <- length(infection_times)
  n_samples <- length(sampling_times)

  if(length(infection_times) == 0) return(rep(0, n_samples))

  titre_matrix <- matrix(0,
                         nrow = length(infection_times),
                         ncol = length(sampling_times))

  for(i in seq_len(n_inf))
  {
    idx             <- which(sampling_times >= infection_times[[i]])
    times_after_inf <- sampling_times[idx]

    sim_titre <- titre_decay_floor(
      par_alpha = perm_rise[[i]],
      par_beta  = temp_rise[[i]],
      par_delta = temp_decay[[i]],
      time      = times_after_inf - infection_times[[i]])

    titre_matrix[i, idx] <- sim_titre
  }

  log_titre_vals <- colSums(titre_matrix)


  log_titre_vals
}


#' Simulate antibody titre trajectories from birth
#'
#' @inheritParams simulate_titre_trajectory
#' @param inf_times_sbs Numeric vector of infection times (in days) relative to
#'  day the individual became susceptible.
#' @param age Numeric. Age of the individual at the time of the first sample
#'  (in years).
#'
#' @returns A data frame containing simulated titre values across the specified
#'  sampling times.
#' @export
#'
#' @examples
#' simulate_titres_seropositive(
#'   inf_times_sbs   = c(1059.646, 2743.543, 4152.091),
#'   sampling_times  = c(41, 74, 122, 157, 290, 468, 787, 1192, 1553),
#'   age             = 10,
#'   subject_id      = 3,
#'   treatment_group = 0,
#'   perm_rise       = 6,
#'   temp_rise       = 2,
#'   meas_sd         = 0)
simulate_titres_seropositive <- function(inf_times_sbs, sampling_times, age,
                                        subject_id, treatment_group,
                                        perm_rise,
                                        temp_rise,
                                        temp_decay,
                                        meas_sd) {
  # We assume that the first blood drawn was taken on the individual's birthday
  sampling_times_rel_to_last_birthday <- sampling_times - min(sampling_times)

  # Duration of protection of maternal antibodies
  d_mAB <- 1

  # Sampling times relative to the day the individual became susceptible
  sampling_times_sbs <- sampling_times_rel_to_last_birthday  +
    365 * (age - d_mAB)

  sim_titre_df <- simulate_titre_trajectory(
    sampling_times  = sampling_times_sbs,
    treatment_group = treatment_group,
    subject_id      = subject_id,
    infection_times = inf_times_sbs,
    perm_rise       = perm_rise,
    temp_rise       = temp_rise,
    temp_decay      = temp_decay,
    meas_sd         = meas_sd)

  sim_titre_df$time <- sampling_times

  sim_titre_df
}

simulate_titres_from_enrolment <- function(infection_times, sampling_times,
                                           perm_rise, temp_rise, temp_decay,
                                           baseline)
{
  n_samples      <- length(sampling_times)
  log_titre_vals <- rep(baseline, n_samples)

  event_times <- c(sampling_times[[1]],
                   infection_times)

  n_events <- length(event_times)

  serostatus <- ifelse(baseline > 0, 1, 0)

  if(n_events > 1)
  {
    for(j in 1:n_events)
    {
      if(j == 1) last_titre_prev_inf <- baseline

      titre_obj <- calculate_interval_titres(
        j                   = j,
        event_times         = event_times,
        n_events            = n_events,
        sampling_times      = sampling_times,
        n_samples           = n_samples,
        perm_rise           = perm_rise,
        temp_rise           = temp_rise,
        temp_decay          = temp_decay,
        serostatus          = serostatus,
        baseline            = baseline,
        last_titre_prev_inf = last_titre_prev_inf)

      interval_idx                 <- titre_obj$interval_idx
      log_titre_vals[interval_idx] <- titre_obj$sim_titres
      last_titre_prev_inf          <- titre_obj$last_titre_prev_inf
      # Any event (vac or inf) will render an individual seropositive
      serostatus                   <- 1
    }
  }

  log_titre_vals
}

calculate_interval_titres <- function(j, event_times, n_events, sampling_times,
                                      n_samples, perm_rise, temp_rise,
                                      temp_decay, serostatus, baseline = NULL,
                                      last_titre_prev_inf)
{
  start_interval <- event_times[[j]]

  end_interval <- ifelse(j == n_events,
                         sampling_times[[n_samples]] + 1, # Extra day
                         event_times[[j + 1]])

  times_within_interval <- sampling_times[sampling_times >= start_interval &
                                            sampling_times <= end_interval]

  # We add `end_interval` to the estimate the last titre right before the
  #  next event
  interval_times  <- c(times_within_interval, end_interval)

  time_since_event <- interval_times - start_interval

  temp_rise_seroneg <- temp_rise # Temporary rise for seronegative individuals

  decay_rate <- temp_decay

  if(j == 1)
  {
    temp_rise_this_inf <- max(0, baseline - perm_rise)
  } else {
    delta_seropos <- ifelse(serostatus == 1,
                            last_titre_prev_inf - perm_rise,
                            0)
    temp_rise_this_inf <- temp_rise_seroneg + delta_seropos
  }

  sim_titres <- titre_decay_floor(par_alpha = perm_rise,
                                  par_beta  = temp_rise_this_inf,
                                  par_delta = decay_rate,
                                  time      = time_since_event)

  interval_idx                 <- which(sampling_times %in%
                                          times_within_interval)
  # Any event (vac or inf) will render an individual seropositive
  serostatus <- 1

  list(interval_idx        = interval_idx,
       sim_titres          = sim_titres[-length(sim_titres)],
       # Last titre of the previous infection
       last_titre_prev_inf = sim_titres[length(sim_titres)],
       serostatus          = serostatus)

}
