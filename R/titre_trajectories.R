#' Simulate Antibody Titre Trajectories for an Individual
#'
#' Generates a simulated time course of true and measured antibody titres
#' for an individual, given sampling times, and infection/vaccination events
#' (exposures). The function applies a rise-and-decay model
#' with measurement error and a limit of detection.
#'
#' @param sampling_times Numeric vector of sampling times, in days relative to
#'   the start of the study.
#' @param subject_id Subject identifier.
#' @param peaks Numeric vector of peak titres immediately after each exposure.
#'   Must have the same length as `exposure_times`.
#' @param perm_rises Numeric vector of long-term plateau titre values following
#'   each exposure. Must have the same length as `exposure_times`.
#' @param decays Numeric vector of decay rates. Must have the same length as
#'   `exposure_times`.
#' @param exposure_times Numeric vector of exposure times relative to enrolment.
#'   Exposures refer to infections and vaccinations. Defaults to `NULL` if no exposures.
#' @param LOD Numeric scalar limit of detection.
#' @param baseline Numeric scalar baseline titre level added to all time points.
#' @param meas_sd Non-negative numeric scalar. Standard deviation of the normal
#'   measurement model for log titres.
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
#' sampling_times <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#' simulate_titre_trajectory(
#'   sampling_times  = sampling_times,
#'   treatment_group = 0,
#'   subject_id      = 1,
#'   exposure_times  = c(300, 900),
#'   peaks           = c(8, 8),
#'   perm_rises      = c(6, 6),
#'   decays          = c(0.003, 0.003),
#'   baseline        = 0,
#'   meas_sd         = 0.3
#' )
simulate_titre_trajectory <- function(sampling_times,
                                      subject_id,
                                      peaks      = NULL,
                                      perm_rises = NULL,
                                      decays     = NULL,
                                      baseline,
                                      meas_sd,
                                      exposure_times = NULL,
                                      LOD            = 1)
{
  if (length(peaks) != length(exposure_times))
  {
    stop("Vectors 'peaks' and 'exposure_times' must have the same length.")
  }

  if (length(perm_rises) != length(exposure_times))
  {
    stop("Vectors 'perm_rises' and 'exposure_times' must have the same length.")
  }

  if (length(decays) != length(exposure_times))
  {
    stop("Vectors 'decays' and 'exposure_times' must have the same length.")
  }

  df <- data.frame(subject_id = subject_id,
                   time       = sampling_times)

  log_titre_vals <- simulate_true_titre_DENV(
    sampling_times  = sampling_times,
    peaks           = peaks,
    perm_rises      = perm_rises,
    decays          = decays,
    baseline        = baseline,
    exposure_times  = exposure_times)

  df$true <- log_titre_vals

  df$meas <- simulate_observed_titres(log_titre_vals,
                                      measurement_error = meas_sd,
                                      LOD = LOD)

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

#' Simulate true log antibody titres for DENV
#'
#' This function generates simulated **true (noise-free)** log-transformed
#' antibody titre values at specified sampling times, given infection events
#' and immune response parameters.
#'
#' It uses the same parameterisation as \code{\link{simulate_titre_trajectory}},
#' but returns only the underlying true trajectory (no measurement error,
#' no limit of detection).
#'
#' Each infection contributes a rise that starts at a peak value and decays
#' exponentially toward a permanent level:
#'
#' \deqn{titre(t) = perm\_rise + (peak - perm\_rise) * exp(-decay * t)}
#'
#' The total titre at each time point is the sum over all infections plus
#' a baseline level.
#'
#' @inheritParams simulate_titre_trajectory
#'
#' @return A numeric vector of length equal to \code{sampling_times},
#' containing simulated **true log antibody titre values**.
#'
#' @export
#'
#' @examples
#' sampling_times <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#'
#' simulate_true_titre_DENV(
#'   sampling_times  = sampling_times,
#'   peaks           = c(8, 7),
#'   perm_rises      = c(4, 3.5),
#'   decays          = c(0.003, 0.002),
#'   baseline        = 2,
#'   treatment_group = 0,
#'   exposure_times = c(300, 900)
#' )
simulate_true_titre_DENV <- function(sampling_times,
                                     peaks,
                                     perm_rises,
                                     decays,
                                     baseline,
                                     exposure_times = NULL)
{
  n_expsr   <- length(exposure_times)
  n_samples <- length(sampling_times)

  if(length(exposure_times) == 0) return(rep(baseline, n_samples))

  titre_matrix <- matrix(0,
                         nrow = n_expsr,
                         ncol = n_samples)

  for(i in seq_len(n_expsr))
  {
    idx             <- which(sampling_times >= exposure_times[[i]])
    times_after_inf <- sampling_times[idx]

    sim_titre <- titre_decay_floor(
      peak      = peaks[[i]],
      perm_rise = perm_rises[[i]],
      decay     = decays[[i]],
      time      = times_after_inf - exposure_times[[i]])

    titre_matrix[i, idx] <- sim_titre
  }

  baseline + colSums(titre_matrix)
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
                                        subject_id,
                                        peaks,
                                        perm_rises,
                                        decays,
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
    subject_id      = subject_id,
    exposure_times  = inf_times_sbs,
    peaks           = peaks,
    perm_rises      = perm_rises,
    decays          = decays,
    baseline        = 0,
    meas_sd         = meas_sd)

  sim_titre_df$time <- sampling_times

  sim_titre_df
}
