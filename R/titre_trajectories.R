#' Simulate Antibody Titre Trajectories for an Individual
#'
#' Generates simulated time courses of true and measured antibody titres
#' for an individual, given sampling times and infection/vaccination events
#' (exposures). The function applies a rise-and-decay model with measurement
#' error and a limit of detection. Multiple antibody markers can be simulated,
#' with homologous and heterologous responses determined by the infecting
#' pathogen and exposure-specific titre multipliers.
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
#' @param meas_sd Numeric vector of standard deviations of the normal
#'   measurement model for log titres. Must contain one value per marker.
#' @param exposure_times Numeric vector of exposure times relative to enrolment.
#'   Exposures refer to infections and vaccinations. Defaults to `NULL` if
#'   there are no exposures.
#' @param LOD Numeric scalar limit of detection.
#' @param n_markers Positive integer giving the number of antibody markers to
#'   simulate. Defaults to `1`.
#' @param infecting_pathogens Integer vector identifying the infecting pathogen
#'   for each exposure. Pathogen indices are zero-based. Must have the same
#'   length as `exposure_times` when `n_markers > 1`. For a single-marker model,
#'   defaults to pathogen `0` for all exposures.
#' @param titre_multipliers List of numeric vectors specifying the
#'   heterologous titre multipliers for each exposure. When `n_markers > 1`,
#'   each exposure must have `n_markers - 1` multipliers, ordered by marker
#'   index while excluding the infecting pathogen. Not used when
#'   `n_markers = 1`.
#'
#' @returns A `data.frame` with one row per marker and sampling time, with
#'   columns:
#' \describe{
#'   \item{subject_id}{Subject identifier.}
#'   \item{time}{Sampling time.}
#'   \item{marker}{Zero-based antibody marker index.}
#'   \item{true}{Simulated true log-titre at each time.}
#'   \item{meas}{Observed log-titre with measurement error and LOD applied.}
#' }
#'
#' @export
#'
#' @examples
#' sampling_times <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#'
#' simulate_titre_trajectory(
#'   sampling_times = sampling_times,
#'   subject_id     = 1,
#'   exposure_times = c(300, 900),
#'   peaks          = c(8, 8),
#'   perm_rises     = c(6, 6),
#'   decays         = c(0.003, 0.003),
#'   meas_sd        = 0.3)
simulate_titre_trajectory <- function(sampling_times,
                                      subject_id,
                                      peaks          = NULL,
                                      perm_rises     = NULL,
                                      decays         = NULL,
                                      meas_sd,
                                      exposure_times = NULL,
                                      LOD            = 1,
                                      n_markers      = 1,
                                      infecting_pathogens = NULL,
                                      titre_multipliers = list())
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

  if(length(meas_sd) != n_markers) stop("meas_sd must have one value per marker")

  df <- data.frame(subject_id = subject_id,
                   time       = sampling_times)

  true_titres <- simulate_true_titre_DENV(
    sampling_times       = sampling_times,
    peaks                = peaks,
    perm_rises           = perm_rises,
    decays               = decays,
    exposure_times       = exposure_times,
    n_markers            = n_markers,
    infecting_pathogens  = infecting_pathogens,
    titre_multipliers    = titre_multipliers)

  df_list <- vector("list", n_markers)

  for(marker_idx in 0:(n_markers - 1))
  {
    true_titre <- true_titres[marker_idx + 1, ]

    df_list[[marker_idx + 1]] <- data.frame(
      subject_id = subject_id,
      time       = sampling_times,
      marker     = marker_idx,
      true       = true_titre,
      meas       = simulate_observed_titres(
        true_titre,
        measurement_error = meas_sd[[marker_idx + 1]],
        LOD = LOD)
    )
  }

  do.call(rbind, df_list)
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
#' simulate_observed_titres(true_titre = true,
#'                          measurement_error = 1.5,
#'                          LOD = 1)
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

#' Simulate True Log Antibody Titres for DENV
#'
#' Generates simulated **true (noise-free)** log-transformed antibody titre
#' values at specified sampling times, given exposure events and immune
#' response parameters.
#'
#' Each exposure generates a homologous response in the marker corresponding
#' to the infecting pathogen. When multiple markers are simulated, the same
#' exposure can also generate heterologous responses in the other markers,
#' scaled by exposure-specific titre multipliers.
#'
#' Each exposure contributes a rise that starts at a peak value and decays
#' exponentially toward a permanent level:
#'
#' \deqn{titre(t) = perm\_rise + (peak - perm\_rise) * exp(-decay * t)}
#'
#' Contributions from multiple exposures are summed within each marker.
#' Titres are zero before any contributing exposure.
#'
#' @inheritParams simulate_titre_trajectory
#'
#' @return A numeric matrix with \code{n_markers} rows and
#'   \code{length(sampling_times)} columns. Rows correspond to antibody
#'   markers and columns correspond to sampling times. Marker indices supplied
#'   through \code{infecting_pathogens} are zero-based.
#'
#' @export
#'
#' @examples
#' sampling_times <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#'
#' simulate_true_titre_DENV(
#'   sampling_times = sampling_times,
#'   peaks          = c(8, 7),
#'   perm_rises     = c(4, 3.5),
#'   decays         = c(0.003, 0.002),
#'   exposure_times = c(300, 900)
#' )
simulate_true_titre_DENV <- function(sampling_times,
                                     peaks,
                                     perm_rises,
                                     decays,
                                     exposure_times = NULL,
                                     n_markers = 1,
                                     infecting_pathogens = NULL,
                                     titre_multipliers = list())
{
  n_expsr   <- length(exposure_times)
  n_samples <- length(sampling_times)

  if(length(peaks) != n_expsr ||
     length(perm_rises) != n_expsr ||
     length(decays) != n_expsr)
  {
    stop("peaks, perm_rises, decays, and exposure_times must have the same length")
  }

  titre_matrix <- matrix(0, nrow = n_markers, ncol = n_samples)

  if(n_expsr == 0) return(titre_matrix)

  if(is.null(infecting_pathogens))
  {
    if(n_markers == 1)
    {
      infecting_pathogens <- rep(0, n_expsr)
    }
    else
    {
      stop("infecting_pathogens must be provided when n_markers > 1")
    }
  }

  if(length(infecting_pathogens) != n_expsr)
  {
    stop("infecting_pathogens and exposure_times must have the same length")
  }

  if(n_markers > 1)
  {
    if(length(titre_multipliers) != n_expsr)
    {
      stop("titre_multipliers and exposure_times must have the same length")
    }

    if(any(lengths(titre_multipliers) != n_markers - 1))
    {
      stop("Each titre_multipliers element must contain n_markers - 1 values")
    }
  }

  for(i in seq_len(n_expsr))
  {
    idx <- which(sampling_times >= exposure_times[[i]])

    if(length(idx) == 0) next

    times_after_exposure <- sampling_times[idx] - exposure_times[[i]]

    contribution <- titre_decay_floor(
      peak      = peaks[[i]],
      perm_rise = perm_rises[[i]],
      decay     = decays[[i]],
      time      = times_after_exposure)

    pathogen <- infecting_pathogens[[i]]


    # Homologous response
    # The plus one is because infecting_pathogen is 0-based
    titre_matrix[pathogen + 1, idx] <-
      titre_matrix[pathogen + 1, idx] + contribution

    # Heterologous responses
    if(n_markers > 1)
    {
      counter <- 1

      for(marker_idx in 0:(n_markers - 1))
      {
        if(marker_idx != pathogen)
        {
          titre_matrix[marker_idx + 1, idx] <-
            titre_matrix[marker_idx + 1, idx] +
            contribution * titre_multipliers[[i]][[counter]]

          counter <- counter + 1

        }
      }
    }
  }

  titre_matrix
}


#' Simulate Antibody Titre Trajectories from Birth
#'
#' Simulates true and measured antibody titre trajectories for a seropositive
#' individual whose infection times are defined relative to the day they became
#' susceptible. Sampling times are converted to the same time scale before
#' simulating the antibody trajectories.
#'
#' Multiple antibody markers can be simulated, with homologous and
#' heterologous responses determined by the infecting pathogen and
#' exposure-specific titre multipliers.
#'
#' @inheritParams simulate_titre_trajectory
#' @param inf_times_sbs Numeric vector of infection times, in days relative to
#'   the day the individual became susceptible.
#' @param age Numeric. Age of the individual at the time of the first sample,
#'   in years.
#'
#' @returns A `data.frame` with one row per marker and sampling time, containing
#'   the subject identifier, sampling time, marker, simulated true log-titre,
#'   and simulated measured log-titre.
#'
#' @export
#'
#' @examples
#' simulate_titres_seropositive(
#'   inf_times_sbs  = c(1059.646, 2743.543, 4152.091),
#'   sampling_times = c(41, 74, 122, 157, 290, 468, 787, 1192, 1553),
#'   age             = 10,
#'   subject_id      = 3,
#'   peaks           = c(6, 6, 6),
#'   perm_rises      = c(2, 2, 2),
#'   decays          = c(0.01, 0.01, 0.01),
#'   meas_sd         = 0
#' )
simulate_titres_seropositive <- function(inf_times_sbs,
                                         sampling_times,
                                         age,
                                         subject_id,
                                         peaks,
                                         perm_rises,
                                         decays,
                                         meas_sd,
                                         n_markers = 1,
                                         infecting_pathogens = NULL,
                                         titre_multipliers = list())
{
  # We assume that the first blood drawn was taken on the individual's birthday
  sampling_times_rel_to_last_birthday <- sampling_times - min(sampling_times)

  # Duration of protection of maternal antibodies
  d_mAB <- 1

  # Sampling times relative to the day the individual became susceptible
  sampling_times_sbs <- sampling_times_rel_to_last_birthday  +
    365 * (age - d_mAB)

  sim_titre_df <- simulate_titre_trajectory(
    sampling_times      = sampling_times_sbs,
    subject_id          = subject_id,
    exposure_times      = inf_times_sbs,
    peaks               = peaks,
    perm_rises          = perm_rises,
    decays              = decays,
    meas_sd             = meas_sd,
    n_markers           = n_markers,
    infecting_pathogens = infecting_pathogens,
    titre_multipliers   = titre_multipliers)

  sim_titre_df$time <- sampling_times

  sim_titre_df
}
