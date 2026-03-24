#' Simulate infection times
#'
#' @param lambda A scalar corresponding to the annual force of infection.
#' @param stop_time A numeric scalar specifying the maximum time horizon for the
#'  simulation.
#' @param min_gap A scalar representing the minimum time that must elapse
#'   between two consecutive infections. This interval is due to temporary
#'   immune protection.
#' @param max_inf An integer indicating the maximum number of infections that
#' an individual can experience from a given pathogen.
#'
#' @returns A numeric vector of infection times occurring before \code{stop_time}
#' @export
#'
#' @examples
#' simulate_infection_times(lambda     = 0.3,
#'                          stop_time  = 365 * 3,
#'                          min_gap    = 365,
#'                          max_inf    = 4)
simulate_infection_times <- function(lambda, stop_time, min_gap, max_inf)
{
  if(lambda <= 0) return(numeric(0))

  daily_lambda <- lambda / 365
  n_serotypes  <- max_inf
  inf_times    <- c()

  current_time <- 0

  while(current_time < stop_time && n_serotypes > 0L)
  {
    inf_time    <- current_time +
      stats::rexp(1, rate = n_serotypes * daily_lambda)

    if(inf_time >= stop_time) break

    inf_times <- c(inf_times, inf_time)

    n_serotypes  <- n_serotypes - 1L
    current_time <- inf_time + min_gap
  }

  inf_times
}

#' Simulate infection times since susceptibility
#'
#' Simulates infection times from the end of maternal antibody protection
#' (assumed to occur 1 year after birth) until the individual's last
#' measurement.
#'
#' @inheritParams simulate_infection_times
#' @param age_enrolment Numeric. Age of the individual at study enrolment, in years.
#' @param enrolment_time Numeric. Time from the start of the study to this
#'   individual's enrolment, in days.
#' @param follow_up_time Numeric. Time from the start of the study to the
#'   individual's last measurement, in days.
#'
#' @return Numeric vector of infection times, in days since becoming susceptible
#'   (that is, since the end of maternal antibody protection).
#'
#' @examples
#' simulate_infection_times_since_susceptibility(
#'   lambda = 0.3,
#'   age_enrolment = 8,      # years
#'   enrolment_time = 0,     # days since study start
#'   follow_up_time = 365 * 5,
#'   min_gap = 365,
#'   max_inf = 4
#' )
#'
#' @export
simulate_infection_times_since_susceptibility <- function(lambda,
                                                          age_enrolment,
                                                          enrolment_time,
                                                          follow_up_time,
                                                          min_gap,
                                                          max_inf)
{
  stop_time <- calculate_stop_time(age_enrolment,
                                   enrolment_time,
                                   follow_up_time)

  simulate_infection_times(lambda, stop_time, min_gap, max_inf)
}

calculate_stop_time <- function(age_enrolment, enrolment_time, follow_up_time)
{
  time_protection_mab <- 1 # year

  # Number of days elapsed from the day the individual became susceptible until
  #  the end of the individual's follow-up.

  (age_enrolment - time_protection_mab) * 365 +
    follow_up_time - enrolment_time
}
