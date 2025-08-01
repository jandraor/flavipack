#' Simulate infection times
#'
#' @param lambda A scalar corresponding to the daily force of infection.
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
  time_to_first_infection      <- stats::rexp(1, rate = lambda)

  time_between_infections_rest <- min_gap + stats::rexp(max_inf - 1,
                                                        rate = lambda)

  time_between_infections <- c(time_to_first_infection,
                               time_between_infections_rest)

  infection_times <- cumsum(time_between_infections)

  infection_times <- infection_times[infection_times <= stop_time]

  infection_times
}

#' Simulates infection times from birth
#'
#' @inheritParams simulate_infection_times
#' @param age_enrolment Numeric. Age of the individual at study enrolment (in years).
#' @param offset_start Numeric. Time from start of study to enrolment (in days).
#' @param stop_time Numeric. Duration from the start to the end of the study (in days)..
#'
#' @return Numeric vector of infection times (relative to birth).
#' @examples
#' simulate_infection_times_since_birth(
#' lambda        = 0.3 / 365,
#' age_enrolment = 8, # in years
#' offset_start  = 0,
#' stop_time     = 365 * 5,
#' min_gap       = 365,
#' max_inf       = 4)
#'
#'@export
simulate_infection_times_since_birth <- function(lambda,
                                                 age_enrolment,
                                                 offset_start,
                                                 stop_time,
                                                 min_gap,
                                                 max_inf)
{
  time_protection_mab <- 1 # year
  # Number of days elapsed from the day the individual became susceptible until
  #  the end of the study.
  time_birth_stop <- (age_enrolment - time_protection_mab) * 365 + stop_time -
    offset_start

  simulate_infection_times(lambda, time_birth_stop, min_gap, max_inf)
}
