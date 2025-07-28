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
  time_to_first_infection      <- rexp(1, rate = lambda)

  time_between_infections_rest <- min_gap + rexp(max_inf - 1, rate = lambda)

  time_between_infections <- c(time_to_first_infection,
                               time_between_infections_rest)

  infection_times <- cumsum(time_between_infections)

  infection_times <- infection_times[infection_times <= stop_time]

  infection_times
}
