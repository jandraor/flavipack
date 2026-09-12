#' Simulate infection times since susceptibility
#'
#' Simulates infection times from the end of maternal antibody protection
#' (assumed to occur 1 year after birth) until the individual's last
#' measurement.
#'
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

  rInfections_DENV(lambda, stop_time, min_gap, max_inf)
}

calculate_stop_time <- function(age_enrolment, enrolment_time, follow_up_time)
{
  time_protection_mab <- 1 # year

  # Number of days elapsed from the day the individual became susceptible until
  #  the end of the individual's follow-up.

  (age_enrolment - time_protection_mab) * 365 +
    follow_up_time - enrolment_time
}
