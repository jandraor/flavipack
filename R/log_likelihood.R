#' Log-likelihood of a dengue infection history under a constant FOI model
#'
#' Computes the log-likelihood of an observed infection history assuming a
#' constant force of infection (FOI) per serotype and a decreasing number of
#' susceptible serotype "buckets" following each infection. The number of
#' initially susceptible serotypes is given by \code{empty_buckets} (typically
#' 4 for seronegative individuals and 3 for seropositive individuals, who have
#' already been infected by at least one serotype before enrolment). After each
#' infection, the number of susceptible serotypes decreases by one, to a
#' minimum of zero.
#'
#' The likelihood is composed of:
#' \itemize{
#'   \item The probability of escaping infection between consecutive infection
#'         events.
#'   \item The hazard of infection at each observed infection time.
#' }
#'
#' An optional minimum gap between infections can be specified. During this
#' interval, individuals are assumed not to be at risk of infection.
#'
#' @param FOI Constant force of infection per serotype.
#' @param start_time Start of the observation period.
#' @param end_time End of the observation period.
#' @param infection_history Numeric vector containing infection times in
#'   increasing order. If \code{NULL}, the likelihood corresponds to escaping
#'   infection throughout the observation period.
#' @param min_gap Minimum period after an infection during which no further
#'   infections can occur. Default is \code{0}.
#' @param empty_buckets Initial number of susceptible serotypes (buckets).
#'   Typically 4 for seronegative individuals and 3 for seropositive (who have
#'   already been infected by at least one serotype before enrolment). Decreases
#'   by one after each infection, to a minimum of zero.
#'
#' @return A numeric value giving the log-likelihood of the infection history.
#'
#' @examples
#' # No infections (seronegative: 4 susceptible serotypes)
#' dInfections_DENV(
#'   FOI = 0.01,
#'   start_time = 0,
#'   end_time = 365,
#'   empty_buckets = 4)
#'
#' # Two infections, seropositive (3 susceptible serotypes)
#' dInfections_DENV(
#'   FOI = 0.01,
#'   start_time = 0,
#'   end_time = 365,
#'   infection_history = c(100, 250),
#'   min_gap = 30,
#'   empty_buckets = 3)
#'
#' @export
dInfections_DENV <- function(FOI, start_time, end_time,
                             infection_history = NULL,
                             min_gap = 0,
                             empty_buckets)
{
  if(is.null(infection_history))
  {
    return (-empty_buckets * FOI * (end_time - start_time))
  }

  risk_start <- start_time

  LL <- 0

  for(inf_time in infection_history)
  {
    risk_end <- inf_time
    interval_length <- risk_end - risk_start

    # prob of escaping infection
    LL <- LL - empty_buckets * FOI * interval_length

    # prob of infection
    LL <- LL + log(empty_buckets * FOI)

    risk_start <- risk_end + min_gap
    empty_buckets <- empty_buckets - 1
  }

  LL <- LL - empty_buckets * FOI * (end_time - risk_start)

  LL
}
