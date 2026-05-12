#' Simulate DENV Long-Term Titre Trajectories
#'
#' This function simulates long-term antibody titre trajectories for dengue
#'  virus (DENV) infections, given an individual's infection history. Each
#'  infection produces a peak antibody titre that may be boosted by subsequent
#'  infections. The boost saturates asymptotically at a maximum fraction
#'  \code{phi} relative to the first infection peak, with the rate of saturation
#'  controlled by \code{beta}. Titres then decay exponentially from the peak
#'  using the corresponding decay rate.
#'
#' @param inf_times Integer vector of ages at which infections occur.
#' @param decay_rate_vec Numeric vector of decay rates for each infection. The
#'   i-th element is used for the i-th infection; if more than four infections
#'   occur, the fourth element is reused for subsequent infections.
#' @param log_first_peak Numeric. Log2-scale peak antibody titre of the first
#'  infection.
#' @param phi Numeric. Asymptotic maximum proportional boost for reinfections,
#'   relative to the first infection peak. For example, \code{phi = 0.5} means
#'   that at saturation the peak titre can increase by up to 50% over the first
#'   peak.
#' @param beta Numeric. Saturation rate controlling how quickly the boost
#'   approaches \code{phi} with successive reinfections.
#' @param subject_id Integer or character. Identifier for the simulated subject.
#' @param final_age Integer. Final age/time index to simulate to; titres will
#'   be returned for ages 1:\code{final_age}.
#'
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{subject_id}{Subject identifier.}
#'   \item{age}{Integer ages from 1 to \code{final_age}.}
#'   \item{titre}{Simulated antibody titre at each age.}
#' }
#'
#' @examples
#' inf_time_list <- replicate(c(5, 7, 10, 15), n = 2)
#' simulate_DENV_long_decay_titres(inf_time_list,
#'                                 decay_rate_vec = c(0.08, 0.04, 0.02, 0.01),
#'                                 log_first_peak = 1.3,
#'                                 phi = 4,
#'                                 beta = 1,
#'                                 final_age = 20)
#'
#' @export
simulate_DENV_long_decay_titres <- function(inf_times_list,
                                            decay_rate_vec,
                                            log_first_peak,
                                            phi,
                                            beta,
                                            final_age)
{
  simulate_DENV_long_decay_titres_cpp(inf_times_list,
                                      decay_rate_vec,
                                      log_first_peak,
                                      phi,
                                      beta,
                                      final_age)
}


#' Exponential antibody titre dynamics with a lower bound
#'
#' This function models antibody titre dynamics following infection,
#' assuming an exponential decay from an initial peak towards a
#' long-term plateau (permanent rise) in the log scale.
#'
#' The titre starts at \code{peak} when \code{time = 0} and decays
#' exponentially at rate \code{decay}, asymptotically approaching
#' \code{perm_rise}.
#'
#' @param peak Numeric scalar. Peak titre immediately after infection.
#' @param perm_rise Numeric scalar. Long-term plateau (asymptotic titre level).
#' @param decay Numeric scalar. Decay rate controlling how quickly titres decline.
#' @param time Numeric vector. Time since infection at which to evaluate titres.
#'
#' @returns A numeric vector of titre values corresponding to each time point.
#' @export
#'
#' @examples
#' titre_decay_floor(
#'   peak      = 8,
#'   perm_rise = 6,
#'   decay     = 0.003,
#'   time      = seq(0, 365 * 3)
#' )
titre_decay_floor <- function(peak, perm_rise, decay, time)
{
  perm_rise + (peak - perm_rise) * exp(-decay * time)
}
