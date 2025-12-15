#' Add titres to infection history
#'
#' @param short_options A list of parameters to configure the short-term decay
#'  rate.
#' @param long_options A list of parameters to configure the long-term decay
#'  rate.
#' @param rise_options A list of parameters to configure the rise.
#' @param inf_history_df A data frame
#'
#' @returns A data frame
#' @export
#'
#' @examples
#' inf_history_df <-
#'   data.frame(subject_id = 1,
#'              year_index = 11:23,
#'              age        = 2:14,
#'              sim_y      = 0)
#' inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1
#'
#' add_titre_to_infection_history(
#' inf_history_df,
#' long_options  = list(L = 0.22, a = -0.96, b = 18.6),
#' short_options = list(intercept = 0.02,
#' slope     = 0.003342),
#' rise_options  = list(intercept = 8,
#'                      slope     = -0.8))
add_titre_to_infection_history <- function(
    inf_history_df,
    short_options,
    long_options,
    rise_options)
{
  # Required columns
  required_cols <- c("subject_id", "year_index", "age", "sim_y")

  missing <- setdiff(required_cols, names(inf_history_df))

  if (length(missing) > 0) {
    stop(paste("Missing required columns in `inf_history_df`:",
               paste(missing, collapse = ", ")))
  }


  inf_history_df$cumulative_infections <- cumsum(inf_history_df$sim_y)

  df_list <- split(inf_history_df, inf_history_df$cumulative_infections)

  n_df <- length(df_list)

  baseline_titre <- 0

  for(i in seq_len(n_df))
  {
    df    <- df_list[[i]]
    n_inf <- unique(df$cumulative_infection)
    n_row <- nrow(df)

    if(n_inf == 0) log_titre <- rep(0, n_row)

    if(n_inf > 0)
    {
      rise     <- estimate_rise(baseline_titre, rise_options)
      log_peak <- baseline_titre + rise

      # short-term decay rate
      short_rate <- estimate_short_rate(rise, short_options)

      # Long-term decay rate
      long_rates <- estimate_long_rate(df$age, long_options)

      long_rates <- long_rates / 360

      A0 <- exp(log_peak)

      # n_years   = n_row embeds an extra index to obtain the baseline titre for
      #   a subsequent infection

      titre <- simulate_antibody_titres(n_years   = n_row,
                                        A0        = A0,
                                        par_rho   = 0.98,
                                        par_sigma = short_rate,
                                        gamma_t   = long_rates)

      baseline_titre <- log(titre[length(titre)])

      log_titre    <- log(titre[-length(titre)])
    }

    df$log_titre <- log_titre
    df_list[[i]] <- df
  }

  output_df <- do.call(rbind, df_list)

  rownames(output_df) <- NULL

  output_df
}


#' Simulate bi-exponential decay
#'
#' @param t Time
#' @param A0 Peak
#' @param par_sigma Rate of short-term decay
#' @param par_gamma Rate of long-term decay
#' @param par_rho Proportion of short-term decay
#'
#' @returns A numeric vector
#' @export
#'
#' @examples
#' bi_exponential_decay(365, exp(10), 0.98, 1/720, 1/30)
bi_exponential_decay <- function(t, A0, par_sigma, par_gamma, par_rho)
{
  A0 * (par_rho * exp(-par_sigma * t) +
          (1 - par_rho) * exp(-par_gamma * t))
}

bi_exponential_step <- function(A0, par_sigma, par_gamma, par_rho)
{
  short_dyn <-  A0 * par_rho * exp(-par_sigma * 360)

  long_dyn  <- A0 * (1 - par_rho) * exp(-par_gamma * 360)

  A <- short_dyn + long_dyn

  list(A       = A,
       new_rho = short_dyn / A)
}

estimate_short_rate <- function(rise, short_options) {

  short_options$intercept + short_options$slope * rise
}

# baseline_titre is in log scale
estimate_rise <- function(baseline_titre, rise_options) {
  max(0, rise_options$intercept + rise_options$slope * max(0, baseline_titre))
}

# Assumes the gamma = f(age), where f is a logistic function
estimate_long_rate <- function(age, long_options)
{
  L <- long_options$L
  a <- long_options$a
  b <- long_options$b

  L / (1 + exp(-a * (age - b)))
}

simulate_antibody_titres <- function(n_years, A0, par_rho,
                                     par_sigma,
                                     gamma_t) {
  A <- numeric(n_years + 1)
  A[1] <- A0

  for (j in seq_len(n_years)) {

    par_gamma <- gamma_t[j]

    output <- bi_exponential_step(
      A0         = A[j],
      par_sigma  = par_sigma,
      par_gamma  = par_gamma,
      par_rho    = par_rho)

    A[j + 1] <- output$A
    par_rho  <- output$new_rho
  }

  A
}

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
#' simulate_DENV_long_decay_titres(
#'   inf_times = c(10, 40, 120),
#'   decay_rate_vec = c(0.05, 0.03, 0.02, 0.02),
#'   log_first_peak = 6,
#'   phi = 0.5,
#'   kappa = 0.4,
#'   subject_id = "z1",
#'   final_age = 200
#' )
#'
#' @export
simulate_DENV_long_decay_titres <- function(inf_times,
                                            decay_rate_vec,
                                            log_first_peak,
                                            phi,
                                            beta,
                                            subject_id,
                                            final_age)
{
  titre_vals <- rep(5, final_age)

  n_inf <- length(inf_times)

  if(n_inf > 0)
  {
    for(inf_idx in 1:n_inf)
    {
      inf_age <- inf_times[[inf_idx]]

      if(inf_idx == 1) A0 <- inv_log2_transform(log_first_peak)

      if(inf_idx > 1)
      {
        multiplier <- 1 + phi * (1 - exp(-beta * (inf_idx - 1)))
        A0 <- inv_log2_transform(log_first_peak * multiplier)
      }

      decay_rate <- ifelse(inf_idx < 4,
                           decay_rate_vec[[inf_idx]], decay_rate_vec[[4]])

      titre_vals[inf_age:final_age] <- A0 * exp(-decay_rate * 0:(final_age - inf_age))
    }
  }

  data.frame(subject_id = subject_id,
             age        = 1:final_age,
             titre      = titre_vals)
}

#' Exponential antibody titre dynamics with a lower bound
#'
#' This function models antibody titre dynamics following infection or vaccination,
#' assuming an exponential decay from a peak towards a permanent baseline level.
#' The titre value starts at \code{par_alpha + par_beta} and decays exponentially
#' at rate \code{par_delta}, asymptotically approaching the lower bound \code{par_alpha}.
#'
#' @param par_alpha Numeric scalar. Permanent baseline (lower bound) of the antibody titre.
#' @param par_beta Numeric scalar. Temporary rise above the baseline immediately after infection.
#' @param par_delta Numeric scalar. Decay rate controlling how quickly titres return to baseline.
#' @param time Numeric vector. Time points at which to evaluate titres.
#'
#' @returns A numeric vector of titre values corresponding to each time point.
#' @export
#'
#' @examples
#' # Example: 3-year decay from a peak titre back to baseline
#' titre_decay_floor(
#'   par_alpha = 6,
#'   par_beta  = 2,
#'   par_delta = 0.003,
#'   time      = seq(0, 365 * 3)
#' )
titre_decay_floor <- function(par_alpha, par_beta, par_delta, time)
{
  par_alpha + par_beta * exp(-par_delta * time)
}
