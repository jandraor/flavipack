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

#' Linear antibody titre dynamics with lower bound
#'
#' @param par_alpha A scalar indicating the permanent rise.
#' @param par_beta A scalar indicating the temporary rise.
#' @param par_delta A scalar indicating the decay rate.
#' @param time Numeric vector. Time points at which to evaluate titres.
#'
#' @returns Numeric vector of titres at each time point.
#' @export
#'
#' @examples
#' titre_decay_floor(par_alpha = 6,
#'                   par_beta  = 2,
#'                   par_delta = 0.003,
#'                   time      = seq(0, 365 * 3))
titre_decay_floor <- function(par_alpha, par_beta, par_delta, time)
{
  pmax((par_alpha + par_beta) - par_delta * time, par_alpha)
}
