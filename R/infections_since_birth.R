#' Simulate yearly dengue infections for a cohort from birth
#'
#' This function simulates annual infections with four dengue serotypes for a
#' cohort of individuals followed from birth up to a specified end year.
#' Individuals begin fully susceptible to all serotypes. Each year immunity
#' wanes, serotype-specific infection probabilities are calculated using a
#' shared force of infection, infection events are drawn, collision events
#' (multiple simultaneous serotype infections) are resolved, and immunity is
#' updated. The function returns a data frame of all infection events.
#'
#' @param lambda_serotype Numeric scalar. The force of infection parameter used
#'   for each of the four serotypes. Higher values imply a higher risk of
#'   infection.
#'
#' @param loss_rate Numeric scalar. Annual rate at which immunity to each
#'   serotype wanes. This value is subtracted from immunity each year and
#'   truncated at zero. A value of 0 corresponds to full susceptibility.
#'
#' @param n_individuals Integer. Number of individuals in the simulated cohort.
#'
#' @param final_age Integer. Total number of simulated years (i.e., the maximum
#'   age tracked in the simulation).
#'
#' @return A data frame with one row per infection event and columns:
#'   \itemize{
#'     \item \code{age}: The year of infection (1 to \code{end_year}).
#'     \item \code{infected_ind}: ID of the infected individual.
#'   }
#'
#' @details
#' The simulation follows these steps each year:
#'
#' \enumerate{
#'   \item **Immunity loss**: Each serotype-specific immunity value decreases by
#'     \code{loss_rate}. Values are truncated to zero. Immunity is bounded in
#'     \eqn{[0, 1]}; 0 means fully susceptible.
#'
#'
#'   \item **Collision resolution**: If an individual draws infection from more
#'     than one serotype, only one is kept, chosen uniformly at random.
#' }
#'
#' @examples
#' simulate_yearly_infections_since_birth(
#'   lambda_serotype = 0.1,
#'   loss_rate       = 0.005,
#'   n_individuals   = 5,
#'   final_age       = 5)
#'
#' @export
simulate_DENV_infections_since_birth <- function(lambda_serotype,
                                                   loss_rate,
                                                   n_individuals,
                                                   final_age)
{
  # 0 means full susceptibility
  current_state   <- matrix(0, nrow = n_individuals, ncol = 4)
  list_annual_inf <- vector("list", final_age) # list of data frames

  for(i in 1:final_age) # years
  {
    #infection_matrix_yr <- matrix(0, nrow = n_individuals, ncol = 4)

    #-----Loss of immunity------------------------------------------------------
    current_state <- (current_state - loss_rate)
    current_state[current_state < 0] <- 0
    #---------------------------------------------------------------------------

    sus_lvl  <- pmax(0, 1 - current_state)
    prob_inf <- 1 - exp(-lambda_serotype * sus_lvl)

    infection_matrix_yr <- matrix(
      stats::rbinom(n_individuals * 4, 1, prob_inf),
      nrow = n_individuals, ncol = 4)

    # Resolve "draws"
    ids_more_1_inf <- which(rowSums( infection_matrix_yr) > 1)

    if(length(ids_more_1_inf) > 0)
    {
      for(j in ids_more_1_inf)
      {
        colliding_st <- which(infection_matrix_yr[j, ] == 1)

        infecting_st <- sample(colliding_st, 1)

        infection_matrix_yr[j, ]             <- 0 # resets
        infection_matrix_yr[j, infecting_st] <- 1 # update
      }
    }

    infected_ind <- which(rowSums(infection_matrix_yr) == 1)

    if(length(infected_ind) > 0)
    {
      current_state[infected_ind, ] <- current_state[infected_ind, ] +
        infection_matrix_yr[infected_ind, ]

      current_state[infected_ind, ][current_state[infected_ind, ] > 1] <- 1

      list_annual_inf[[i]] <- data.frame(infected_ind = infected_ind, age = i)
    }
  }

  do.call("rbind", list_annual_inf)
}

#' Simulate dengue infection outcomes for an existing cohort data frame
#'
#' This function assigns simulated dengue infection events to individuals in an
#' existing cohort data frame.
#'
#' @inheritParams simulate_DENV_infections_since_birth
#'
#' @param cohort_df A data frame containing at least two columns:
#'   \itemize{
#'     \item \code{subject_id}: Unique identifier for each individual.
#'     \item \code{age_sample}: Age (in years) corresponding to each row.
#'   }
#'   The function will add a logical column \code{infection} indicating whether
#'   a simulated infection occurred at that age.
#'
#' @return A data frame identical to \code{cohort_df} but with an additional
#'   logical column \code{infection}, which is \code{TRUE} if the individual
#'   experienced an infection at the corresponding age and \code{FALSE}
#'   otherwise.
#'
#'
#' @examples
#' cohort_df <- data.frame(
#'   subject_id = rep(1:5, each = 5),
#'   age_sample = rep(1:5, times = 5))
#'
#' simulate_DENV_infections_cohort(
#'   lambda_serotype = 0.1,
#'   loss_rate       = 0.01,
#'   cohort_df       = cohort_df)
#'
#' @export
simulate_DENV_infections_cohort <- function(lambda_serotype, loss_rate, cohort_df)
{
  final_age     <- max(cohort_df$age_sample)
  subject_ids   <- unique(cohort_df$subject_id)
  n_individuals <- length(subject_ids)

  cohort_df$key <- paste(cohort_df$subject_id, cohort_df$age_sample, sep = "_")

  sim_inf <- simulate_DENV_infections_since_birth(lambda_serotype,
                                                    loss_rate,
                                                    n_individuals,
                                                    final_age)

  sim_inf$key <- paste(sim_inf$infected_ind, sim_inf$age, sep = "_")

  cohort_df$infection <- cohort_df$key %in% sim_inf$key

  cohort_df$key <- NULL

  cohort_df
}
