#' Simulate yearly dengue infections for a cohort from birth
#'
#' This function simulates annual infections with four dengue serotypes for a
#' cohort of individuals followed from birth up to a specified maximum age.
#' Individuals begin fully susceptible to all serotypes. Each year immunity to
#' homotypic reinfection declines wanes, infection events are drawn using a
#' shared force of infection, collision events (multiple serotypes competing to
#' be the infecting serotype) are resolved, and immunity is updated. The
#' function returns a data frame of all infection events.
#'
#' @param lambda_serotype Numeric scalar or numeric vector. If scalar, a
#'   constant force of infection (FOI) is used for all years. If a vector,
#'   \code{lambda_serotype[t]} gives the FOI in calendar year \code{t}.
#'
#' @param loss_rate Numeric scalar. Annual rate at which immunity to each
#'   serotype declines. This value is subtracted from immunity each year and
#'   truncated at zero. A value of 0 corresponds to full susceptibility.
#'
#' @param n_individuals Integer. Number of individuals in the simulated cohort.
#'
#' @param final_age Integer. Total number of simulated years (i.e., the maximum
#'   age tracked in the simulation).
#'
#' @param r_i Numeric vector of length \code{n_individuals} or \code{NULL}.
#'   Individual-level multiplicative risk (frailty) that scales the force of
#'   infection. If \code{NULL}, all individuals have \code{r_i = 1}.
#'
#' @param birth_index Integer vector of length \code{n_individuals} or
#'   \code{NULL}. Birth-year index used to map age to calendar year when
#'   \code{lambda_serotype} is time-varying. For an individual at age
#'   \code{a}, the FOI index is \code{birth_index + a}. Required when
#'   \code{lambda_serotype} has length greater than 1.
#'
#' @return A data frame with one row per infection event and columns:
#'   \itemize{
#'     \item \code{subject_id}: ID of the infected individual.
#'     \item \code{age}: Age (year) of infection (1 to \code{final_age}).
#'     \item \code{serotype}: Infecting serotype (1--4).
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
#'   \item **Infection draw**: For each individual, the probability of infection
#'     with any serotype is computed from the shared force of infection and the
#'     sum of serotype-specific susceptibilities, scaled by \code{r_i}.
#'
#'   \item **Collision resolution**: If more than one serotype is eligible to
#'     cause infection (i.e., an individual has nonzero susceptibility to
#'     multiple serotypes), exactly one serotype is retained.
#'
#'   \item **Serotype assignment**: Conditional on infection, exactly one
#'     serotype is selected with probability proportional to susceptibility.
#' }
#'
#' @examples
#' simulate_DENV_infections_since_birth(
#'   lambda_serotype = 0.1,
#'   loss_rate       = 0.005,
#'   n_individuals   = 5,
#'   final_age       = 5
#' )
#'
#' @export
simulate_DENV_infections_since_birth <- function(lambda_serotype,
                                                 loss_rate,
                                                 n_individuals,
                                                 final_age,
                                                 r_i = NULL,
                                                 birth_index = NULL)
{
  if (is.null(r_i)) {
    r_i <- rep.int(1, n_individuals)
  } else if (length(r_i) != n_individuals) {
    stop("r_i must be NULL or have length n_individuals.")
  }

  time_varying <- length(lambda_serotype) > 1L

  if (time_varying)
  {
    if (is.null(birth_index)) {
      stop("birth_index must be provided (length n_individuals) when lambda is time-varying.")
    }
    if (length(birth_index) != n_individuals) {
      stop("birth_index must have length n_individuals.")
    }

    max_year_needed <- max(birth_index) + final_age

    if (length(lambda_serotype) < max_year_needed) {
      stop("lambda_serotype is too short: need length >= max(birth_index) + final_age.")
    }
  } else {
    if (!is.null(birth_index))
      stop("birth_index must be NULL when lambda_serotype is scalar.")
  }

  # 0 means full susceptibility
  current_state   <- matrix(0, nrow = n_individuals, ncol = 4)

  res_ind <- vector("list", final_age)
  res_age <- vector("list", final_age)
  res_st  <- vector("list", final_age)


  for(i in 1:final_age) # years
  {
    #-----Loss of immunity------------------------------------------------------
    current_state <- current_state - loss_rate
    truncate_at_zero(current_state)
    #---------------------------------------------------------------------------

    sus_lvl   <- 1 - current_state

    if (time_varying) {
      lambda_i <- lambda_serotype[birth_index + i]
    } else {
      lambda_i <- lambda_serotype
    }

    p_any_inf <- 1 - exp(- r_i * lambda_i * rowSums(sus_lvl))

    infected <- stats::runif(n_individuals) < p_any_inf

    idx <- which(infected)

    n_infected <- length(idx)

    if(n_infected > 0)
    {
      w      <- sus_lvl[idx, , drop = FALSE]
      rs     <- rowSums(w)
      norm_w <- w / rs

      u <- stats::runif(nrow(norm_w))

      st <- 1L +
        (u > norm_w[, 1]) +
        (u > norm_w[, 1] + norm_w[, 2]) +
        (u > norm_w[, 1] + norm_w[, 2] + norm_w[, 3])

      coord                <- cbind(idx, st)
      current_state[coord] <- 1

      res_ind[[i]] <- idx
      res_age[[i]] <- rep.int(i, length(idx))
      res_st[[i]]  <- st
    }
  }

  data.frame(
    subject_id = unlist(res_ind, use.names = FALSE),
    age        = unlist(res_age, use.names = FALSE),
    serotype   = unlist(res_st,  use.names = FALSE))
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
simulate_DENV_infections_cohort <- function(lambda_serotype, loss_rate,
                                            cohort_df)
{
  final_age     <- max(cohort_df$age_sample)
  subject_ids   <- unique(cohort_df$subject_id)
  n_individuals <- length(subject_ids)

  cohort_df$key <- paste(cohort_df$subject_id, cohort_df$age_sample, sep = "_")

  sim_inf <- simulate_DENV_infections_since_birth(lambda_serotype,
                                                  loss_rate,
                                                  n_individuals,
                                                  final_age)

  if(is.character(subject_ids))
  {
    sim_inf$infected_ind <- subject_ids[sim_inf$subject_id]
  }

  sim_inf$key <- paste(sim_inf$subject_id, sim_inf$age, sep = "_")

  cohort_df$infection <- cohort_df$key %in% sim_inf$key

  cohort_df$key <- NULL

  cohort_df
}
