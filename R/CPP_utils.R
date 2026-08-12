#' Create a subject-level input list from multiple data frames
#'
#' This function constructs a list of subject-specific inputs by combining
#' participant metadata, titre measurements, symptom (infection) times, and
#' optional vaccination, guessed-infection and flicked-infection data. Each
#' subject is returned as a list element containing their longitudinal
#' measurements and associated covariates.
#'
#' @param part_df A data frame containing participant-level information.
#'   Must include the columns \code{subject_id}, \code{age_enrolment},
#'   \code{location}, and \code{serostatus}.
#' @param titre_df A data frame containing titre measurements. Must include the
#'   columns \code{subject_id}, \code{time}, \code{meas}, and \code{marker}.
#' @param symp_df A data frame containing symptom (e.g., infection) times with
#'   columns \code{subject_id} and \code{time}. Times after a subject's last
#'   titre observation are dropped.
#' @param vacc_df An optional data frame containing vaccination records with
#'   columns \code{subject_id}, \code{time}, and \code{baseline_titre}.
#'   Defaults to \code{NULL}. Times more than one time unit after a subject's
#'   last titre observation are dropped.
#' @param guess_df An optional data frame of initial guesses for infection
#'   times, with columns \code{subject_id} and \code{time}. Defaults to
#'   \code{NULL}.
#' @param n_markers An integer specifying the number of markers. (Currently
#'   unused but retained for compatibility or future extensions.)
#' @param flick_df An optional data frame of infection times to flick, with
#'   columns \code{subject_id} and \code{time}. Defaults to \code{NULL}. Each
#'   flick time must match one of the subject's retained symptomatic infection
#'   times; an error is raised otherwise, or if a subject has no retained
#'   infections at all.
#'
#' @return A list where each element corresponds to a subject. Each subject is
#'   represented as a list with the following components:
#'   \describe{
#'     \item{subject_id}{Unique subject identifier, coerced to character}
#'     \item{age_enrolment}{Age at enrolment}
#'     \item{location}{Subject location}
#'     \item{serostatus}{Baseline serostatus}
#'     \item{obs_times}{List of unique observation times from \code{titre_df}}
#'     \item{measurements}{Named list with one element per level of
#'       \code{marker}, each a vector of measurements for that marker}
#'     \item{infection_times}{(Optional) Sorted list of infection times, present
#'       only if \code{symp_df} contains times within the subject's follow-up}
#'     \item{vaccination}{(Optional) List of vaccination records, each a list
#'       with elements \code{time} and \code{baseline_titre}. Present only if
#'       \code{vacc_df} is supplied and contains rows for the subject}
#'     \item{guessed_infections}{(Optional) Sorted list of guessed infection
#'       times. Present only if \code{guess_df} is supplied and contains rows
#'       for the subject}
#'     \item{flicks}{(Optional) Sorted list of flicked infection times. Present
#'       only if \code{flick_df} is supplied and contains rows for the subject}
#'   }
#'
#' @examples
#' part_df <- data.frame(
#'   subject_id = c(1, 2),
#'   age_enrolment = c(30, 25),
#'   location = c("A", "B"),
#'   serostatus = c(1, 0)
#' )
#'
#' titre_df <- data.frame(
#'   subject_id = c(1, 1, 2),
#'   time = c(0, 1, 0),
#'   meas = c(10, 20, 15),
#'   marker = "mean"
#' )
#'
#' symp_df <- data.frame(
#'   subject_id = c(1),
#'   time = c(0.5)
#' )
#'
#' create_input_list(part_df, titre_df, symp_df, n_markers = 1)
#'
#' @export
create_input_list <- function(part_df, titre_df, symp_df,
                              vacc_df = NULL,
                              guess_df = NULL,
                              n_markers,
                              flick_df = NULL)
{
  #----checks-------------------------------------------------------------------
  required_cols <- c("subject_id", "age_enrolment", "location", "serostatus")
  check_df_colnames(part_df, required_cols)

  required_cols <- c("subject_id", "time", "meas", "marker")
  check_df_colnames(titre_df, required_cols)
  #-----------------------------------------------------------------------------
  subject_ids <- unique(part_df$subject_id)

  lapply(subject_ids, \(s_id) {

    subject_titre <- titre_df[titre_df$subject_id == s_id, ]

    subject_info <- part_df[part_df$subject_id == s_id, ]

    subject_symp <- symp_df[symp_df$subject_id == s_id, ]

    obs_times_vctr <- unique(subject_titre$time)

    obs_times <- I(as.list(obs_times_vctr))

    measurements <- create_measurement_object(subject_titre)

    obj <- list(
      subject_id      = as.character(s_id),
      age_enrolment   = subject_info$age_enrolment,
      location        = subject_info$location,
      serostatus      = subject_info$serostatus,
      obs_times       = obs_times,
      measurements    = measurements)

    if(nrow(subject_symp) > 0)
    {
      dropout_time <- max(subject_titre$time)

      subject_symp <- subject_symp[subject_symp$time <= dropout_time, ,
                                   drop = FALSE]

      if (nrow(subject_symp) > 0)
      {
        obj$infection_times <- I(as.list(sort(subject_symp$time)))
      }
    }

    if(!is.null(vacc_df))
    {
      subject_vacc <- vacc_df[vacc_df$subject_id == s_id, , drop = FALSE]

      if (nrow(subject_vacc) > 0)
      {
        dropout_time <- max(subject_titre$time) + 1

        subject_vacc <- subject_vacc[
          subject_vacc$time <= dropout_time,
          ,
          drop = FALSE]

        if (nrow(subject_vacc) > 0)
        {
          obj$vaccination <- lapply(seq_len(nrow(subject_vacc)), \(i) {
            list(time           = subject_vacc$time[i],
                 baseline_titre = subject_vacc$baseline_titre[i])

          })
        }
      }
    }

    if(!is.null(guess_df)) obj <- add_guessed_infections(obj, guess_df, s_id)

    if(!is.null(flick_df))
    {
      subject_flick_df <- flick_df[flick_df$subject_id == s_id, , drop = FALSE]

      if(nrow(subject_flick_df) > 0)
      {
        if(!nrow(subject_symp) > 0)
        {
          stop(paste0("Subject id: ", s_id,
                      ". There are no infections to flick"), call. = FALSE)
        }

        obj <- add_flicked_infections(obj, flick_df, s_id, subject_symp)
      }
    }


    obj
  })
}

check_df_colnames <- function(df, required_cols)
{
  missing_cols  <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from the data frame: ",
      paste(missing_cols, collapse = ", "))
  }
}

create_measurement_object <- function(subject_titre_df)
{
  df_list <- split(subject_titre_df, subject_titre_df$marker)

  lapply(df_list, \(df){
    I(df$meas)
  })
}

add_guessed_infections <- function(obj, guess_df, s_id)
{
  subject_guess_df<- guess_df[guess_df$subject_id == s_id, , drop = FALSE]

  if(nrow(subject_guess_df) > 0)
  {
    obj$guessed_infections <- I(as.list(sort(subject_guess_df$time)))
  }

  obj
}

add_flicked_infections <- function(obj, flick_df, s_id, subject_symp)
{
  subject_flick_df <- flick_df[flick_df$subject_id == s_id, , drop = FALSE]

  if(nrow(subject_flick_df) > 0)
  {
    invalid <- setdiff(subject_flick_df$time, subject_symp$time)

    if (length(invalid) > 0)
    {
      stop(paste0("Subject id: ", s_id,
                  ". Infections to flick don't correspond to symp infections: ",
                  paste(invalid, collapse = ", ")),
           call. = FALSE)
    }

    obj$flicks <- I(as.list(sort(subject_flick_df$time)))
  }

  obj
}
