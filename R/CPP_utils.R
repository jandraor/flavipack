#' Create a subject-level input list from multiple data frames
#'
#' This function constructs a list of subject-specific inputs by combining
#' participant metadata, titre measurements, symptom (infection) times,
#' and optional vaccination data. Each subject is returned as a list element
#' containing their longitudinal measurements and associated covariates.
#'
#' @param part_df A data frame containing participant-level information.
#'   Must include the columns \code{subject_id}, \code{age_enrolment},
#'   \code{location}, and \code{serostatus}.
#' @param titre_df A data frame containing titre measurements with at least
#'   \code{subject_id}, \code{time}, and \code{meas} columns.
#' @param symp_df A data frame containing symptom (e.g., infection) times with
#'   columns \code{subject_id} and \code{time}.
#' @param vacc_df An optional data frame containing vaccination times with
#'   columns \code{subject_id} and \code{time}. Defaults to \code{NULL}.
#' @param n_markers An integer specifying the number of markers. (Currently
#'   unused but retained for compatibility or future extensions.)
#'
#' @return A list where each element corresponds to a subject. Each subject is
#'   represented as a list with the following components:
#'   \describe{
#'     \item{subject_id}{Unique subject identifier}
#'     \item{age_enrolment}{Age at enrolment}
#'     \item{location}{Subject location}
#'     \item{serostatus}{Baseline serostatus}
#'     \item{obs_times}{Vector of observation times from \code{titre_df}}
#'     \item{measurements}{List containing measurement vectors}
#'     \item{infection_times}{(Optional) Sorted vector of infection times}
#'     \item{vaccination_times}{(Optional) Sorted vector of vaccination times}
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
#'   meas = c(10, 20, 15)
#' )
#'
#' symp_df <- data.frame(
#'   subject_id = c(1),
#'   time = c(0.5)
#' )
#'
#' create_input_list(part_df, titre_df, symp_df)
#'
#' @export
create_input_list <- function(part_df, titre_df, symp_df, vacc_df = NULL,
                              n_markers)
{
  required_cols <- c("subject_id", "age_enrolment", "location", "serostatus")

  missing_cols  <- setdiff(required_cols, names(part_df))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from the data frame: ",
      paste(missing_cols, collapse = ", "))
  }

  subject_ids <- unique(part_df$subject_id)

  lapply(subject_ids, \(s_id) {

    subject_titre <- titre_df[titre_df$subject_id == s_id, ]

    subject_info <- part_df[part_df$subject_id == s_id, ]

    subject_symp <- symp_df[symp_df$subject_id == s_id, ]

    obs_times <- I(as.list(subject_titre$time))

    measurements <- lapply(subject_titre$meas, \(x) I(list(x)))

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

    if (!is.null(vacc_df))
    {
      subject_vacc <- vacc_df[vacc_df$subject_id == s_id, , drop = FALSE]

      if (nrow(subject_vacc) > 0)
      {
        dropout_time <- max(subject_titre$time)

        subject_vacc <- subject_vacc[
          subject_vacc$time <= dropout_time,
          ,
          drop = FALSE]

        if (nrow(subject_vacc) > 0)
        {
          obj$vaccination_times <- I(as.list(sort(subject_vacc$time)))
        }
      }
    }

    obj
  })
}
