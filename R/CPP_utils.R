create_df_for_CPP <- function(titre_df, part_df, symp_df)
{
  cpp_df <- titre_df |> left_join(part_df) |>
    group_by(subject_id) |>
    arrange(subject_id) |>
    mutate(obs_id = row_number()) |>
    select(subject_id, age_enrolment, country, serostatus,
           obs_id, time, meas) |>
    pivot_wider(
      names_from = obs_id,
      values_from = c(time, meas),
      names_glue = "{.value}_{obs_id}") |>
    ungroup()

  col_order <- c("subject_id", "age_enrolment", "country", "serostatus")
  n_pairs <- (ncol(cpp_df) - 4) / 2

  for (i in seq_len(n_pairs)) {
    col_order <- c(col_order, paste0("time_", i), paste0("meas_", i))
  }

  cpp_df <- cpp_df[, col_order]

  symp_dates_df <- format_symp_infections(symp_df)

  # for (i in seq_len(n_pairs)) {
  #   col_order <- c(col_order, paste0("date_", i), paste0("meas_", i))
  # }

  cpp_df <- cpp_df |>
    left_join(symp_dates_df , by = "subject_id") |>
    mutate(across(everything(), ~ replace_na(., -1)))
}

format_symp_infections <- function(symp_df)
{
  required_cols <- c("subject_id", "infection_time")
  missing_cols  <- setdiff(required_cols, names(symp_df))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from the data frame: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  symp_df$obs_id <- paste0(
    "symp_",
    ave(symp_df$subject_id, symp_df$subject_id, FUN = seq_along)
  )

  symp_df$code <- 2

  wide_df <- tidyr::pivot_wider(symp_df,
                         id_cols = "subject_id",
                         names_from = "obs_id",
                         values_from = c("infection_time", "code"),
                         names_sep = "_")

  wide_df <- wide_df[, c("subject_id",
                      "infection_time_symp_1",
                      "code_symp_1",
                      "infection_time_symp_2",
                      "code_symp_2")]

  as.data.frame(wide_df)
}

#' Create a wide input data frame for model fitting
#'
#' Combines participant-level data, titre measurements, and symptomatic infection
#' data into a single wide-format data frame. For each participant, titre
#' measurements are stored as alternating time and measurement columns, followed
#' by symptomatic infection times and their corresponding infection code.
#'
#' Missing measurement and infection entries are padded with `-1`.
#'
#' @param part_df A data frame containing participant-level information. Must
#'   include the columns `"subject_id"`, `"age_enrolment"`, `"country"`, and
#'   `"serostatus"`.
#' @param titre_df A data frame containing titre measurements for each
#'   participant. Must include `"subject_id"`, `"time"`, and `"meas"` columns.
#' @param symp_df A data frame containing symptomatic infection records for each
#'   participant. Must include `"subject_id"` and `"time"` columns.
#' @param max_n_meas Integer. The maximum number of titre measurements to
#'   allocate per participant.
#' @param max_n_symp Integer. The maximum number of symptomatic infections to
#'   allocate per participant.
#'
#' @return A wide-format data frame with:
#' \describe{
#'   \item{subject_id, age_enrolment, country, serostatus}{Participant-level
#'   covariates from `part_df`.}
#'   \item{time_1, meas_1, ..., time_n, meas_n}{Titre measurement times and
#'   values, up to `max_n_meas`.}
#'   \item{infection_time_symp_1, code_symp_1, ..., infection_time_symp_n,
#'   code_symp_n}{Symptomatic infection times and codes, up to
#'   `max_n_symp`. Infection codes are set to `2`.}
#' }
#'
#' Empty entries are filled with `-1`.
#'
#' @details
#' The output matrix has `2 * (max_n_meas + max_n_symp)` columns:
#' two columns per titre measurement (`time`, `meas`) and two columns per
#' symptomatic infection (`infection_time`, `code`).
#'
#' Titre measurements are written first for each participant, followed by
#' symptomatic infections.
#'
#' @examples
#' part_df <- data.frame(
#'   subject_id = c(1, 2),
#'   age_enrolment = c(34, 29),
#'   country = c("UK", "FR"),
#'   serostatus = c(1, 0)
#' )
#'
#' titre_df <- data.frame(
#'   subject_id = c(1, 1, 2),
#'   time = c(0, 30, 0),
#'   meas = c(10, 20, 15)
#' )
#'
#' symp_df <- data.frame(
#'   subject_id = c(1, 2),
#'   time = c(45, 60)
#' )
#'
#' create_input_df(
#'   part_df = part_df,
#'   titre_df = titre_df,
#'   symp_df = symp_df,
#'   max_n_meas = 2,
#'   max_n_symp = 1
#' )
#'
#' @export
create_input_df <- function(part_df, titre_df, symp_df, max_n_meas, max_n_symp)
{
  required_cols <- c("subject_id", "age_enrolment", "country", "serostatus")

  missing_cols  <- setdiff(required_cols, names(part_df))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from the data frame: ",
      paste(missing_cols, collapse = ", "))
  }

  n_individuals <- nrow(part_df)

  df1 <- part_df[, required_cols]

  m <- matrix(NA, nrow = n_individuals, ncol = 2 * (max_n_meas + max_n_symp))

  subject_ids <- df1$subject_id

  for(i in seq_along(subject_ids))
  {
    s_id <- subject_ids[i]

    subject_titre_df <- titre_df[titre_df$subject_id == s_id, ]

    n_meas <- nrow(subject_titre_df)

    index <- 1

    for(j in seq_len(max_n_meas))
    {
      m[i, index] <- subject_titre_df$time[j]
      index       <- index + 1

      m[i, index] <- subject_titre_df$meas[j]
      index       <- index + 1
    }

    subject_symp_df <- symp_df[symp_df$subject_id == s_id, ]

    n_symp_inf <- nrow(subject_symp_df)

    if(n_symp_inf > 0)
    {
      for(k in seq_len(n_symp_inf))
      {
        m[i, index] <- subject_symp_df$time[k]
        index       <- index + 1

        m[i, index] <- 2
        index       <- index + 1
      }
    }
  }

  m[is.na(m)] <- -1

  df2 <- as.data.frame(m)

  names_titres   <- paste0(c("time_", "meas_"), rep(1:max_n_meas, each = 2))
  names_symp_inf <- paste0(c("infection_time_symp_", "code_symp_"),
                           rep(1:max_n_symp, each = 2))

  colnames(df2) <- c(names_titres, names_symp_inf)

  df <- cbind(df1, df2)
}
