check_columns <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Data frame is missing required columns: ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }
}
