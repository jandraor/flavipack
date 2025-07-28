#' Scaled log2-transform of antibody titres
#'
#' Applies a log2 transformation to antibody titres relative to a baseline of 10,
#'  returning 1 for a titre of 10 and incrementing by 1 for each doubling.
#' @param titre Numeric vector of titres.
#'
#' @return Numeric vector of log2-transformed titres.
#' @export
#'
#' @examples
#' log2_transform(10)   # 1
#' log2_transform(20)   # 2
log2_transform <- function(titre)
{
  1 + log(titre / 10) / log(2)
}
