#' Scaled log2-transform of antibody titres
#'
#' Applies a log2 transformation to antibody titres relative to a baseline of 10.
#' A titre of 10 returns 1, and each doubling of the titre increments the value by 1.
#'
#' @param titre Numeric vector of original antibody titres.
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

#' Inverse of the scaled log2-transform of antibody titres
#'
#' Converts log2-transformed titres back to the original scale.
#' This is the inverse of `log2_transform()`. For example, a log2-transformed
#' titre of 1 corresponds to an original titre of 10, and each increment of 1
#' doubles the original titre.
#'
#' @param log2_titre Numeric vector of log2-transformed titres.
#'
#' @return Numeric vector of original antibody titres.
#' @export
#'
#' @examples
#' inv_log2_transform(1)  # 10
#' inv_log2_transform(2)  # 20
inv_log2_transform <- function(log2_titre)
{
  exp((log2_titre - 1) * log(2) + log(10))
}
