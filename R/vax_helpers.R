#' Check multi-vaccination inputs
#'
#' @description
#' A helper function to check whether a list of vaccination objects has been
#' passed as input.
#'
#' @param x A list to be checked as a list of `<daedalus_vaccination>`.
#' Must have a minimum length of 2 elements.
#'
#' @return
#' Either `x`, if it is a list of `<daedalus_vaccination>`, or an error
#' side-effect if not.
#'
#' @keywords internal
validate_vax_list_input <- function(x) {
  is_good_npi_list <- checkmate::test_list(
    x,
    c("daedalus_vaccination", "NULL"),
    min.len = 2
  )

  if (!is_good_npi_list) {
    cli::cli_abort(
      "{.code x} may only be a list of >= 2 {.cls daedalus_vaccination}, \\
      but other classes were found, or the list has only one element."
    )
  }

  x
}
