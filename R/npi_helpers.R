#' Check multi-NPI inputs
#'
#' @description
#' A helper function to check whether a list of NPI objects has been passed as
#' input. This could live in \pkg{daedalus} similar to
#' `daedalus:::validate_infection_list_input()`, but lives here instead.
#'
#' @param x A list to be checked as a list of `<daedalus_npi>`. Must have a
#' minimum length of 2 elements.
#'
#' @param timed_only A boolean defaulting to `TRUE` for whether only
#' time-limited NPIs are allowed. This is typically the case for a real-time
#' modelling exercise.
#'
#' @return
#' Either `x`, if it is a list of `<daedalus_npi>`, or an error side-effect if
#' not.
#'
#' @keywords internal
validate_npi_list_input <- function(x, timed_only = TRUE) {
  is_good_npi_list <- checkmate::test_list(
    x,
    c("daedalus_npi", "NULL"),
    min.len = 2
  )

  if (!is_good_npi_list) {
    cli::cli_abort(
      "{.code x} may only be a list of >= 2 {.cls daedalus_npi}, \\
      but other classes were found, or the list has only one element."
    )
  }

  if (timed_only) {
    is_timed_npi_list <- all(
      vapply(
        Filter(Negate(is.null), x),
        function(z) {
          z$identifier == "custom_timed"
        },
        FUN.VALUE = logical(1L)
      )
    )

    if (!is_timed_npi_list) {
      cli::cli_abort(
        c(
          "{.code x} is a list of {.cls daedalus_npi}, but some or all \\
          NPIs are reactive to model state; expected only time-limited NPIs!",
          i = "Use {.fn daedalus::daedalus_timed_npi} to create timed-limited \\
          NPIs."
        )
      )
    }
  }

  x
}
