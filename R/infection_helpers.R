#' Generate multiple `<daedalus_infection>`s from parameter distributions
#'
#' @param name An infection name from among `daedalus::epidemic_names`.
#'
#' @param samples The number of samples to generate.
#'
#' @param param_distributions A named list of `<distribution>` class
#' objects provided by \pkg{distributional}, with names corresponding to
#' `daedalus::infection_parameter_names`. These are used to generate `samples`
#' draws from each distribution for the corresponding infection parameters.
#' Arguments which vary by age are supported, with the drawn and scaled value
#' treated as the median of the profile vector.
#' See **Examples**.
#'
#' @param param_ranges An optional named list of two-element vectors, giving the
#' ranges to which samples drawn using `param_distributions` should be rescaled.
#' This allows easy rescaling of values drawn from bounded distributions (e.g.
#' the Beta distribution).
#' If passed, the list must have names corresponding to `param_distribution`
#' names. When empty, the generated values are used without scaling.
#' See **Examples**.
#'
#' @return A list of `<daedalus_infection>` objects.
#'
#' @export
#'
#' @examples
#' # make 10 infection objects varying in R0, with a skewed distribution
#' # scaled between 1.0 and 2.0
#' make_infection_samples(
#'   "influenza_2009",
#'   samples = 3,
#'   list(
#'     r0 = distributional::dist_beta(2, 5)
#'   ),
#'   list(
#'     r0 = c(0.1, 0.2)
#'   )
#' )
make_infection_samples <- function(name,
                                   param_distributions,
                                   param_ranges = NULL,
                                   samples = 100) {
  # check inputs
  checkmate::assert_subset(
    name, daedalus::epidemic_names,
    empty.ok = FALSE
  )
  checkmate::assert_count(samples, positive = TRUE)

  # check parameter distributions
  checkmate::assert_list(
    param_distributions, "distribution",
    any.missing = FALSE
  )
  param_names <- names(param_distributions)
  checkmate::assert_subset(
    param_names,
    daedalus::infection_parameter_names
  )

  # check parameter rescale ranges
  if (!is.null(param_ranges)) {
    checkmate::assert_list(param_ranges, "numeric")
    invisible(
      lapply(
        param_ranges, checkmate::assert_numeric,
        lower = 0, finite = TRUE, any.missing = FALSE
      )
    )
    param_rescale_names <- names(param_ranges)
    checkmate::assert_subset(param_rescale_names, param_names)
  }

  # create default values for use
  infection_default <- daedalus::daedalus_infection(name)

  # generate draws from distributions
  param_distributions_list <- Reduce(param_distributions, f = `c`)
  params_list <- distributional::generate(param_distributions_list, samples)
  names(params_list) <- param_names

  # scale params to range if provided
  # NOTE: ORDER OF RANGE VECTOR MATTERS, reverse scaling possible
  if (length(param_ranges) > 0) {
    params_list[names(param_ranges)] <- Map(
      params_list[names(param_ranges)], param_ranges,
      f = scales::rescale
    )
  }

  are_good_params <- vapply(
    params_list, checkmate::test_numeric, logical(1),
    lower = 0, finite = TRUE, any.missing = FALSE
  )
  no_negs_generated <- all(are_good_params)
  which_has_negs <- which(!are_good_params) # nolint lintr cannot parse glue xpr

  if (!no_negs_generated) {
    cli::cli_abort(
      "Negative values generated while sampling these parameters:
        {.str {param_names[which_has_negs]}}.
        Check distribution parameterisation."
    )
  }

  # handle numeric vectors if any using profile of default and mean
  # NOTE: consider standardising as a function, see scales::rescale_mid()
  if (any(param_names %in% NAMES_VECTOR_INF_PARAMS)) {
    # get subset of names provided
    vector_params <- intersect(param_names, NAMES_VECTOR_INF_PARAMS)
    params_list[vector_params] <- lapply(
      vector_params, function(x) {
        lapply(params_list[[x]], function(vl) {
          value_profile <- daedalus::get_data(infection_default, x)
          value_profile <- value_profile / mean(value_profile)
          vl * value_profile
        })
      }
    )
  }

  params_list <- purrr::list_transpose(params_list)

  infection_list <- lapply(
    params_list, function(pl) {
      do.call(
        daedalus::daedalus_infection,
        c(
          list(name = name),
          pl
        )
      )
    }
  )

  infection_list
}
