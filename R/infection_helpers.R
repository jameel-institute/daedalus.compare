#' Generate multiple infections with parameter distributions
#'
#' @param name An infection name from among `daedalus::epidemic_names`.
#' @param samples The number of samples to generate.
#'
#' @param param_distributions A named list of `<distribution>` class
#' objects provided by \pkg{distributional}, with names corresponding to
#' `daedalus::infection_parameter_names`. These are used to generate `samples`
#' draws from each distribution for the corresponding infection parameters.
#' Arguments which vary by age are supported, with the drawn and scaled value
#' treated as the median of the profile vector.
#'
#' @param param_ranges An optional named list of two-element vectors, giving the
#' ranges to which samples drawn using `param_distributions` should be rescaled.
#' This allows easy rescaling of values drawn from bounded distributions (e.g.
#' the Beta distribution).
#' If passed, the list must have names corresponding to `param_distribution`
#' names. When empty, the generated values are used without scaling.
#'
#' @return A list of `<daedalus_infection>` objects.
#' @export
make_infection_samples <- function(name, samples = 100,
                                   param_distributions = list(
                                     r0 = distributional::dist_beta(2, 5)
                                   ),
                                   param_ranges = list(
                                     r0 = c(1.0, 2.0)
                                   )) {
  # check inputs
  checkmate::assert_subset(
    name, daedalus::epidemic_names,
    empty.ok = FALSE
  )
  checkmate::assert_count(samples, positive = TRUE)

  # create default values for use
  infection_default <- daedalus::daedalus_infection(name)

  param_names <- names(param_distributions)
  checkmate::assert_subset(
    param_names,
    # setdiff(
    daedalus::infection_parameter_names,
    #   NAMES_VECTOR_INF_PARAMS
    # )
  )

  param_rescale_names <- names(param_ranges)
  checkmate::assert_subset(param_rescale_names, param_names)

  checkmate::assert_list(
    param_distributions, "distribution",
    any.missing = FALSE
  )
  # TODO: more checks needed for range
  checkmate::assert_list(param_ranges, "numeric")

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
  which_has_negs <- which(!are_good_params)

  if (!no_negs_generated) {
    cli::cli_abort(
      "Negative values generated while sampling these parameters:
        {.str {param_names[which_has_negs]}}.
        Check distribution parameterisation."
    )
  }

  # handle numeric vectors if any using profile of default and mean
  # TODO: consider standardising as a function, see scales::rescale_mid()
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
