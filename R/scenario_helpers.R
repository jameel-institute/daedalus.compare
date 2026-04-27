#' Calculate an interval
#'
#' @param x A numeric vector.
#' @param level The interval level as a numeric vector  with values in the range
#' \eqn{[0, 100]}. Defaults to `95` for the 95% interval.
#'
#' @return A numeric vector of the same length as `level`.
#' @export
ci <- function(x, level = 95) {
  checkmate::assert_numeric(x, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(level, 0, 100, finite = TRUE, any.missing = FALSE)

  z <- 1.0 - (1.0 - (level / 100)) / 2

  stats::qnorm(z) * (mean(x) / sqrt(length(x)))
}

#' Get pandemic costs from a list of model outputs
#'
#' @description See `daedalus::get_costs()` for more information.
#'
#' @param l A list of `<daedalus_output>` objects; each object will be passed to
#' `daedalus::get_costs()`.
#'
#' @param names A character vector, intended to apply to elements of `l`.
#'
#' @return A `<data.table>`.
#'
#' @keywords internal output_helpers
get_costs_list <- function(l, names) {
  # no input checks on internal functions
  costs <- lapply(l, daedalus::get_costs, summarise_as = "domain")

  costs <- do.call(rbind, costs)
  cost_dt <- data.table::as.data.table(costs)
  cost_dt$tag <- names

  data.table::melt(
    cost_dt,
    id.vars = "tag",
    variable.name = "domain",
    value.name = "cost"
  )
}

#' Get epidemic summary data from a list of model outputs
#'
#' @inheritParams get_costs_list
#'
#' @param ... Additional arguments passed to [daedalus::get_epidemic_summary()].
#'
#' @return A `<data.table>`.
#'
#' @keywords internal output_helpers
get_summary_list <- function(l, names, ...) {
  # no input checks on internal functions
  summary_dt <- Map(
    l,
    names,
    f = function(x, n) {
      df <- daedalus::get_epidemic_summary(x, ...)
      df$tag <- n

      df
    }
  )

  data.table::rbindlist(summary_dt)
}

#' Get incidence data from a list of model outputs
#'
#' @inheritParams get_costs_list
#'
#' @return A `<data.table>`.
#'
#' @keywords internal output_helpers
get_epidata_list <- function(l, names) {
  # no input checks on internal functions
  compartment <- NULL
  value <- NULL
  df_list <- lapply(l, daedalus::get_incidence)

  total_hosp_list <- lapply(l, function(x) {
    z <- daedalus::get_data(x)
    data.table::setDT(z)

    # NOTE: pass option of summarising by age group
    z <- z[
      data.table::like(compartment, "hospitalised"),
      list(value = sum(value)),
      by = "time"
    ]
    z$measure <- "total_hosp"

    z
  })

  hosp_dt <- Map(
    df_list,
    total_hosp_list,
    names,
    f = function(incidence_df, hosp_df, n) {
      dt <- data.table::rbindlist(
        list(incidence_df, hosp_df)
      )
      dt$tag <- n

      dt
    }
  )

  data.table::rbindlist(hosp_dt)
}

#' Get details of economic costs from model outputs
#'
#' @param l A list of `<daedalus.output>` objects.
#'
#' @return A `<data.table>`.
#'
#' @keywords internal output_helpers
get_econ_costs_list <- function(l) {
  # no input checks on internal functions
  econ_costs_list <- lapply(l, function(x) {
    z <- daedalus::get_costs(x)[["economic_costs"]]
    as.data.frame(z[c("economic_cost_closures", "economic_cost_absences")])
  })

  data.table::rbindlist(econ_costs_list)
}

#' Run multiple DAEDALUS scenarios
#'
#' @inheritParams daedalus::daedalus_multi_infection
#'
#' @param response_strategy A time-limited response strategy specified as a
#' single `<daedalus_npi>` created using [daedalus::daedalus_timed_npi()],
#' or a list of such objects, or `NULL` for no response. Lists may include
#' `NULL`, which is useful when comparing response scenarios against the
#' no-response counterfactual.
#' If a list is passed, all elements must be named, or the list may have no
#' names, in which case synthetic names will be assigned and a message printed
#' to screen. Lists with some elements named and some unnamed are not accepted.
#'
#' @param vaccination_strategy A vaccination strategy specified as a
#' single `<daedalus_vaccination>` created using
#' [daedalus::daedalus_vaccination()], or a list of such objects, or `NULL`
#' for no response. Lists may include `NULL`, which is useful when comparing
#' vaccination scenarios against the no-vaccination counterfactual.
#' If a list is passed, all elements must be named, or the list may have no
#' names, in which case synthetic names will be assigned and a message printed
#' to screen. Lists with some elements named and some unnamed are not accepted.
#'
#' @param time_end A vector of integer-ish numbers giving the durations over
#' which to run scenarios. Each scenario is run for each `time_end`. The
#' intention is to be able to run response scenarios for different durations,
#' to be able to generate a time-series of how different costs accumulate.
#'
#' @return A `<data.table>`, with data held in list-columns.
#'
#' @export
run_scenarios <- function(
  country,
  infection,
  response_strategy = NULL,
  vaccination_strategy = NULL,
  time_end = 100,
  initial_state_manual = NULL
) {
  # daedalus::daedalus_* should bubble up input errors

  # handle response objects
  # this function allows passing a list of responses
  is_good_response <- checkmate::test_multi_class(
    response_strategy,
    c("list", "daedalus_npi"),
    null.ok = TRUE
  )
  if (!is_good_response) {
    cli::cli_abort(
      "Got {.code response_strategy} of class \\
      {.cls {class(response_strategy)}}, but only \\
      {.cls daedalus_npi}, {.cls list}, or {.code NULL} are allowed."
    )
  }

  # handle possible inputs and generate names
  if (daedalus::is_daedalus_npi(response_strategy)) {
    is_timed_npi <- response_strategy$identifier == "custom_timed"
    if (!is_timed_npi) {
      cli::cli_abort(
        c(
          "{.code x} is a {.cls daedalus_npi}, but it is \\
          reactive to model state; expected a time-limited NPI!",
          i = "Use {.fn daedalus::daedalus_timed_npi} to create a \\
          timed-limited NPI."
        )
      )
    }
    response_strategy <- list(response_strategy)
    response_names <- "custom_timed"
  } else if (is.list(response_strategy)) {
    response_strategy <- validate_npi_list_input(response_strategy)

    has_good_names <- checkmate::test_named(response_strategy, "unique")
    if (has_good_names) {
      response_names <- names(response_strategy)
    } else {
      cli::cli_inform(
        c(
          "Got {.code response_strategy} as a list without names! Assigning \\
          manually constructed names of the form {.code response_*}.",
          i = "Assigning names to strategies helps with identifying scenarios!"
        )
      )
      response_names <- sprintf("response_%i", seq_along(response_strategy))
    }
  } else {
    response_strategy <- list(response_strategy)
    response_names <- "none"
  }

  # handle vaccination objects
  # this function allows passing a list of vaccinations
  is_good_vax <- checkmate::test_multi_class(
    vaccination_strategy,
    c("list", "daedalus_vaccination"),
    null.ok = TRUE
  )
  if (!is_good_vax) {
    cli::cli_abort(
      "Got {.code vaccination_strategy} of class \\
      {.cls {class(vaccination_strategy)}}, but only \\
      {.cls daedalus_vaccination}, {.cls list}, or {.code NULL} are allowed."
    )
  }

  # handle possible inputs and generate names
  if (daedalus::is_daedalus_vaccination(vaccination_strategy)) {
    vaccination_strategy <- list(vaccination_strategy)
    vaccination_names <- "custom_vaccination"
  } else if (is.list(vaccination_strategy)) {
    vaccination_strategy <- validate_vax_list_input(vaccination_strategy)

    has_good_names <- checkmate::test_named(vaccination_strategy, "unique")
    if (has_good_names) {
      vaccination_names <- names(vaccination_strategy)
    } else {
      cli::cli_inform(
        c(
          "Got {.code vaccination_strategy} as a list without names! \\
          Assigning manually constructed names of the form \\
          {.code vaccination_*}.",
          i = "Assigning names to vaccination strategies helps with \\
          identifying scenarios!"
        )
      )
      vaccination_names <- sprintf(
        "vaccination_%i",
        seq_along(vaccination_strategy)
      )
    }
  } else {
    vaccination_strategy <- list(vaccination_strategy)
    vaccination_names <- "none"
  }

  scenarios <- data.table::CJ(
    response = response_strategy,
    vaccination = vaccination_strategy,
    time_end = time_end,
    sorted = FALSE
  )

  # get name combinations grid
  scenario_names <- data.table::CJ(
    response = response_names,
    vaccination = vaccination_names,
    time_end = time_end,
    sorted = FALSE
  )

  # get range names from disease_x or synthetic names
  # assuming 4-digit list lengths
  disease_tags <- names(infection)
  if (is.null(disease_tags)) {
    disease_tags <- formatC(
      seq_along(infection),
      flag = "0",
      width = floor(log10(length(infection))) + 1
    )
  }

  scenarios$output <- Map(
    scenarios$response,
    scenarios$vaccination,
    scenarios$time_end,
    f = function(resp, vaccination, t_end) {
      daedalus::daedalus_multi_infection(
        country,
        infection,
        response_strategy = resp,
        vaccine_investment = vaccination,
        time_end = t_end,
        initial_state_manual = initial_state_manual
      )
    }
  )

  # replace raw response strings and numerics with names
  # and apply a conversion to character vector from list
  scenarios$response <- scenario_names$response
  scenarios$vaccination <- scenario_names$vaccination

  # returned as a `<data.table>` to allow list-columns
  scenarios
}

#' Get summary data from DAEDALUS scenarios
#'
#' @param dt A `<data.table>` resulting from `run_scenarios()`.
#'
#' @param disease_tags A character vector giving names for replicates within
#' each scenario.
#'
#' @param format A string for whether the data should be returned in `"long"` or
#' `"wide"` format. Returning a 'wide' data.frame is only recommended when there
#' is a small number of `disease_tags`, typically a triplet of "low", "medium",
#' and "high" risk scenarios.
#'
#' @param ... Additional arguments passed to [daedalus::get_epidemic_summary()].
#' If passing `groups`, only `format = "long"` is supported.
#'
#' @return A `<data.frame>` summarising epidemiological outcomes: cumulative
#' infections, deaths, and hospitalisations over the timeframe of the modelled
#' epidemics.
#'
#' @keywords output_helpers
#'
#' @export
get_summary_data <- function(
  dt,
  disease_tags,
  format = c("long", "wide"),
  ...
) {
  # input checks
  checkmate::assert_data_table(dt, any.missing = FALSE)
  checkmate::assert_character(disease_tags)

  epi_summary <- NULL
  dt$epi_summary <- lapply(
    dt$output,
    get_summary_list,
    disease_tags,
    ...
  )

  format <- rlang::arg_match(format)
  # nolint start
  # prevent linting data.table syntax for column selection
  cols_to_keep <- setdiff(colnames(dt), "output")
  dt <- dt[, cols_to_keep, with = FALSE]
  # nolint end

  dt <- dt[,
    unlist(epi_summary, recursive = FALSE),
    by = c("response", "vaccination", "time_end")
  ]

  dt <- switch(
    format,
    long = {
      dt
    },
    wide = {
      data.table::dcast(
        dt,
        response + time_end + measure ~ tag,
        value.var = "value"
      )
    }
  )

  data.table::setDF(dt)

  dt
}

#' Get cost data from DAEDALUS scenarios
#'
#' @inheritParams get_summary_data
#'
#' @return A `<data.frame>` summarising epidemic costs over the timeframe of the
#' modelled epidemics.
#'
#' @return A `<data.frame>` of the costs for each model scenario.
#'
#' @keywords output_helpers
#'
#' @export
get_cost_data <- function(
  dt,
  disease_tags = "default",
  format = c("long", "wide")
) {
  # NOTE: add input checks
  checkmate::assert_data_table(dt)
  checkmate::assert_character(disease_tags)

  format <- rlang::arg_match(format)
  # NOTE: assume names taken from infection list
  # NOTE: dt must be a data.table for list-columns

  costs <- NULL
  dt$costs <- lapply(dt$output, get_costs_list, disease_tags)

  dt <- dt[,
    unlist(costs, recursive = FALSE),
    by = c("response", "vaccination", "time_end")
  ]

  dt <- switch(
    format,
    long = {
      dt
    },
    wide = {
      data.table::dcast(
        dt,
        response + vaccination + domain + time_end ~ tag,
        value.var = "cost"
      )
    }
  )

  data.table::setDF(dt)

  dt
}

#' Get economic cost details from an output data.table
#'
#' @inheritParams get_summary_data
#'
#' @return A `<data.frame>` summarising economic costs over the timeframe of the
#' modelled epidemics, broken down into costs due to restrictions and illness-
#' related absences.
#'
#' @keywords output_helpers
#'
#' @export
get_econ_cost_data <- function(dt) {
  # NOTE: add input checks
  # NOTE: dt must be a data.table for list-columns
  # NOTE: no tracking of infection tags, no option to return wide format

  econ_costs <- NULL
  dt$econ_costs <- lapply(dt$output, get_econ_costs_list)

  dt <- dt[,
    unlist(econ_costs, recursive = FALSE),
    by = c("response", "vaccination", "time_end")
  ]

  # switch to long format
  dt <- data.table::melt(
    dt,
    id.vars = c("response", "vaccination", "time_end"),
    value.name = "cost",
    variable.name = "cost_type"
  )

  data.table::setDF(dt)

  dt
}

#' Get epidemiological curves from model data
#'
#' @inheritParams get_summary_data
#'
#' @return A `<data.frame>` of individuals in each epidemiological compartments
#' over the timeframe of the modelled epidemics.
#'
#' @keywords output_helpers
#'
#' @export
get_epicurve_data <- function(
  dt,
  disease_tags = "default",
  format = c("long", "wide")
) {
  # NOTE: add input checks
  # NOTE: assume names taken from infection list
  # NOTE: dt must be a data.table for list-columns

  format <- rlang::arg_match(format)
  epidata <- NULL
  dt$epidata <- lapply(dt$output, get_epidata_list, names = disease_tags)

  dt <- dt[,
    unlist(epidata, recursive = FALSE),
    by = c("response", "vaccination", "time_end")
  ]

  dt <- switch(format, long = dt, wide = {
    data.table::dcast(
      dt,
      time + response + vaccination + measure + time_end ~ tag,
      value.var = "value"
    )
  })

  data.table::setDF(dt)

  dt
}
