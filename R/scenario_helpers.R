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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
get_epidata_list <- function(l, names) {
  # no input checks on internal functions
  compartment <- NULL
  value <- NULL
  df_list <- lapply(l, daedalus::get_incidence)

  total_hosp_list <- lapply(l, function(x) {
    z <- daedalus::get_data(x)
    data.table::setDT(z)

    # NOTE: pass option of summarising by age group
    z <- z[compartment == "hospitalised", list(value = sum(value)), by = "time"]
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
#' @keywords internal
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
#' @param time_end A vector of integer-ish numbers giving the durations over
#' which to run scenarios. Each scenario is run for each `time_end`. The
#' intention is to be able to run response scenarios for different durations,
#' to be able to generate a time-series of how different costs accumulate.
#'
#' @return A `<data.table>`, potentially with data held in list-columns.
#'
#' @export
run_scenarios <- function(
  country,
  infection,
  response_strategy = "none",
  response_time = 30,
  response_duration = 365,
  time_end = 100,
  initial_state_manual = NULL
) {
  # daedalus::daedalus_* should bubble up input errors

  # handle custom and pre-defined response scenarios
  # this function allows passing a list of responses
  if (is.list(response_strategy)) {
    resp_predef <- unlist(Filter(is.character, response_strategy))
    checkmate::assert_subset(resp_predef, names(daedalus.data::closure_data))
    custom_responses <- Filter(is.numeric, response_strategy)
    resp_names_custom <- names(custom_responses)

    if (is.null(resp_names_custom)) {
      resp_names_custom <- glue::glue(
        "custom_response_{seq_along(custom_responses)}"
      )
    }
  }
  scenarios <- data.table::CJ(
    response = response_strategy,
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
    scenarios$time_end,
    f = function(resp, t_end) {
      daedalus::daedalus_multi_infection(
        country,
        infection,
        response_time = response_time,
        response_duration = response_duration,
        response_strategy = resp,
        time_end = t_end,
        initial_state_manual = initial_state_manual
      )
    }
  )

  # replace raw response strings and numerics with names
  if (is.list(response_strategy)) {
    scenarios$response <- c(resp_predef, resp_names_custom)
  }

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
    by = c("response", "time_end")
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

  dt <- dt[, unlist(costs, recursive = FALSE), by = c("response", "time_end")]

  dt <- switch(
    format,
    long = {
      dt
    },
    wide = {
      data.table::dcast(
        dt,
        response + domain + time_end ~ tag,
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
#' @export
get_econ_cost_data <- function(dt) {
  # NOTE: add input checks
  # NOTE: dt must be a data.table for list-columns
  # NOTE: no tracking of infection tags, no option to return wide format

  econ_costs <- NULL
  dt$econ_costs <- lapply(dt$output, get_econ_costs_list)

  dt <- dt[,
    unlist(econ_costs, recursive = FALSE),
    by = c("response", "time_end")
  ]

  # switch to long format
  dt <- data.table::melt(
    dt,
    id.vars = c("response", "time_end"),
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

  dt <- dt[, unlist(epidata, recursive = FALSE), by = c("response", "time_end")]

  dt <- switch(format, long = dt, wide = {
    data.table::dcast(
      dt,
      time + response + measure + time_end ~ tag,
      value.var = "value"
    )
  })

  data.table::setDF(dt)

  dt
}
