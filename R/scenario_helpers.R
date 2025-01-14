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
#' @param l A list of `<daedalus_output>` objects; each object will be passed to
#' `daedalus::get_costs()`.
#' @param names A character vector, intended to apply to elements of `l`.
#' @param summarise_as A string specifying whether costs are split by domain.
#' See `daedalus::get_costs()` for more information.
#'
#' @keywords internal
get_costs_list <- function(l,
                           names,
                           summarise_as = c("domain", "total")) {
  costs <- lapply(l, daedalus::get_costs, summarise_as)

  costs <- switch(summarise_as,
    domain = {
      costs <- do.call(rbind, costs)
      cost_df <- data.table::as.data.table(costs)
      cost_df$tag <- names

      dt <- data.table::melt(
        cost_df,
        id.vars = "tag",
        variable.name = "domain",
        value.name = "cost"
      )

      data.table::setDF(dt)
    },
    total = {
      data.frame(
        cost = unlist(costs),
        domain = "total",
        tag = names,
        stringsAsFactors = FALSE
      )
    }
  )

  costs
}

#' Get epidemic summary data from a list of model outputs
#'
#' @inheritParams get_costs_list
#'
#' @keywords internal
get_summary_list <- function(l,
                             names = c("lower", "mean", "upper")) {
  summary_data <- Map(
    l, names,
    f = function(x, n) {
      df <- daedalus::get_epidemic_summary(x)
      df$tag <- n

      df
    }
  )

  dt <- data.table::rbindlist(summary_data)

  data.table::setDF(dt)
}

#' Get incidence data from a list of model outputs
#'
#' @inheritParams get_costs_list
#'
#' @keywords internal
get_epidata_list <- function(l, names) {
  compartment <- NULL
  value <- NULL
  df_list <- lapply(l, daedalus::get_incidence)

  total_hosp_list <- lapply(l, function(x) {
    z <- daedalus::get_data(x)
    data.table::setDT(z)

    # TODO: pass option of summarising by age group
    z <- z[compartment == "hospitalised",
      list(value = sum(value)),
      by = "time"
    ]
    z$measure <- "total_hosp"

    z
  })

  dt_list <- Map(
    df_list, total_hosp_list, names,
    f = function(incidence_df, hosp_df, n) {
      dt <- data.table::rbindlist(
        list(incidence_df, hosp_df)
      )
      dt$tag <- n

      dt
    }
  )

  dt <- data.table::rbindlist(dt_list)

  dt
}

#' Run multiple DAEDALUS scenarios
#'
#' @inheritParams daedalus::daedalus_rtm
#' @param duration A vector of integer-ish numbers giving the durations over
#' which to run scenarios. Each scenario is run for each duration.
#'
#' @export
run_scenarios <- function(country,
                          disease_x,
                          response = c("none", "elimination"),
                          response_time_start = 0,
                          response_time_end = 0,
                          duration = 100) {
  # TODO: add input checking

  # handle custom and pre-defined response scenarios
  if (is.list(response)) {
    resp_predef <- unlist(Filter(is.character, response))
    checkmate::assert_subset(resp_predef, names(daedalus::closure_data))
    custom_responses <- Filter(is.numeric, response)
    resp_names_custom <- names(custom_responses)

    if (is.null(resp_names_custom)) {
      resp_names_custom <- glue::glue(
        "custom_response_{seq_along(custom_responses)}"
      )
    }
  }
  scenarios <- data.table::CJ(
    response = response, duration = duration, sorted = FALSE
  )

  # get range names from disease_x or synthetic names
  # assuming 4-digit list lengths
  disease_tags <- names(disease_x)
  if (is.null(disease_tags)) {
    disease_tags <- sprintf("replicate_%04i", seq_along(disease_x))
  }

  scenarios$output <- Map(
    scenarios$response, scenarios$duration,
    f = function(resp, dur) {
      daedalus::daedalus_rtm(
        country, disease_x,
        time_end = dur,
        response_time_start = response_time_start,
        response_time_end = response_time_end,
        response_strategy = resp
      )
    }
  )

  # replace raw response strings and numerics with names
  if (is.list(response)) {
    scenarios$response <- c(resp_predef, resp_names_custom)
  }

  # returned as a `<data.table>` to allow list-columns
  scenarios
}

#' Get summary data from DAEDALUS scenarios
#'
#' @param dt A `<data.table>` resulting from `run_scenarios()`.
#' @param disease_tags A character vector giving names for replicates within
#' each scenario.
#' @param format A string for whether the data should be returned in `"long"` or
#' `"wide"` format. Returning a 'wide' data.frame is only recommended when there
#' is a small number of `disease_tags`, typically a triplet of "low", "medium",
#' and "high" risk scenarios.
#'
#' @return A `<data.frame>` summarising epidemiological outcomes: cumulative
#' infections, deaths, and hospitalisations over the timeframe of the modelled
#' epidemics.
#'
#' @export
get_summary_data <- function(dt, disease_tags,
                             format = c("long", "wide")) {
  # TODO: add input checks, dt must be a data.table
  epi_summary <- NULL
  dt$epi_summary <- lapply(
    dt$output, get_summary_list, disease_tags
  )

  format <- rlang::arg_match(format)
  # nolint start
  # prevent linting data.table syntax for column selection
  cols_to_keep <- setdiff(colnames(dt), "output")
  dt <- dt[, ..cols_to_keep]
  # nolint end

  dt <- dt[, unlist(epi_summary, recursive = FALSE),
    by = c("response", "duration")
  ]

  dt <- switch(format,
    long = {
      dt
    },
    wide = {
      data.table::dcast(
        dt,
        response + duration + measure ~ tag,
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
#' @export
get_cost_data <- function(dt, disease_tags = "default",
                          format = c("long", "wide")) {
  # TODO: add input checks
  # NOTE: assume names taken from infection list
  # NOTE: dt must be a data.table for list-columns

  costs <- NULL
  dt$costs <- lapply(dt$output, get_costs_list, disease_tags, "domain")

  dt <- dt[, unlist(costs, recursive = FALSE),
    by = c("response", "duration")
  ]

  dt <- switch(format,
    long = {
      dt
    },
    wide = {
      data.table::dcast(
        dt,
        time + response + measure + duration ~ tag,
        value.var = "value"
      )
    }
  )

  data.table::setDF(dt)

  dt
}

#' Get details of economic costs from model outputs
#'
#' @param l A list of `<daedalus.output>` objects.
#'
#' @keywords internal
get_econ_costs_list <- function(l) {
  econ_costs_list <- lapply(l, function(x) {
    z <- daedalus::get_costs(x)
    z[["economic_costs"]][
      c("economic_cost_closures", "economic_cost_absences")
    ]
  })

  econ_costs_list <- lapply(econ_costs_list, as.data.frame)

  data.table::rbindlist(econ_costs_list)
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
  # TODO: add input checks
  # NOTE: dt must be a data.table for list-columns
  # NOTE: no tracking of infection tags, no option to return wide format

  econ_costs <- NULL
  dt$econ_costs <- lapply(dt$output, get_econ_costs_list)

  dt <- dt[, unlist(econ_costs, recursive = FALSE),
    by = c("response", "duration")
  ]

  # switch to long format
  dt <- data.table::melt(
    dt,
    id.vars = c("response", "duration"),
    value.name = "cost", variable.name = "cost_type"
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
get_epicurve_data <- function(dt, disease_tags = "default",
                              format = c("long", "wide")) {
  # TODO: add input checks
  # NOTE: assume names taken from infection list
  # NOTE: dt must be a data.table for list-columns

  format <- rlang::arg_match(format)
  epidata <- NULL
  dt$epidata <- lapply(dt$output, get_epidata_list, names = disease_tags)

  dt <- dt[, unlist(epidata, recursive = FALSE),
    by = c("response", "duration")
  ]

  dt <- switch(format,
    long = dt,
    wide = {
      data.table::dcast(
        dt,
        time + response + measure + duration ~ tag,
        value.var = "value"
      )
    }
  )

  data.table::setDF(dt)

  dt
}
