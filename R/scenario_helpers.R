#' Calculate a confidence interval
#'
#' @param x
#' @param level
#'
#' @export
ci <- function(x, level = 95) {
  z <- 1.0 - (1.0 - (level / 100)) / 2

  stats::qnorm(z) * (mean(x) / sqrt(length(x)))
}

#' Get costs from a list of model outputs
#'
#' @param l
#' @param names
#' @param summarise_as
#'
#' @export
get_costs_list <- function(l,
                           names = c("lower", "mean", "upper"),
                           summarise_as = "domain") {
  checkmate::assert_character(
    names,
    len = length(l)
  )
  checkmate::assert_choice(
    summarise_as, c("domain", "total")
  )
  costs <- lapply(l, get_costs, summarise_as)

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

#' Get epidemic summary from a list of model outputs
#'
#' @param l
#' @param names
#' @param ...
#'
#' @export
get_summary_list <- function(l,
                             names = c("lower", "mean", "upper"),
                             ...) {
  summary_data <- Map(
    l, names,
    f = function(x, n) {
      df <- get_epidemic_summary(x, ...)
      df$tag <- n

      df
    }
  )

  dt <- data.table::rbindlist(summary_data)

  data.table::setDF(dt)
}

# #' Get incidence from a list of model outputs
# #'
# #' @export
# get_epidata_list <- function(l,
#                              names = c("lower", "mean", "upper"),
#                              ...) {
#   df_list <- lapply(l, f = function(x) {
#     get_incidence(x, ...)
#   })

#   #' Title
#   #'
#   #' @param x
#   #'
#   #' @return
#   #' @export
#   #'
#   #' @examples
#   total_hosp_list <- lapply(l, f = function(x) {
#     z <- get_data(x)
#     z <- dplyr::group_by(
#       z, "time"
#     ) %>%
#       dplyr::filter(compartment == "hospitalised") %>%
#       dplyr::summarise(
#         total_hosp = sum(value)
#       )

#     z
#   })

#   # df_list = Map(df_list, total_hosp_list, names,
#   #               f = function(x, y, n) {
#   #
#   #               })

#   dplyr::bind_rows(df_list)
# }

#' Run multiple DAEDALUS scenarios
#'
#' @param country
#' @param disease_x
#' @param response
#' @param response_time_start
#' @param response_time_end
#' @param duration
#' @param ...
#'
#' @export
run_scenarios <- function(country,
                          disease_x,
                          response = c("none", "elimination"),
                          response_time_start = 0,
                          response_time_end = 0,
                          duration = 100,
                          ...) {
  # TODO: add input checking

  scenarios <- data.table::CJ(
    response = response, duration = duration
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
      daedalus_rtm(
        country, disease_x,
        time_end = dur,
        response_time_start = response_time_start,
        response_time_end = response_time_end,
        response_strategy = resp,
        ...
      )
    }
  )

  # returned as a `<data.table>` to allow list-columns
  scenarios
}

#' Get summary data from DAEDALUS scenarios
#'
#' @param df
#' @param disease_tags
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_summary_data <- function(dt, disease_tags = "default", ...) {
  # TODO: add input checks, dt must be a data.table
  dt$epi_summary <- lapply(
    dt$output, get_summary_list, disease_tags, ...
  )

  # nolint start
  # prevent linting data.table syntax for column selection
  cols_to_keep <- setdiff(colnames(dt), "output")
  dt <- dt[, ..cols_to_keep]
  # nolint end

  dt <- dt[, unlist(dt$epi_summary, recursive = FALSE),
    by = c("response", "duration")
  ]

  data.table::dcast(
    dt,
    response + duration + measure ~ tag,
    value.var = "value"
  )
}

#' Get cost data from DAEDALUS scenarios
#'
#' @param dt
#' @param disease_tags
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_cost_data <- function(dt, disease_tags = "default") {
  # TODO: add input checks
  # NOTE: assume names taken from infection list
  # NOTE: dt must be a data.table for list-columns

  dt$costs <- lapply(dt$output, get_costs_list, disease_tags, "domain")

  dt <- dt[, unlist(dt$costs, recursive = FALSE),
    by = c("response", "duration")
  ]

  data.table::dcast(
    dt,
    response + duration + domain ~ tag,
    value.var = "cost"
  )
}
