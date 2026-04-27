# test basic scenario runners
test_that("daedalus.compare: Basic usage of `run_scenarios()`", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list
    )
  )

  # expect output is a data table with list columns
  output <- run_scenarios(
    "GBR",
    infection_list
  )
  checkmate::expect_data_table(
    output
  )
  checkmate::expect_list(output$output, "list")
})

test_that("Fn `get_epicurve_data()` works", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  disease_tags <- sprintf("disease_%i", seq_along(infection_list))

  cty <- "GBR"
  # single custom response
  response <- daedalus::daedalus_timed_npi(
    10,
    20,
    list(rep(0.5, daedalus:::N_ECON_SECTORS)),
    cty
  )

  output <- run_scenarios(
    cty,
    infection_list,
    list(custom = response, custom2 = response)
  )

  expect_no_condition(
    get_epicurve_data(output)
  )
  checkmate::expect_data_frame(
    get_epicurve_data(output)
  )
  expect_no_condition(
    get_epicurve_data(output, disease_tags, "wide")
  )

  epicurve_data <- get_epicurve_data(output)
  checkmate::expect_subset(
    unique(epicurve_data[["measure"]]),
    c(
      "daily_deaths",
      "daily_infections",
      "daily_hospitalisations",
      "total_hosp"
    )
  )
})

test_that("Fn `get_summary_data()` works", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  cty <- "GBR"
  # single custom response
  response <- daedalus::daedalus_timed_npi(
    10,
    20,
    list(rep(0.5, daedalus:::N_ECON_SECTORS)),
    cty
  )

  output <- run_scenarios(
    cty,
    infection_list,
    list(custom = response, custom2 = response)
  )
  disease_tags <- sprintf("tag_%i", seq_along(infection_list))

  expect_no_condition(
    get_summary_data(output, disease_tags)
  )
  checkmate::expect_data_frame(
    get_summary_data(output, disease_tags)
  )

  expect_no_condition(
    get_summary_data(output, disease_tags, "wide")
  )

  expect_no_condition(
    get_summary_data(output, disease_tags, "wide", measures = "deaths"),
  )

  expect_no_condition(
    get_summary_data(
      output,
      disease_tags,
      "long",
      measures = c("deaths", "infections"),
      groups = "age_group"
    ),
  )

  group_wanted <- "age_group"
  measures_expected <- c("total_deaths", "epidemic_size") # different from args
  summary_data <- get_summary_data(
    output,
    disease_tags,
    "long",
    measures = c("deaths", "infections"),
    groups = group_wanted
  )
  checkmate::expect_subset(
    group_wanted,
    colnames(summary_data)
  )
  checkmate::expect_subset(
    measures_expected,
    summary_data$measure
  )

  expect_error(
    get_summary_data(
      output,
      disease_tags,
      "long",
      measures = c("deaths", "infections"),
      groups = "dummy_group"
    ),
    "Expected `groups` to be either `NULL` or a character vector"
  )

  expect_error(
    get_summary_data(
      output,
      disease_tags,
      "long",
      measures = "dummy_measure",
      groups = "age_group"
    ),
    "`measures` must be one of"
  )
})

test_that("Fn `get_cost_data()` works", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  cty <- "GBR"
  # single custom response
  response <- daedalus::daedalus_timed_npi(
    10,
    20,
    list(rep(0.5, daedalus:::N_ECON_SECTORS)),
    cty
  )

  output <- run_scenarios(
    cty,
    infection_list,
    list(custom = response, custom2 = response)
  )
  disease_tags <- sprintf("tag_%i", seq_along(infection_list))

  # TODO: NOT WORKING
  expect_no_condition(
    get_cost_data(output, disease_tags)
  )
  checkmate::expect_data_frame(
    get_cost_data(output, disease_tags)
  )

  checkmate::expect_data_frame(
    get_cost_data(output, disease_tags, "wide")
  )
})

test_that("Fn `get_econ_cost_data()` works", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  cty <- "GBR"
  # single custom response
  response <- daedalus::daedalus_timed_npi(
    10,
    20,
    list(rep(0.5, daedalus:::N_ECON_SECTORS)),
    cty
  )

  output <- run_scenarios(
    cty,
    infection_list,
    list(custom = response, custom2 = response)
  )
  disease_tags <- sprintf("tag_%i", seq_along(infection_list))

  # TODO: NOT WORKING
  expect_no_condition(
    get_econ_cost_data(output)
  )
  checkmate::expect_data_frame(
    get_econ_cost_data(output)
  )
})

test_that("Scenarios runner: errors and messages", {
  # fails with single list infection
  infection_list <- list(daedalus::daedalus_infection("sars_cov_1"))

  expect_error(
    run_scenarios(
      "GBR",
      infection_list
    ),
    "must be a list of >= 2"
  )

  # errors on country
  expect_error(
    run_scenarios(
      "XYZ",
      infection_list
    ),
    "`code` must be one of"
  )
  expect_error(
    run_scenarios(
      c("GBR", "CAN"),
      infection_list
    ),
    "`code` must be one of"
  )
})
