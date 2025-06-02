# test basic scenario runners
test_that("Running multiple infection and response scenarios", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      response_strategy = "none"
    )
  )
  # expect output is a data table with list columns
  output <- run_scenarios(
    "GBR",
    infection_list,
    response_strategy = "none"
  )
  checkmate::expect_data_table(
    output
  )
  checkmate::expect_list(output$output, "list")

  # expect running scenarios works for multiple pre-defined responses
  responses <- c("none", "elimination")
  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      response_strategy = responses
    )
  )
  output <- run_scenarios(
    "GBR",
    infection_list,
    response_strategy = responses
  )
  checkmate::expect_data_table(
    output,
    nrows = length(responses)
  )
  checkmate::expect_list(output$output, "list")
})

test_that("Running custom response scenarios", {
  # single custom response
  response <- rep(0.5, daedalus:::N_ECON_SECTORS)
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  # NOTE: daedalus_multi_infection only works on lists w/ len > 1
  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      list(custom = response)
    )
  )
  output <- run_scenarios(
    "GBR",
    infection_list,
    list(custom = response)
  )
  checkmate::expect_data_table(
    output
  )

  # NOTE: check that each output is a list of daedalus_output
  # this is because we no longer allow single infection use of
  # daedalus_multi_infection
  checkmate::expect_list(
    output$output[[1]],
    "daedalus_output"
  )

  # multiple custom responses
  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      list(custom = response, custom2 = response)
    )
  )
  output <- run_scenarios(
    "GBR",
    infection_list,
    list(response, response)
  )
  checkmate::expect_data_table(
    output
  )
  checkmate::expect_list(
    output$output[[1]],
    "daedalus_output"
  )

  # custom responses and multiple infections
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      list(custom = response, custom2 = response)
    )
  )

  # custom and pre-defined responses
  # TODO: examine how names and openness coefficients are returned
  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      list(response, "elimination")
    )
  )
})

test_that("Get epi curve data", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  disease_tags <- sprintf("disease_%i", seq_along(infection_list))
  response <- rep(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR",
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
})

test_that("Get epi summary data", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  response <- rep(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR",
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
})

test_that("Get costs data", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  response <- rep(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR",
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

test_that("Get econ costs data", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  response <- rep(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR",
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

  infection_list <- list(
    daedalus::daedalus_infection("sars_cov_1"),
    daedalus::daedalus_infection("sars_cov_1")
  )
  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      "dummy_response"
    ),
    "Got an unexpected value for `response_strategy`"
  )
  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      c("elimination", "dummy")
    ),
    "Got an unexpected value for `response_strategy`"
  )

  # response time can be 0 if the strategy is 'none' or NULL
  expect_no_condition(
    run_scenarios(
      "GBR",
      infection_list,
      response_time = 0
    )
  )
  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      "elimination",
      response_time = 0
    ),
    "Expected `response_time` to be between 1 and 100"
  )

  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      response_strategy = "elimination",
      response_duration = -1
    ),
    "(Expected)*(single positive integer-like)"
  )
  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      response_strategy = "elimination",
      response_duration = c(100, 365)
    ),
    "(Expected)*(single positive integer-like)"
  )

  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      time_end = -1
    ),
    "Expected `time_end` to be a single positive integer-like"
  )
})
