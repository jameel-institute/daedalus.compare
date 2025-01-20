# test basic scenario runners
test_that("Running multiple infection and response scenarios", {
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  expect_no_condition(
    run_scenarios(
      "GBR", infection_list,
      response = "none"
    )
  )
  # expect output is a data table with list columns
  output <- run_scenarios(
    "GBR", infection_list,
    response = "none"
  )
  checkmate::expect_data_table(
    output
  )
  checkmate::expect_list(output$output, "list")

  # expect running scenarios works for multiple pre-defined responses
  responses <- c("none", "elimination")
  expect_no_condition(
    run_scenarios(
      "GBR", infection_list,
      response = responses
    )
  )
  output <- run_scenarios(
    "GBR", infection_list,
    response = responses
  )
  checkmate::expect_data_table(
    output,
    nrows = length(responses)
  )
  checkmate::expect_list(output$output, "list")
})

test_that("Running custom response scenarios", {
  # single custom response
  response <- seq(0.5, daedalus:::N_ECON_SECTORS)
  expect_no_condition(
    run_scenarios(
      "GBR", "influenza_2009", list(custom = response)
    )
  )
  output <- run_scenarios(
    "GBR", "influenza_2009", list(custom = response)
  )
  checkmate::expect_data_table(
    output
  )
  checkmate::expect_list(
    output$output, "daedalus_output"
  )

  # multiple custom responses
  expect_no_condition(
    run_scenarios(
      "GBR", "influenza_2009",
      list(custom = response, custom2 = response)
    )
  )
  output <- run_scenarios(
    "GBR", "influenza_2009", list(response, response)
  )
  checkmate::expect_data_table(
    output
  )
  checkmate::expect_list(
    output$output, "daedalus_output"
  )

  # custom responses and multiple infections
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  expect_no_condition(
    run_scenarios(
      "GBR", infection_list,
      list(custom = response, custom2 = response)
    )
  )

  # custom and pre-defined responses
  # TODO: examine how names and openness coefficients are returned
  expect_no_condition(
    run_scenarios(
      "GBR", "influenza_2009", list(response, "elimination")
    )
  )
})

test_that("Get epi curve data", {
  infection_list <- make_infection_samples(
    "influenza_2009", list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  disease_tags <- sprintf("disease_%i", seq_along(infection_list))
  response <- seq(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR", infection_list,
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
    "influenza_2009", list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  response <- seq(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR", infection_list,
    list(custom = response, custom2 = response)
  )
  disease_tags <- sprintf("tag_%i", seq_along(infection_list))

  # TODO: NOT WORKING
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
    "influenza_2009", list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  response <- seq(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR", infection_list,
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
    "influenza_2009", list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )
  response <- seq(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR", infection_list,
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
