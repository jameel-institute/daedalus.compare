# test basic scenario runners
test_that("Generating multiple infection samples", {
  expect_no_condition(
    make_infection_samples("influenza_1918")
  )
  checkmate::expect_list(
    make_infection_samples("influenza_1918"),
    "daedalus_infection"
  )
})

test_that("Running multiple infection and response scenarios", {
  infection_list <- make_infection_samples("influenza_2009", samples = 10)

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
  expect_data_table(
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
  expect_data_table(
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
  expect_data_table(
    output
  )
  expect_list(
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
  expect_data_table(
    output
  )
  expect_list(
    output$output, "daedalus_output"
  )

  # custom responses and multiple infections
  infection_list <- make_infection_samples("influenza_2009", samples = 10)
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
  infection_list <- make_infection_samples("influenza_2009", samples = 10)
  response <- seq(0.5, daedalus:::N_ECON_SECTORS)

  output <- run_scenarios(
    "GBR", infection_list,
    list(custom = response, custom2 = response)
  )

  # TODO: NOT WORKING
  expect_no_condition(
    get_epicurve_data(output)
  )
})
