test_that("daedalus.compare: run multiple infections and responses", {
  cty <- "GBR"

  # single custom response
  response <- daedalus::daedalus_timed_npi(
    10,
    20,
    list(rep(0.5, daedalus:::N_ECON_SECTORS)),
    cty
  )

  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 10
  )

  # NOTE: daedalus_multi_infection only works on lists w/ len > 1
  expect_no_condition(
    run_scenarios(
      cty,
      infection_list,
      response
    )
  )
  output <- run_scenarios(
    cty,
    infection_list,
    response
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
      cty,
      infection_list,
      list(custom = response, custom2 = response)
    )
  )
  output <- run_scenarios(
    cty,
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
})

test_that("Multi-response comparisons: errors and messages", {
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
    "(Got `response_strategy` of class)*(character)"
  )

  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      list(
        NULL,
        daedalus::daedalus_timed_npi(
          10,
          20,
          list(rep(0.5, 45)),
          "GBR"
        ),
        "dummy"
      )
    ),
    "(other)*(classes were found)"
  )

  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      list(
        NULL,
        daedalus::daedalus_npi(
          "elimination",
          "GBR",
          "sars_cov_1"
        )
      )
    ),
    "some or all NPIs are reactive to model state"
  )

  expect_error(
    run_scenarios(
      "GBR",
      infection_list,
      daedalus::daedalus_npi(
        "elimination",
        "GBR",
        "sars_cov_1"
      )
    ),
    "it is reactive to model state"
  )
})
