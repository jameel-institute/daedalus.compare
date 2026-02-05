test_that("Basic multi-vaccination scenario runs", {
  country <- "GB"
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 2
  )

  # single vaccination
  late_vaccination <- daedalus::daedalus_vaccination("low", country)
  expect_no_condition(
    run_scenarios(
      country,
      infection_list,
      vaccination_strategy = late_vaccination,
      time_end = 100
    )
  )

  expected_vax_name <- "custom_vaccination"
  output <- run_scenarios(
    country,
    infection_list,
    vaccination_strategy = late_vaccination,
    time_end = 100
  )
  expect_identical(
    output$vaccination,
    expected_vax_name
  )

  # multiple vaccinations
  no_vaccination <- NULL
  early_vaccination <- daedalus::daedalus_vaccination("high", country)

  expect_no_condition(
    run_scenarios(
      country,
      infection_list,
      vaccination_strategy = list(
        none = no_vaccination,
        late_vaccination = late_vaccination,
        early_vaccination = early_vaccination
      ),
      time_end = 100
    )
  )

  output <- run_scenarios(
    country,
    infection_list,
    vaccination_strategy = list(
      none = no_vaccination,
      late_vaccination = late_vaccination,
      early_vaccination = early_vaccination
    ),
    time_end = 100
  )

  expect_s3_class(output, "data.table")

  expected_names <- c("response", "vaccination", "time_end", "output")
  checkmate::expect_names(
    names(output),
    must.include = expected_names
  )

  # unnamed list of vaccinations
  expect_message(
    run_scenarios(
      country,
      infection_list,
      vaccination_strategy = list(
        no_vaccination,
        late_vaccination,
        early_vaccination
      ),
      time_end = 100
    ),
    "Assigning manually constructed names"
  )

  vax_strategy_list <- list(
    no_vaccination,
    late_vaccination,
    early_vaccination
  )
  output <- run_scenarios(
    country,
    infection_list,
    vaccination_strategy = vax_strategy_list,
    time_end = 100
  )
  expected_names <- sprintf("vaccination_%i", seq_along(vax_strategy_list))
  expect_identical(
    output$vaccination,
    expected_names
  )
})

test_that("Multi-response, multi-vaccination scenario runs", {
  country <- "GB"
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 2
  )

  no_vaccination <- NULL
  late_vaccination <- daedalus::daedalus_vaccination("low", country)
  early_vaccination <- daedalus::daedalus_vaccination("high", country)

  # single custom response
  response <- daedalus::daedalus_timed_npi(
    10,
    20,
    list(rep(0.5, daedalus:::N_ECON_SECTORS)),
    country
  )

  expect_no_condition(
    run_scenarios(
      country,
      infection_list,
      response_strategy = list(
        none = NULL,
        custom_npi = response
      ),
      vaccination_strategy = list(
        none = no_vaccination,
        late_vaccination = late_vaccination,
        early_vaccination = early_vaccination
      ),
      time_end = 100
    )
  )

  output <- run_scenarios(
    country,
    infection_list,
    response_strategy = list(
      none = NULL,
      custom_npi = response
    ),
    vaccination_strategy = list(
      none = no_vaccination,
      late_vaccination = late_vaccination,
      early_vaccination = early_vaccination
    ),
    time_end = 100
  )

  # expect 6 unique entries
  expected_unique <- 6L
  expect_identical(
    data.table::uniqueN(
      output,
      by = c("response", "vaccination")
    ),
    expected_unique
  )
})

test_that("Multi-vaccination comparisons: errors and messages", {
  country <- "GB"
  infection_list <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = 2
  )

  no_vaccination <- NULL
  late_vaccination <- daedalus::daedalus_vaccination("low", country)
  early_vaccination <- daedalus::daedalus_vaccination("high", country)

  expect_error(
    run_scenarios(
      country,
      infection_list,
      vaccination_strategy = list(
        none = no_vaccination
      ),
      time_end = 100
    ),
    "other classes were found, or the list has only one element."
  )

  expect_error(
    run_scenarios(
      country,
      infection_list,
      vaccination_strategy = "dummy_vax",
      time_end = 100
    ),
    "Got `vaccination_strategy` of class"
  )
})
