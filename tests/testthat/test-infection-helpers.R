test_that("Making infection samples works", {
  n_samples <- 10L
  expect_no_condition(
    make_infection_samples(
      "influenza_2009",
      list(r0 = distributional::dist_beta(2, 5)),
      samples = n_samples
    )
  )
  samples <- make_infection_samples(
    "influenza_2009",
    list(r0 = distributional::dist_beta(2, 5)),
    samples = n_samples
  )
  checkmate::expect_list(
    samples,
    "daedalus_infection",
    any.missing = FALSE,
    len = n_samples
  )

  # passing parameter distributions
  expect_no_condition(
    make_infection_samples(
      "influenza_2009",
      samples = n_samples,
      list(
        sigma = distributional::dist_beta(2, 5)
      )
    )
  )

  expect_no_condition(
    make_infection_samples(
      "influenza_2009",
      samples = n_samples,
      list(
        sigma = distributional::dist_beta(2, 5)
      ),
      list(
        sigma = c(0.1, 0.2)
      )
    )
  )

  # passing age-varying parameters
  expect_no_condition(
    make_infection_samples(
      "influenza_2009",
      samples = n_samples,
      list(
        eta = distributional::dist_beta(2, 5)
      ),
      list(
        eta = c(0.1, 0.2)
      )
    )
  )
})

test_that("Infection samples errors", {
  expect_error(
    make_infection_samples("dummy")
  )
  expect_error(
    make_infection_samples("dummy", -10)
  )
  expect_error(
    make_infection_samples("dummy", Inf)
  )

  # negatives generated
  dist_normal <- distributional::dist_normal()
  expect_error(
    make_infection_samples(
      "influenza_2009",
      list(r0 = dist_normal)
    ),
    regexp = "Negative values generated"
  )
})
