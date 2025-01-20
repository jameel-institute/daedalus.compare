test_that("Miscellaneous helper functions", {
  # test `ci()`
  expect_no_condition(
    ci(1:10)
  )
  expect_no_condition(
    ci(1:10, c(90, 91))
  )
  expect_error(
    ci(1:10, -1)
  )
  expect_error(
    ci(1:10, 101)
  )
  expect_error(
    ci(1:10, Inf)
  )
  expect_error(
    ci(1:10, NULL)
  )
})
