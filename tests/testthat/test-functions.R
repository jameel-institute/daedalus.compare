test_that("function works", {
  x <- c(2, 3)
  expect_equal(
    get_manual_mean(x),
    mean(x)
  )
})
