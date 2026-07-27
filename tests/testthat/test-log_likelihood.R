test_that("dInfections_DENV() works with no infections",
{
  actual <- dInfections_DENV(0.025, 0, 100, empty_buckets = 4)

  expected <- -10

  expect_equal(actual, expected)
})

test_that("dInfections_DENV() works with one infection",
{
  infection_history <- c(40)

  actual <- dInfections_DENV(0.025, 0, 100,
                             infection_history,
                             min_gap = 30,
                             empty_buckets = 4)

  expected <- -0.1 * 40.0 +
    log(0.1) -
    0.075 * 30.0

  expect_equal(actual, expected)
})
