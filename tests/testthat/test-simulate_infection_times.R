test_that("simulate_infection_times work()",
{
  set.seed(1305)

  actual <- simulate_infection_times(daily_lambda = 0.3 / 365,
                                     stop_time    = 365 * 5,
                                     min_gap      = 0,
                                     max_inf      = 4)

  expected <- c(340.2998, 445.2786, 539.4458)

  expect_equal(actual, expected, tolerance = 1e-4)

  set.seed(1305)

  actual <- simulate_infection_times(daily_lambda = 0.3 / 365,
                                     stop_time    = 365 * 3,
                                     min_gap      = 365,
                                     max_inf      = 4)

  expected <- c(340.2998, 445.2786 + 365)

  expect_equal(actual, expected, tolerance = 1e-4)
})

