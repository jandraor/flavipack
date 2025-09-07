test_that("calculate_log_lik() works",
{
  measurements <- c(6.325498, 6.345967, 6.298308, 6.511790, 6.427388, 7.273151,
                    6.309557, 7.001809, 6.448628)

  sampling_times <- sampling_times <- c(537, 571, 627, 655, 808, 999,
                                        1340, 1749, 2138)

  infection_times <- c(1656.6)
  actual <- calculate_log_lik(infection_times,
                              sampling_times ,
                              perm_rise  = 6,
                              temp_rise  = 2,
                              temp_decay = 0.5,
                              measurements,
                              meas_sd = 0.1)

  expected <- -166.292066

  expect_equal(actual, expected)
})
