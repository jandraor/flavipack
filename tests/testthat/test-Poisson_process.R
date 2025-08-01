test_that("simulate_infection_times() works", {

  sim_times <- simulate_infection_times(lambda    = 0.3 / 365, # Daily FOI
                                        stop_time = 365 * 3,
                                        min_gap   = 365,
                                        max_inf   = 4)

  expect_equal(length(sim_times) <= 4, TRUE)
})

test_that("simulate_infection_times_since_birth() works", {

  actual <- simulate_infection_times_since_birth(
    lambda        = 0.3 / 365,
    age_enrolment = 8, # in years
    offset_start  = 0,
    stop_time     = 365 * 5,
    min_gap       = 365,
    max_inf       = 4)

  expect_true(is.numeric(actual))
})
