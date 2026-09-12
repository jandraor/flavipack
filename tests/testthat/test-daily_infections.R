#--------simulate_infection_times_since_susceptibility--------------------------

test_that("wrapper matches simulate_infection_times with derived stop time", {

  age_enrolment  <- 5
  enrolment_time <- 20
  follow_up_time <- 365 * 5
  lambda         <- 0.1
  min_gap        <- 365
  max_inf        <- 4

  set.seed(999)
  x1 <- simulate_infection_times_since_susceptibility(
    lambda         = lambda,
    age_enrolment  = age_enrolment,
    enrolment_time = enrolment_time,
    follow_up_time = follow_up_time,
    min_gap        = min_gap,
    max_inf        = max_inf)

  stop_time <- calculate_stop_time(
    age_enrolment  = age_enrolment,
    enrolment_time = enrolment_time,
    follow_up_time = follow_up_time)

  set.seed(999)
  x2 <- rInfections_DENV(
    FOI           = lambda,
    end_time      = stop_time,
    min_gap       = min_gap,
    empty_buckets = max_inf)

  expect_equal(x1, x2)
})

test_that("calculate_stop_time() works",
          {
            actual <- calculate_stop_time(10,
                                          follow_up_time = 2000,
                                          enrolment_time = 100)

            expect_equal(actual, 5185)
          })

