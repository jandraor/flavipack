test_that("simulate_infection_times() returns empty when lambda <= 0",
{
  expect_equal(simulate_infection_times(0, 2000, 365, 4), numeric(0))
  expect_equal(simulate_infection_times(-1, 2000, 365, 4), numeric(0))
})

test_that("simulate_infection_times() never exceeds max_inf", {
  x <- simulate_infection_times(0.8, 2000, 365, 4)
  expect_lte(length(x), 4)
})

test_that("simulate_infection_times() infection times are increasing", {
  x <- simulate_infection_times(0.1, 2000, 365, 4)
  expect_false(is.unsorted(x))
})

test_that("simulate_infection_times() infection times are before stop_time", {
  x <- simulate_infection_times(0.1, 2000, 365, 4)
  expect_true(all(x < 2000))
})

test_that("simulate_infection_times() minimum gap is respected", {
  x <- simulate_infection_times(0.1, 2000, 365, 4)
  expect_true(length(x) < 2 || all(diff(x) >= 365))
})

test_that("mean times to first and second infections match theory", {
  set.seed(123)

  lambda    <- 0.1
  min_gap   <- 365
  max_inf   <- 4
  n         <- 20000

  sims <- replicate(n, {
    simulate_infection_times(
      lambda    = lambda,
      stop_time = 1e6,
      min_gap   = min_gap,
      max_inf   = max_inf
    )
  }, simplify = FALSE)

  first_times  <- vapply(sims, `[`, numeric(1), 1)
  second_times <- vapply(sims, `[`, numeric(1), 2)

  daily_lambda <- lambda / 365

  theory_first  <- 1 / (4 * daily_lambda)
  theory_second <- 1 / (4 * daily_lambda) + min_gap + 1 / (3 * daily_lambda)

  mcse_first  <- sd(first_times) / sqrt(n)
  mcse_second <- sd(second_times) / sqrt(n)

  expect_lt(abs(mean(first_times)  - theory_first),  3 * mcse_first)
  expect_lt(abs(mean(second_times) - theory_second), 3 * mcse_second)
})

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
  x2 <- simulate_infection_times(
    lambda    = lambda,
    stop_time = stop_time,
    min_gap   = min_gap,
    max_inf   = max_inf)

  expect_equal(x1, x2)
})

test_that("calculate_stop_time() works",
{
  actual <- calculate_stop_time(10,
                                follow_up_time = 2000,
                                enrolment_time = 100)

  expect_equal(actual, 5185)
})
