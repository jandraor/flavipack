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

test_that("rInfections_DENV() returns empty when lambda <= 0",
{
  expect_equal(rInfections_DENV(0, 2000, 365, 4), numeric(0))
  expect_equal(rInfections_DENV(-1, 2000, 365, 4), numeric(0))
})

test_that("rInfections_DENV() never exceeds max_inf", {
  x <- rInfections_DENV(0.8, 2000, 365, 4)
  expect_lte(length(x), 4)
})

test_that("rInfections_DENV() infection times are increasing", {
  x <- rInfections_DENV(0.1, 2000, 365, 4)
  expect_false(is.unsorted(x))
})

test_that("rInfections_DENV() infection times are before stop_time", {
  x <- rInfections_DENV(0.1, 2000, 365, 4)
  expect_true(all(x < 2000))
})

test_that("rInfections_DENV() minimum gap is respected", {
  x <- rInfections_DENV(0.1, 2000, 365, 4)
  expect_true(length(x) < 2 || all(diff(x) >= 365))
})

test_that("mean times to first and second infections match theory", {
  set.seed(123)

  lambda    <- 0.1
  min_gap   <- 365
  max_inf   <- 4
  n         <- 20000

  sims <- replicate(n, {
    rInfections_DENV(
      FOI           = lambda,
      end_time      = 1e6,
      min_gap       = min_gap,
      empty_buckets = max_inf
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

