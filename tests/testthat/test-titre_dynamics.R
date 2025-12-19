test_that("simulate_antibody_titres() works", {

  actual <- simulate_antibody_titres(n_years   = 2,
                                     A0        = exp(8),
                                     par_rho   = 0.98,
                                     par_sigma = 0.04,
                                     gamma_t   = c(1/720, 1/720))

  expected <- bi_exponential_decay(c(0, 360, 720),
                                   A0 =  exp(8),
                                   par_sigma = 0.04,
                                   par_gamma = 1/720,
                                   par_rho   = 0.98)


  expect_equal(actual, expected)

  actual <- simulate_antibody_titres(n_years   = 2,
                                     A0        = exp(8),
                                     par_rho   = 0.98,
                                     par_sigma = 0.04,
                                     gamma_t   = c(1/720, 1/1000))

  first_two_years <- bi_exponential_decay(
    t         = c(0, 360),
    A0        =  exp(8),
    par_sigma = 0.04,
    par_gamma = 1/720,
    par_rho   = 0.98)

  new_rho <- exp(8) * 0.98 * exp(-360 * 0.04) / first_two_years[[2]]


  last_year <- bi_exponential_decay(
    t         = 360,
    A0        =  first_two_years[[2]],
    par_sigma = 0.04,
    par_gamma = 1/1000,
    par_rho   = new_rho)

  expected <- c(first_two_years, last_year)

  expect_equal(actual, expected)
})

test_that("titre_decay_floor() works", {

  actual <- titre_decay_floor(par_alpha = 1, par_beta = 2, par_delta = 0.02,
                              time = 100)

  expected <- 1 + 2 * exp(-0.02 * 100)

  expect_equal(actual, expected)
})

test_that("simulate_DENV_long_decay_titres() works",
{
  inf_times      <- c(5, 7)

  actual <- simulate_DENV_long_decay_titres(
    inf_times      = list(inf_times),
    decay_rate_vec = c(0.1, 0.05, 0.025, 0.0115),
    log_first_peak = 1,
    phi            = 3,
    beta           = 1,
    final_age      = 10)

  titre_vals <- c(rep(5, 4),
                  10 * exp(-0.1 * (0:1)),
                  24.020082691 * exp(-0.05 *(0:3)))

  expected <- matrix(titre_vals, ncol = 10)

  expect_equal(actual, expected)
})
