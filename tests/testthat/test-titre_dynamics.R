test_that("add_titre_to_infection_history() works",
          {
            inf_history_df <-
              data.frame(subject_id = 1,
                         year_index = 11:23,
                         age        = 2:14,
                         sim_y      = 0)

            inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1

            actual <- add_titre_to_infection_history(
              inf_history_df,
              long_options  = list(L = 0.44, a = 0, b = 0),
              short_options = list(intercept = 0.02,
                                   slope     = 0.003342),
              rise_options  = list(intercept = 8,
                                   slope     = -0.8))

            expected <- inf_history_df

            first_inf <- bi_exponential_decay(
              c(0, 360, 720), exp(8),
              par_sigma = 0.02 + 0.003342*8,
              par_gamma = 0.22 / 360,
              par_rho = 0.98) |> log()

            rise <- 8 + first_inf[3] * -0.8

            second_inf <- bi_exponential_decay(
              360 * 0:3,
              exp(first_inf[3] + rise),
              0.02 + 0.003342 * rise,
              par_gamma = 0.22 / 360,
              par_rho = 0.98) |> log()

            rise <- 8 + second_inf[4] * -0.8

            third_inf <- bi_exponential_decay(
              360 * 0:1,
              A0 = exp(second_inf[4] + rise),
              0.02 + 0.003342 * rise,
              par_gamma = 0.22 / 360,
              par_rho = 0.98) |>
              log()

            rise <- 8 + third_inf[2] * -0.8

            fourth_inf <- bi_exponential_decay(
              360 * 0:1,
              A0 = exp(third_inf[2] + rise),
              0.02 + 0.003342 * rise,
              par_gamma = 0.22/ 360,
              par_rho = 0.98) |>
              log()

            expected$cumulative_infections <- cumsum(expected$sim_y)

            expected$log_titre <- c(rep(0, 5), first_inf[1:2],
                                    second_inf[1:3],
                                    third_inf[1],
                                    fourth_inf)

            expect_equal(actual, expected)
          })

test_that("add_titre_to_infection_history() validates input",
          {
            expect_error(add_titre_to_infection_history(
              data.frame(subject_id = 1, year_index = 10, current_age = 10, sim_y = 0),
              short_options = list(),
              long_options = list(),
              rise_options = list()))
          })

test_that("estimate_long_rate() works",
          {
            long_options <- list(L = 0.22, a = -0.96, b = 18.6)

            age <- 7:9

            actual <- estimate_long_rate(age = 7:9, long_options)

            expected <- 0.22 / (1 + exp(0.96 * (7:9 - 18.6)))

            expect_equal(actual, expected)
          })

test_that("estimate_rise() works",
          {
            rise_options <- list(intercept = 8,
                                 slope     = -0.8)

            actual <- estimate_rise(baseline_titre = 0, rise_options)

            expected <- 8

            expect_equal(actual, expected)

            actual <- estimate_rise(baseline_titre = 1000, rise_options)

            expected <- 0

            expect_equal(actual, expected)

            actual <- estimate_rise(baseline_titre = -1000, rise_options)

            expected <- 8

            expect_equal(actual, expected)
          })

test_that("estimate_short_rate() works", {

  short_options <- list(intercept = 0.02,
                        slope     = 0.003342)

  actual <- estimate_short_rate(rise = 0, short_options)


  expected <- 0.02

  expect_equal(actual, expected)
})

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
              inf_times      = inf_times,
              decay_rate_vec = c(0.1, 0.05, 0.025, 0.0115),
              log_first_peak = 1,
              phi            = 2,
              beta           = 1,
              subject_id     = 1,
              final_age      = 10)

            expected <- data.frame(subject_id = 1,
                                   age        = 1:10,
                                   titre = c(rep(5, 4),
                                             10 * exp(-0.1 * (0:1)),
                                             24.020082691 * exp(-0.05 *(0:3))))

            expect_equal(actual, expected)
          })
