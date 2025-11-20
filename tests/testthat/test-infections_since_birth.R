test_that("simulate_DENV_infections_since_birth() works",
{
  set.seed(123)

  actual <- simulate_DENV_infections_since_birth(lambda_serotype = 0.1,
                                                   loss_rate       = 0.005,
                                                   n_individuals   = 5,
                                                   final_age       = 5)

  expected <- data.frame(infected_ind = c(1, 5, 3, 5, 1),
                         age          = c(1, 1, 2, 2, 5))

  expect_equal(actual, expected)
})

test_that("simulate_DENV_infections_cohort() works",
{
  cohort_df <- data.frame(subject_id = c(rep(1, 3), rep(2, 2)),
                          age_sample = c(5:7, 2:3))

  set.seed(1226)

  actual <- simulate_DENV_infections_cohort(lambda_serotype = 0.1,
                                       loss_rate       = 0.005,
                                       cohort_df       = cohort_df)

  expected <- cohort_df

  expected$infection <- c(FALSE, FALSE, FALSE, TRUE, FALSE)

  expect_equal(actual, expected)
})
