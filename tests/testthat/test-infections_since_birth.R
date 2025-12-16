test_that("simulate_DENV_infections_since_birth() works",
{
  set.seed(123)

  sim_df <- simulate_DENV_infections_since_birth(lambda_serotype = 0.1,
                                                   loss_rate     = 0.005,
                                                   n_individuals = 5,
                                                   final_age     = 5)

  expect_true(is.data.frame(sim_df))

  expect_equal(
    colnames(sim_df),
    c("infected_ind", "age", "serotype")
  )
})

test_that("simulate_DENV_infections_cohort() handles character ids",
{
  cohort_df <- data.frame(subject_id = c(rep("A", 3), rep("B", 2)),
                          age_sample = c(5:7, 2:3))

  set.seed(1226)

  sim_df <- simulate_DENV_infections_cohort(lambda_serotype = 0.1,
                                            loss_rate       = 0.005,
                                            cohort_df       = cohort_df)

  expect_true(all(sim_df$subject_id %in% c("A", "B")))
})
