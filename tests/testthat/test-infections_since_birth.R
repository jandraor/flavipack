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
    c("subject_id", "age", "serotype")
  )
})

test_that("simulate_DENV_infections_since_birth() returns the expected value",
{
  final_age     <- 20
  n_indiv       <- 1e5
  n_runs        <- 25
  result_matrix <- matrix(NA, nrow = final_age, ncol = n_runs)

  for(i in 1:n_runs)
  {
    sim_df     <- simulate_DENV_infections_since_birth(0.1, 0.005,   n_indiv ,
                                                       final_age)
    pct_df     <- as.data.frame(table(sim_df$age))
    pct_df$pct <- pct_df$Freq / n_indiv

    result_matrix[, i ] <- pct_df$pct
  }

  actual <- rowMeans(result_matrix)

  expected <- c(0.32971560, 0.30647732, 0.28466016, 0.26359112, 0.24434680,
                0.22586132, 0.20883036, 0.19299904, 0.17829556,
                0.16498116, 0.15247208, 0.14114520, 0.13092424, 0.12174656,
                0.11316416, 0.10574648, 0.09891168, 0.09284392,
                0.08750268, 0.08261564)

  expect_equal(actual, expected, tolerance = 1e-3)
})

test_that("simulate_DENV_infections_since_birth() handles heterogeneity",
{
  set.seed(123)

  sim_df <- simulate_DENV_infections_since_birth(
    lambda_serotype = 1,
    loss_rate       = 0,
    n_individuals   = 5,
    final_age       = 5,
    r_i             = c(1, 1, 1, 1, 0))

  expect_equal(!(5 %in% sim_df$subject_id), TRUE)
})

#-------------------------------------------------------------------------------
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
