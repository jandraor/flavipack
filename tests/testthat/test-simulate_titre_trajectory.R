test_that("simulate_titre_trajectory() works", {

  subject_id      <- 1
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
  serostatus      <- 0 # Seronegative
  treatment_group <- 0 # Placebo

  set.seed(1652)

  actual <- simulate_titre_trajectory(sampling_times,
                                      serostatus,
                                      treatment_group,
                                      subject_id = subject_id)

  expected <- data.frame(subject_id = subject_id,
                         time       = sampling_times,
                         true_titre = c(rep(0, length(sampling_times) - 2),
                                        7.10032337, 6.000000))

  expect_equal(nrow(actual), length(sampling_times))

  expect_s3_class(actual, "data.frame")

  expect_equal(actual, expected)
})
