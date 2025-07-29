test_that("simulate_titre_trajectory() works", {

  subject_id      <- 1
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
  serostatus      <- 0 # Seronegative
  treatment_group <- 0 # Placebo

  actual <- simulate_titre_trajectory(sampling_times,
                                      serostatus,
                                      treatment_group,
                                      subject_id      = subject_id,
                                      infection_times = 1047.108,
                                      meas_sd         = 0.3)

  expected <- data.frame(subject_id = subject_id,
                         time       = sampling_times,
                         true       = c(rep(0, length(sampling_times) - 2),
                                        7.10032400, 6.000000))

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)

  # Two infections

  t_inf1 <- 300

  t_inf2 <- 900

  actual <- simulate_titre_trajectory(sampling_times,
                                      serostatus,
                                      treatment_group,
                                      subject_id      = subject_id,
                                      infection_times = c(t_inf1, t_inf2),
                                      meas_sd         = 0.3)

  first_inf <- 8 - 0.003 * (c(sampling_times[5:6], t_inf2) - t_inf1)

  second_inf <- pmax(first_inf[3] + 2  - 0.003 * (sampling_times[7:9] - t_inf2),
                    6)


  expected <- data.frame(subject_id = subject_id,
                         time       = sampling_times,
                         true       = c(rep(0, 4),
                                        first_inf[1:2],
                                        second_inf))

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)
})
