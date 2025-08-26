test_that("simulate_titre_trajectory() works", {

  subject_id      <- 1
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
  treatment_group <- 0 # Placebo

  actual <- simulate_titre_trajectory(sampling_times,
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

test_that("simulate_titre_trajectory() works for reconstructing titres since
          birth",
{
  sampling_times <- c(3285, 3318, 3366, 3401, 3534, 3712, 4031, 4436, 4797)

  inf_times <- c(1059.646, 2743.543, 4152.091)

  actual <- simulate_titre_trajectory(
    sampling_times  = sampling_times,
    treatment_group = 0, # Placebo
    subject_id      = 3,
    infection_times = inf_times,
    meas_sd         = 0.3)

  first_infection_last_value <- max(8 - 0.003 * (inf_times[[2]] - inf_times[[1]]), 6)

  second_inf <- pmax(8 - 0.003 * (c(sampling_times[1:7], inf_times[[3]]) - inf_times[[2]]), 6)

  third_inf <- titre_decay_floor(6, 2, 0.003,
                                 sampling_times[8:9] - inf_times[[3]])

  expected <- data.frame(subject_id = 3,
                         time       = sampling_times,
                         true       = c(second_inf[1:7], third_inf))

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)
})



test_that("simulate_titre_trajectory() handles no infections",
{
  sampling_times <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)

  actual <- simulate_titre_trajectory(
    sampling_times  = sampling_times ,
    treatment_group = 0,
    subject_id      = 1,
    infection_times = numeric(0),
    meas_sd         = 0.3)

  expected <- data.frame(subject_id = 1,
                         time       = sampling_times,
                         true       = 0,
                         meas       = 0)

  expect_equal(actual, expected)
})

#simulate_titres_seropositive---------------------------------------------------

test_that("simulate_titres_seropositive() works",
{
  age <- 10
  # Since becoming susceptible (sbs)
  inf_times_sbs <- c(1059.646, 2743.543, 4152.091)

  inf_times <- c(1059.646, 2743.543, 4152.091)

  # relative to study start
  sampling_times <- c(41, 74, 122, 157, 290, 468, 787, 1192, 1553)

  sampling_times_rel_to_last_birthday <- sampling_times - min(sampling_times)

  sampling_times_sbs <- sampling_times_rel_to_last_birthday  + 365 * (age - 1)

  actual <- simulate_titres_seropositive(inf_times_sbs   = inf_times_sbs,
                                        sampling_times  = sampling_times,
                                        age             = age,
                                        subject_id      = 3,
                                        treatment_group = 0,
                                        meas_sd         = 0.3)

  first_infection_last_value <- max(8 - 0.003 * (inf_times_sbs[[2]] - inf_times_sbs[[1]]), 6)

  second_inf <- titre_decay_floor(6, 2, 0.003,
                                  c(sampling_times_sbs[1:7],
                                    inf_times_sbs[[3]]) - inf_times_sbs[[2]])

  third_inf <- titre_decay_floor(6, 2, 0.003,
                                 sampling_times_sbs[8:9] - inf_times_sbs[[3]])

  expected <-  data.frame(subject_id = 3,
                          time       = sampling_times,
                          true       = c(second_inf[1:7], third_inf))

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)
})
