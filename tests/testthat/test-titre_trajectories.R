test_that("simulate_titre_trajectory() works", {

  subject_id      <- 1
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
  treatment_group <- 0 # Placebo

  actual <- simulate_titre_trajectory(sampling_times,
                                      treatment_group,
                                      subject_id      = subject_id,
                                      infection_times = 1047.108,
                                      perm_rise       = 6,
                                      temp_rise       = 2,
                                      temp_decay      = 0.003,
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
                                      perm_rise       = 6,
                                      temp_rise       = 2,
                                      temp_decay      = 0.003,
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
    perm_rise       = 6,
    temp_rise       = 2,
    temp_decay      = 0.003,
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

#simulate_true_titres-----------------------------------------------------------


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
                                        perm_rise       = 6,
                                        temp_rise       = 2,
                                        temp_decay      = 0.003,
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

test_that("simulate_titres_seropositive() works",
{
  perm_rise <- 7.449828
  temp_rise <- 2.286833

  max_val <- perm_rise + 4 * temp_rise

  actual <- simulate_titres_seropositive(
    inf_times_sbs   = c(558.673, 1596.572, 2489.329, 3590.192),
    sampling_times  = c(40, 75, 143, 173, 299, 480, 851),
    age             = 15,
    subject_id      = 23,
    treatment_group = 0, # Placebo
    perm_rise       = perm_rise,
    temp_rise       = temp_rise,
    temp_decay      = 0.0001959646,
    meas_sd         = 0.1)

  expect_equal(all(actual$true <= max_val), TRUE)
})

test_that("simulate_titres_from_enrolment() works",
{
  infection_times <- NULL

  sampling_times <- c(537, 571, 627, 655, 808, 999, 1340, 1749, 2138)

  perm_rise  <- 6
  temp_rise  <- 2
  temp_decay <- 0.003

  actual <- simulate_titres_from_enrolment(infection_times, sampling_times,
                                             perm_rise, temp_rise, temp_decay,
                                             baseline = 6)

  expected <- rep(6, 9)

  expect_equal(actual, expected)

  infection_times <- 700

  actual <- simulate_titres_from_enrolment(infection_times, sampling_times,
                                           perm_rise, temp_rise, temp_decay,
                                           baseline = 6)

  after_1st_inf <- titre_decay_floor(perm_rise, temp_rise, temp_decay,
                                     sampling_times[5:9] - infection_times)


  expected <- c(rep(6, 4), after_1st_inf)

  expect_equal(actual, expected)

  #=============================================================================

  infection_times <- c(700, 1400)

  actual <- simulate_titres_from_enrolment(infection_times, sampling_times,
                                           perm_rise, temp_rise, temp_decay,
                                           baseline = 7)

  before_first_inf <- titre_decay_floor(
    par_alpha = perm_rise,
    par_beta  = 7 - perm_rise,
    par_delta = temp_decay,
    time      = c(sampling_times[1:4], infection_times[[1]]) -
      sampling_times[[1]])

  after_1st_inf <- titre_decay_floor(
    par_alpha = perm_rise,
    par_beta  = before_first_inf[[5]] - perm_rise + temp_rise,
    par_delta = temp_decay,
    time      = c(sampling_times[5:7], infection_times[[2]]) -
      infection_times[[1]])

  after_2nd_inf <- titre_decay_floor(
    par_alpha = perm_rise,
    par_beta  = after_1st_inf[[4]] - perm_rise + temp_rise,
    par_delta = temp_decay,
    time      = sampling_times[8:9] - infection_times[[2]])

  expected <- c(before_first_inf[1:4],
                after_1st_inf[1:3],
                after_2nd_inf)

  expect_equal(actual, expected)
})
