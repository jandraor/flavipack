test_that("simulate_titre_trajectory() works", {

  subject_id      <- 1
  inf_times       <- 1047.108
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
  treatment_group <- 0 # Placebo

  actual <- simulate_titre_trajectory(sampling_times,
                                      treatment_group,
                                      subject_id      = subject_id,
                                      infection_times = inf_times,
                                      perm_rise   = 1,
                                      temp_rise   = 2,
                                      temp_decay  = 0.003,
                                      meas_sd     = 0.3)

  expected <- data.frame(subject_id = subject_id,
                         time       = sampling_times,
                         true       = c(rep(0, length(sampling_times) - 2),
                                        1.81340282, 1.25019146))

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)

  # Two infections

  inf_times <- c(300, 900)

  perm_rise_pri <- 1
  temp_rise_pri <- 2
  perm_rise_sec <- 2
  temp_rise_sec <- 3

  temp_decay <- 0.003

  actual <- simulate_titre_trajectory(sampling_times,
                                      treatment_group,
                                      subject_id      = subject_id,
                                      infection_times = inf_times,
                                      perm_rise       = c(perm_rise_pri,
                                                          perm_rise_sec),
                                      temp_rise       = c(temp_rise_pri,
                                                          temp_rise_sec),
                                      temp_decay      = c(temp_decay,
                                                          temp_decay),
                                      meas_sd         = 0.3)

  titre_matrix <- matrix(0,
                         nrow = length(inf_times),
                         ncol = length(sampling_times))

  idx_1st   <- which(sampling_times > inf_times[[1]])
  times_1st <- sampling_times[idx_1st]

  titre_matrix[1, idx_1st] <- titre_decay_floor(
    par_alpha = perm_rise_pri,
    par_beta  = temp_rise_pri,
    par_delta = temp_decay,
    time      = times_1st - inf_times[[1]])

  idx_2nd   <- which(sampling_times > inf_times[[2]])
  times_2nd <- sampling_times[idx_2nd]

  titre_matrix[2, idx_2nd] <- titre_decay_floor(
    par_alpha = perm_rise_sec,
    par_beta  = temp_rise_sec,
    par_delta = temp_decay,
    time      = times_2nd - inf_times[[2]])

  true_titre <- colSums(titre_matrix)

  expected <- data.frame(subject_id = subject_id,
                         time       = sampling_times,
                         true       = true_titre)

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

#simulate_true_titres_DENV------------------------------------------------------

test_that("simulate_true_titre() works",
{
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)

  inf_times <-c(300, 900)

  perm_rise_pri <- 1
  temp_rise_pri <- 2

  perm_rise_sec <- 2
  temp_rise_sec <- 3

  temp_decay    <- 0.003

  actual <- simulate_true_titre_DENV(
    sampling_times   = sampling_times,
    perm_rise       = c(perm_rise_pri,
                        perm_rise_sec),
    temp_rise       = c(temp_rise_pri,
                        temp_rise_sec),
    temp_decay       = c(temp_decay,
                         temp_decay),
    treatment_group  = 0, # Placebo
    infection_times  = inf_times)

  titre_matrix <- matrix(0,
                         nrow = length(inf_times),
                         ncol = length(sampling_times))

  idx_1st   <- which(sampling_times > inf_times[[1]])
  times_1st <- sampling_times[idx_1st]

  titre_matrix[1, idx_1st] <- titre_decay_floor(
    par_alpha = perm_rise_pri,
    par_beta  = temp_rise_pri,
    par_delta = temp_decay,
    time      = times_1st - inf_times[[1]])

  idx_2nd   <- which(sampling_times > inf_times[[2]])
  times_2nd <- sampling_times[idx_2nd]

  titre_matrix[2, idx_2nd] <- titre_decay_floor(
    par_alpha = perm_rise_sec,
    par_beta  = temp_rise_sec,
    par_delta = temp_decay,
    time      = times_2nd - inf_times[[2]])

  expected <- colSums(titre_matrix)

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

  perm_rise_pri <- 1
  temp_rise_pri <- 2
  perm_rise_sec <- 2
  temp_rise_sec <- 3

  temp_decay <- 0.003

  actual <- simulate_titres_seropositive(inf_times_sbs   = inf_times_sbs,
                                        sampling_times  = sampling_times,
                                        age             = age,
                                        subject_id      = 3,
                                        treatment_group = 0,
                                        perm_rise       = c(perm_rise_pri,
                                                            perm_rise_sec,
                                                            perm_rise_sec),
                                        temp_rise       = c(temp_rise_pri,
                                                            temp_rise_sec,
                                                            temp_rise_sec),
                                        temp_decay      = rep(0.003, 3),
                                        meas_sd         = 0.3)

  titre_matrix <- matrix(0, nrow = length(inf_times),
                         ncol = length(sampling_times))

  idx_1st   <- which(sampling_times_sbs > inf_times[[1]])
  times_1st <- sampling_times_sbs[idx_1st]

  titre_matrix[1, idx_1st] <- titre_decay_floor(
    par_alpha = perm_rise_pri,
    par_beta  = temp_rise_pri,
    par_delta = temp_decay,
    time      = times_1st - inf_times[[1]])

  idx_2nd   <- which(sampling_times_sbs > inf_times[[2]])
  times_2nd <- sampling_times_sbs[idx_2nd]

  titre_matrix[2, idx_2nd] <- titre_decay_floor(
    par_alpha = perm_rise_sec,
    par_beta  = temp_rise_sec,
    par_delta = temp_decay,
    time      = times_2nd - inf_times[[2]])

  idx_3rd   <- which(sampling_times_sbs > inf_times[[3]])
  times_3rd <- sampling_times_sbs[idx_3rd]

  titre_matrix[3, idx_3rd] <- titre_decay_floor(
    par_alpha = perm_rise_sec,
    par_beta  = temp_rise_sec,
    par_delta = temp_decay,
    time      = times_3rd - inf_times[[3]])

  true_titre <- colSums(titre_matrix)

  expected <-  data.frame(subject_id = 3,
                          time       = sampling_times,
                          true       = true_titre)

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)
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
