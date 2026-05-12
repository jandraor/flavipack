test_that("simulate_titre_trajectory() works", {

  subject_id      <- 1
  inf_times       <- 1047.108
  sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)

  actual <- simulate_titre_trajectory(sampling_times,
                                      subject_id      = subject_id,
                                      exposure_times  = inf_times,
                                      peaks           = 3,
                                      perm_rises      = 1,
                                      decays          = 0.003,
                                      baseline        = 0,
                                      meas_sd         = 0.3)

  expected <- data.frame(subject_id = subject_id,
                         time       = sampling_times,
                         true       = c(rep(0, length(sampling_times) - 2),
                                        1.81340282, 1.25019146))

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)

  # Two infections

  inf_times <- c(300, 900)

  peak_pri      <- 3
  perm_rise_pri <- 1

  peak_sec      <- 5
  perm_rise_sec <- 2

  temp_decay <- 0.003

  actual <- simulate_titre_trajectory(sampling_times,
                                      subject_id      = subject_id,
                                      exposure_times  = inf_times,
                                      peaks           = c(peak_pri,
                                                          peak_sec),
                                      perm_rises      = c(perm_rise_pri,
                                                          perm_rise_sec),
                                      decays          = c(temp_decay,
                                                          temp_decay),
                                      baseline        = 0,
                                      meas_sd         = 0.3)
  titre_matrix <- matrix(0,
                         nrow = length(inf_times),
                         ncol = length(sampling_times))

  idx_1st   <- which(sampling_times > inf_times[[1]])
  times_1st <- sampling_times[idx_1st]

  titre_matrix[1, idx_1st] <- titre_decay_floor(
    peak      = peak_pri,
    perm_rise = perm_rise_pri,
    decay     = temp_decay,
    time      = times_1st - inf_times[[1]])

  idx_2nd   <- which(sampling_times > inf_times[[2]])
  times_2nd <- sampling_times[idx_2nd]

  titre_matrix[2, idx_2nd] <- titre_decay_floor(
    peak      = peak_sec,
    perm_rise = perm_rise_sec,
    decay     = temp_decay,
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
    subject_id      = 1,
    baseline        = 0,
    exposure_times  = numeric(0),
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

  peak_pri      <- 3
  perm_rise_pri <- 1

  peak_sec <- 5
  perm_rise_sec <- 2


  temp_decay    <- 0.003

  actual <- simulate_true_titre_DENV(
    sampling_times   = sampling_times,
    peaks            = c(peak_pri,
                         peak_sec),
    perm_rises       = c(perm_rise_pri,
                         perm_rise_sec),
    decays           = c(temp_decay,
                         temp_decay),
    baseline         = 0,
    exposure_times   = inf_times)

  titre_matrix <- matrix(0,
                         nrow = length(inf_times),
                         ncol = length(sampling_times))

  idx_1st   <- which(sampling_times > inf_times[[1]])
  times_1st <- sampling_times[idx_1st]

  titre_matrix[1, idx_1st] <- titre_decay_floor(
    perm_rise = perm_rise_pri,
    peak      = peak_pri,
    decay     = temp_decay,
    time      = times_1st - inf_times[[1]])

  idx_2nd   <- which(sampling_times > inf_times[[2]])
  times_2nd <- sampling_times[idx_2nd]

  titre_matrix[2, idx_2nd] <- titre_decay_floor(
    perm_rise = perm_rise_sec,
    peak      = peak_sec,
    decay     = temp_decay,
    time      = times_2nd - inf_times[[2]])

  expected <- colSums(titre_matrix)

  expect_equal(actual, expected)
})

test_that("simulate_true_titre_DENV() works",
{
  sampling_times <- c(537, 571, 627, 655, 808, 999, 1340, 1749, 2138)

  actual <- simulate_true_titre_DENV(
    sampling_times = sampling_times,
    perm_rises = 5,
    peaks      = 5.9,
    decays     = 0.001,
    baseline   = 0,
    exposure_times = 537)

  sim_times <- sampling_times - min(sampling_times)

  expected <- 5 + 0.9 * exp(-0.001 * sim_times)

  expect_equal(actual, expected)
})

test_that("simulate_true_titre_DENV() adds the baseline",
{
  sampling_times <- c(100, 200, 300)

  actual <- simulate_true_titre_DENV(sampling_times  = sampling_times,
                                     peaks           = 6,
                                     perm_rises      = 2,
                                     decays          = 0.01,
                                     baseline        = 3,
                                     exposure_times  = 150)

  expected <- c(3, 5 + 4 * exp(-0.01 * c(50, 150)))

  expect_equal(actual, expected)
})

test_that("simulate_true_titre_DENV() returns the baseline with no infections",
{
  sampling_times <- c(100, 200, 300)

  actual <- simulate_true_titre_DENV(sampling_times = sampling_times,
                                     peaks           = numeric(0),
                                     perm_rises      = numeric(0),
                                     decays          = numeric(0),
                                     baseline        = 3,
                                     exposure_times = numeric(0))

  expected <- c(3, 3, 3)

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

  peak_pri      <- 3
  perm_rise_pri <- 1

  peak_sec      <- 5
  perm_rise_sec <- 2

  decay <- 0.003

  actual <- simulate_titres_seropositive(inf_times_sbs   = inf_times_sbs,
                                        sampling_times  = sampling_times,
                                        age             = age,
                                        subject_id      = 3,
                                        peaks           = c(peak_pri,
                                                            peak_sec,
                                                            peak_sec),
                                        perm_rises      = c(perm_rise_pri,
                                                            perm_rise_sec,
                                                            perm_rise_sec),
                                        decays          = rep(decay, 3),
                                        meas_sd         = 0.3)

  titre_matrix <- matrix(0, nrow = length(inf_times),
                         ncol = length(sampling_times))

  idx_1st   <- which(sampling_times_sbs > inf_times[[1]])
  times_1st <- sampling_times_sbs[idx_1st]

  titre_matrix[1, idx_1st] <- titre_decay_floor(
    peak      = peak_pri,
    perm_rise = perm_rise_pri,
    decay     = decay,
    time      = times_1st - inf_times[[1]])

  idx_2nd   <- which(sampling_times_sbs > inf_times[[2]])
  times_2nd <- sampling_times_sbs[idx_2nd]

  titre_matrix[2, idx_2nd] <- titre_decay_floor(
    peak      = peak_sec,
    perm_rise = perm_rise_sec,
    decay     = decay,
    time      = times_2nd - inf_times[[2]])

  idx_3rd   <- which(sampling_times_sbs > inf_times[[3]])
  times_3rd <- sampling_times_sbs[idx_3rd]

  titre_matrix[3, idx_3rd] <- titre_decay_floor(
    peak      = peak_sec,
    perm_rise = perm_rise_sec,
    decay     = decay,
    time      = times_3rd - inf_times[[3]])

  true_titre <- colSums(titre_matrix)

  expected <-  data.frame(subject_id = 3,
                          time       = sampling_times,
                          true       = true_titre)

  expect_equal(actual[, 1:3], expected)

  expect_equal("meas" %in% colnames(actual), TRUE)
})
