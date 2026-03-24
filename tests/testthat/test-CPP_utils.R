test_that("format_symp_infections() works",
{
  symp_df <- data.frame(subject_id      = c(1, 1),
                        infection_time  = c(10, 20))

  actual <- format_symp_infections(symp_df)

  expected <- data.frame(subject_id            = 1,
                         infection_time_symp_1 = 10,
                         code_symp_1           = 2,
                         infection_time_symp_2 = 20,
                         code_symp_2           = 2)

  expect_equal(actual, expected)

  symp_df <- data.frame(subject_id      = c(1, 1, 2),
                        infection_time  = c(10, 20, 15))

  actual <- format_symp_infections(symp_df)

  expected <- data.frame(subject_id            = c(1, 2),
                         infection_time_symp_1 = c(10, 15),
                         code_symp_1           = c(2, 2),
                         infection_time_symp_2 = c(20, NA),
                         code_symp_2           = c(2, NA))

  expect_equal(actual, expected)
})

test_that("create_input_df() works",
{
  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        country       = "COL",
                        serostatus    = 1)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(10, 20),
                         meas = c(2, 4))

  symp_df <- data.frame(subject_id = 1, time = c(11, 15))


  actual <- create_input_df(part_df, titre_df, symp_df,
                            max_n_meas = 2,
                            max_n_symp = 2)

  expected <- data.frame(
    subject_id            = 1,
    age_enrolment         = 10,
    country               = "COL",
    serostatus            = 1,
    time_1                = 10,
    meas_1                = 2,
    time_2                = 20,
    meas_2                = 4,
    infection_time_symp_1 = 11,
    code_symp_1           = 2,
    infection_time_symp_2 = 15,
    code_symp_2           = 2)

  expect_equal(actual, expected)
})

test_that("create_input_df() works with no symp infections",
{
  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        country       = "COL",
                        serostatus    = 1)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(10, 20),
                         meas = c(2, 4))

  symp_df <- data.frame()


  actual <- create_input_df(part_df, titre_df, symp_df,
                            max_n_meas = 2,
                            max_n_symp = 2)

  expected <- data.frame(
    subject_id            = 1,
    age_enrolment         = 10,
    country               = "COL",
    serostatus            = 1,
    time_1                = 10,
    meas_1                = 2,
    time_2                = 20,
    meas_2                = 4,
    infection_time_symp_1 = -1,
    code_symp_1           = -1,
    infection_time_symp_2 = -1,
    code_symp_2           = -1)

  expect_equal(actual, expected)
})

test_that("create_input_df() handles less titres meas than expected",
{
  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        country       = "COL",
                        serostatus    = 1)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(10),
                         meas       = c(2))

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))


  actual <- create_input_df(part_df, titre_df, symp_df,
                            max_n_meas = 2,
                            max_n_symp = 2)

  expected <- data.frame(
    subject_id            = 1,
    age_enrolment         = 10,
    country               = "COL",
    serostatus            = 1,
    time_1                = 10,
    meas_1                = 2,
    time_2                = -1,
    meas_2                = -1,
    infection_time_symp_1 = 11,
    code_symp_1           = 2,
    infection_time_symp_2 = 15,
    code_symp_2           = 2)

  expect_equal(actual, expected)

})
