test_that("create_input_list() works",
{
  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        location      = "COL",
                        serostatus    = TRUE)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(20),
                         meas       = c(2),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  actual <- create_input_list(part_df, titre_df, symp_df, n_markers = 1)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     = list(
          "mean" = I(c(2))
         ),
         infection_times  = I(list(11, 15))))

  expect_equal(actual, expected)
})

test_that("create_input_list() works with multiple individuals",
{
  part_df <- data.frame(subject_id    = c(1, 10),
                        age_enrolment = c(10, 8),
                        location      = c("COL", "BRA"),
                        serostatus    = c(1, 0))

  titre_df <- data.frame(subject_id = c(1, 1, 10, 10),
                         time       = c(5, 20, 4, 8),
                         meas       = c(5, 6, 7, 7),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  actual <- create_input_list(part_df, titre_df, symp_df, n_markers = 1)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = 1,
         obs_times        = I(list(5, 20)),
         measurements     = list("mean" = I(c(5, 6))),
         infection_times  = I(list(11, 15))
         ),
    list(subject_id       = "10",
         age_enrolment    = 8,
         location         = "BRA",
         serostatus       = 0,
         obs_times        = I(list(4, 8)),
         measurements     = list("mean" = I(c(7, 7)))
         )
    )

  expect_equal(actual, expected)
})

test_that("create_input_list() works with vaccination",
{
  part_df <- data.frame(subject_id    = c(1, 10, 20),
                        age_enrolment = c(10, 8, 5),
                        location      = c("COL", "BRA", "THA"),
                        serostatus    = c(1, 0, 1))

  titre_df <- data.frame(subject_id = c(1, 1, 10, 10, 20, 20),
                         time       = c(5, 20, 4, 8, 4, 8),
                         meas       = c(5, 6, 7, 7, 10, 15),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  vacc_df <- data.frame(subject_id     = 20,
                        time           = 4.5,
                        baseline_titre = 10)

  actual <- create_input_list(part_df, titre_df, symp_df, vacc_df,
                              n_markers = 1)

  expected <- list(
    list(subject_id        = "1",
         age_enrolment     = 10,
         location          = "COL",
         serostatus        = 1,
         obs_times         = I(list(5, 20)),
         measurements      = list("mean" = I(c(5, 6))),
         infection_times   = I(list(11, 15))
         ),
    list(subject_id        = "10",
         age_enrolment     = 8,
         location          = "BRA",
         serostatus        = 0,
         obs_times         = I(list(4, 8)),
         measurements      = list("mean" = I(c(7, 7)))
         ),
    list(subject_id        = "20",
         age_enrolment     = 5,
         location          = "THA",
         serostatus        = 1,
         obs_times         = I(list(4, 8)),
         measurements      = list("mean" = I(c(10, 15))),
         vaccination       = list(
           list(time           = 4.5,
                baseline_titre = 10))
         )
    )

  expect_equal(actual, expected)
})

test_that("create_input_list() works with multiple vaccinations",
{
  part_df <- data.frame(subject_id    = c(30),
                        age_enrolment = c(9),
                        location      = c("DOM"),
                        serostatus    = c(1))

  titre_df <- data.frame(subject_id = c(30),
                         time       = c(5, 20, 30, 40),
                         meas       = c(1, 2, 3, 4),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 30,
                        time = c(25, 35))

  vacc_df <- data.frame(subject_id     = c(30, 30),
                        time           = c(5.5, 20.5),
                        baseline_titre = c(1, 2))

  actual <- create_input_list(part_df, titre_df, symp_df, vacc_df,
                              n_markers = 1)

  expected <- list(
    list(subject_id        = "30",
         age_enrolment     = 9,
         location          = "DOM",
         serostatus        = 1,
         obs_times         = I(list(5, 20, 30, 40)),
         measurements      = list("mean" = I(c(1, 2, 3, 4))),
         infection_times   = I(list(25, 35)),
         vaccination       = list(
           list(time = 5.5,
                baseline_titre = 1),
           list(time           = 20.5,
                baseline_titre = 2))
         )
    )

  expect_equal(actual, expected)
})

test_that("create_input_list() excludes exposures after dropout",
{
  part_df <- data.frame(subject_id    = c(12),
                        age_enrolment = c(10),
                        location      = c("COL"),
                        serostatus    = c(1))

  titre_df <- data.frame(subject_id = c(12),
                         time       = c(5, 10, 15),
                         meas       = c(5, 5, 5),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 12,
                        time       = 16)

  vacc_df <- data.frame(subject_id     = 12,
                        time           = 5.5,
                        baseline_titre = 5)

  actual <- create_input_list(part_df, titre_df, symp_df, vacc_df,
                              n_markers = 1)

  expected <- list(
    list(subject_id    = "12",
         age_enrolment = 10,
         location      = "COL",
         serostatus    = 1,
         obs_times     = I(list(5, 10, 15)),
         measurements  = list("mean" = I(c(5, 5, 5))),
         vaccination  = list(
           list(time           = 5.5,
                baseline_titre = 5))
    )
  )

  expect_equal(actual, expected)
})

test_that("create_input_list() handles individuals dropping out before a
          follow-up measurement after vaccination",
{
  part_df <- data.frame(subject_id    = c(177),
                        age_enrolment = c(8),
                        location      = c("COL"),
                        serostatus    = c(1))

  titre_df <- data.frame(subject_id = c(177),
                         time       = c(394),
                         meas       = c(7),
                         marker     = "mean")

  symp_df <- data.frame()

  vacc_df <- data.frame(subject_id     = 177,
                        time           = 394.5,
                        baseline_titre = 7)

  actual <- create_input_list(part_df, titre_df, symp_df, vacc_df,
                              n_markers = 1)

  expected <- list(
    list(subject_id    = "177",
         age_enrolment = 8,
         location      = "COL",
         serostatus    = 1,
         obs_times     = I(list(394)),
         measurements  = list("mean" = I(c(7))),
         vaccination  = list(
           list(time           = 394.5,
                baseline_titre = 7))
    )
  )

  expect_equal(actual, expected)
})

test_that("create_input_list() handles multiple markers",
{
  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        location      = "COL",
                        serostatus    = TRUE)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(20),
                         meas       = c(2, 2.5, 1.5, 2),
                         marker     = c("D1", "D2", "D3", "D4"))

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  actual <- create_input_list(part_df, titre_df, symp_df, n_markers = 4)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     =
           list(
             "D1" = I(c(2)),
             "D2" = I(c(2.5)),
             "D3" = I(c(1.5)),
             "D4" = I(c(2))
           ),
         infection_times  = I(list(11, 15))))

  expect_equal(actual, expected)
})

test_that("create_input_list() handle starting infections", {

  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        location      = "COL",
                        serostatus    = TRUE)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(20),
                         meas       = c(2),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  guess_df <- data.frame(subject_id = 1,
                         time       = 1)

  actual <- create_input_list(part_df  = part_df,
                              titre_df = titre_df,
                              symp_df  = symp_df,
                              n_markers = 1,
                              guess_df = guess_df)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     = list(
           "mean" = I(c(2))
         ),
         infection_times    = I(list(11, 15)),
         guessed_infections = I(list(1))
    )
  )

  expect_equal(actual, expected)
})

test_that("create_input_list() handle symptomatic infections with flicks", {

  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        location      = "COL",
                        serostatus    = TRUE)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(20),
                         meas       = c(2),
                         marker     = "mean")

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  guess_df <- data.frame(subject_id = 1,
                         time       = 1)

  flick_df <- data.frame(subject_id = 1,
                         time       = 15)

  actual <- create_input_list(part_df  = part_df,
                              titre_df = titre_df,
                              symp_df  = symp_df,
                              n_markers = 1,
                              guess_df = guess_df,
                              flick_df = flick_df)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     = list(
           "mean" = I(c(2))
         ),
         infection_times    = I(list(11, 15)),
         guessed_infections = I(list(1)),
         flicks             = I(list(15))
    )
  )

  expect_equal(actual, expected)

  flick_df <- data.frame(subject_id = 1,
                         time       = 16)

  expect_error(
    create_input_list(
      part_df   = part_df,
      titre_df  = titre_df,
      symp_df   = symp_df,
      n_markers = 1,
      guess_df  = guess_df,
      flick_df  = flick_df)
  )
})

test_that("create_input_list() handle symptomatic infections with flicks", {

  part_df <- data.frame(subject_id    = c(1, 10),
                        age_enrolment = c(5, 10),
                        location      = c("COL", "LKA"),
                        serostatus    = c(TRUE, TRUE))

  titre_df <- data.frame(subject_id = c(1, 10),
                         time       = c(20, 20),
                         meas       = c(2, 2),
                         marker     = c("mean", "mean"))

  symp_df <- data.frame(subject_id = 10,
                        time       = 15)

  flick_df <- data.frame(subject_id = 10,
                         time       = 15)

  actual <- create_input_list(part_df  = part_df,
                              titre_df = titre_df,
                              symp_df  = symp_df,
                              n_markers = 1,
                              flick_df = flick_df)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 5,
         location         = "COL",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     = list(
           "mean" = I(c(2))
         )),
    list(subject_id       = "10",
         age_enrolment    = 10,
         location         = "LKA",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     = list(
           "mean" = I(c(2))
         ),
         infection_times    = I(list(15)),
         flicks             = I(list(15))
    )
  )

  expect_equal(actual, expected)
})


test_that("create_measurement_object() works",
{
  subject_titre_df <- data.frame(subject_id = 1,
                                 time       = c(20),
                                 meas       = c(2, 2.5, 1.5, 2),
                                 marker     = c("D1", "D2", "D3", "D4"))

  actual <- create_measurement_object(subject_titre_df)

  expected <- list(
    "D1" = I(c(2)),
    "D2" = I(c(2.5)),
    "D3" = I(c(1.5)),
    "D4" = I(c(2))
  )

  expect_equal(actual, expected)
})

