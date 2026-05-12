test_that("create_input_list() works",
{
  part_df <- data.frame(subject_id    = 1,
                        age_enrolment = 10,
                        location      = "COL",
                        serostatus    = TRUE)

  titre_df <- data.frame(subject_id = 1,
                         time       = c(20),
                         meas       = c(2))

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  actual <- create_input_list(part_df, titre_df, symp_df, n_markers = 1)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = TRUE,
         obs_times        = I(list(20)),
         measurements     = list(I(list(2))),
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
                         meas       = c(5, 6, 7, 7))

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  actual <- create_input_list(part_df, titre_df, symp_df, n_markers = 1)

  expected <- list(
    list(subject_id       = "1",
         age_enrolment    = 10,
         location         = "COL",
         serostatus       = 1,
         obs_times        = I(list(5, 20)),
         measurements     = list(I(list(5)), I(list(6))),
         infection_times  = I(list(11, 15))
         ),
    list(subject_id       = "10",
         age_enrolment    = 8,
         location         = "BRA",
         serostatus       = 0,
         obs_times        = I(list(4, 8)),
         measurements     = list(I(list(7)), I(list(7)))
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
                         meas       = c(5, 6, 7, 7, 10, 15))

  symp_df <- data.frame(subject_id = 1,
                        time = c(11, 15))

  vacc_df <- data.frame(subject_id = 20,
                        time       = 5)

  actual <- create_input_list(part_df, titre_df, symp_df, vacc_df,
                              n_markers = 1)

  expected <- list(
    list(subject_id        = "1",
         age_enrolment     = 10,
         location          = "COL",
         serostatus        = 1,
         obs_times         = I(list(5, 20)),
         measurements      = list(I(list(5)), I(list(6))),
         infection_times   = I(list(11, 15))
         ),
    list(subject_id        = "10",
         age_enrolment     = 8,
         location          = "BRA",
         serostatus        = 0,
         obs_times         = I(list(4, 8)),
         measurements      = list(I(list(7)), I(list(7)))
         ),
    list(subject_id        = "20",
         age_enrolment     = 5,
         location          = "THA",
         serostatus        = 1,
         obs_times         = I(list(4, 8)),
         measurements      = list(I(list(10)), I(list(15))),
         vaccination_times = I(list(5))
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
                         meas       = c(5, 5, 5))

  symp_df <- data.frame(subject_id = 12,
                        time       = 16)

  vacc_df <- data.frame(subject_id = 12,
                        time       = 6)

  actual <- create_input_list(part_df, titre_df, symp_df, vacc_df,
                              n_markers = 1)

  expected <- list(
    list(subject_id        = "12",
         age_enrolment     = 10,
         location          = "COL",
         serostatus        = 1,
         obs_times         = I(list(5, 10, 15)),
         measurements      = list(I(list(5)),
                                  I(list(5)),
                                  I(list(5))),
         vaccination_times = I(list(6))
    )
  )

  expect_equal(actual, expected)
})
