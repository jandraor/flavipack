test_that("log2_transform() works",
{
  expect_equal(log2_transform(10), 1)

  expect_equal(log2_transform(5), 0)
})
