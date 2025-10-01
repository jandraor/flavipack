test_that("log2_transform() works",
{
  expect_equal(log2_transform(10), 1)

  expect_equal(log2_transform(5), 0)
})

test_that("inv_log2_transform() works",
{
  expect_equal(inv_log2_transform(1), 10)

  expect_equal(inv_log2_transform(0), 5)
})
