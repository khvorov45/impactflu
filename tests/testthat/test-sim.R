# Tests of simulation functions
# Arseniy Khvorov
# Created 2019/12/24
# Last edit 2019/12/24

test_that("count generation works", {
  test_counts <- generate_counts(
    init_pop_size = 1e6, n_timepoints = 304,
    overall_prop = 0.55, mean = 100, sd = 50
  )
  expect_equal(length(test_counts), 304)
})
