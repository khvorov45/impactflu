# Tests of simulation functions
# Arseniy Khvorov
# Created 2019/12/24
# Last edit 2019/12/24

test_that("count generation works", {
  set.seed(1)
  test_counts <- generate_counts(
    init_pop_size = 1e6, n_timepoints = 304,
    overall_prop = 0.55, mean = 100, sd = 50
  )
  expect_equal(length(test_counts), 304)
  expect_equal(sum(test_counts) / 1e6, 0.55, tol = 0.01)
})

test_that("simulation works", {
  pop <- sim_ideal(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    ve = 0.48,
    lag = 14L,
    seed = 1L
  )
  expect_equal(attr(pop, "seed"), 1L)
  expect_equal(attr(pop, "init_pop_size"), 1e6L)
  expect_equal(attr(pop, "lag"), 14L)
  expect_equal(pop$timepoint, 1L:304L)
})
