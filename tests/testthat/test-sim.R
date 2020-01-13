# Tests of simulation functions
# Arseniy Khvorov
# Created 2019/12/24
# Last edit 2019/12/24

library(testthat)

test_that("count generation works", {
  set.seed(1)
  test_counts <- generate_counts(
    init_pop_size = 1e6, n_timepoints = 304,
    overall_prop = 0.55, mean = 100, sd = 50
  )
  expect_equal(length(test_counts), 304)
  expect_equal(sum(test_counts) / 1e6, 0.55, tol = 0.01)
})

test_that("date generation works", {
  timepoints <- 1:5
  start <- lubridate::ymd("2000-01-01")
  expect_equal(
    generate_dates(timepoints, start, "day"),
    c(
      lubridate::ymd("2000-01-01"), lubridate::ymd("2000-01-02"),
      lubridate::ymd("2000-01-03"), lubridate::ymd("2000-01-04"),
      lubridate::ymd("2000-01-05")
    )
  )
  expect_equal(
    generate_dates(timepoints, start, "month"),
    c(
      lubridate::ymd("2000-01-01"), lubridate::ymd("2000-02-01"),
      lubridate::ymd("2000-03-01"), lubridate::ymd("2000-04-01"),
      lubridate::ymd("2000-05-01")
    )
  )
  expect_equal(
    generate_dates(timepoints, start, "year"),
    c(
      lubridate::ymd("2000-01-01"), lubridate::ymd("2001-01-01"),
      lubridate::ymd("2002-01-01"), lubridate::ymd("2003-01-01"),
      lubridate::ymd("2004-01-01")
    )
  )
  expect_error(
    generate_dates(timepoints, start, "wrong"),
    "unrecognised unit 'wrong'"
  )
})

test_that("simulation works with no lag, no deaths and 0 dur", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    infections_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    deaths_novac = rep(0L, 304L),
    ve = 0.48,
    lag = 0L,
    dur = 0L,
    seed = 1L,
    deterministic = TRUE
  )
  expect_named(
    pop, c(
      "timepoint", "vaccinations", "infections_novac", "deaths_novac",
      "ve", "pflu", "popn", "pvac", "b", "b_og", "A", "C", "D", "e", "e_og",
      "f", "f_og", "I", "J", "infections", "avert"
  ))
  expect_equal(attr(pop, "seed"), 1L)
  expect_equal(attr(pop, "init_pop_size"), 1e6L)
  expect_equal(attr(pop, "lag"), 0L)
  expect_equal(attr(pop, "deterministic"), TRUE)
  with(pop, {
    expect_equal(timepoint, 1L:304L)
    expect_equal(pflu, infections_novac / dplyr::lag(popn, default = 1e6L))
    expect_equal(
      infections,
      as.integer(round(pflu * dplyr::lag(A, default = 1e6L), 0)) +
      as.integer(round(pflu * dplyr::lag(C, default = 0L), 0))
    )
    expect_equal(popn, dplyr::lag(popn, default = 1e6L) - infections_novac)
    expect_equal(avert, infections_novac - infections)
    expect_equal(
      pvac, vaccinations /
        (dplyr::lag(A, default = 1e6L) + dplyr::lag(I, default = 0L))
    )
    expect_equal(b, as.integer(round(pvac * dplyr::lag(A, default = 1e6L), 0)))
    expect_equal(b_og, b) # No time to get infected
    expect_equal(
      A, dplyr::lag(A, default = 1e6L) -
        as.integer(round(dplyr::lag(A, default = 1e6L) * pflu), 0) - b
    )
    expect_equal(
      C, dplyr::lag(C, default = 0L) -
        as.integer(round(dplyr::lag(C, default = 0L) * pflu, 0)) +
        as.integer(round(b * (1 - ve), 0))
    )
    expect_equal(
      D, dplyr::lag(D, default = 0L) + b - as.integer(round(b * (1 - ve), 0))
    )
    expect_equal(
      e, as.integer(round(dplyr::lag(A, default = 0L) * pflu, 0))
    )
    expect_equal(e, e_og)
    expect_equal(
      f, as.integer(round(dplyr::lag(C, default = 0L) * pflu, 0))
    )
    expect_equal(f, f_og)
    expect_equal(
      I, dplyr::lag(I, default = 0L) + e -
        as.integer(round(dplyr::lag(I, default = 0L) * pvac, 0))
    )
    expect_equal(
      J, dplyr::lag(J, default = 0L) + f +
        as.integer(round(dplyr::lag(I, default = 0L) * pvac, 0))
    )
  })
})

test_that("simulation works with lag", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 150, 50),
    infections_novac = generate_counts(1e6L, 304L, 0.35, 150, 35),
    deaths_novac = rep(0L, 304L),
    ve = 0.48,
    lag = 3L,
    dur = 0L,
    deterministic = TRUE
  )
  with(pop, {
    b_1 <- b_og -
      as.integer(round(b_og * dplyr::lead(pflu, n = 1L, default = 0L)))
    b_2 <- b_1 -
      as.integer(round(b_1 * dplyr::lead(pflu, n = 2L, default = 0L)))
    expect_equal(
      b, b_2 -
        as.integer(round(b_2 * dplyr::lead(pflu, n = 3L, default = 0L)))
    )
    expect_equal(
      C, dplyr::lag(C, n = 1L, default = 0L) +
        as.integer(round(dplyr::lag(b, n = 3L, default = 0L) * (1 - 0.48))) -
        as.integer(round(dplyr::lag(C, n = 1L, default = 0L) * pflu))
    )
    expect_equal(
      D, dplyr::lag(D, n = 1L, default = 0L) +
        dplyr::lag(b, n = 3L, default = 0L) -
        as.integer(round(dplyr::lag(b, n = 3L, default = 0L) * (1 - 0.48)))
    )
  })
})

test_that("simulation works with non-0 dur", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 150, 50),
    infections_novac = generate_counts(1e6L, 304L, 0.35, 150, 35),
    deaths_novac = rep(0L, 304L),
    ve = 0.48,
    lag = 0L,
    dur = 3L,
    deterministic = TRUE
  )
  with(pop, {
    expect_equal(
      I, dplyr::lag(I, default = 0L) + dplyr::lag(e, n = 3L, default = 0L) -
        as.integer(round(dplyr::lag(I, default = 0L) * pvac))
    )
    expect_equal(
      J, dplyr::lag(J, default = 0L) + dplyr::lag(f, n = 3L, default = 0L) +
        as.integer(round(dplyr::lag(I, default = 0L) * pvac))
    )
  })
})

test_that("random simulation works", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    infections_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    deaths_novac = rep(0L, 304L),
    ve = 0.48,
    lag = 0L,
    dur = 14L,
    seed = 1L,
    deterministic = FALSE
  )
  sum1 <- sum(pop$avert)
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    infections_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    deaths_novac = rep(0L, 304L),
    ve = 0.48,
    lag = 0L,
    dur = 14L,
    seed = 1L,
    deterministic = FALSE
  )
  sum2 <- sum(pop$avert)
  expect_equal(sum1, sum2)
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    infections_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    deaths_novac = rep(0L, 304L),
    ve = 0.48,
    lag = 0L,
    dur = 14L,
    seed = 2L,
    deterministic = FALSE
  )
  sum3 <- sum(pop$avert)
  expect_true(sum1 != sum3)
})
