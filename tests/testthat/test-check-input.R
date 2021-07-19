test_that("check_input_methods works", {
  expect_error(
    check_input_methods(
      init_pop_size = -1,
      vaccinations = 0,
      cases = 0,
      ve = 0.48
    ),
    "init_pop_size must be greater than 0"
  )
  expect_error(
    check_input_methods(
      init_pop_size = 100,
      vaccinations = 0,
      cases = 0,
      ve = 0.48
    ),
    "length of cases should be greater than 1"
  )
  expect_error(
    check_input_methods(
      init_pop_size = 100,
      vaccinations = 0,
      cases = generate_counts(1e6L, 304L, 0.12, 190, 35),
      ve = 0.48
    ),
    paste0(
      "length of cases \\(304\\) should match ",
      "length of vaccinations \\(1\\)"
    )
  )
  expect_error(
    check_input_methods(
      init_pop_size = 100,
      vaccinations = generate_counts(1e6L, 304L, 0.5, 50, 35),
      cases = generate_counts(1e6L, 304L, 0.12, 190, 35),
      ve = c(0.48, 0.52)
    ),
    paste0(
      "length of ve \\(2\\) should be either 1 or the same as ",
      "the length of vaccinations and cases \\(304\\)"
    )
  )
  expect_error(
    check_input_methods(
      init_pop_size = 100,
      vaccinations = generate_counts(1e6L, 304L, 0.5, 50, 35),
      cases = generate_counts(1e6L, 304L, 0.12, 190, 35),
      ve = c(0.48)
    ),
    NA
  )
})
