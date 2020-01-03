# Tests of methods
# Arseniy Khvorov
# Created 2020/01/03
# Last edit 2020/01/03

library(dplyr)

test_that("method1 works", {
  pop <- sim_ideal(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    ve = 0.48,
    lag = 0L,
    seed = 1L,
    deterministic = TRUE
  )
  pop_monthly <- pop %>%
    mutate(
      dates = generate_dates(timepoint, lubridate::ymd("2017/08/01"), "day"),
      month = lubridate::month(dates),
      year = lubridate::year(dates)
    ) %>%
    group_by(year, month) %>%
    summarise(
      vaccinations = sum(vaccinations),
      cases = sum(cases),
      ve = mean(ve)
    ) %>%
    ungroup()
  m1 <- method1(
    attr(pop, "init_pop_size"),
    pop_monthly$vaccinations, pop_monthly$cases, pop_monthly$ve
  )
  with(m1, {
    expect_equal(vaccinations, pop_monthly$vaccinations)
    expect_equal(cases, pop_monthly$cases)
    expect_equal(ve, pop_monthly$ve)
    expect_equal(vc_lag, (pvac + lag(pvac, default = 0)) / 2)
    expect_equal(
      pops,
      as.integer((lag(pops, default = attr(pop, "init_pop_size")) -
        lag(cases, default = 0L)) * (1 - vc_lag * ve))
    )
    expect_equal(pflu, cases / pops)
    expect_equal(
      popn, lag(popn, default = attr(pop, "init_pop_size")) -
        lag(cases_novac, default = 0L)
    )
    expect_equal(
      avert, cases_novac - cases
    )
  })
})
