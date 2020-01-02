# Methods from Tokars (2018)
# Arseniy Khvorov
# Created 2020/01/03
# Last edit 2020/01/03

#' Methods from Tokars (2018)
#'
#' Method 1 was described as current. Method 3 was determined to be least
#' biased.
#'
#' @param init_pop_size Integer initial population size
#' @param vaccinations Integer vector counts of vaccinations
#' @param cases Integer vector counts of cases
#' @param ve Vector vaccine effectiveness. If length 1, assumed to not vary with
#'   time.
#'
#' @return
#' A \link[tibble]{tibble} with the following columns (method-dependent):
#'   \item{cases}{Observed cases}
#'   \item{vaccinations}{Observed vaccinations}
#'   \item{ve}{Assumed vaccine effectiveness}
#'   \item{vc_lag}{Vaccine coverage lagged}
#'   \item{pops}{Susceptible population}
#'   \item{pflu}{Infection risk}
#'   \item{popn}{Non-cases is absence of vaccination}
#'   \item{cases_novac}{Cases in absense of vaccination}
#'   \item{avert}{Expected number of vaccinations}
#'
#' @references Tokars JI, Rolfes MA, Foppa IM, Reed C. An evaluation and update
#'   of methods for estimating the number of influenza cases averted by
#'   vaccination in the United States. Vaccine. 2018;36(48):7331–7337.
#'   doi:10.1016/j.vaccine.2018.10.026
#'
#' @importFrom tibble as_tibble
#'
#' @export
method1 <- function(init_pop_size, vaccinations, cases, ve) {
  check_counts(vaccinations, cases, ve)
  as_tibble(method1_cpp(init_pop_size, vaccinations, cases, ve))
}
