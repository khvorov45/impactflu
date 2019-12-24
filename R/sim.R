# Simulation functions
# Arseniy Khvorov
# Created 2019/12/24
# Last edit 2019/12/24

#' Simulate an ideal population
#'
#' Simulates an ideal population using the reference model from Tokars (2018).
#'
#' @param init_pop_size Inital population size
#' @param vaccinations Number of vaccinations at every timepoint
#' @param cases_novac Number of cases at every timepoint
#' @param ve Vaccine effectiveness
#' @param lag Lag period measured in timepoints
#' @param seed Integer seed to use
#'
#' @return A \link[tibble]{tibble} with the following columns:
#'   \item{timepoint}{Index of timepoint}
#'   \item{vaccinations}{Expected number of vaccinations}
#'   \item{cases_novac}{Expected number of cases in absence of vaccination}
#'   \item{ve}{Expected vaccine effectiveness}
#'   \item{pflu}{Flu incidence}
#'   \item{cases}{Actual number of cases}
#'   \item{popn}{Non-cases in absence of vaccination}
#'   \item{pvac}{Proportion of starting population vaccinated}
#'   \item{b}{Number vaccinated at that time}
#'   \item{A_to_E}{Number moved from A to E at that time}
#'   \item{A}{Non-vaccinated non-cases}
#'   \item{B}{Vaccinated non-cases lagging}
#'   \item{E}{Non-vaccinated cases}
#'
#' @references Tokars JI, Rolfes MA, Foppa IM, Reed C. An evaluation and update
#'   of methods for estimating the number of influenza cases averted by
#'   vaccination in the United States. Vaccine. 2018;36(48):7331–7337.
#'   doi:10.1016/j.vaccine.2018.10.026
#'
#' @importFrom tibble as_tibble
#'
#' @export
sim_ideal <- function(init_pop_size,
                      vaccinations,
                      cases_novac,
                      ve,
                      lag,
                      seed = sample.int(.Machine$integer.max, 1)) {
  ideal_pop <- sim_ideal_cpp(init_pop_size, vaccinations, cases_novac, ve, lag)
  attr(ideal_pop, "seed") <- seed
  attr(ideal_pop, "init_pop_size") <- init_pop_size
  attr(ideal_pop, "lag") <- lag
  as_tibble(ideal_pop)
}

#' Generate normal counts
#'
#' Generates counts from a normal distribution density function.
#'
#' @param init_pop_size Inital population size
#' @param n_timepoints Number of timepoints
#' @param overall_prop Overall proportion of the population to be included in
#'   the counts over all the timepoints
#' @param mean Mean of the normal distribution
#' @param sd Standard deviation of the normal distribution
#'
#' @return An integer vector of counts of length \code{n_timepoints}
#'
#' @importFrom stats dnorm
#'
#' @export
#'
#' @examples
#' # Tokars (2018) vaccinations
#' generate_counts(1e6, 304, 0.55, 100, 50)
#' # Tokars (2018) cases
#' generate_counts(1e6, 304, 0.12, 190, 35)
generate_counts <- function(init_pop_size, n_timepoints,
                            overall_prop, mean, sd) {
  densities <- dnorm(1:n_timepoints, mean, sd)
  density_coef <- overall_prop * init_pop_size / sum(densities)
  counts <- densities * density_coef
  as.integer(counts)
}
