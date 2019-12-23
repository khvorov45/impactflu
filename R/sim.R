# Simulation functions
# Arseniy Khvorov
# Created 2019/12/24
# Last edit 2019/12/24

#' Simulate an ideal population
#'
#' @param init_pop_size Inital population size
#' @param vaccinations Number of vaccinations at every timepoint
#' @param cases_novac Number of cases at every timepoint
#' @param ve Vaccine effectiveness
#' @param lag Lag period measured in timepoints
#'
#' @return A \link[tibble]{tibble} with the following columns:
#' \item{timepoint}{Index of timepoint}
#'
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#'
#' @export
sim_ideal <- function(init_pop_size,
                      vaccinations,
                      cases_novac,
                      ve,
                      lag) {
  sim_ideal_cpp(init_pop_size, vaccinations, cases_novac, ve, lag) %>%
    as_tibble()
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
