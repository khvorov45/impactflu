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
