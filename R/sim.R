#' Simulate an ideal population
#'
#' Simulates an ideal population using the reference model from Tokars (2018).
#'
#' @param init_pop_size Integer initial population size
#' @param vaccinations Integer vector number of vaccinations at every timepoint
#' @param infections_novac Integer vector number of infections at every
#'   timepoint
#' @param ve Vaccine effectiveness (proportion)
#' @param lag Integer lag period measured in timepoints
#'
#' @return A \link[tibble]{tibble} with the following columns:
#'
#'   \item{timepoint}{Index of timepoint}
#'
#'   \item{vaccinations}{Expected number of vaccinations}
#'
#'   \item{infections_novac}{Expected number of infections in absence of
#'   vaccination}
#'
#'   \item{ve}{Expected vaccine effectiveness}
#'
#'   \item{pflu}{Flu incidence}
#'
#'   \item{infections}{Actual number of infections}
#'
#'   \item{popn}{Non-cases in absence of vaccination}
#'
#'   \item{pvac}{Proportion of starting population vaccinated}
#'
#'   \item{b}{Number vaccinated at that time who didn't get infected later}
#'
#'   \item{b_og}{Number vaccinated at that time}
#'
#'   \item{A}{Non-vaccinated non-cases}
#'
#'   \item{C}{Vaccinated susceptible}
#'
#'   \item{D}{Vaccinated immune}
#'
#'   \item{E}{Non-vaccinated infections cumulative total}
#'
#'   \item{F}{Vaccinated infections cumulative total}
#'
#' @references Tokars JI, Rolfes MA, Foppa IM, Reed C. An evaluation and update
#'   of methods for estimating the number of influenza cases averted by
#'   vaccination in the United States. Vaccine. 2018;36(48):7331–7337.
#'   doi:10.1016/j.vaccine.2018.10.026
#'
#' @importFrom tibble as_tibble
#'
#' @export
#'
#' @examples
#' # Population from Tokars (2018)
#' nsam <- 1e6L
#' ndays <- 304L
#' pop_tok <- sim_reference(
#'   init_pop_size = nsam,
#'   vaccinations = generate_counts(nsam, ndays, 0.55, mean = 100, sd = 50),
#'   infections_novac = generate_counts(nsam, ndays, 0.12, mean = 190, sd = 35),
#'   ve = 0.48,
#'   lag = 14
#' )
#' head(pop_tok)
#' sum(pop_tok$avert)
sim_reference <- function(init_pop_size,
                          vaccinations,
                          infections_novac,
                          ve,
                          lag) {
  if (length(ve) == 1) ve <- rep(ve, length(vaccinations))
  ideal_pop <- sim_reference_cpp(
    init_pop_size, vaccinations, infections_novac, ve, lag
  )
  attr(ideal_pop, "init_pop_size") <- init_pop_size
  attr(ideal_pop, "lag") <- lag
  as_tibble(ideal_pop)
}

#' Generate normal counts
#'
#' Generates counts from a normal distribution density function.
#'
#' @param init_pop_size Initial population size
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
#' vacs_tok <- generate_counts(1e6, 304, 0.55, 100, 50)
#' # Tokars (2018) cases
#' casen_tok <- generate_counts(1e6, 304, 0.12, 190, 35)
generate_counts <- function(init_pop_size, n_timepoints,
                            overall_prop, mean, sd) {
  densities <- dnorm(1:n_timepoints, mean, sd)
  density_coef <- overall_prop * init_pop_size / sum(densities)
  counts <- densities * density_coef
  as.integer(counts)
}

#' Generate dates
#'
#' Generate dates given timepoint indices, start date and step unit
#'
#' @param timepoints Integer vector timepoint indices
#' @param start Date of index 1
#' @param unit "year" "month" or "day"
#'
#' @return A vector of dates the same length as \code{timepoints}
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom lubridate day<- month<- year<- day month year ymd
#'
#' @export
#'
#' @examples
#' # Dates from Tokars (2018)
#' timepoints <- 1L:304L
#' dates <- generate_dates(timepoints, lubridate::ymd("2017-08-01"), "day")
generate_dates <- function(timepoints, start, unit) {
  if (unit == "day") {
    day(start) <- day(start) + timepoints - 1
  } else if (unit == "month") {
    month(start) <- month(start) + timepoints - 1
  } else if (unit == "year") {
    year(start) <- year(start) + timepoints - 1
  } else {
    abort(glue("unrecognised unit '{unit}'"))
  }
  start
}
