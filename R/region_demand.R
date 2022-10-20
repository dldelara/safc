#' Calculate demand for each region
#'
#' @param destinations Dataframe of the destinations of each point of demand
#' @param origins Dataframe of the origins of each point of supply
#' @param demand Demand associated with each destination region
#' @param population Population counts in each destination subregion
#' @param n_iter Number of iterations to run the bootstrapping algorithm
#' @return A list containing a dataframe with the access amounts and a dataframe with each
#' iteration's demand.
#' @export
region_demand = function(destinations, origins, demand, population, n_iter) {
  #. = NAME = V1 = V2 = V3 = V4 = dem_met = demand_met = fips = geometry = pop_65up = shp = NULL
  start      = Sys.time()
  orgs       = as.matrix(origins[, 1:2])
  dests      = as.matrix(destinations[, 1:2])
  supply     = origins$tot_benes
  catch_fips = tibble(fips = destinations$fips) %>%
    group_by(fips) %>%
    group_indices(fips) %>%
    `-`(1)
  catchment = NULL
  tot_dem = matrix(0, nrow = length(unique(catch_fips)), n_iter)
  for (i in 1:n_iter) {
    if (i == 1) cat("Generating catchment estimates...\n")
    demand = destinations %>% group_by(fips) %>%
      transmute(demand = rmultinom(1, demand, pop_65up)) %>%
      .$demand
    catchment = as_tibble(.catchment_areas(orgs, dests, supply, demand, catch_fips)) %>%
      rename(fips = V1, facility = V2, demand = V3, cumdemand = V4) %>%
      mutate(fips = unique(destinations$fips)[fips + 1])
    tot_dem[, i] = catchment %>%
      group_by(fips) %>%
      summarize(tot_dem = sum(demand)) %>%
      right_join(tibble(fips = unique(destinations$fips)), by = "fips") %>%
      arrange(fips) %>%
      transmute(tot_dem = ifelse(is.na(tot_dem), 0, tot_dem)) %>% as_vector()
    cat("Iteration", i, "/", n_iter, "\r")
    if (i == n_iter) {
      cat("\n")
      print(round(Sys.time() - start, 1))
    }
  }
  access = tibble(fips = unique(destinations$fips), dem_met = apply(tot_dem, 1, median)) %>%
    right_join(shp, by = "fips") %>%
    mutate(demand_met = dem_met / demand) %>%
    select(fips, dem_met, demand, demand_met, name, geometry)
  if (n_iter == 1) return(list(access = access, catchment = catchment))
  return(list(tot_dem = tot_dem, access = access))
}
