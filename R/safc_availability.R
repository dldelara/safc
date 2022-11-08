#' Calculate demand for each region
#'
#' @param supply Tibble with origin matrix and supply for each origin
#' @param demand Tibble with destination matrix and demand for each destination. May contain
#' optional bootstrapping information.
#' @param region Tibble with information to aggregate from destination to region level
#' @param distmat Optional distance matrix. If \code{distmat = NULL}, a distance matrix will be
#' generated automatically using the origin/destination information
#' @param max_dist A double describing the cutoff distance in the distance matrix calculation
#' @param n_boot Number of iterations to run the bootstrapping algorithm
#' @return A list containing model information and the estimated demand met for each region
#' @export
safc_availability = function(supply, demand, region, distmat = NULL, max_dist = Inf, n_boot = NULL) {
  params = list(
    n_origin = nrow(supply),
    n_dest   = nrow(demand),
    n_region = nrow(region)
  )
  if (is.null(distmat)) {
    if (is.infinite(max_dist)) warning("'max_dist' not specified. Distance matrix calculations may take long time")
    distmat = get_distmat(supply$origins, demand$destinations, max_dist, params$n_origin)
    if (!is.infinite(max_dist)) params$max_dist = max_dist
  }
  if (!is.null(n_boot)) {
    params$n_boot = n_boot
    demand$demand = get_demboot(
      as.numeric(factor(demand$region)) - 1,
      region$tot_demand,
      demand$weight,
      n_boot,
      params$n_dest,
      params$n_region
    )
  }
  availability = mutate(
    region,
    availability = get_availability(
      supply$supply,
      demand$demand,
      distmat,
      as.numeric(factor(demand$region)) - 1,
      params$n_region,
      params$n_origin
    )[, 1]
  )
  mod = list(
    availability = availability,
    supply = supply,
    demand = demand,
    distmat = distmat
  )
  mod
}
