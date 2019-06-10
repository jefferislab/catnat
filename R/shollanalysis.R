#' Sholl analysis on neuron skeletons
#'
#' @description Functions for Sholl analysis of neuronal skeletons
#'
#' @param neuron a neuron object
#' @param start the origin from which spheres are grown for the Sholl analysis
#' @param starting.radius the radius of the first sphere. Defaults to the radius step
#' @param ending.radius the radius of the last sphere. If NULL the distance to the furthest dendritic point from the start point is taken
#' @param radius.step the change in radius between successive spheres. Defaults to one 100th of the radius of the ending sphere
#' @return a data.frame of spheres radii and the number of dendritic intersections at each radius
#' @export
#' @rdname sholl_data
sholl_data <- function(neuron, start = colMeans(xyzmatrix(neuron)), starting.radius = radius.step, ending.radius = NULL, radius.step = ending.radius/100){
  message("Reminder: Neuron should be appropriately resampled and dendrites marked (Label = 3)")
  unit.vector <- function(x) {x / sqrt(sum(x^2))}
  dend = neuron$d
  dend$dists = nabor::knn(data = matrix(start,ncol=3), query = nat::xyzmatrix(neuron),k=1)$nn.dists
  if(is.null(ending.radius)){
    ending.radius = max(dend$dists)
  }
  radii = seq(from = starting.radius, to = ending.radius, by = radius.step)
  sholl = data.frame(radii = radii, intersections = 0)
  for(n in 1:length(radii)){
    r = radii[n]
    segments = neuron$SegList
    for(segment in segments){
      p = dend[segment,]
      dists = (nabor::knn(data = matrix(start,ncol=3), query = nat::xyzmatrix(p),k=1)$nn.dists - r) >= 0
      sholl[n,]$intersections = sholl[n,]$intersections + lengths(regmatches(paste(dists,collapse=""), gregexpr("TRUEFALSE|FALSETRUE", paste(dists,collapse=""))))
    }
  }
  sholl
}

