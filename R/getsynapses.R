#' Get 3D coordinates of synapse (not connector) positions
#'
#' @param x a neuronlist or neuron object
#' @param target whether post or presynapse are ot be returned, or both
#' @param polypre whether to consider the number of presynapses as a multiple of the numbers of connections each makes
#' @param ... additional arguments passed to methods.
#'
#' @return Anatomically accurate synapse position (i.e. not just connector positions) as 3D coordinates
#' @export
#' @seealso \code{\link{clusterbysynapses}} \code{\link{flow.centrality}}
get.synapses <-function(x, target = c("BOTH", "PRE", "POST"), polypre = T, ...) UseMethod("get.synapses")

#' @export
#' @rdname get.synapses
get.synapses.neuron <- function (x, target = c("BOTH", "PRE", "POST"), polypre = T, ...){
  if (target%in%c("POST","BOTH")) {
    syns.in = x$connectors[x$connectors[,3]==1,][,1]
    point.no = rownames(x$d)[x$d[,"PointNo"]%in%syns.in]
    points = nat::xyzmatrix(x$d[point.no,])
    points = cbind(points, prepost = 1)
  }
  if (target%in%c("PRE","BOTH")) {
    pres = x$connectors[x$connectors[,3]==0,][,2]
    syns.out = x$connectors[,1][match(pres, x$connectors[,2])]
    point.no = rownames(x$d)[x$d[,"PointNo"]%in%syns.out]
    points.out = nat::xyzmatrix(x$d[point.no,])
    points.out = cbind(points.out, prepost = 0)
    if(target == "BOTH"){
      points = rbind(points, points.out)
    }else{ points = points.out}
  }
  points
}

#' @export
#' @rdname get.synapses
get.synapses.neuronlist <- function (x, target = c("BOTH", "PRE", "POST"), polypre = T, ...){
  points = nat::nlapply(x, get.synapses.neuron, target = target, polypre = polypre)
  points = do.call(rbind, points)
  points
}
