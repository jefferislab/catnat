#' Assign cell body side based on neuron name
#'
#' @description Assigns cell body side based in neuron name.
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param ... additional arguments passed to methods
#'
#' @return Someneuronlist with cell sidedness in the metadata
#' @export
assignside <- function(someneuronlist, ...){
  sdf=as.data.frame(someneuronlist)
  sdf=transform(sdf, side=factor(ifelse(grepl("right|Right|_r$|R$|r$|left|Left|_l$|L$|l$", name),ifelse(grepl("right|Right|_r|R$|r$", name),"R","L"), "NA")))
  attr(someneuronlist,'df')=sdf
  return(someneuronlist)
}

# And old and now redundnt function
#convert <- function(someneuronlist, factor = 1/1e3){
#for (neuron in 1:length(someneuronlist)){
#  if (length(someneuronlist[[neuron]]$d) == 7){
#    someneuronlist[[neuron]]$d$X <- someneuronlist[[neuron]]$d$X*factor
#    someneuronlist[[neuron]]$d$Y <- someneuronlist[[neuron]]$d$Y*factor
#    someneuronlist[[neuron]]$d$Z <- someneuronlist[[neuron]]$d$Z*factor
#  }
#  if (length(someneuronlist[[neuron]]$connectors) == 6){
#    someneuronlist[[neuron]]$connectors$x <- someneuronlist[[neuron]]$connectors$x*factor
#    someneuronlist[[neuron]]$connectors$y <- someneuronlist[[neuron]]$connectors$y*factor
#    someneuronlist[[neuron]]$connectors$z <- someneuronlist[[neuron]]$connectors$z*factor
#  }
#}
#return (someneuronlist)
#}

# an old and now redundant function
#join.neuronlists <- function(...){
#  argg <- list(...)
#  skids = c()
#  for (item in 1:length(argg)){
#    skids = c(skids, c(as.data.frame(argg[[item]])$skid))
#  }
#  neurons = subset(db, as.data.frame(db)$skid%in%skids)
#  return(assignside(neurons))
#}

#primary.neurite.first <- function(someneuron, k = 100){      # Find the first 100 points of the primary neurite
 # som = soma.neuron(someneuron)
#  if (is.na(som[1])){
#    if (length(someneuron[[1]]$tags$soma[[1]])>0){
#      som = matrix(xyzmatrix(someneuron)[someneuron[[1]]$d$PointNo%in%someneuron[[1]]$tags$soma,], ncol = 3)
 #   }else{
  #    som = matrix(xyzmatrix(someneuron)[someneuron[[1]]$StartPoint,], ncol = 3)
#    }
#  }
#  p = nat::xyzmatrix(someneuron)
#  n = nabor::knn(p, som, ifelse(nrow(p)>k,k,nrow(p)))
#  m = p[c(n$nn.idx),]
#}

#' Returns the primary neurite of a neuron
#'
#' @description Returns the primary neurite of a neuron, defined as the cable between soma and first branch point
#'
#' @param x a neuron or neuronlist object
#' @param resample The newspacing with which to evenly resample each neuron. Can be set to F to prevent resampling.
#' @param k the number of nodes from the soma to include
#' @param kee.pnt whether to keep just the primary neurite tract of the skeleton, or remove it and keep the neuron's arbour instead
#' @param ... additional arguments passed to methods
#'
#' @return A neuron pruned to its primary dendrite
#' @export
primary.neurite<-function(x,resample = 1, keep.pnt = T, ...) UseMethod("primary.neurite")

#' @export
#' @rdname primary.neurite
primary.neurite.neuron <- function(x, resample = 1, keep.pnt = T,...){
  som <- if (is.null(x$tags$soma)){
    warning("No soma found, using startpoint")
    som = x$StartPoint
  } else somapos(x)
  sl=as.seglist(x, flatten=T, all=T)
  som_seg=which(sapply(sl, function(x) som%in%x))
  if(length(som_seg)>1){
    warning("more than one segment contains root, choosing first!")
    som_seg=som_seg[1]
  } else if(length(som_seg)==0) stop("no segment found for soma!")

  pn=nat::prune_vertices(x, verticestoprune = sl[[som_seg]], invert = keep.pnt)
  if(is.numeric(resample)){
    resample(pn, stepsize=resample)
  }else{
    pn
  }
}

#' @export
#' @rdname primary.neurite
primary.neurite.neuronlist <- function(x, resample = 1, keep.pnt = T, ...){
  nlapply(x, primary.neurite.neuron, OmitFailures = T, resample = 1, keep.pnt = T,...)
}

#' @export
#' @rdname primary.neurite
starting.neurite<-function(x, k, ...) UseMethod("starting.neurite")

#' @export
#' @rdname primary.neurite
starting.neurite.neuron <- function(x, k = 100, ...){      # Find the first 100 points of the primary neurite
  som = soma.neuron(x)
  if (is.na(som[1])){
    if (length(x[[1]]$tags$soma[[1]])>0){
      som = matrix(xyzmatrix(x)[x[[1]]$d$PointNo%in%x[[1]]$tags$soma,], ncol = 3)
    }else{
      som = matrix(xyzmatrix(x)[x[[1]]$StartPoint,], ncol = 3)
    }
  }
  p = nat::xyzmatrix(x)
  n = nabor::knn(p, som, ifelse(nrow(p)>k,k,nrow(p)))
  m = p[c(n$nn.idx),]
}

#' @export
#' @rdname primary.neurite
starting.neurite.neuronlist  <- function(x, k = 100, ...){
  nlapply(x, starting.neurite.neuron, OmitFailures = T)
}

# @export
# @rdname primary.neurite
longestpathfromsoma = function (x, UseStartPoint = TRUE, SpatialWeights = TRUE, invert = FALSE,
                        rval = c("neuron", "length", "ids"), model = NULL)
{
  ng <- as.ngraph(n, weights = SpatialWeights)
  rval = match.arg(rval)
  if (invert && rval == "length")
    stop("invert=TRUE is not implemented for rval='length'")
  if (UseStartPoint) {
    lps = igraph::shortest.paths(graph = ng, n$StartPoint, to = n$EndPoints,
                         mode = "all")
    if (rval == "length")
      return(max(lps))
    to = n$EndPoints[which.max(lps)]
    longestpath = igraph::get.shortest.paths(ng, from = n$StartPoint,
                                     to = to, mode = "all")$vpath[[1]]
  }
  else {
    if (rval == "length") {
      return(igraph::diameter(ng, directed = FALSE))
    }
    else {
      longestpath = igraph::get.diameter(ng, directed = FALSE)
    }
  }
  if (rval == "ids") {
    if (invert) {
      ie = igraph::difference(igraph::E(ng), igraph::E(ng,
                                                       path = longestpath))
      edgemat = igraph::ends(ng, ie, names = FALSE)
      return(unique(as.integer(t(edgemat))))
    }
    else return(as.integer(longestpath))
  }
  prune_edges(ng, edges = longestpath, invert = !invert)
}






