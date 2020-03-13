#' Calculate the cable length inside of different neuropils in a segmented .surf brain
#'
#' @description Calculates the cable length supplied by neurons 'x' to different brain regions defined in 'brain'.
#'
#' @param x a set of neurons or, for points_in_neuropil a mxn matrix of 3D points
#' @param neuropil the code for the neuropil volume to be searched.
#' @param brain the .surf brainspace in which the neurons are registered, must be segmented into neuropils
#' @param method whether to calculate cable inside a given neuropil volume, or synapse number (PRE or POST)
#' @param min.endpoints the minimum number of endpoints a neuron must have in a neuropil to be counted as included in it
#' @param alpha the alpha given to the ashape3d() function to generate neuropil objects by which to calculate point inclusion.
#' If \code{NULL} then \code{nat::pointsinside} is used.
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname inside_neuropils
inside_neuropils<-function(x, brain, method = c("cable","PRE","POST"), min.endpoints = 1, alpha=NULL, ...) UseMethod("inside_neuropils")

#' @export
#' @rdname inside_neuropils
inside_neuropils.neuron <- function(x, brain, method = c("cable","PRE","POST"), min.endpoints = 1, alpha=NULL, ...){
  method = match.arg(method)
  sapply(brain$RegionList, function(n) in_neuropil(x=x,method = method, brain = brain,neuropil=n,min.endpoints=min.endpoints,alpha=alpha))
}

#' @export
#' @rdname inside_neuropils
inside_neuropils.neuronlist <- function(x, brain, method = c("cable","PRE","POST"), min.endpoints = 1, alpha=NULL, ...){
  method = match.arg(method)
  nat::nlapply(x, inside_neuropils.neuron, brain=brain, method=method, min.endpoints = min.endpoints, alpha = alpha, ...)
}

#' @export
#' @rdname inside_neuropils
points_in_neuropil <- function(x, brain, alpha = NULL, ...){
  method = match.arg(method)
  nps = brain$RegionList
  df = cbind(as.data.frame(x),neuropil=0)
  for (n in nps){
    neuropil = subset(brain, n)
    if(!is.null(alpha)){
      neuropil = alphashape3d::ashape3d(xyzmatrix(neuropil),alpha=alpha)
      a = alphashape3d::inashape3d(points=nat::xyzmatrix(x),as3d=neuropil,indexAlpha = "ALL")
    }else{
      a = nat::pointsinside(x=nat::xyzmatrix(x),surf=neuropil)
    }
    df$neuropil[which(a==T)] = n
  }
  df
}

#' @export
#' @rdname inside_neuropils
in_neuropil<-function(x,
                      brain,
                      neuropil,
                      method = c("cable","PRE","POST"),
                      min.endpoints =1,alpha=NULL, ...) UseMethod("in_neuropil")

#' @export
#' @rdname inside_neuropils
in_neuropil.neuron <- function(x,
                               brain,
                               neuropil,
                               method = c("cable","PRE","POST"),
                               min.endpoints =1,
                               alpha=NULL,
                               ...){
  method = match.arg(method)
  neuropil = subset(brain,neuropil)
  if(method == "cable"){
    points = x$d
  }else if (method == "PRE"){
    points = nat::xyzmatrix(x$connectors[x$connectors$prepost==0,])
  }else if (method == "POST"){
    points = nat::xyzmatrix(x$connectors[x$connectors$prepost==1,])
  }
  if(!is.null(alpha)){
    neuropil = alphashape3d::ashape3d(nat::xyzmatrix(neuropil),alpha=alpha)
    go = sum(alphashape3d::inashape3d(points=endings(x),as3d=neuropil, indexAlpha = "ALL"))>min.endpoints
    if(go){
      points.in = points[alphashape3d::inashape3d(points=nat::xyzmatrix(points),as3d=neuropil, indexAlpha = "ALL"),]
    }
  }else{
    go = sum(nat::pointsinside(x=endings(x),surf=neuropil))>min.endpoints
    if(go){
      points.in = points[nat::pointsinside(x=nat::xyzmatrix(points),surf=neuropil),]
    }
  }
  if(go & method == "cable"){
    v = rownames(points.in)
    pruned = tryCatch(nat::prune_vertices(x,verticestoprune=v,invert=TRUE),error = function(e)NULL)
    summary(pruned)$cable.length
  }else if (go){
    nrow(points.in)
  }else{
    0
  }
}


#' @export
#' @rdname inside_neuropils
in_neuropil.neuronlist <- function(x,
                                   brain = nat.flybrains::FCWBNP.surf,
                                   neuropil = "LH_R",
                                   method = c("cable","PRE","POST"),
                                   min.endpoints =1,
                                   alpha=NULL, ...){
  method = match.arg(method)
  nat::nlapply(x, in_neuropil.neuron, brain=brain, neuropil=neuropil,min.endpoints=min.endpoints,alpha=alpha,...)
}

# hidden
endings <- function(x){
  points=x$d[nat::endpoints(x$d)[which(nat::endpoints(x$d)!=nat::rootpoints(x))],]
  EndNo = nat::endpoints(x$d)[which(nat::endpoints(x$d)!=nat::rootpoints(x))]
  points = subset(x$d,PointNo%in%EndNo)
  nat::xyzmatrix(points)
}
