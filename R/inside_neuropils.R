#' Calculate the cable length inside of different neuropils in a segmented .surf brain
#'
#' @description Calculates the cable length supplied by neurons 'x' to different brain regions defined in 'brain'.
#'
#' @param x a set of neurons or, for points_in_neuropil a mxn matrix of 3D points
#' @param brain the .surf brainspace in which the neurons are registered, must be segmented into neuropils
#' @param method whether to calculate cable inside a given neuropil volume, or synapse number (PRE or POST)
#' @param min.endpoints the minimum number of endpoints a neuron must have in a neuropil to be counted as included in it
#' @param alpha the alpha given to the ashape3d() function to generate neuropil objects by which to calculate point inclusion
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname inside_neuropils
inside_neuropils<-function(x, brain = nat.flybrains::FCWBNP.surf, method = c("Cable","PRE","POST"), min.endpoints = 1,alpha=30, ...) UseMethod("inside_neuropils")

#' @export
#' @rdname inside_neuropils
inside_neuropils.neuron <- function(x, brain = nat.flybrains::FCWBNP.surf, method = c("Cable","PRE","POST"), min.endpoints = 1,alpha=30, ...){
  if(!requireNamespace('nat.flybrains', quietly = TRUE))
    stop("You must install suggested package nat.flybrains to use this function!")
  method = method[1]
  sapply(brain$RegionList, function(n) in_neuropil(x=x,method = method, brain = brain,neuropil=n,min.endpoints=min.endpoints,alpha=alpha))
}

#' @export
#' @rdname inside_neuropils
inside_neuropils.neuronlist <- function(x, brain = nat.flybrains::FCWBNP.surf, method = c("Cable","PRE","POST"), min.endpoints = 1,alpha=30, ...){
  nat::nlapply(x, inside_neuropils.neuron, brain=brain, method=method)
}

#' @export
#' @rdname inside_neuropils
points_in_neuropil <- function(x, brain, alpha = 30, ...){
  nps = brain$RegionList
  df = cbind(as.data.frame(x),neuropil=0)
  for (n in nps){
    neuropil = subset(brain, n)
    neuropil = alphashape3d::ashape3d(xyzmatrix(neuropil),alpha=alpha)
    a = alphashape3d::inashape3d(points=nat::xyzmatrix(x),as3d=neuropil,indexAlpha = "ALL")
    df$neuropil[which(a==T)] = n
  }
  df
}

#' @export
#' @rdname inside_neuropils
in_neuropil<-function(x,method = c("Cable","PRE","POST"),brain = nat.flybrains::FCWBNP.surf,neuropil = "LH_R",min.endpoints =1,alpha=alpha, ...) UseMethod("in_neuropil")

#' @export
#' @rdname inside_neuropils
in_neuropil.neuron <- function(x,method = c("Cable","PRE","POST"),brain = nat.flybrains::FCWBNP.surf,neuropil = "LH_R",min.endpoints =1,alpha=alpha){
  neuropil = subset(brain,neuropil)
  neuropil = alphashape3d::ashape3d(nat::xyzmatrix(neuropil),alpha=alpha)
  endings <- function(x){
    points=x$d[nat::endpoints(x$d)[which(nat::endpoints(x$d)!=nat::rootpoints(x))],]
    EndNo = nat::endpoints(x$d)[which(nat::endpoints(x$d)!=nat::rootpoints(x))]
    points = subset(x$d,PointNo%in%EndNo)
    nat::xyzmatrix(points)
  }
  if(sum(alphashape3d::inashape3d(points=endings(x),as3d=neuropil, indexAlpha = "ALL"))>min.endpoints){
    points = x$d
    points.in = points[alphashape3d::inashape3d(points=nat::xyzmatrix(points),as3d=neuropil, indexAlpha = "ALL"),]
    v = rownames(points.in)
    if("catmaidneuron"%in%class(x)){
      pruned = tryCatch(prune_vertices.catmaidneuron(x,verticestoprune=v,invert=TRUE),error = function(e)NULL)
      class(pruned) = c("catmaidneuron","neuron")
    }else{
      pruned = tryCatch(nat::prune_vertices(x,verticestoprune=v,invert=TRUE),error = function(e)NULL)
    }
    if(!is.null(pruned)&method=="Cable"){
      summary(pruned)$cable.length
    }else if(!is.null(pruned)&method%in%c("PRE","POST")){
      nrow(get.synapses(pruned,target=method))
    }else{0}
  }else{0}
}


#' @export
#' @rdname inside_neuropils
in_neuropil.neuronlist <- function(x,method = c("Cable","PRE","POST"),brain = nat.flybrains::FCWBNP.surf,neuropil = "LH_R",min.endpoints =1,alpha=alpha){
  nat::nlapply(x, in_neuropil.neuron, brain=brain, neuropil=neuropil,min.endpoints=min.endpoints,alpha=alpha)
}
