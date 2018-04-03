#' Calculate the cable length inside of different neuropils in a segmented .surf brain
#'
#' @description Calculates the cable length supplied by neurons 'x' to different brain regions defined in 'brain'.
#'
#' @param x a set of neurons or, for points_in_neuropil a mxn matrix of 3D points
#' @param brain the .surf brainspace in which the neurons are registered, must be segmented into neuropils
#' @param method whether to use the neurons' axons or dendrites, or both
#' @param min.endpoints the minimum number of endpoints a neuron must have in a neuropil to be counted as included in it
#' @param alpha the alpha given to the ashape3d() function to generate neuropil objects by which to calculate point inclusion
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname cable_inside_neuropils
cable_inside_neuropils<-function(x, brain = nat.flybrains::FCWBNP.surf, method = c("neurites","axons","dendrites"), min.endpoints = 1,alpha=30, ...) UseMethod("cable_inside_neuropils")

#' @export
#' @rdname cable_inside_neuropils
cable_inside_neuropils.neuron <- function(x, brain = nat.flybrains::FCWBNP.surf, method = c("neurites","axons","dendrites"), min.endpoints = 1,alpha=30, ...){
  if(!requireNamespace('nat.flybrains', quietly = TRUE))
    stop("You must install suggested package nat.flybrains to use this function!")
  method = method[1]
  targets = c(-100:100)
  if (method=="axons"){targets = c(-2,2)}
  if (method=="dendrites"){targets = c(-3,3)}
  endings <- function(x){
    points=x$d[nat::endpoints(x)[which(endpoints(x)!=nat::rootpoints(x))],]
    nat::xyzmatrix(points[points$Label%in%targets,])
  }
  sapply(brain$RegionList, function(n) in_neuropil(x=x,brain = brain,neuropil=n,min.endpoints=min.endpoints,alpha=alpha))
}

#' @export
#' @rdname cable_inside_neuropils
cable_inside_neuropils.neuronlist <- function(x, brain = nat.flybrains::FCWBNP.surf, method = c("neurites","axons","dendrites"), min.endpoints = 1,alpha=30, ...){
  nat::nlapply(x, cable_inside_neuropils.neuron, brain=brain, method=method)
}

#' @export
#' @rdname cable_inside_neuropils
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
#' @rdname cable_inside_neuropils
in_neuropil<-function(x,brain = nat.flybrains::FCWBNP.surf,neuropil = "LH_R",min.endpoints =1,alpha=alpha, ...) UseMethod("in_neuropil")

#' @export
#' @rdname cable_inside_neuropils
in_neuropil.neuron <- function(x,brain = nat.flybrains::FCWBNP.surf,neuropil = "LH_R",min.endpoints =1,alpha=alpha){
  neuropil = subset(brain,n)
  neuropil = alphashape3d::ashape3d(nat::xyzmatrix(neuropil),alpha=alpha)
  if(sum(alphashape3d::inashape3d(points=endings(x),as3d=neuropil, indexAlpha = "ALL"))>min.endpoints){
    points = x$d[x$d$Label%in%targets,]
    points.in = points[alphashape3d::inashape3d(points=nat::xyzmatrix(points),as3d=neuropil, indexAlpha = "ALL"),]
    v = rownames(points.in)
    pruned = tryCatch(nat::prune_vertices(x,verticestoprune=v,invert = TRUE),error = function(e)NULL)
    if(!is.null(pruned)){
      summary(pruned)$cable.length
    }else{0}
  }else{0}
}

#' @export
#' @rdname cable_inside_neuropils
in_neuropil.neuronlist <- function(x,brain = nat.flybrains::FCWBNP.surf,neuropil = "LH_R",min.endpoints =1,alpha=alpha){
  nat::nlapply(x, in_neuropil.neuron, brain=brain, neuropil=neuropil,min.endpoints=min.endpoints,alpha=alpha)
}
