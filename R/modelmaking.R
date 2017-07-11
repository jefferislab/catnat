#' Generate a neuroanatomicla alpha shape from connector and/or tree node data
#'
#' @description implementation of the algorithm for clustering neurons by synapse location from Schlegel et al. (2016). Assumes neurons are scaled to microns.
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param substrate whether to make the model based off of connectors, neuron cable or both
#' @param maxdistance for automated cluster identification. Maximum distance at which nodes can be part of a cluster
#' @param groupsize an integer number of nearest neighbours to find using nabor::knn()
#' @param selection whether or not to interactively select values for maxdistance and groupsize.
#' @param chosen.points whether to feed the function pre-chosen points. A matrix for 3D points
#' @param alpha a single value or vector of values for α, fed to alpshaped3d::ashape3d(). Selection is subsequently interactive
#' @param auto.selection whether to try and remove points based on interactively chosen values for 'groupsize' and 'maxdistance'
#' @param ... additional arguments passed to methods
#'
#' @return An alphashape object
#' @export
make.anatomical.model <- function(someneuronlist, substrate = c("connectors","cable", "both"), maxdistance = 10000, groupsize = 100, alpha = 3000, auto.selection = T, chosen.points = NULL)
{
  if (substrate=="connectors"){synapse.points = nat::xyzmatrix(catmaid::connectors(someneuronlist))
  }else if(substrate =="cable"){synapse.points = nat::xyzmatrix(someneuronlist)
  }else if (substrate == "both"){synapse.points = rbind(nat::xyzmatrix(someneuronlist), nat::xyzmatrix(catmaid::connectors(someneuronlist)))}
  rgl::open3d()
  # Generate a good density of points to define neuropil
  if (auto.selection == TRUE){
    progress = "n"
    while (progress == "n"){
      groupsize <- as.numeric (readline(prompt="Select a value for the cluster groupsize  "))
      maxdistance <- as.numeric (readline(prompt="Select a value for maximum distance between points  "))
      neighbours = nabor::knn(synapse.points, synapse.points, k = groupsize)
      loose <- apply(neighbours$nn.dists, 1, function(x) {(any(as.numeric(x[1:ncol(neighbours$nn.dists)]) > maxdistance))})
      keep = c(neighbours$nn.idx[,1][!loose])
      close.points = synapse.points[keep,]
      rgl::clear3d();rgl::points3d(close.points, cl = 'black'); rgl::points3d(synapse.points, col = 'red')
      progress = readline(prompt="Continue? y/n  ")
    }
  }
  else{
    neighbours = nabor::knn(synapse.points, synapse.points, k = groupsize)
    loose <- apply(neighbours$nn.dists, 1, function(x) {(any(as.numeric(x[1:ncol(neighbours$nn.dists)]) > maxdistance))})
    keep = c(neighbours$nn.idx[,1][!loose])
    close.points = synapse.points[keep,]
    rgl::clear3d();rgl::points3d(close.points, cl = 'black'); rgl::points3d(synapse.points, col = 'red')
  }
  # Manual point deselection
  selected.points = unique(close.points)
  if (!is.null(chosen.points)){ selected.points = chosen.points}
  progress = readline(prompt="Remove (r), add (a) or save (s) points?  ")
  while (progress != "s"){
    if (progress == 'r'){
      remove.points <- select3d()
      removed.points <- remove.points(selected.points)
      selected.points = subset(selected.points, !removed.points)
      rgl::clear3d(); rgl::points3d(selected.points); rgl::points3d(synapse.points, col = 'red')
    }
    if (progress == 'a'){
      add.points <- select3d()
      added.points = subset(synapse.points, add.points(synapse.points))
      selected.points = rbind(selected.points, added.points)
      rgl::clear3d(); rgl::points3d(selected.points); rgl::points3d(synapse.points, col = 'red')
    }
    progress = readline(prompt="Remove (r), add (a) or save (s) points?  ")
  }
  progress = "n"
  while (progress == "n"){
    alpha <- as.numeric (readline(prompt="Select a value for alpha  "))
    alphashape = alphashape3d::ashape3d(unique(selected.points), alpha = alpha)
    plot(alphashape)
    progress = readline(prompt="Continue? y/n  ")
  }
  return (alphashape)
}

# Old and now redundant function
make.3d.tract=function (someneuronlist, savefile, maxdistance = 10, groupsize = 10, alpha = 300, selection = FALSE, scale = NULL, originalneurons = NULL,  chosen.points = NULL)
{
  neuron.points = xyzmatrix(someneuronlist)
  if (is.null(scale) == F) { neuron.points = xyzmatrix(someneuronlist)*scale}
  open3d()
  # Generate a good density of points to define tract
  if (selection == TRUE){
    progress = "n"
    while (progress == "n"){
      groupsize <- as.numeric (readline(prompt="Select a value for the groupsize  "))
      maxdistance <- as.numeric (readline(prompt="Select a value for maximum distance between points  "))
      neighbours = nabor::knn(neuron.points, neuron.points, k = groupsize)
      loose <- apply(neighbours$nn.dists, 1, function(x) {(any(as.numeric(x[1:ncol(neighbours$nn.dists)]) > maxdistance))})
      keep = c(neighbours$nn.idx[,1][!loose])
      close.points = neuron.points[keep,]
      clear3d();points3d(close.points, cl = 'black'); points3d(neuron.points, col = 'red')
      progress = readline(prompt="Continue? y/n  ")
    }
  }
  else{
    neighbours = nabor::knn(neuron.points, neuron.points, k = groupsize)
    loose <- apply(neighbours$nn.dists, 1, function(x) {(any(as.numeric(x[1:ncol(neighbours$nn.dists)]) > maxdistance))})
    keep = c(neighbours$nn.idx[,1][!loose])
    close.points = neuron.points[keep,]
    clear3d();points3d(close.points, cl = 'black'); points3d(neuron.points, col = 'red')
  }
  # Manual point deselection
  selected.points = unique(close.points)
  if (is.null(chosen.points) == F){ selected.points = chosen.points}
  progress = readline(prompt="Remove (r), add (a) or save (s) points?  ")
  while (progress != "s"){
    if (progress == 'r'){
      remove.points <- select3d()
      removed.points <- remove.points(selected.points)
      selected.points = subset(selected.points, !removed.points)
      clear3d(); points3d(selected.points); points3d(neuron.points, col = 'red')
    }
    if (progress == 'a'){
      add.points <- select3d()
      added.points = subset(neuron.points, add.points(neuron.points))
      selected.points = rbind(selected.points, added.points)
      clear3d(); points3d(selected.points); points3d(neuron.points, col = 'red')
    }
    progress = readline(prompt="Remove (r), add (a) or save (s) points?  ")
  }
  if (is.null(originalneurons) == F){
    synapse.points = find.neuropil.points(originalneurons, savefile, maxdistance, groupsize, alpha, selection = TRUE, scale = scale, tract = selected.points)
    selected.points = rbind(selected.points, synapse.points)
  }
  saveRDS(selected.points, file = savefile)
  cat("selected tract points saved as", savefile)
  progress = "n"
  while (progress == "n"){
    alpha <- as.numeric (readline(prompt="Select a value for alpha  "))
    alphashape = ashape3d(unique(selected.points), alpha = alpha)
    plot(alphashape)
    progress = readline(prompt="Continue? y/n  ")
  }
  return (alphashape)
}

#Convert alpha to mesh3d object

find.neuropil.points=function (someneuronlist, name, maxdistance = 1000, groupsize = 10, selection = FALSE, scale = TRUE, tract = NULL)
{
  synapse.points = get.synapses(someneuronlist)
  if (is.null(scale) == F) { synapse.points = get.synapses(someneuronlist)*scale}
  open3d()
  if (is.null(tract) == F) { points3d(tract, col = 'blue')}
  # Generate a good density of points to define neuropil
  if (selection == TRUE){
    progress = "n"
    while (progress == "n"){
      groupsize <- as.numeric (readline(prompt="Select a value for the groupsize  "))
      maxdistance <- as.numeric (readline(prompt="Select a value for maximum distance between points  "))
      neighbours = knn(synapse.points, synapse.points, k = groupsize)
      loose <- apply(neighbours$nn.dists, 1, function(x) {(any(as.numeric(x[1:ncol(neighbours$nn.dists)]) > maxdistance))})
      keep = c(neighbours$nn.idx[,1][!loose])
      close.points = synapse.points[keep,]
      clear3d();points3d(close.points, cl = 'black'); points3d(synapse.points, col = 'red')
      progress = readline(prompt="Continue? y/n  ")
    }
  }
  else{
    neighbours = knn(synapse.points, synapse.points, k = groupsize)
    loose <- apply(neighbours$nn.dists, 1, function(x) {(any(as.numeric(x[1:ncol(neighbours$nn.dists)]) > maxdistance))})
    keep = c(neighbours$nn.idx[,1][!loose])
    close.points = synapse.points[keep,]
    clear3d();points3d(close.points, cl = 'black'); points3d(synapse.points, col = 'red')
  }
  # Manual point deselection
  selected.points = unique(close.points)
  progress = readline(prompt="Remove points? y/n  ")
  while (progress == "y"){
    remove.points <- select3d()
    removed.points <- remove.points(selected.points)
    selected.points = subset(selected.points, !removed.points)
    clear3d(); points3d(selected.points); points3d(synapse.points, col = 'red')
    progress = readline(prompt="Remove more points? y/n  ")
  }
  return (selected.points)
}


#' @importFrom Morpho crossProduct
spinny <- function(object, target){
  rotate = 'go!'
  first = 0
  while (rotate != 'e'){
    if (first == 0){
      v1 = select.points(object)
      rgl::pop3d();v2 = select.points(object)
      rgl::pop3d();v3 = select.points(object)
    }
    if (rotate == 1){v1 = select.points(object)}
    if (rotate == 2){v2 = select.points(object)}
    if (rotate == 3){v3 = select.points(object)}
    # Check all is okay by plotting plane
    normal <- Morpho::crossProduct(v2-v1,v3-v1)
    zeroPro <- Morpho::points2plane(rep(0,3),v1,normal)
    ## get sign of normal displacement from zero
    sig <- sign(crossprod(-zeroPro,normal))
    d <- sig*norm(zeroPro,"2")
    rgl::planes3d(normal[1],normal[2],normal[3],d=d, alpha = 0.4)
    m = Morpho::mirror2plane(object, v1 = c(v1), v2 = c(v2), v3 = c(v3))
    rgl::points3d(object, col = 'blue')
    rgl::points3d(target, col = 'cyan')
    rgl::points3d(m, col = 'grey')
    rgl::points3d(rbind(v1,v2,v3), col = c(1,2,3), size = 40)
    first = first + 1
    rotate = readline(prompt="Change V1, V2, V3 or exit (e)?  ")
  }
  return(m)
}

#' Find neurites inside of an mesh3d object
#'
#' @param x 3D points
#' @param surf A 3D shape, typically a mesh3d object
#' @param rval Whetehr to return a logical, distances or a mesh3d object
#' @param ... additional arguments passed to methods
#'
#' @return A neuronlist
#' @export
#' @seealso \code{\link{neurons.inside}}
pointsinsidemesh <- function (x, surf, rval = c("logical", "distance", "mesh3d"),...)
{
  if (!requireNamespace("Rvcg", quietly = TRUE))
    stop("Please install suggested library Rvcg to use pointsinside")
  rval = match.arg(rval)
  pts = xyzmatrix(x)
  if (inherits(surf, "hxsurf")) {
    surf = as.mesh3d(surf, ...)
  }
  rmesh = Rvcg::vcgClostKD(pts, surf, sign = TRUE)
  switch(rval, logical = is.finite(rmesh$quality) & rmesh$quality >
           0 & rmesh$quality < 1e+12, distance = rmesh$quality,
         mesh3d = rmesh)
}


#' @aliases applyTransform
#' @importFrom Morpho applyTransform
applyTransform.neuron <- function(neuron, trafo, inverse = F){
  xyzmatrix(neuron$d)<-Morpho::applyTransform(xyzmatrix(neuron$d), trafo = trafo, inverse = inverse)
  neuron
}

#' @aliases applyTransform
#' @importFrom Morpho applyTransform
applyTransform.neuronlist <- function(someneuronlist, trafo, inverse = F){
  nlapply(someneuronlist, applyTransform.neuron, trafo, inverse)
}


#' Get the transformation matrix from the Morpho::mirror function
#'
#' @param x k x 3 matrix or mesh3d
#' @param subsample integer: use only a subset for icp matching
#' @param pcAlign if TRUE, the icp will be preceeded by an alignment of the principal axis (only used if icpiter > 0), currently only works for 3D data
#' @param mirroraxis integer: which axis to mirror at
#' @param initPC logical: if TRUE the data will be prealigned by its principal axes
#' @param initCenter logical: if TRUE and initPC=FALSE, x will be translated to its centroid before mirroring
#' @param mc.cores use parallel processing to find best alignment to original shape
#' @param ... additional arguments passed to methods
#'
#' @export
make.mirror.transform <- function (x, icpiter = 50, subsample = NULL, pcAlign = FALSE,
                                   mirroraxis = 1, initPC = TRUE, initCenter = TRUE, mc.cores = 2)
{
  m <- ncol(x)
  if (m == 2) {
    x <- cbind(x, 0)
    pcAlign <- FALSE
  }
  if (initPC) {
    pca <- Morpho::prcompfast(x, scale. = F)
    pca$rotation <- cbind(rbind(pca$rotation, 0), 0)
    pca$rotation[4, 4] <- 1
  }
  else if (initCenter) {
    xcenter <- scale(x, scale = F)
    rotmat <- Morpho::computeTransform(xcenter, x)
    pca <- list(x = xcenter, rotation = rotmat, center = rep(0,
                                                             3))
  }
  else {
    pca <- list(x = x, rotation = diag(m + 1), center = rep(0,
                                                            3))
  }
  mirmat <- diag(c(1, 1, 1))
  mirmat[mirroraxis, mirroraxis] <- -1
  out <- pca$x %*% t(mirmat)
  if (pcAlign)
    out <- Morpho::pcAlign(out, pca$x, iterations = icpiter, subsample = subsample,
                   mc.cores = mc.cores)
  else if (icpiter > 0)
    out <- Morpho::icpmat(out, pca$x, icpiter, subsample = subsample)
  out <- Morpho::applyTransform(out, pca$rotation, inverse = !initPC)
  out <- t(t(out) + pca$center)
  if (m == 2)
    out <- out[, 1:2]
  return(Morpho::computeTransform(x,out,type="tps"))
}




