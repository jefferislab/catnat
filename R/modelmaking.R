# Model construction related functions

make.neuropil=function (someneuronlist, savefile, maxdistance = 1000, groupsize = 10, alpha = 3000, selection = FALSE, nano2microns = F, chosen.points = NULL)
{
  #synapse.points = xyzmatrix(connectors(someneuronlist))
  if (nano2microns == T) { synapse.points = get.connectors(someneuronlist)/1e3}
  else {synapse.points = get.connectors(someneuronlist)}
  open3d()
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
  if (is.null(chosen.points) == F){ selected.points = chosen.points}
  progress = readline(prompt="Remove (r), add (a) or save (s) points?  ")
  while (progress != "s"){
    if (progress == 'r'){
      remove.points <- select3d()
      removed.points <- remove.points(selected.points)
      selected.points = subset(selected.points, !removed.points)
      clear3d(); points3d(selected.points); points3d(synapse.points, col = 'red')
    }
    if (progress == 'a'){
      add.points <- select3d()
      added.points = subset(synapse.points, add.points(synapse.points))
      selected.points = rbind(selected.points, added.points)
      clear3d(); points3d(selected.points); points3d(synapse.points, col = 'red')
    }
    progress = readline(prompt="Remove (r), add (a) or save (s) points?  ")
  }
  saveRDS(selected.points, file = savefile)
  cat("selected points saved as", savefile)
  progress = "n"
  while (progress == "n"){
    alpha <- as.numeric (readline(prompt="Select a value for alpha  "))
    alphashape = ashape3d(unique(selected.points), alpha = alpha)
    plot(alphashape)
    progress = readline(prompt="Continue? y/n  ")
  }
  return (alphashape)
}


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
    neighbours = knn(neuron.points, neuron.points, k = groupsize)
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
      clear3d(); points3d(selected.points); points3d(neuron.points, col = 'red'); points3d(nc82, col = 'grey')
    }
    if (progress == 'a'){
      add.points <- select3d()
      added.points = subset(neuron.points, add.points(neuron.points))
      selected.points = rbind(selected.points, added.points)
      clear3d(); points3d(selected.points); points3d(neuron.points, col = 'red'); points3d(nc82, col = 'grey')
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


find.neuropil.points=function (someneuronlist, name, maxdistance = 1000, groupsize = 10, selection = FALSE, scale = TRUE, tract = NULL)
{
  synapse.points = get.connectors(someneuronlist)
  if (is.null(scale) == F) { synapse.points = get.connectors(someneuronlist)*scale}
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


transform.vtk = function (vtk, transformations){
  positions = xyzmatrix(ReadVTKLandmarks(vtk, item = "points"))
  if (is.list(transformations) == F){
    cat("Single transformation")
    positions <- xform(positions, transformations)
  }
  if (is.list(transformations) == T){
    for (transformation in transformations){
      positions <- xform(positions, transformation)
    }
  }
  return (positions)
}

combine.alphashape = function (ashapelist)
{
  #positions = matrix(ncol = 3)
  #triangles = matrix(ncol = 3)
  # colnames(triangles) = c("tr1", "tr2", "tr2", "on.ch", "attached", "rhoT", "muT", "MuT", "fc:")
  initial = ashapelist[[1]]
  for (a in 2:length(ashapelist))
  {
    ashape = ashapelist[[a]]
    count = nrow(initial$x)
    initial$x <- rbind(initial$x, ashape$x)
    ashape$tetra[,1:4] <- ashape$tetra[,1:4] + count
    initial$tetra <- rbind(initial$tetra, ashape$tetra)
    ashape$triang[,1:3] <- ashape$triang[,1:3] + count
    initial$triang <- rbind(initial$triang, ashape$triang)
    ashape$edge[,1:2] <- ashape$edge[,1:2] + count
    initial$edge <- rbind(initial$edge, ashape$edge)
    ashape$vertex[,1] <- ashape$vertex[,1] + count
    initial$vertex <- rbind(initial$vertex, ashape$vertex)
  }
  return (initial)
}


spinny <- function(object, target){
  rotate = 'go!'
  first = 0
  while (rotate != 'e'){
    if (first == 0){
      v1 = select.points(object)
      pop3d();v2 = select.points(object)
      pop3d();v3 = select.points(object)
    }
    if (rotate == 1){v1 = select.points(object)}
    if (rotate == 2){v2 = select.points(object)}
    if (rotate == 3){v3 = select.points(object)}
    # Check all is okay by plotting plane
    normal <- crossProduct(v2-v1,v3-v1)
    zeroPro <- points2plane(rep(0,3),v1,normal)
    ## get sign of normal displacement from zero
    sig <- sign(crossprod(-zeroPro,normal))
    d <- sig*norm(zeroPro,"2")
    planes3d(normal[1],normal[2],normal[3],d=d, alpha = 0.4)
    m = mirror2plane(object, v1 = c(v1), v2 = c(v2), v3 = c(v3))
    points3d(object, col = 'blue')
    points3d(target, col = 'cyan')
    points3d(m, col = 'grey')
    points3d(rbind(v1,v2,v3), col = c(1,2,3), size = 40)
    first = first + 1
    rotate = readline(prompt="Change V1, V2, V3 or exit (e)?  ")
  }
  return(m)
}

