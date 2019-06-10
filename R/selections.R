# Functions for manipulating neuronlist objects
get.connectors=function (someneuronlist, target = c("BOTH", "PRE", "POST")){
  dot.points = c()
  for (neuron in 1:length(someneuronlist)){
    conn = someneuronlist[[neuron]]$connectors
    if (target == "POST") { conn = subset(conn, conn$prepost == 0) }
    if (target == "PRE") { conn = subset(conn, conn$prepost == 1) }
    if (length(conn) > 0 && nrow(conn) > 0){
      p = xyzmatrix(conn)
      dot.points = rbind(dot.points, p)
    }
  }
  return (dot.points)
}


#' Select  neurons in space
#'
#' @description Select and deselect neurons in space interactively. Black = selected, red = deselected
#'
#' @param someneuronlist a neuronlist object
#' @param ... additional arguments passed to methods, i.e. plot3d()
#'
#' @return Selected neurons
#' @export
select_neurons <- function (someneuronlist, ...){
  thechosen = someneuronlist
  progress = 'y'
  rgl::open3d(); rgl::plot3d(thechosen)
  progress = readline(prompt="Add (a) or remove (r) neurons, or exit (e)?  ")
  while (progress != "e"){
    if (progress == "a"){
      keeps = find.neuron(rval = "neuronlist", db = someneuronlist)
      if (length(keeps) > 0){ if (length(thechosen)>0){
        c = duplicated(c(names(thechosen), names(keeps)))
        thechosen = c(thechosen, subset(keeps, !as.data.frame(keeps)$skid%in%c(names(thechosen), names(keeps))[c]))
      }
       else { thechosen = keeps }
      }
      clear3d()
      if (length(thechosen) > 0) {plot3d(thechosen)}
      plot3d(someneuronlist, col = 'grey')    }
    if (progress == "r"){
      removals = c()
      remove <- find.neuron(rval = "names", db = someneuronlist)
      for (neuron in 1:length(thechosen)){
        for (name in remove){
          if (name == names(thechosen[neuron])){
            removals = c(removals, neuron)
          }
       }
     }
      if (length(removals > 0)){ thechosen = thechosen[-removals] }
      clear3d()
      if (length(thechosen) > 0) {plot3d(thechosen)}
      plot3d(someneuronlist, col = 'grey')
    }
    progress = readline(prompt="Add (a) or remove (r) neurons, or exit (e)?  ")
  }
  return (thechosen)
}

#' @export
#' @rdname select_neurons
deselect_neurons =function (someneuronlist,...){
  select_neurons(someneuronlist,...)
}

neurons.inside.alphashape <- function(alpha, db, synapse = "BOTH", degree = NULL){
  selection = c()
  for (neuron in 1:length(db)){
    neuron = db[neuron]
    xyz = get.connectors(neuron, target = synapse)
    no = nrow(xyz)
    if (!is.null(no)&&no > 0){
      p = alphashape3d::inashape3d(alpha, indexAlpha = 1, xyz)
      inside = sum(p)
      if(!is.null(degree)){
        if (inside/no > degree){selection = c(selection, TRUE)
        }else{selection = c(selection, FALSE)}
      }else{
        if (inside > 0){selection = c(selection, TRUE)
        }else{selection = c(selection, FALSE)}
      }
    }else{
      selection = c(selection, FALSE)
    }
  }
  return(db[selection])
}

connectors.inside <- function(skids, alpha, direction = "BOTH", degree = NULL){
  selection = c()
  for (neuron in 1:length(db)){
    neuron = db[neuron]
    xyz = get.connectors(neuron, target = synapse)
    no = nrow(xyz)
    if (!is.null(no)&&no > 0){
      p = alphashape3d::inashape3d(alpha, indexAlpha = 1, xyz)
      inside = sum(p)
      if(!is.null(degree)){
        if (inside/no > degree){selection = c(selection, TRUE)
        }else{selection = c(selection, FALSE)}
      }else{
        if (inside > 0){selection = c(selection, TRUE)
        }else{selection = c(selection, FALSE)}
      }
    }else{
      selection = c(selection, FALSE)
    }
  }
  db[selection]
}

## Moved to nat ###
#' Select 3D points in space
#'
#' @description Select and deselect 3D points in space interactively. Black = selected, red = deselected
#'
#' @param points a mxn matrix
#' @param plot3d Object to pass to the rgl::plot3d function
#'
#' @return A mxn matrix
#' @export
select_points <- function (points, plot3d = NULL){
  plot3d(plot3d,col="grey")
  points = nat::xyzmatrix(points)
  selected.points = unique(points)
  points3d(selected.points)
  progress = readline(prompt="Add (a) or remove (r) points, or exit (e)?  ")
  while (progress != "e"){
    if (progress == "a"){
      keeps = rgl::select3d()
      keep.points <- keeps(unique(points))
      keep.points = subset(unique(points), keep.points)
      selected.points = rbind(selected.points, keep.points)
      rgl::clear3d();rgl::plot3d(plot3d); rgl::points3d(selected.points); rgl::points3d(unique(points), col = 'red')
    }
    if (progress == "r"){
      remove.points <- select3d()
      removed.points <- remove.points(selected.points)
      selected.points = subset(selected.points, !removed.points)
    }
    clear3d();plot3d(plot3d,col="grey")
    if (length(selected.points) > 0) {points3d(selected.points)}
    points3d(unique(points), col = 'red')
    progress = readline(prompt="Add (a) or remove (r) neurons, or exit (e)?  ")
  }
  return (selected.points)
}
