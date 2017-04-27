## These are flycircuit functions

#' Fetch Flycircuit neuron skeletons from the Tawian Flycircuit server
#'
#' @description Assigns cell body side based in neuron name.
#'
#' @param fcneurons a vectors of valid FlyCircuit neuron ids
#' @param xform_version the version of nat.templatebrains::xform_brain to use
#' @param female whether the neurons to be transformed are female (TRUE) or male (FALSE)
#' @param x some neuron or neuronlist
#' @param ... additional arguments passed to methods
#'
#' @return A neuronlist of FlyCirucit neurons reigstered in the intersex FCWB brain space
#' @export
get.skeleton.from.flycircuit <- function(fcneurons, xform_version=1, ...){
  ids = c()
  fcns = neuronlist()
  for (n in 1:length(fcneurons)){
    swc=sprintf("http://flycircuit.tw/download/swc/%s.swc", fcneurons[n])
    ofcn=tryCatch(nat::read.neuron(swc), error = function(e) NULL)
    if(!is.null(ofcn)) {
      fcns = c(fcns, as.neuronlist(ofcn))
      ids = c(ids, fcneurons[n])
    }
  }
  names(fcns) = ids
  fcns=Chiang2FCWB(fcns)
  fcns = nat::nlapply(fcns,reroot.flycircuit.neuron)
  fcns
}

#' @rdname get.skeleton.from.flycircuit
reroot.flycircuit.neuron <- function(x){
  x =as.neuron(as.ngraph(x), origin = which(x$d$Label==4))
  x$d$Label = 0
  x$d$Label[x$StartPoint] = 1
  x
}

#' @rdname get.skeleton.from.flycircuit
Chiang2FCWB <- function(x, female = grepl("-F-",x), xform_version=1) {
  template_to_use=ifelse(female, "chiangf","chiangm")
  if(xform_version>1) template_to_use=paste0(template_to_use, xform_version)
  nat::nmapply(xform_brain, x, sample=template_to_use, MoreArgs = list(reference=nat.flybrains::FCWB))
}


#' Apply a transform to a neruon/neuronlist
#'
#' @description Apply a transform to a neruon/neuronlist object using Morpho::applyTransform
#'
#' @param x a neuron/neuronlist object
#' @param trafo A valid transformation for Morpho:applyTransform
#' @param inverse whether to calculate the inverse of trafo
#' @param ... additional arguments passed to methods
#'
#' @return Transformed neuron/neuronlist object
#' @export
napplyTransform<-function(x, trafo, inverse = F, ...) UseMethod("napplyTransform")

#' @export
#' @rdname napplyTransform
napplyTransform.neuron <- function(x, trafo, inverse = F,...){
  xyzmatrix(x$d)<-Morpho::applyTransform(xyzmatrix(x$d), trafo = trafo, inverse = inverse)
  if (!is.null(x$connectors)){
    xyzmatrix(x$connectors)<-Morpho::applyTransform(xyzmatrix(x$connectors), trafo = trafo, inverse = inverse)
  }
  x
}

#' @export
#' @rdname napplyTransform
napplyTransform.neuronlist <- function(x, trafo, inverse = F,...){
  nlapply(x, napplyTransform.neuron, trafo, inverse)
}

#' Assign axon/dendrite split to skeletons
#'
#' @description manually cycle through and assign a label to points in a neuron to mark out axonic, dednritic and mixed/uncertain cable
#'
#' @param someneuronlist a neuron/neuronlist object
#' @param resample the unit lengths to which neurons should be resampled
#' @param ... additional arguments passed to methods
#'
#' @return Neuronlist with polarity assignantion marked in the neuron$d dataframe of each neuron object within that neuronlist
#' @export
#' @importFrom grDevices colorRampPalette
assign.cable.polarity <- function(someneuronlist,resample=1,...){
  jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  for (neuron in 350:length(someneuronlist)){
    rgl::clear3d();rgl::plot3d(nat.flybrains::FCWB)
    print(neuron)
    print(names(someneuronlist[neuron]))
    rgl::plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F)
    cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
    points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
    progress = readline(prompt="Change? y/n   ")
    while (progress=="y"){
      message("Select axonic nodes")
      continue = readline(prompt="Select? y/n   ")
      while(continue=="y"){
        s = select.points(xyzmatrix(someneuronlist[neuron][[1]]),plot3d=someneuronlist[neuron][1])
        if (nrow(s)>0){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 2}
        rgl::clear3d();rgl::plot3d(nat.flybrains::FCWB);rgl::plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F);
        cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
        rgl::points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
        continue = readline(prompt="Select again? y/n   ")
        if (nrow(s)>0){if (continue=="y"){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 0}}
      }
      message("Select dendritic nodes")
      continue = readline(prompt="Select? y/n   ")
      while(continue=="y"){
        s = select.points(xyzmatrix(someneuronlist[neuron][[1]]),plot3d=someneuronlist[neuron][1])
        if (nrow(s)>0){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 3}
        rgl::clear3d();rgl::plot3d(nat.templatebrains::FCWB);rgl::plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F);
        cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
        rgl::points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
        continue = readline(prompt="Select again? y/n   ")
        if (nrow(s)>0){if (continue=="y"){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 0}}
      }
      message("Select mixed/uncertain nodes")
      continue = readline(prompt="Select? y/n   ")
      while(continue=="y"){
        s = select.points(xyzmatrix(someneuronlist[neuron][[1]]),plot3d = someneuronlist[neuron][1])
        if (nrow(s)>0){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 8}
        rgl::clear3d();rgl::plot3d(nat.flybrains::FCWB);rgl::plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F);
        cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
        rgl::points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
        continue = readline(prompt="Select again? y/n   ")
        if (nrow(s)>0){if (continue=="y"){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 0}}
      }
      someneuronlist[neuron][[1]]$d[1,]$Label = 1
      print(neuron)
      progress = readline(prompt="Review? y/n   ")
    }
    someneuronlist = nlapply(someneuronlist,resample,stepsize=resample)
    correct.labels.neuron<-function(neuron){
      neuron$d$Label[!neuron$d$Label%in%c(0:10)] = 0
      neuron
    }
    someneuronlist = nlapply(someneuronlist,correct.labels.neuron)
    rgl::clear3d()
  }
  someneuronlist
}


#' Extract axonic/dendritic points from a neuron/neuronlist
#'
#' @description Extract axonic/dendritic points/endpoints from a neuron/neuronlist object
#'
#' @param x a neuron/neuronlist object, which has its axons/dendrites labelled in swc format in its neuron$d dataframes
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname extract.cable
axonic.points<-function(x, ...) UseMethod("axonic.points")
#' @export
#' @rdname extract.cable
dendritic.points<-function(x, ...) UseMethod("dendritic.points")
#' @export
#' @rdname extract.cable
mixed.points<-function(x, ...) UseMethod("mixed.points")
#' @rdname extract.cable
axonic.points.neuron <- function(x){
  points=x$d
  xyzmatrix(points[points$Label%in%c(-2,2),])
}
#' @rdname extract.cable
dendritic.points.neuron <- function(x){
  points=x$d
  xyzmatrix(points[points$Label%in%c(-3,3),])
}
#' @rdname extract.cable
mixed.points.neuron <- function(x){ # Mised also means that I do not know
  points=x$d
  xyzmatrix(points[points$Label%in%c(8),])
}
#' @rdname extract.cable
dendritic.points.neuronlist <- function(x){
  do.call(rbind,nlapply(x,dendritic.points.neuron))
}
#' @rdname extract.cable
axonic.points.neuronlist <- function(x){
  do.call(rbind,nlapply(x,axonic.points.neuron))
}
#' @rdname extract.cable
mixed.points.neuronlist <- function(x){
  do.call(rbind,nlapply(x,mixed.points.neuron))
}
#' @export
#' @rdname extract.cable
axonal.endings <- function(x){
  points=x$d[nat::endpoints(x)[which(endpoints(x)!=rootpoints(x))],]
  xyzmatrix(points[points$Label%in%c(-2,2),])
}
#' @export
#' @rdname extract.cable
dendritic.endings <- function(x){
  points=x$d[nat::endpoints(x)[which(endpoints(x)!=rootpoints(x))],]
  xyzmatrix(points[points$Label%in%c(-3,3),])
}


#' Re-root neurons to their soma
#'
#' @description Cycle through and manually re-root neurons to their soma
#'
#' @param someneuronlist a neuron/neuronlist object
#' @param brain a brain to plot from the nat.templatebrains package e.g. FCWB
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
correctsoma <- function(someneuronlist, brain = NULL,...){
  correctedsomas = neuronlist()
  nopen3d()
  for (n in 1:length(someneuronlist)){
    print(n)
    w = someneuronlist[[n]]
    print(names(someneuronlist[n]))
    if(!is.null(brain)){rgl::plot3d(brain)}
    rgl::plot3d(w, soma = T)
    eps.xyz=w$d[nat::endpoints(w),]
    progress =F
    while(progress == F){
      cat ("Rotate brain and then hit [enter] to continue")
      line <- readline()
      message("Select new root from highlighted endpoints")
      selected.point <- select3d()
      selected.point <- selected.point(nat::xyzmatrix(eps.xyz))
      selected.point <- eps.xyz$PointNo[selected.point]
      if (length(selected.point)>1|length(selected.point)==0){
        message("Multiple end points selected, try again")
      }else{
        corrected =as.neuron(as.ngraph(w), origin = selected.point)
        rgl::plot3d(corrected, soma = T, col = "blue")
        progress = readline(prompt="Good enough? T/F  ")
      }
    }
    rgl::clear3d()
    correctedsomas = c(correctedsomas, as.neuronlist(corrected))
  }
  attr(correctedsomas, "df") = attr(someneuronlist, "df")
  correctedsomas
}

#' Generate a connectivity matrix based on euclidean distance between points
#'
#' @description Generates an 'overlap matrix' of overlap scores between neurons in the 'neurons' and 'targets' pools.
#' For every point in a given neuron in 'neurons', a distance score is calculated to every point in a neuron in 'targets'.
#' The sum of this score is added to the final output matrix. The score is calculated as e(-d^2/2δ^2), where d is the euclidean distance between the two points,
#' and δ is the expected distance in um that is considered 'close'.
#'
#' @param neurons first set of neurons
#' @param targets second set of neurons
#' @param neuropil an as3d object of the neuropil in which to consider connectivity. Defaults to whole brain.
#' @param delta the distance (in um) at which a synapse might occur
#' @param split with a CATMAID neuron, whether or not to split the neuron using flow centrality
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
overlap.connectivity.matrix <- function(neurons,targets,neuropil = NULL,delta =1){
  if(length(neuropil)>0){
    points = rbind(axonic.points.neuronlist(neurons),mixed.points.neuronlist(neurons))
    points = points[inashape3d(points=points,as3d=neuropil),]
    neurons = nlapply(neurons, nat::prune, target = points, keep = 'near', maxdist = 1,OmitFailures=T)
    points = rbind(dendritic.points.neuronlist(targets),mixed.points.neuronlist(targets))
    points = points[inashape3d(points=points,as3d=neuropil),]
    targets = nlapply(targets, nat::prune, target = points, keep = 'near', maxdist = 1,OmitFailures=T)
    correct.d <- function(neuron){
      neuron$d$Label=3
      neuron
    }
    targets = nlapply(targets, correct.d) # Correct label attribution
  }
  score.matrix = matrix(0,nrow = length(neurons),ncol = length(targets))
  rownames(score.matrix) = names(neurons)
  colnames(score.matrix) = names(targets)
  for (n in 1:length(neurons)){
    message("Working on neuron ", n,"/",length(neurons))
    a = axonic.points(neurons[[n]])
    targets.d = nlapply(targets, dendritic.points)
    s = sapply(targets.d, function(x)sum(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)))) # Score similar to that in Schlegel et al. 2015
    score.matrix[n,] = s
  }
  score.matrix
}

#' @rdname overlap.connectivity.matrix
overlap.connectivity.matrix.catmaid <- function(neurons,targets,neuropil = NULL,delta =1, split = TRUE){
  if (split==TRUE){
    message("Calculating flow centrality, splitting neurons")
    neurons.f = flow.centrality(neurons)
    neurons = neurites(neurons.f,fragment="axons")
    targets.f = flow.centrality(targets)
    message("Calculating flow centrality, splitting targets")
    targets = neurites(targets.f,fragment="dendrites")
  }
  if(length(neuropil)>0){
    points = xyzmatrix(neurons)[inashape3d(points=xyzmatrix(neurons),as3d=neuropil),]
    neurons = nlapply(neurons, nat::prune, target = points, keep = 'near', maxdist = 1,OmitFailures=T)
    points = xyzmatrix(neurons)[inashape3d(points=xyzmatrix(neurons),as3d=neuropil),]
    targets = nlapply(targets, nat::prune, target = points, keep = 'near', maxdist = 1,OmitFailures=T)
  }
  score.matrix = matrix(0,nrow = length(neurons),ncol = length(targets))
  rownames(score.matrix) = names(neurons)
  colnames(score.matrix) = names(targets)
  for (n in 1:length(neurons)){
    message("Working on neuron ", n,"/",length(neurons))
    a = xyzmatrix(neurons[[n]])
    targets.d = nlapply(targets, xyzmatrix)
    s = sapply(targets.d, function(x)sum(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)))) # Score similar to that in Schlegel et al. 2015
    score.matrix[n,] = s
  }
  score.matrix
}


polaritycluster <- function(someneuronlist, sigma = 1, omega = 1, symmetric = T){
  m = matrix(nrow = length(someneuronlist), ncol = length(someneuronlist))
  colnames(m) = rownames(m) = names(someneuronlist)
  for (neuron in 1:length(someneuronlist)){
    g = as.data.frame(axon.points(someneuronlist[neuron][[1]]))
    if(plyr::empty(g)){ m[neuron,] = 0}
    for (neuron2 in 1:length(someneuronlist)){
      if(plyr::empty(g)){ m[neuron,] = break}
      t = as.data.frame(axon.points(someneuronlist[neuron2][[1]]))
      scores = c()
      for (syn in 1:nrow(g)){
        if(plyr::empty(t)){scores = c(scores, 0); break}
        n = nabor::knn(t, nat::xyzmatrix(g[syn,]), k =1)
        close = t[n$nn.idx,]
        gn = nabor::knn(nat::xyzmatrix(g[syn,]), g, k =1)
        gn = sum(gn$nn.dists<omega) - 1
        tn = nabor::knn(nat::xyzmatrix(close), t, k =1)
        tn = sum(tn$nn.dists<omega) - 1
        multiplier = exp(-abs(gn-tn)/(tn+gn))
        if(is.infinite(multiplier)|is.na(multiplier)|is.nan(multiplier)){ multiplier = 0}
        score = exp((-n$nn.dist^2)/2*(sigma)^2)*multiplier
        scores = c(scores, score)
      }
      m[neuron,neuron2] = sum(scores, na.rm = T)
    }
  }
  if (symmetric == T){
    #pmean <- function(x,y) (x+y)/2 take the mean, or take the lowest score as below.
    m[] <- pmin(m, matrix(m, nrow(m), byrow=TRUE))
  }
  return(m)
}

#' Calculate the cable length inside of different neuropils in a segmented .surf brain
#'
#' @description Calculates the cable length supplied by neurons 'x' to different brain regions defined in 'brain'.
#'
#' @param x a set of neurons
#' @param brain the .surf brainspace in which the neurons are registered, must be segmented into neuropils
#' @param method whether to use the neurons' axons or dendrites, or both
#' @param stepsize the unit to which neurons should be resampled for counting
#' @param min.endpoints the minimum number of endpoints a neuron must have in a neuropil to be counted as included in it
#' @param alpha the alpha given to the ashape3d() function to generate neuropil objects by which to calculate point inclusion
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
cable.inside.neuropils<-function(x, brain = nat.flybrains::FCWBNP.surf, method = c("neurites","axons","dendrites"), stepsize = 0.1, min.endpoints = 2,alpha=30, ...) UseMethod("cable.inside.neuropils")

#' @rdname cable.inside.neuropils
cable.inside.neuropils.neuron <- function(x, brain = nat.flybrains::FCWBNP.surf, method = c("neurites","axons","dendrites"), stepsize = 0.1, min.endpoints = 2,alpha=30){
  if(!requireNamespace('nat.flybrains', quietly = TRUE))
    stop("You must install suggested package nat.flybrains to use this function!")
  neuron = resample(x, stepsize = stepsize)
  targets = c(0,2,3,8)
  if (method=="axons"){targets = c(2)}
  if (method=="dendrites"){targets = c(3)}
  endings <- function(neuron){
    points=neuron$d[nat::endpoints(neuron)[which(endpoints(neuron)!=rootpoints(neuron))],]
    xyzmatrix(points[points$Label%in%targets,])
  }
  in.neuropil <- function(neuron,neuropil,stepsize =0.1, min.endpoints =2,alpha=30){
    neuropil = alphashape3d::ashape3d(xyzmatrix(neuropil),alpha=alpha)
    if(sum(alphashape3d::inashape3d(points=endings(neuron),as3d=neuropil))>min.endpoints){
      points = neuron$d[neuron$d$Label%in%targets,]
      sum(alphashape3d::inashape3d(points=nat::xyzmatrix(points),as3d=neuropil))
    }else{0}
  }
  sapply(brain$RegionList, function(x) in.neuropil(neuron,neuropil=subset(brain, x),stepsize =stepsize, min.endpoints=min.endpoints,alpha=alpha))*stepsize
}

#' @rdname cable.inside.neuropils
cable.inside.neuropils.neuronlist <- function(x, brain = nat.flybrains::FCWBNP.surf, method = c("neurites","axons","dendrites"), stepsize = 0.1,min.endpoints = 2,alpha=30){
  nlapply(x, cable.inside.neuropils.neuron, brain=brain, method=method,stepsize=stepsize)
}


#' Write ordered swc
#'
#' @description  Write separate, ordered swc files for the MATLAB SPIn package.
#' One can also write separate files for a neuron's axon, dendrite and other cable
#'
#' @param neuron a neuron object
#' @param file file path to which to write output
#' @param ... additional arguments passed to methods
#'
#' @export
#' @rdname write.spin.swc
#' @importFrom stats complete.cases
write.spin.split.swc<-function(neuron,file){
  dend = neuron$d[neuron$d$Label==3,]
  axon = neuron$d[neuron$d$Label==2,]
  rest = neuron$d[!neuron$d$Label%in%c(2,3),]
  dendrites = prune_vertices(neuron,verticestoprune = c(axon$PointNo,rest$PointNo))
  dendrites.parents = match(dendrites$d$Parent,dendrites$d$PointNo)
  dendrites.parents[is.na(dendrites.parents)] = -1
  dendrites$d$PointNo = 1:nrow(dendrites$d)
  dendrites$d$Parent = dendrites.parents
  axons = prune_vertices(neuron,c(dend$PointNo,rest$PointNo))
  axons.parents = match(axons$d$Parent,axons$d$PointNo)
  axons.parents[is.na(axons.parents)] = -1
  axons$d$PointNo = 1:nrow(axons$d)
  axons$d$Parent = axons.parents
  i = 1
  if (axons$nTrees>1){
    s = axons$SubTrees
  }else{
    s = list(axons$SegList)
  }
  for (n in 1:length(s)){
    if (length(unique(unlist(s[n])))<2){
      neuron$d$Label[neuron$d$PointNo%in%unique(unlist(s[n]))] = 0
    }else{
      axon = axons
      axon$d$PointNo = match(axon$d$PointNo,unique(unlist(s[n])))
      axon$d =axon$d[complete.cases(axon$d),]
      axon$d$Parent = ifelse(!is.na(match(axon$d$Parent,unique(unlist(s[n])))),match(axon$d$Parent,unique(unlist(s[n]))),-1)
      axon$d = axon$d[order(axon$d$PointNo),]
      write.neuron(axon,file=gsub(".swc",paste("_",i,".swc",sep=""),file),Force=T)
      i = i+1
    }
  }
  if (dendrites$nTrees>1){
    s = dendrites$SubTrees
  }else{
    s = list(dendrites$SegList)
  }
  for (n in 1:length(s)){
    if (length(unique(unlist(s[n])))<2){
      neuron$d$Label[neuron$d$PointNo%in%unique(unlist(s[n]))] = 0
    }else{
      dendrite = dendrites
      dendrite$d$PointNo = match(dendrite$d$PointNo,unique(unlist(s[n])))
      dendrite$d =dendrite$d[complete.cases(dendrite$d),]
      dendrite$d$Parent = ifelse(!is.na(match(dendrite$d$Parent,unique(unlist(s[n])))),match(dendrite$d$Parent,unique(unlist(s[n]))),-1)
      dendrite$d = dendrite$d[order(dendrite$d$PointNo),]
      write.neuron(dendrite,file=gsub(".swc",paste("_",i,".swc",sep=""),file),Force=T)
      i = i+1
    }
  }
  s = unique(unlist(as.seglist(neuron)))
  if (neuron$nTrees>1){
    s = unlist(neuron$SubTrees)
  }
  neuron$d$PointNo = match(neuron$d$PointNo,s) # Re-name points according to SegList
  neuron$d$Parent = ifelse(!is.na(match(neuron$d$Parent,s)),match(neuron$d$Parent,s),-1) # Given the above, correct Parent attribution
  neuron$d = neuron$d[order(neuron$d$PointNo),] # Order the swc matrix
  write.neuron(neuron,file = file,format="swc",Force=T) # Write file
}

#' @export
#' @rdname write.spin.swc
write.spin.swc <- function(neuron, file){
  s = unique(unlist(as.seglist(neuron)))
  if (neuron$nTrees>1){
    s = unlist(neuron$SubTrees)
  }
  neuron$d$PointNo = match(neuron$d$PointNo,s) # Re-name points according to SegList
  neuron$d$Parent = ifelse(!is.na(match(neuron$d$Parent,s)),match(neuron$d$Parent,s),-1) # Given the above, correct Parent attribution
  neuron$d = neuron$d[order(neuron$d$PointNo),] # Order the swc matrix
  neuron$d$Label = 0 # Labels are all wrong when downloaded from Flycircuit, make them 0 for 'unknwon'
  neuron$d$Label[1] = 1 # Add soma label
  write.neuron(neuron,file = file,format="swc",Force=T) # Write file
}



#' Plot a neuron with different coloured synapses depending on partners
#'
#' @description  Plot a neuron with different coloured synapses depending on partners.
#' Points indicate inputs synapses, and asterixes indicate output synapses.
#'
#' @param neuron a neuron to plot
#' @param skids the skeleton ids of neurons with which synapses are to be shown. Defaults to all partners.
#' @param col the colour of the neuron skeleton
#' @param inputs whether or not to show input synapses are to be shown
#' @param outputs whether or not output synapses are to be shown
#' @param printout whether or not to plot a basic legend to indicate what colours mean which neuron
#' @param ... additional arguments passed to methods
#'
#' @export
#' @importFrom graphics legend plot plot.new
synapsecolours.neuron <-function(neuron, skids = NULL, col = "black", inputs = T, outputs = T,printout=F){
  if(is.neuronlist(neuron)){neuron = neuron[[1]]}
  rgl::plot3d(neuron,WithNodes=F,soma=T,col=col,lwd=2)
  graphics::plot.new()
  if (inputs)
    inputs = neuron$connectors[neuron$connectors$prepost==1,]
  c = subset(catmaid_get_connectors(inputs$connector_id),post==neuron$skid)[,-3]
  inputs = merge(inputs,c, all.x=F,all.y=F)
  if(!is.null(skids)){inputs = subset(inputs,pre%in%skids)}
  colours = data.frame(pre =unique(inputs$pre),col =rainbow(length(unique(inputs$pre))))
  inputs = merge(inputs,colours, all.x=T,all.y=F)
  if (nrow(inputs)>0){
    points3d(nat::xyzmatrix(inputs),col = inputs$col)
    if(printout)
      legend("left",legend=catmaid_get_neuronnames(colours$pre),fill=colours$col,cex=2/nrow(colours))
  }
  if (outputs)
    outputs = neuron$connectors[neuron$connectors$prepost==0,]
  c = subset(catmaid_get_connectors(outputs$connector_id),pre==neuron$skid)[,-2]
  outputs = merge(outputs,c, all.x=F,all.y=F)
  if(!is.null(skids)){outputs = subset(outputs,post%in%skids)}
  colours = data.frame(post =unique(outputs$post),col =rainbow(length(unique(outputs$post))))
  outputs = merge(outputs,colours, all.x=T,all.y=F)
  if (nrow(outputs)>0){
    text3d(nat::xyzmatrix(outputs),text = "*",col = outputs$col,cex=2)
    if(printout)
      legend("right",legend=catmaid_get_neuronnames(colours$post),col=colours$col,cex=2/nrow(colours))
  }
}


#' Average non-branching tracts
#'
#' @description  Function returns an 'average tract' when given a neuronlist of non-branching cable.
#' The function uses NBlast to find the the neuron within the group that has the highest average similarity score with the rest of
#' the given cable
#'
#' @param cable somneuronlist object
#' @param sigma smoothing parameter given to nat::smooth_neuron, i.e. the standard deviation of the Gaussian smoothing kernel (which has the same spatial units as the object being smoothed)
#' @param stepsize the spacing to which cable should be resampled
#' @param mode Either 1 or 2. Method 1 simply takes the best NBlast match and, for each of its nodes, finds the mean xyz coordinate for that node index in the cable group (soma has index 1), all resampled to the stepsize value.
#' Method 2 is similar, but takes the closest node in each of the other neurons in the cable neuronlist rather than the nodes of the same index.
#' @param ... additional arguments passed to methods
#'
#' @return a single neuron object, the averaged tract, and containing values for standard deviation for each 3D point.
#' @export
#' @importFrom nabor knn
average.tracts <- function(cable, sigma = 6, mode = c(1,2),stepsize = 1,...){
  if (length(cable)>1){
    colSD <- function(m){apply(m, 2, sd)}
    roots = t(sapply(cable,function(x)nat::xyzmatrix(x)[nat::rootpoints(x),]))
    root.sd = apply(roots, 2, sd)
    root = colMeans(roots)
    cable = nat::nlapply(cable,resample,stepsize =stepsize)
    c.dps = nat::dotprops(cable, OmitFailures = T)
    csmat=nat.nblast::nblast_allbyall(c.dps)
    best = cable[names(which.max(colMeans(csmat)))]
    best.points = nat::xyzmatrix(best)
    points = sapply(cable, function(x) nat::xyzmatrix(x))
    end.points = t(sapply(points, function(p) p[nrow(p),]))
    points.averaged = matrix(root,ncol = 3)
    points.sd = matrix(root.sd, ncol = 3)
    points.used = matrix(ncol=3)
    if (mode==1){
      for(b in 2:nrow(best.points)){
        points.keep = colMeans(do.call(rbind,lapply(points, function(x) tryCatch(x[b,], error = function(e) NULL))))
        sd.keep = colSD(do.call(rbind,lapply(points, function(x) tryCatch(x[b,], error = function(e) NULL))))
        points.averaged = rbind(points.averaged,points.keep)
        points.sd = rbind(points.sd,sd.keep)
      }
    }
    if (mode==2){
      for(b in 2:nrow(best.points)){
        selected.points = do.call(rbind,lapply(points, function(x) x[nabor::knn(data = x,query =matrix(best.points[b,],ncol=3),k=1)$nn.idx[1],]))
        reused = matrix(selected.points[apply(selected.points,1,function(x) sum(x%in%points.used)>=3),],ncol=3)
        ends = apply(reused,1,function(x) sum(x%in%end.points)>=3)
        remove = names(ends)[ends]
        selected.points = selected.points[!rownames(selected.points)%in%ends,]
        points = points[!names(points)%in%remove]
        points.keep = colMeans(selected.points)
        sd.keep = colSD(selected.points)
        points.averaged = rbind(points.averaged,points.keep)
        points.sd = rbind(points.sd,sd.keep)
        points.used = rbind(points.used, selected.points)
      }
    }
    rownames(points.sd) = rownames(points.averaged) = NULL
    best[[1]]$d[c("X","Y","Z")] = points.averaged
    colnames(points.sd) = c("sdX","sdY","sdZ")
    best[[1]]$d = cbind(best[[1]]$d,points.sd,n=length(cable))
    as.neuron(nat::smooth_neuron(best[[1]],sigma=sigma))
  }else {cable[[1]]}
}

#' Assign identitiy to lateral horn neurons
#'
#' @description  Given a neuron, this function assigns it to a lateral horn primary neurite tract, anatomy group and cell type as described by Fretcher et al. 2017.
#' Note should be taken with what side of the brain the sample group is on. Note that the function will force an assignation even if the neurons in someneuronlist belong to
#' unidentiified tracts/anatomy group/cell types. Relies on NBlast. Skeletons must have their somas as their root. This fucntion provides a good guess, however if it gets a neurons assignment
#' incorrect, that does not necessarily mean the neuron in question is not covered by the naming system.
#'
#' @param someneuronlist a neuronlist object
#' @param most.lhns a dataset containing example neurons of different primary neurite tracts, anatomy groups and cell types. Defaults to the Flycircuit neurons and dye-fills used in Fretcher et al. 2017, as do the two arguments below
#' @param most.lhns.dps a dotprops object of the above
#' @param most.lhns.pnts.dps the primary neurite tracts to search against as a dotrpops object
#' @param brain the brainspace of the someneuronlist. If left NULL, assumes the space is FCWB
#' @param someassignedneuronlist a neuronlist that has been passed through assign_lh_neuron
#' @param ... additional arguments passed to methods
#'
#' @return a neuronlist with the best guess for primary neurite tract, anatomy group and cell type listed in its metadata. Quality of this estimation depends on quality of the skeleton objects used in this function and their registration ot FCWB space
#' @export
#' @importFrom nat.templatebrains xform_brain
assign_lh_neuron <- function(someneuronlist, most.lhns = catnat::most.lhns, most.lhns.dps = catnat::most.lhns.dps, most.lhns.pnts.dps = catnat::most.lhns.pnts.dps, brain = NULL, strahler=1){
  most.lhns = subset(most.lhns, pnt!="notLHproper")
  if (!is.null(brain)){ most.lhns = xform_brain(most.lhns, sample = nat.flybrains::FCWB, reference = brain)}
  message("Generating primary neurites across the LHNs")
  someneuronlist.pnts = suppressWarnings(primary.neurite.neuronlist(someneuronlist,resample=1))
  message("Generating dotprops objects")
  most.lhns.dps = subset(most.lhns.dps, pnt!="notLHproper")
  most.lhns.pnts.dps = most.lhns.pnts.dps[most.lhns.pnts.dps[,"good.trace"]==T]
  someneuronlist.dps = rescue.dps(someneuronlist, resample = 1,k=5)
  someneuronlist.pnts.dps = rescue.dps(someneuronlist.pnts, resample = 1, k=5)
  # Now try to find the tract a neuron fits into
  message("Assigning primary neurites")
  if(length(someneuronlist.pnts.dps)!=length(someneuronlist)){warning("Neurons dropped!")}
  if (length(someneuronlist)==1){
    results1 = as.matrix(nat.nblast::nblast(query = someneuronlist.pnts.dps, target = most.lhns.pnts.dps, UseAlpha = T))
    results2 = as.matrix(nat.nblast::nblast(target = someneuronlist.pnts.dps, query = most.lhns.pnts.dps, UseAlpha = T))
    results = (results1+results2)/2
  }else{
    results1 = nat.nblast::nblast(query = someneuronlist.pnts.dps, target = most.lhns.pnts.dps, UseAlpha = T,.parallel=TRUE)
    results2 = nat.nblast::nblast(target = someneuronlist.pnts.dps, query = most.lhns.pnts.dps, UseAlpha = T,.parallel=TRUE)
    results = (results1+t(results2))/2
  }
  pct = most.lhns.pnts.dps[,"pnt"][apply(results,2,which.max)]
  attr(someneuronlist, "df") = cbind(attr(someneuronlist, "df"), pnt = pct)
  message("Assigning anatomy group and cell types")
  anatomy.group = cell.type = c()
  for (n in 1:length(someneuronlist)){
    cluster.dps = subset(most.lhns.dps, pnt==pct[n])
    results1 = as.matrix(nat.nblast::nblast(query = someneuronlist.dps[n], target = cluster.dps, UseAlpha = T,.parallel=TRUE))
    results2 = as.matrix(nat.nblast::nblast(target = someneuronlist.dps[n], query = cluster.dps, UseAlpha = T,.parallel=TRUE))
    results = (results1+results2)/2
    anatomy.group = c(anatomy.group,cluster.dps[,"anatomy.group"][apply(results,2,which.max)])
    cell.type = c(cell.type,cluster.dps[,"cell.type"][apply(results,2,which.max)])
  }
  attr(someneuronlist, "df") = cbind(attr(someneuronlist, "df"), anatomy.group = anatomy.group,cell.type = cell.type)
  return(someneuronlist)
}

#' @export
#' @rdname assign_lh_neuron
scan_lh_matches <-function(someassignedneuronlist){
  load(system.file("data/most.lhns.rda", package = 'catnat'))
  rgl::open3d(userMatrix=structure(c(0.999696135520935, -0.0139423999935389, 0.0203206837177277,
                              0, -0.0140735022723675, -0.999880969524384, 0.00632292777299881,
                              0, 0.0202301144599915, -0.00660702586174011, -0.999773621559143,
                              0, 0, 0, 0, 1), .Dim = c(4L, 4L)),zoom=0.746215641498566)
  message("Cell type: Green, Other Anatomy Group members: Red")
  for(n in 1:length(em.lhns.assigned)){
    rgl::plot3d(nat.flybrains::FCWB)
    rgl::plot3d(em.lhns.assigned[n],soma=T,col="black",lwd=3)
    cluster = as.character(em.lhns.assigned[,"anatomy.group"][n])
    type = as.character(em.lhns.assigned[,"cell.type"][n])
    pt = as.character(em.lhns.assigned[,"pnt"][n])
    rgl::plot3d(subset(most.lhns,cell.type==type),col="green",soma=T,lwd=2)
    rgl::plot3d(subset(most.lhns,anatomy.group==cluster),col="red",soma=T,lwd=1)
    #rgl::plot3d(subset(most.lhns,pnt==pt),col="grey",soma=T)
    message("Press ENTER to continue")
    progress = readline()
    rgl::clear3d()
  }
}

#' Generate a dps object without dropping small neurons
#'
#' @description  Generates a dotprops object, changing the resample size in order to make sure all neurons are kept
#'
#' @param someneuronlist a neuronlist object
#' @param resample.unit the desired resampling. Smaller neurons will have a lower resampling interval than this.
#' @param ... additional arguments passed to methods
#'
#' @return a dotprops object
#' @export
rescue.dps <- function(someneuronlist,resample.unit=1,...){
  someneuronlist.dps = nat::dotprops(someneuronlist, resample = resample.unit, OmitFailures = T)
  no.points = summary(someneuronlist)$cable.length
  tooshort = someneuronlist[!names(someneuronlist)%in%names(someneuronlist.dps)]
  tooshort.dps = nat::nlapply(tooshort, function(x) nat::dotprops(x, resample = resample.unit/min(no.points), OmitFailures = F))
  someneuronlist.dps = c(someneuronlist.dps,tooshort.dps)[names(someneuronlist)]
  someneuronlist.dps
}


#' NBlast two different sets of neurons forwards and backwards
#'
#' @description  NBlast two different sets of neurons forwards and backwards, and take the mean of these results
#'
#' @param group1 a neuronlist object
#' @param group2 a neuronlist object, defaults to group1
#' @param sd Standard deviation to use in distance dependence of nblast v1 algorithm. Ignored when version=2
#' @param version the version of the algorithm to use (the default, 2, is the latest)
#' @param UseAlpha whether to consider local directions in the similarity calculation (default: FALSE)
#' @param normalised whether to divide scores by the self-match score of the query
#' @param OmitFailures Whether to omit neurons for which FUN gives an error. The default value (NA) will result in nblast stopping with an error message the moment there is an eror
#' @param smat the score matrix to be used. See ?nlast.
#' @param ... additional arguments passed to methods
#'
#' @return a dotprops object
#' @export
nblast_bothways<-function(group1,group2=group1,smat = NULL,
                          sd = 3, version = c(2, 1), normalised = FALSE, UseAlpha = FALSE,
                          OmitFailures = NA){
  nblast.forward = nat.nblast::nblast(query=group1,target=group2,UseAlpha=UseAlpha,normalised=normalised)
  nblast.backward = nat.nblast::nblast(query=group2,target=group1,UseAlpha=UseAlpha,normalised=normalised)
  (nblast.forward+t(nblast.backward))/2
}

