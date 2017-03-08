## These are flycircuit functions

#' Fetch Flycircuit neuron skeletons from the Tawian Flycircuit server
#'
#' @description Assigns cell body side based in neuron name.
#'
#' @param fcneurons a vectors of valid FlyCircuit neuron ids
#' @param ... additional arguments passed to methods
#'
#' @return A neuronlist of FlyCirucit neurons reigstered in the intersex FCWB brain space
#' @export
#' @rdname get.skeleton.from.flycircuit
get.skeleton.from.flycircuit <- function(fcneurons, ...){
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
  fcns = c(Chiang2FCWB(fcns[grepl("-F-",names(fcns))],sex="F"), Chiang2FCWB(fcns[grepl("-M-",names(fcns))],sex="M"))
  names(fcns) = ids
  fcns = nat::nlapply(fcns,reroot.flycircuit.neuron)
  fcns
}

#' @rdname get.skeleton.from.flycircuit
reroot.flycircuit.neuron <- function(neuron){
  neuron =as.neuron(as.ngraph(neuron), origin = which(neuron$d$Label==4))
  neuron$d$Label = 0
  neuron$d$Label[neuron$StartPoint] = 1
  neuron
}

#' @rdname get.skeleton.from.flycircuit
Chiang2FCWB <- function(x, sex = 'F'){
  if (sex == "M"){
    system.file("extdata/CMTKreg/", package = 'catnat')
    affinetransform.m = readRDS(system.file("extdata/CMTKreg/InitialAffine/initialiseCMTKreg_ChiangMaleTowardsFCWB.rds", package = 'catnat'))
    affinetransform.m2 = readRDS(system.file("extdata/CMTKreg/InitialAffine/finalaffine_ChiangMaleTowardsFCWB.rds", package = 'catnat'))
    x = napplyTransform.neuronlist(x, affinetransform.m)
    x = nat::xform(x,reg=system.file("extdata/CMTKreg/Registration/warp/FCWB_typicalbrainmale_01_warp_m0g80c8e1e-1x26r4.list/", package = 'catnat'))
    x = napplyTransform.neuronlist(x, affinetransform.m2)
  }
  if (sex == "F"){
    affinetransform.f = readRDS(system.file("extdata/CMTKreg/InitialAffine/initialiseCMTKreg_ChiangMaleTowardsFCWB.rds", package = 'catnat'))
    affinetransform.f2 = readRDS(system.file("extdata/CMTKreg/InitialAffine/finalaffine_ChiangFemaleTowardsFCWB.rds", package = 'catnat'))
    x = napplyTransform.neuronlist(x, affinetransform.f)
    x = nat::xform(x,reg=system.file("extdata/CMTKreg/Registration/warp/FCWB_typicalbrainfemale_01_warp_m0g80c8e1e-1x26r4.list/", package = 'catnat'))
    x = napplyTransform.neuronlist(x, affinetransform.f2)
  }
  x
}

#' Apply a transform to a neruon/neuronlist
#'
#' @description Apply a transform to a neruon/neuronlist object using Morpho::applyTransform
#'
#' @param someneuronlist a neuron/neuronlist object
#' @param trafo A valid transformation for Morpho:applyTransform
#' @param inverse whether to calculate the inverse of trafo
#' @param ... additional arguments passed to methods
#'
#' @return Transformed neuron/neuronlist object
#' @export
#' @rdname napplyTransform
napplyTransform<-function(someneuronlist, trafo, inverse = F, ...) UseMethod("napplyTransform")

#' @export
#' @rdname primary.neurite
napplyTransform.neuron <- function(neuron, trafo, inverse = F,...){
  xyzmatrix(neuron$d)<-Morpho::applyTransform(xyzmatrix(neuron$d), trafo = trafo, inverse = inverse)
  if (!is.null(neuron$connectors)){
    xyzmatrix(neuron$connectors)<-Morpho::applyTransform(xyzmatrix(neuron$connectors), trafo = trafo, inverse = inverse)
  }
  neuron
}

#' @export
#' @rdname napplyTransform
napplyTransform.neuronlist <- function(someneuronlist, trafo, inverse = F,...){
  nlapply(someneuronlist, napplyTransform.neuron, trafo, inverse)
}

#' Assign axon/dendrite split to skeletons
#'
#' @description manually cycle through and assign a label to points in a neuron to mark out axonic, dednritic and mixed/uncertain cable
#'
#' @param someneuronlist a neuron/neuronlist object
#' @param ... additional arguments passed to methods
#'
#' @return Neuronlist with polarity assignantion marked in the neuron$d dataframe of each neuron object within that neuronlist
#' @export
#' @rdname assign.cable.polarity
assign.cable.polarity <- function(someneuronlist,...){
  jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  for (neuron in 1:length(someneuronlist)){
    clear3d();plot3d(FCWB)
    print(neuron)
    print(names(someneuronlist[neuron]))
    plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F)
    cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
    points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
    progress = readline(prompt="Change? y/n   ")
    while (progress=="y"){
      message("Select axonic nodes")
      continue = readline(prompt="Select? y/n   ")
      while(continue=="y"){
        s = select.points(xyzmatrix(someneuronlist[neuron][[1]]))
        if (nrow(s)>0){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 2}
        clear3d();plot3d(FCWB);plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F);
        cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
        points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
        continue = readline(prompt="Select again? y/n   ")
        if (nrow(s)>0){if (continue=="y"){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 0}}
      }
      message("Select dendritic nodes")
      continue = readline(prompt="Select? y/n   ")
      while(continue=="y"){
        s = select.points(xyzmatrix(someneuronlist[neuron][[1]]))
        if (nrow(s)>0){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 3}
        clear3d();plot3d(FCWB);plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F);
        cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
        points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
        continue = readline(prompt="Select again? y/n   ")
        if (nrow(s)>0){if (continue=="y"){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 0}}
      }
      message("Select mixed/uncertain nodes")
      continue = readline(prompt="Select? y/n   ")
      while(continue=="y"){
        s = select.points(xyzmatrix(someneuronlist[neuron][[1]]))
        if (nrow(s)>0){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 8}
        clear3d();plot3d(FCWB);plot3d(someneuronlist[neuron][[1]],col="black",soma=T,WithNodes=F);
        cols = ifelse(someneuronlist[neuron][[1]]$d$Label<0,(someneuronlist[neuron][[1]]$d$Label*-1)+2,someneuronlist[neuron][[1]]$d$Label)+1
        points3d(xyzmatrix(someneuronlist[neuron][[1]]),col = jet.colors(7)[cols])
        continue = readline(prompt="Select again? y/n   ")
        if (nrow(s)>0){if (continue=="y"){someneuronlist[neuron][[1]]$d[xyzmatrix(someneuronlist[neuron][[1]])%in%s[,1],]$Label = 0}}
      }
      someneuronlist[neuron][[1]]$d[1,]$Label = 1
      print(neuron)
      progress = readline(prompt="Review? y/n   ")
    }
    clear3d()
  }
  someneuronlist
}


#' Extract axonic/dendritic points from a neuron/neuronlist
#'
#' @description Extract axonic/dendritic points/endpoints from a neuron/neuronlist object
#'
#' @param someneuronlist a neuron/neuronlist object, which has its axons/dendrites labelled in swc format in its neuron$d dataframes
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname extract.cable
axonic.points<-function(someneuronlist, ...) UseMethod("axonic.points")
#' @export
#' @rdname extract.cable
dendritic.points<-function(someneuronlist, ...) UseMethod("dendritic.points")
#' @export
#' @rdname extract.cable
mixed.points<-function(someneuronlist, ...) UseMethod("mixec.points")
#' @rdname extract.cable
axonic.points.neuron <- function(neuron){
  points=neuron$d
  xyzmatrix(points[points$Label%in%c(-2,2),])
}
#' @rdname extract.cable
dendritic.points.neuron <- function(neuron){
  points=neuron$d
  xyzmatrix(points[points$Label%in%c(-3,3),])
}
#' @rdname extract.cable
mixed.points.neuron <- function(neuron){ # Mised also means that I do not know
  points=neuron$d
  xyzmatrix(points[points$Label%in%c(8),])
}
#' @rdname extract.cable
dendritic.points.neuronlist <- function(someneuronlist){
  do.call(rbind,nlapply(someneuronlist,dendritic.points.neuron))
}
#' @rdname extract.cable
axonic.points.neuronlist <- function(someneuronlist){
  do.call(rbind,nlapply(someneuronlist,axonic.points.neuron))
}
#' @rdname extract.cable
mixed.points.neuronlist <- function(someneuronlist){
  do.call(rbind,nlapply(someneuronlist,mixed.points.neuron))
}
#' @export
#' @rdname extract.cable
axonal.endings <- function(neuron){
  points=neuron$d[nat::endpoints(neuron)[which(endpoints(neuron)!=rootpoints(neuron))],]
  xyzmatrix(points[points$Label%in%c(-2,2),])
}
#' @export
#' @rdname extract.cable
dendritic.endings <- function(neuron){
  points=neuron$d[nat::endpoints(neuron)[which(endpoints(neuron)!=rootpoints(neuron))],]
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
#' @rdname correctsoma
correctsoma <- function(someneuronlist, brain = NULL,...){
  correctedsomas = neuronlist()
  nopen3d()
  for (n in 1:length(someneuronlist)){
    print(n)
    w = someneuronlist[[n]]
    print(names(someneuronlist[n]))
    if(!is.null(brain)){plot3d(brain)}
    plot3d(w, soma = T)
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
        plot3d(corrected, soma = T, col = "blue")
        progress = readline(prompt="Good enough? T/F  ")
      }
    }
    clear3d()
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
#' @param targets secodn set of neurons
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname overlap.connectivity.matrix
overlap.connectivity.matrix <- function(neurons,targets,neuropil = NULL,delta =1){
  if(length(neuropil)>0){
    points = rbind(dendritic.points.neuronlist(neurons),mixed.points.neuronlist(neurons))
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

#' Generate a connectivity matrix based on euclidean distance between points
#'
#' @description Generates an 'overlap matrix' of overlap scores between neurons in the 'neurons' and 'targets' pools.
#' For every point in a given neuron in 'neurons', a distance score is calculated to every point in a neuron in 'targets'.
#' The sum of this score is added to the final output matrix. The score is calculated as e(-d^2/2δ^2), where d is the euclidean distance between the two points,
#' and δ is the expected distance in um that is considered 'close'.
#'
#' @param neurons first set of neurons
#' @param targets second set of neurons
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname cable.inside.neuropils
cable.inside.neuropils<-function(neuron, brain = FCWBNP.surf, method = c("neurites","axons","dendrites"), stepsize = 0.1, min.endpoints = 2,alpha=30, ...) UseMethod("cable.inside.neuropils")

#' @rdname cable.inside.neuropils
cable.inside.neuropils.neuron <- function(neuron, brain = FCWBNP.surf, method = c("neurites","axons","dendrites"), stepsize = 0.1, min.endpoints = 2,alpha=30){
  require(nat.flybrains)
  require(alphashape3d)
  neuron = resample(neuron, stepsize = stepsize)
  targets = c(0,2,3,8)
  if (method=="axons"){targets = c(2)}
  if (method=="dendrites"){targets = c(3)}
  endings <- function(neuron){
    points=neuron$d[nat::endpoints(neuron)[which(endpoints(neuron)!=rootpoints(neuron))],]
    xyzmatrix(points[points$Label%in%targets,])
  }
  in.neuropil <- function(neuron,neuropil,stepsize =0.1, min.endpoints =2,alpha=30){
    neuropil = ashape3d(xyzmatrix(neuropil),alpha=alpha)
    if(sum(alphashape3d::inashape3d(points=endings(neuron),as3d=neuropil))>min.endpoints){
      points = neuron$d[neuron$d$Label%in%targets,]
      sum(alphashape3d::inashape3d(points=nat::xyzmatrix(points),as3d=neuropil))
    }else{0}
  }
  n = sapply(brain$RegionList, function(x) in.neuropil(neuron,neuropil=subset(brain, x),stepsize =stepsize, min.endpoints=min.endpoints,alpha=alpha))*stepsize
}

#' @rdname cable.inside.neuropils
cable.inside.neuropils.neuronlist <- function(someneuronlist, brain = FCWBNP.surf, method = c("neurites","axons","dendrites"), stepsize = 0.1,min.endpoints = 2,alpha=30){
  nlapply(someneuronlist, cable.inside.neuropils.neuron, brain=brain, method=method,stepsize=stepsize)
}


#' Write ordered swc
#'
#' @description  Write separate, ordered swc files for the MATLAB SPIn package.
#' One can also write separate files for a neuron's axon, dendrite and other cable
#'
#' @param neuron a neuron object
#' @param targets secodn set of neurons
#' @param ... additional arguments passed to methods
#'
#' @export
#' @rdname write.spin.swc
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
#' @param inputs whether or not to show input synapses are to be shown
#' @param outputs whether or not output synapses are to be shown
#' @param printout whether or not to plot a basic legend to indicate what colours mean which neuron
#' @param ... additional arguments passed to methods
#'
#' @export
#' @rdname synapsecolours.neuron
synapsecolours.neuron <-function(neuron, skids = NULL, col = "black", inputs = T, outputs = T,printout=F){
  if(is.neuronlist(neuron)){neuron = neuron[[1]]}
  plot3d(neuron,WithNodes=F,soma=T,col=col,lwd=2)
  plot.new()
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
