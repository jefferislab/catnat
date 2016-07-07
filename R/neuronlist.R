# Some extra neuronlist functions
convert <- function(someneuronlist, factor = 1/1e3){
  for (neuron in 1:length(someneuronlist)){
    if (length(someneuronlist[[neuron]]$d) == 7){
      someneuronlist[[neuron]]$d$X <- someneuronlist[[neuron]]$d$X*factor
      someneuronlist[[neuron]]$d$Y <- someneuronlist[[neuron]]$d$Y*factor
      someneuronlist[[neuron]]$d$Z <- someneuronlist[[neuron]]$d$Z*factor
    }
    if (length(someneuronlist[[neuron]]$connectors) == 6){
      someneuronlist[[neuron]]$connectors$x <- someneuronlist[[neuron]]$connectors$x*factor
      someneuronlist[[neuron]]$connectors$y <- someneuronlist[[neuron]]$connectors$y*factor
      someneuronlist[[neuron]]$connectors$z <- someneuronlist[[neuron]]$connectors$z*factor
    }
  }
  return (someneuronlist)
}


assignside <- function(someneuronlist){
  sdf=as.data.frame(someneuronlist)
  sdf=transform(sdf, side=factor(ifelse(grepl("right|Right|_r$|R$|r$|left|Left|_l$|L$|l$", name),ifelse(grepl("right|Right|_r|R$|r$", name),"R","L"), "NA")))
  attr(someneuronlist,'df')=sdf
  return(someneuronlist)
}


join.neuronlists <- function(...){
  argg <- list(...)
  skids = c()
  for (item in 1:length(argg)){
    skids = c(skids, c(as.data.frame(argg[[item]])$skid))
  }
  neurons = subset(db, as.data.frame(db)$skid%in%skids)
  return(assignside(neurons))
}


primary.neurite <- function(someneuron, k = 100){      # Find the first 100 points of the primary neurite
  som = soma.neuron(someneuron)
  if (is.na(som[1])){
    if (length(someneuron[[1]]$tags$soma[[1]])>0){
      som = matrix(xyzmatrix(someneuron)[someneuron[[1]]$d$PointNo%in%someneuron[[1]]$tags$soma,], ncol = 3)
    }else{
      som = matrix(xyzmatrix(someneuron)[someneuron[[1]]$StartPoint,], ncol = 3)
    }
  }
  p = nat::xyzmatrix(someneuron)
  n = nabor::knn(p, som, ifelse(nrow(p)>k,k,nrow(p)))
  m = p[c(n$nn.idx),]
}

