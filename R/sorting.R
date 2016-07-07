# Sorting neurons
assignsides.l1 <- function(db, nan2microns = F){
  scale = 1000
  if (nan2microns == F){
    scale = 1
    cat("Note that query neurons should be scaled to microns")
  }
  r.l1.alpha = readRDS("right_half_l1mesh") # Fix this at some point, add it to package?
  fb = db
  #good_soma=sapply(db, function(x) !is.null(x$tags$soma))
  #fb = db[good_soma]
  #if (length(db) != length(fb)){
  #  warning(paste("Dropping", length(db) - length(fb), "neurons without tagged somas", sep = " "))
  #}
  fdf=as.data.frame(fb)
  fdf=cbind(fdf, side = matrix(0,ncol=1, nrow=nrow(fdf)))
  for (neuron in 1:length(fb)){
    k = 0
    rep = 'yes'
    while (rep == 'yes'){
      k = k + 100
      som = soma.neuron(fb[[neuron]])
      if (is.na(som[1])){
        if (length(fb[[neuron]]$tags$soma[[1]])>0){
          som = matrix(xyzmatrix(fb[[neuron]])[fb[[neuron]]$d$PointNo%in%fb[[neuron]]$tags$soma,], ncol = 3)
        }else{
          som = matrix(xyzmatrix(fb[[neuron]])[fb[[neuron]]$StartPoint,], ncol = 3)
        }
      }
      p = nat::xyzmatrix(fb[[neuron]])
      n = nabor::knn(p, som, ifelse(nrow(p)>k,k,nrow(p)))
      m = p[c(n$nn.idx),]
      r = alphashape3d::inashape3d(r.l1.alpha, points = m/scale) # scale to microns
      rn = sum(r, na.rm=TRUE)
      nn = ifelse(nrow(m) >0, nrow(m),1)
      e = rn/nn
      if (e > 0.6){
        fdf[neuron,]$side = "right"
        rep = "no"
      } else if (e < 0.4){
        fdf[neuron,]$side = "left"
        rep = "no"
      } else{
        rep = "yes"
        if (nrow(p) < k | k > 500){
          r = alphashape3d::inashape3d(r.l1.alpha, points = som/scale) # scale to microns
          pos = ifelse(r, 'right','left')
          fdf[neuron,]$side = pos
          rep = "no"
        }
      }
    }
  }
  # attr(fb,'df')=fdf
  return(fdf)
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


assignside.names <- function(someneuronlist){
  sdf=as.data.frame(someneuronlist)
  sdf=transform(sdf, side=factor(ifelse(grepl("right|Right|_r$|R$|r$|left|Left|_l$|L$|l$", name),ifelse(grepl("right|Right|_r|R$|r$", name),"R","L"), "NA")))
  attr(someneuronlist,'df')=sdf
  return(someneuronlist)
}


# Sorting neurons
assignsynapses <- function(fb, thresh = 0.7, nan2microns = F){
  scale = 1000
  if (nan2microns == F){
    scale = 1
    cat("Note that query neurons should be scaled to microns")
  }
  good_neurites = sapply(fb, function(x) ifelse(!is.list(x),is.list(x),length(x$connectors$x) > 3))
  fb = fb[good_neurites]
  if (sum(good_neurites) < length(good_neurites)){
    warning(paste(length(good_neurites)-sum(good_neurites), "fragment(s) dropped", sep = " "))
  }
  r.l1.alpha = readRDS("~/projects/alexanalysis/leftrightfinal/right_half_l1mesh") # Fix this at some point, add it to package?
  fdf=as.data.frame(fb)
  fdf=cbind(fdf, main.input.side = matrix(0,ncol=1, nrow=nrow(fdf)), main.output.side = matrix(0,ncol=1, nrow=nrow(fdf)))
  for (neuron in 1:length(fb)){
    pre.points = get.connectors(fb[neuron], target = "PRE")
    post.points = get.connectors(fb[neuron], target = "POST")
    pre = alphashape3d::inashape3d(r.l1.alpha, points = pre.points/scale) # scale to microns
    post = alphashape3d::inashape3d(r.l1.alpha, points = post.points/scale) # scale to microns
    pre = sum(pre, na.rm=TRUE)/length(pre)
    post = sum(post, na.rm=TRUE)/length(post)
    if (is.na(pre)|is.null(pre)){
      fdf[neuron,]$main.input.side = NA
    } else if (pre > thresh){
      fdf[neuron,]$main.input.side = "right"
    } else if (pre < (1-thresh)){
      fdf[neuron,]$main.input.side = "left"
    } else {
      fdf[neuron,]$main.input.side = "both"
    }
    if (is.na(post)|is.null(post)){
      fdf[neuron,]$main.output.side = NA
    } else if (post > thresh){
      fdf[neuron,]$main.output.side = "right"
    } else if (post < (1-thresh)){
      fdf[neuron,]$main.output.side = "left"
    } else {
      fdf[neuron,]$main.output.side = "both"
    }
  }
  # attr(fb,'df')=fdf
  return(fdf)
}



find.symmetric <- function(someneuronlist, brain, sym = NULL){
  res = c()
  good_neurites = sapply(someneuronlist, function(x) ifelse(!is.list(x),is.list(x),length(x$d$X) > 3))
  someneuronlist = someneuronlist[good_neurites]
  if (sum(good_neurites) < length(good_neurites)){
    warning(paste(length(good_neurites)-sum(good_neurites), "fragment(s) dropped", sep = " "))
  }
  for (someneuron in 1:length(someneuronlist)){
    print(someneuron)
    someneuron = someneuronlist[someneuron]
    if (!is.null(sym)){
      someneuron = applyTransform.neuronlist(someneuron, symtransformation)
      xyzmatrix(brain) <- Morpho::applyTransform(xyzmatrix(brain), symtransformation)
    }
    if(!nat::is.dotprops(someneuron)){
      someneuron = nat::dotprops(someneuron, k=5, .progress = 'text', resample = 1, OmitFailures = F)
    }
    someneuron.flipped = nat::mirror(someneuron, mirrorAxisSize = nat::boundingbox(brain))
    names(someneuron.flipped) = paste(names(someneuron.flipped), "flp", sep = "")
    result = nat.nblast::nblast(someneuron, target = c(someneuron, someneuron.flipped), normalised = T)[2,]
    names(result) = names(someneuron)
    res = c(res, result)
  }
  return(res)
}

