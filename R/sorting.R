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
