# Alphashape3d related functions

combine.alphashape = function (ashapelist)
{
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


transform3dalphashape = function (ashape, transformations){
  positions = ashape$x
  if (is.list(transformations) == F){
    cat("Single transformation")
    positions <- xform(positions, transformations)
  }
  if (is.list(transformations) == T){
    for (transformation in transformations){
      positions <- xform(positions, transformation)
    }
  }
  shape = ashape3d(positions, alpha = ashape$alpha)
  return (shape)
}

WriteVTKalphashape <-function(ashape, filename, title = filename, datatype=c("float","double")){
  d = ashape$x
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())

  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)

  write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)

  data = ashape$triang
  keeps <- apply(data, 1, function(x) {( any(as.numeric(x[9]) > 1))} ) # Removes rows for triangles not included in the alphashape, for the chosen alpha. Includes other simplexes: interior, regular and singular.
  mx = data.matrix(data[keeps,][,1:3]-1) # VTK files are 0 indexed
  numpoints = rep(3, nrow(mx))
  mx = cbind(numpoints, mx)
  cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
  write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
}


neurites.inside <- function(alpha, db, degree = NULL){
  selection = c()
  for (neuron in 1:length(db)){
    neuron = db[neuron]
    xyz = xyzmatrix(neuron)
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

