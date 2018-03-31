#' Create a skeleton from a skeletonised neuron saved as a .Nrrd file
#'
#' @description Skeletonise a neuron in Fiji and then use this function to retrieve it as a skeleton
#'
#' @param files paths to saved .Nrrd files that have been skeletonsied in Fiji
#' @param connection.distance maximal connection distance between points
#' @param k value used in nearest neighbours search to identify close points and draw lines between them
#' @param distance.steps each round of the algorithm connects the cosest points to the start point / leaf node, that are within the search range, starting at distance.steps and increasing by this amount to the maximal connection.distance
#' @param ... additional arguments passed to methods.
#'
#' @return A neuronlist object
#' @export create_skeleton
create_skeleton_from_nrrd <- function(files, connection.distance = 25, k = "all", distance.steps = 0.5, ...){
  nl = nat::neuronlist()
  message("Note: Argument 'files' should be the paths neurons skeletonised in Fiji and saved as a .Nrrd files")
  for(file in files){
    message(paste0("Reading Nrrd file: ",file, "     ",match(file,files),"/",length(files)))
    img = nat::dotprops(file)$points # Read in point information
    message("Generating draft neuron skeleton in SWC format...")
    swc = data.frame(PointNo = 1:nrow(img),Label = 0,img,W = -2, Parent = -1, tree = NA) # Create SWC data.frame
    if(k=="all"){
      k = nrow(swc)
    }
    nears = nabor::knn(data = nat::xyzmatrix(swc), query = nat::xyzmatrix(swc), k = k) # Find the closest k nodes
    # Is k big enough?
    if(k!="all"){
      undershot = sum(nears[[2]][,k]<=connection.distance)
      if(undershot>0){
        message(paste0("k may be too small. The kth nearest node is less than the connection.distance for ",undershot,"/",nrow(swc)," nodes"))
      }
    }
    skel = 1 # We will grow our tree from, as start, the root for the first skeleton as the first point
    used = c()
    d = min(nears[[2]][,-1]) # Start with the smallest connection distance, then grow the connection distance to its maximum. dd is the lower bound.
    parents = nears[[1]][,-1][which.min(nears[[2]][,-1])] # Start with a faux 'root'
    while(sum(is.na(swc$tree))>0){
      if(all(parents%in%used)){ # If we run out of nodes in the connection range
        used = c()
        skel = skel + 1; leaves = c()
        dists = nears[[2]][-1*used,-1]
        indices = nears[[1]][-1*used,-1]
        parents = indices[which.min(dists)]
        parents = subset(swc,is.na(tree))$PointNo[[1]]
      }
      search = nears[[1]]; search[!nears[[2]]<=d] = NA # Get the closest matches in space to each non-zero point in the image stack
      used = c(used,parents) # Parent nodes cannot be assigned to matches already within the tree
      children = search[parents,]# Find the right match row
      children = unique(children[!is.na(children)]) # Get children, i.e. candidate child nodes, within the connection distance
      children = children[!children%in%used] # Remove the nodes we have already added to our skeleton
      swc[unique(c(parents,children)),"tree"] = skel # So now we'll assign them to a skeleton
      if(length(children[!is.na(children)])>0){ # If there are children
        if(length(children)>1){
          closest.parents = apply(search[children,],1, function(x) x[x%in%parents][1])
          if(sum(is.na(closest.parents))>0){ # The the closest parent is not within the top k matches, look the other way
            lost = children[is.na(closest.parents)]
            closest.parents[is.na(closest.parents)] = sapply(lost,function(x) rep(parents,k)[which.max(search[parents,]==x)[1]])
          }
        }else{
          closest.parents = search[children,][search[children,]%in%parents][1]
          if(sum(is.na(closest.parents))>0){ # The the closest parent is not within the top k matches, look the other way
            closest.parents = rep(parents,k)[which.max(search[parents,]==children)][1]
          }
        }
        swc[children,"Parent"] = closest.parents # Assign to the closest parent
        parents = children # Now grow tree from the old children
      }
    }
    # Now we have created a bunch of skeletons at the minimal connection distance, d.
    # We need to grow our search aread to d.max and so merge these different skeletons
    dd = d # Minimal
    d = d + distance.steps # Maximal
    while(d<=connection.distance){
      # Get the indices within the search distance
      search = nears[[1]]
      search[!nears[[2]]<=d] = NA
      search[!nears[[2]]>dd] = NA
      # Get the distances
      dists = nears[[2]]
      m = reshape2::melt(dists)
      m = subset(m,value<=d)
      m = subset(m,value>dd)
      # Order indices to grow by shortest distance to next node
      m = m[order(m$value),]
      for(parent in m$Var1){
        t = swc[parent,"tree"]
        used = subset(swc,tree==t)$PointNo
        children = search[parent,]# Find the right match row
        children = unique(children[!is.na(children)]) # Get children, i.e. candidate child nodes, within the connection distance
        children = children[!children%in%used] # Remove the nodes we have already added to our tree
        children = children[!duplicated(swc[children,"tree"])] # Children must be in different trees from one another
        if(length(children)>0){
          for(child in children){ # Now combine the trees
            child.tree = swc[as.character(child),"tree"]
            n = subset(swc,tree==child.tree)
            if(nrow(n)>1){ # Connect another tree
              n = nat::as.neuron(nat::as.ngraph(n), origin = child)$d
              n$tree = t
              n[n$PointNo==child,"Parent"] = parent
              swc[n$PointNo,] = n
            }else if (nrow(n)==1){ # Connect isolated node
              n$tree = t
              n[n$PointNo==child,"Parent"] = parent
              swc[n$PointNo,] = n
            }
          }
        }
      }
      # Increase the search distance
      dd = d # New min
      d = d + distance.steps # New max
    }
    #swc = swc[,1:7]
    isolated = sapply(swc$PointNo,function(x) (!x%in%swc$Parent)&(swc[x,"Parent"]==-1))
    swc = swc[!isolated,] # Remove isolated nodes
    nn = nat::as.neuron(swc)
    nl = c(nl,nat::as.neuronlist(nn))
    names(nl) = c(names(nl),file)
  }
  attr(nl,"df") = data.frame(file=names(nl))
  nl
}
