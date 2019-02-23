#' Upload neuron(s) to CATMAID
#'
#' @description  Uploads neurons to CATMAID, names them and annotates them.
#' Please use with caution, as you could be heavily adding to a live tracing environment.
#' @param swc a neuron or neuronlist object, or (a) local file path(s) to your saved .swc file(s).
#' @param name whatever you want to name your uploaded neurons. If a single character, then it will be added to all uploaded neurons. Else, can be a character vector the same length as swc.
#' @param annotations a character vector of annotations, to be added to all ofthe uploaded neurons
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_upload_neurons
catmaid_upload_neurons <- function (swc, name ="neuron SWC upload", annotations = "SWC upload",
                                    include.tags = TRUE,include.connectors = TRUE,
                                    search.range.nm = 100,
                                    pid = 1, conn = NULL, max.upload = 10, ...) {
  if(is.neuronlist(swc)|is.neuron(swc)){
    temp.files = paste0(names(swc),".swc")
    nat::write.neurons(swc, format= "swc",dir = tempdir(), files = temp.files, ...)
    swc = paste0(tempdir(),temp.files)
  }
  if (length(swc)>max.upload){
    stop(paste0('You are uploading a large number of SWC files.
                Are you sure you want to do this?
                If so change the value of the max.upload argument'))
  }
  if(length(name)==1&length(swc)>1){
    name = paste0(name,"_",seq_along(swc))
  } else if (length(name)>length(swc)|length(swc)>length(name)){
    stop(paste0('The names argument mut be a vector or length 1, to be applied to all neurons
                or a vector the same length as swc'))
  }
  skids = c()
  for(file in 1:length(swc)){
    post_data = list()
    post_data["name[1]"] = as.list(name[file])
    post_data["file[1]"] = list(upload_file(swc[file]))
    path = sprintf("/%d/skeletons/import", pid)
    res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                 simplifyVector = T, conn = conn, ...)
    invisible(catmaid:::catmaid_error_check(res))
    if(include.tags&is.neuron(swc[[file]])){
      fafbseg_transfer_tags(new.neuron = res$skeleton_id, old.neuron = swc[[file]], search.range.nm = search.range.nm, pid = pid, conn = conn, ...)
    }
    if(include.connectors&is.neuron(swc[[file]])){
      fafbseg_transfer_connectors(new.neuron = res$skeleton_id, old.neuron = swc[[file]], links = TRUE, search.range.nm = search.range.nm, pid = pid, conn = conn, ...)
    }
    skids = c(skids,res$skeleton_id)
  }
  if(!is.null(annotations)){
    sapply(seq_along(skids), function(x) catmaid_rename_neuron(skids = skids[x], names = name[x], pid = pid, conn = conn))
    catmaid::catmaid_set_annotations_for_skeletons(skids = as.numeric(skids),
                                                   annotations = annotations,
                                                   pid = pid, conn = conn, ...)
    message(paste0("Annotations ", annotations," set to new skeletons: ", res$skeleton_id))
  }
}

# fafbseg.skid = 210077691
# skid = 2428780

#' Transfer information between v14 and v14-seg neurons
#'
#' @description  Functions to enable the mapping of connectors and tags between cognate FAFB v14 and v14-seg neurons, to help enable neuronal reconstruction using the partial auto-segmentation in v14-seg space
#' @param new.neuron a neuron object, or the skeleton ID of a neuron in v14 space
#' @param fafbseg.skid a neuron object, or the skeleton ID of a neuron in v14-seg space
#' @param search.range.nm the maimum distance between points in the v14 and v14-seg skeleton, that is acceptable for the transfer ot tag/connector information. Not set to 0 be default to accept small variations in node position that might occur between the two tracing environments, due to tracer edit.s
#' @param links whether or not to link transferred connectors  to post-pre synaptic sites on the neuron identified with the new.neuron argument
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_set_labels
#' @export
#' @rdname fafbseg_transfer
fafbseg_transfer_tags<- function(new.neuron, old.neuron, search.range.nm = 100, pid = 1, conn = NULL, ...){
  if(!is.neuron(new.neuron)){
    new.neuron = catmaid::read.neuron.catmaid(new.neuron)
  }
  if(!is.neuron(old.neuron)){
    old.neuron = read.neurons.fafbseg(old.neuron)[[1]]
  }
  labels = unlist(sapply(1:length(old.neuron$tags),function(t) rep(names(old.neuron$tags[t]),length(old.neuron$tags[[t]]))))
  tagged.points = unlist(old.neuron$tags)
  tagged.points = tagged.points[tagged.points%in%old.neuron$d$PointNo]
  tagged.points = subset(old.neuron$d,PointNo%in%tagged.points)
  tagged.points$tag = labels
  near = nabor::knn(query = tagged.points[,c('X','Y', 'Z')], data =nat::xyzmatrix(new.neuron), k=1, radius=search.range.nm)$nn.idx
  near = near[near!=0]
  labels = tagged.points$tag[near!=0]
  nodes = new.neuron$d$PointNo[near]
  if(length(nodes)!=length(labels)){
    stop("error with assigning tag information")
  }
  if(length(nodes)>0){
    for(i in 1:length(nodes)){
      catmaid::catmaid_set_labels(node = nodes[i] , labels = labels[i], type= "treenode", pid = pid, conn = conn, ...)
    }
  }else{
    NULL
  }
}
#' @export
#' @rdname fafbseg_transfer
fafbseg_transfer_connectors<- function(new.neuron, old.neuron, search.range.nm = 100, links = TRUE, pid = 1, conn = NULL, ...){
  if(!is.neuron(new.neuron)){
    new.neuron = catmaid::read.neuron.catmaid(new.neuron)
  }
  if(!is.neuron(old.neuron)){
    old.neuron = read.neurons.fafbseg(old.neuron)[[1]]
  }
  c.df = old.neuron$connector
  near = nabor::knn(query = nat::xyzmatrix(c.df), data =nat::xyzmatrix(new.neuron), k=1, radius=search.range.nm)$nn.idx
  new.c.df = c.df[near!=0,]
  new.c.df$near_id = new.neuron$d[near,"PointNo"]
  if(nrow(new.c.df)>0){
    for(i in 1:nrow(new.c.df)){
      post_data = list()
      post_data["x"] = new.c.df[i,"x"]
      post_data["y"] = new.c.df[i,"y"]
      post_data["z"] = new.c.df[i,"z"]
      path = sprintf("/%d/connector/create", pid)
      res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                   simplifyVector = T, conn = conn, ...)
      invisible(catmaid:::catmaid_error_check(res))
      if(links){
        from_id = ifelse(new.c.df[i,"prepost"]==0,new.c.df[i,"near_id"],res$connector_id)
        to_id = ifelse(new.c.df[i,"prepost"]!=0,new.c.df[i,"near_id"],res$connector_id)
        link_type = ifelse(new.c.df[i,"prepost"]==0,"presynaptic_to","postsynaptic_to")
        state = sprintf('{"%s":"%s"}, {"%s":"%s"}',from_id, res$connector_edition_time, to_id, res$connector_edition_time)
        post_data = list()
        post_data["from_id"] = from_id
        post_data["to_id"] = to_id
        post_data["link_type"] = link_type
        post_data["state"] = state
        path = sprintf("/%d/link/create", pid)
        res2 = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                     simplifyVector = T, conn = conn, ...)
        invisible(catmaid:::catmaid_error_check(res2))
      }
    }
  }
}

#' Find the location of specified tags for a CATMAID neuron
#'
#' @description  Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
#' @param x a neuron or neuronlist object
#' @param tag a single character specifying which tag to look for. Defaults to TODO.
#' @param url if TRUE (defualt) a list of URLs pertaining to specified tag locations are returned. If FALSE, a data.frame subsetted from x$d is returned, reporting treenode ID and X,Y,Z positions for specified tags
#' @param server the CATMAID server for which URLs should be created if url == TRUE
#' @export
#' @rdname catmaid_get_tag
catmaid_get_tag<-function(x, tag = "TODO", url = TRUE, server = "https://neuropil.janelia.org/tracing/fafb/v14/") UseMethod("catmaid_get_TODO")
catmaid_get_tag.neuron <- function(x, url = TRUE, server = "https://neuropil.janelia.org/tracing/fafb/v14/"){
  TODO = unique(unlist(x$tags[[tag]]))
  if(is.null(TODO)){
    NULL
  }else{
    df = subset(x$d,PointNo%in%TODO)
    if(url){
      catmaid_url = paste0(server, "?pid=1")
      catmaid_url = paste0(catmaid_url, "&zp=", df[["Z"]])
      catmaid_url = paste0(catmaid_url, "&yp=", df[["Y"]])
      catmaid_url = paste0(catmaid_url, "&xp=", df[["X"]])
      catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
      catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
      invisible(catmaid_url)
    }
    else{
      df
    }
  }
}
catmaid_get_tag.neuronlist <- function(x, server = "https://neuropil.janelia.org/tracing/fafb/v14/"){
  if(url){
    unlist(lapply(x,catmaid_get_tag.neuron, url=url,server = server))
  }else{
    do.call(rbind,lapply(x,catmaid_get_tag.neuron, url=url,server = server))
  }
}

#' Find the TODO tagged merge sites between a given neuron and the CATMAID database at large
#'
#' @description  Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
#' @param TODO a data.frame, as returned by catmaid_get_tag with tag = "TODO" and url = FALSE
#' @param skid the skeleton ID of a neuron for which merges are to be found
#' @param fafbseg whether or not to use fafbseg::read_brainmaps_meshes on the TODO tag locations and restrict possible merges to neurons within the search range and within the cognate auto-segmented volume
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @export
#' @rdname catmaid_find_likely_merge
catmaid_find_likely_merge <- function(TODO, skid, pid=1, fafbseg = FALSE, min_nodes = 2, search.range.nm = 100){
  if(!is.neuron(x)){
    stop("x is not a neuron object")
  }
  possible.merges = data.frame()
  for(i in 1:nrow(TODO)){
    todo = df[i,]
    seg = tryCatch(fafbseg::brainmaps_xyz2id(todo[,c('X','Y', 'Z')]),error=function(e)NULL)
    if(is.null(seg)){
      warning("cannot find auto-segmented volume containing TODO tag")
    }else{
      bbx = rbind(todo[,c('X','Y', 'Z')]-search.range.nm,todo[,c('X','Y', 'Z')]+search.range.nm)
      skids.bbx = catmaid_skeletons_in_bbx(bbx=bbx,pid=pid,min_nodes=min_nodes)
      neuron.bbx = catmaid::read.neurons.catmaid(skids.bbx, OmitFailures = TRUE)
      neuron.bbx = neuron.bbx[setdiff(names(neuron.bbx),skid)]
      if(fafbseg){
        require(fafbseg)
        s = fafbseg::find_merged_segments(seg)
        vol = tryCatch(fafbseg::read_brainmaps_meshes(s), error = function(e) "read error")
        in.vol = sapply(neuron.bbx,function(n) sum(nat::pointsinside(nat::xyzmatrix(n),surf=vol))>0)
        neuron.bbx = neuron.bbx[in.vol]
        message(length(neuron.bbx), " fragments found in fafb auto-traced segment corresponding with TODO tag")
      }else{
        message(length(neuron.bbx), " fragments found within ", search.range.nm," nm of TODO tag")
      }
      neuron.bbx.d = do.call(rbind,lapply(1:length(neuron.bbx),function(n) cbind(neuron.bbx[[n]]$d,skid=names(neuron.bbx[n]))))
      near = nabor::knn(query = todo[,c('X','Y', 'Z')], data =nat::xyzmatrix(neuron.bbx.d),k=10, radius=search.range.nm)$nn.idx
      near = near[near!=0]
      near = neuron.bbx.d[near,]
      near$TODO = todo$PointNo
      near$merger = skid
      possible.merges = rbind(possible.merges,near)
    }
  }
  possible.merges
}

#' Interactively choose to join neurons in CATMAID via rgl in R
#'
#' @description Interactively choose to join neurons in CATMAID via rgl in R. Use with caution in live tracing environments.
#' @param possible.merges a data frame of possible mergers between tree nodes for neurons in a CATMAID database, as returned via catmaid_find_likely_merge
#' @param brain the brain to plot while visualising potential mergers using rgl. Defaults to NULL, no brain plotted.
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @export
#' @rdname catmaid_interactive_join
catmaid_interactive_join <- function(possible.merges, brain = NULL, pid = 1, conn = conn, ...){
  nat::nopen3d()
  neurons = catmaid::read.neurons.catmaid(skids=unique(possible.merges$skid))
  for(todo in unique(possible.merges$TODO)){
    TODO.possible = subset(possible.merges,TODO==todo)
    continue = FALSE
    i = 1
    while(!progress){
      if(!is.null(brain)){
        rgl::plot3d(brain,alpha=0.1,col="grey")
      }
      skid = possible.merges[i,"skid"]
      neuron = neurons[skid]
      merger = catmaid_get_treenodes_detail(tnids = possible.merges[i,"merger"], pid = pid, conn = conn, ...)
      merger.neuron = neurons[merger$skid]
      rgl::plot3d(neuron,col="black",lwd=2,soma=T)
      rgl::plot3d(merger.neuron,col="red",lwd=2,soma=T)
      rgl::spheres3d(nat::xymatrix(merger), col="orange")
      rgl::spheres3d(nat::xymatrix(possible.merges[i,]), col="grey")
      progress = readLines("Make merge? y = yes, n = no, c = cycle  :")
      if(progress=="y"){
        nat::npop3d()
        rgl::plot3d(merger.neuron,col="green",lwd=2,soma=T)
        sure = readLines("Sure? y = yes, n = no  :")
        if(sure=="y"){
          catmaid_join_skeletons(from_treenode_id = possible.merges[i,"merger"], to_treenode_id = possible.merges[i,"TODO"], pid = pid, conn = conn, ...)
          continue=TRUE
        }else{
          i = ifelse(i==nrow(TODO.possible),i,i+1)
        }
      }else if(progress=="n"){
        TODO.possible = TODO.possible[-i,]
        i = ifelse(i==nrow(TODO.possible),i,i+1)
      }else{
        i = ifelse(i==nrow(TODO.possible),i,i+1)
      }
    }
  }
}


# to_treenode_id = 1254672855
# from_treenode_id = 676710716

#' Programmatically join CATMAID skeletons
#'
#' @description  Programmatically join CATMAID skeletons. Use with caution in live tracing environments.
#' @param from_treenode_id the treenode ID of the downstream neuron
#' @param from_treenode_id the treenode ID of the upstream neuron
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_join_skeletons
catmaid_join_skeletons <- function(from_treenode_id, to_treenode_id, pid = 1, conn = NULL, ...){
  post_data = list()
  post_data["from_id"] = to_treenode_id
  post_data["to_id"] = from_treenode_id
  path = sprintf("/%d/skeleton/join", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                               simplifyVector = T, conn = conn, ...)
  invisible(catmaid:::catmaid_error_check(res))
  message(res$message)
}

#' Search for CATMAID skeletons within a volume
#'
#' @description  Programmatically search for skeleton IDs pertaining to neurons within a search volume defined by a bounding box.
#' @param bbx the bounding box (a matrix of 2 rows and 3 columns) describing a search volume
#' @param min_nodes the minimum number of nodes a neuron in the search area must have (includes ndoes outside search area)
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_skeletons_in_bbx
catmaid_skeletons_in_bbx <- function(bbx, min_nodes = 2, pid = 1, conn = NULL, ...){
  post_data = list()
  post_data["minx"] = bbx[1,1]
  post_data["miny"] = bbx[1,2]
  post_data["minz"] = bbx[1,3]
  post_data["maxx"] = bbx[2,1]
  post_data["maxy"] = bbx[2,2]
  post_data["maxz"] = bbx[2,3]
  post_data["min_nodes"] = min_nodes
  path = sprintf("/%d/skeletons/in-bounding-box", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                               simplifyVector = T, conn = conn, ...)
  res
}

#' Try to detect whether a neuron your are about to import into CATMAID might cause a duplication in that CATMAID instance
#'
#' @description  Programmatically search for skeleton IDs pertaining to neurons within a search volume defined by a bounding box.
#' @param neuron the neuron object you are thinking about uploading
#' @param tolerance how many potentially duplictated nodes  you will tolerate
#' @param search.range.nm determines the size of the bounding box around each node in neuron to search for a duplicated. Defaults to 1 nm
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_skeletons_in_bbx
catmaid_duplicated <- function(neuron, tolerance = 0.5, search.range.nm = 1, pid = 1, conn = NULL){
  duplicated = c()
  for(i in 1:nrow(neuron$d)){
    bbx = rbind(nat::xymatrix(neuron$d[i,])-search.range.nm, nat::xymatrix(neuron$d[i,])+search.range.nm)
    skids.bbx = catmaid_skeletons_in_bbx(bbx=bbx,pid=pid, min_nodes = 2, pid = pid, conn = conn)
    duplicated = c(duplicated,length(skids.bbx)>0)
  }
  if(sum(duplicated)/nrow(neuron$d)>tolerance){
    TRUE
  }else{
    FALSE
  }
}
