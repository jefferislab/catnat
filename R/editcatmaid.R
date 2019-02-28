#' Upload neuron(s) to CATMAID
#'
#' @description  Uploads neurons to CATMAID, names them and annotates them.
#' Please use with caution, as you could be heavily adding to a live tracing environment.
#' @param swc a neuron or neuronlist object, or (a) local file path(s) to your saved .swc file(s).
#' @param name whatever you want to name your uploaded neurons. If a single character, then it will be added to all uploaded neurons. Else, can be a character vector the same length as swc.
#' @param annotations a character vector of annotations, to be added to all of the uploaded neurons
#' @param return.new.skids if TRUE, the new skids created in the CATMAID instance specificed by conn, are returned
#' @param include.tags whether of not, if swc is a CATMAID neuron/neuronlist, to transfer its tags to its newly uploaded cognate
#' @param include.connectors whether of not, if swc is a CATMAID neuron/neuronlist, to transfer its connectors to its newly uploaded cognate
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_upload_neurons
catmaid_upload_neurons <- function (swc = NULL, name ="neuron SWC upload", annotations = NULL,
                                    include.tags = TRUE,include.connectors = TRUE,
                                    search.range.nm = 0, return.new.skids = FALSE,
                                    pid = 1, conn = NULL, max.upload = 10, ...) {
  annotations = c(annotations,"neuron upload")
  if(nat::is.neuronlist(swc)|nat::is.neuron(swc)){
    neurons = nat::as.neuronlist(swc)
    if(nat::is.neuron(swc)){
     temp.files = "1.swc"
     nat::write.neuron(swc, format= "swc",dir = tempdir(), file = temp.files, Force = TRUE)
    }else{
     temp.files = paste0(names(swc),".swc")
     nat::write.neurons(swc, format= "swc",dir = tempdir(), files = temp.files, Force = TRUE)
    }
    swc = paste0(tempdir(),"/",temp.files)
  }else{
    neurons = NULL
  }
  if (length(swc)>max.upload){
    stop(paste0('You are uploading a large number of SWC files.
                Are you sure you want to do this?
                If so change the value of the max.upload argument'))
  }
  if(!is.null(name)){
    if(length(name)==1&length(swc)>1){
      name = paste0(name,"_",seq_along(swc))
    } else if (length(name)>length(swc)|length(swc)>length(name)){
      stop(paste0('The names argument mut be a vector or length 1, to be applied to all neurons
                  or a vector the same length as swc'))
    }
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
    catmaid::catmaid_set_annotations_for_skeletons(skids = res$skeleton_id,
                                                   annotations = annotations,
                                                   pid = pid, conn = conn, ...)
    tryCatch(catmaid::catmaid_rename_neuron(skids = res$skeleton_id, names = name[file], pid = pid, conn = conn, ...),error = function(e) warning("Could not name new neuron ",res$skeleton_id))
    message(paste0("Annotations and names set to new skeleton: ", res$skeleton_id))
    new.neuron = catmaid::read.neuron.catmaid(res$skeleton_id, pid = pid, conn = conn, ...)
    if(include.tags&is.neuron(neurons[[file]])){
      fafbseg_transfer_tags(new.neuron = new.neuron, old.neuron = neurons[[file]], search.range.nm = search.range.nm, pid = pid, conn = conn, ...)
    }
    if(include.connectors&is.neuron(neurons[[file]])){
      fafbseg_transfer_connectors(new.neuron = new.neuron, old.neuron = neurons[[file]], links = TRUE, search.range.nm = search.range.nm, pid = pid, conn = conn, ...)
    }
    skids = c(skids,res$skeleton_id)
  }
  #sapply(seq_along(skids), function(x) catmaid_rename_neuron(skids = skids[x], names = name[x], pid = pid, conn = conn, ...))
  if(return.new.skids){
    return(skids)
  }
}

#' Transfer information between v14 and v14-seg neurons
#'
#' @description  Functions to enable the mapping of connectors and tags between cognate FAFB v14 and v14-seg neurons, to help enable neuronal reconstruction using the partial auto-segmentation in v14-seg space
#' @param new.neuron a neuron object, or the skeleton ID of a neuron in v14 space
#' @param old.neuron a neuron object, or the skeleton ID of a neuron in v14-seg space
#' @param offset whether or not to expect a displacement in space between the new and old neurons. If so, search.range.nm is used to find the nearest nodes to which to assign tags/connectors, within the radius it specifies
#' @param search.range.nm the maimum distance between points in the v14 and v14-seg skeleton, that is acceptable for the transfer ot tag/connector information. Not set to 0 be default to accept small variations in node position that might occur between the two tracing environments, due to tracer edit.s
#' @param links whether or not to link transferred connectors  to post-pre synaptic sites on the neuron identified with the new.neuron argument
#' @param pid project id for conn. Defaults to 1
#' @param pid project id for conn2. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param conn2 CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_set_labels
#' @export
#' @rdname fafbseg_transfer
fafbseg_transfer_tags<- function(new.neuron, old.neuron, offset = TRUE, search.range.nm = 1000, pid = 1, conn = NULL, pid2 = 1, conn2 = fafb_seg_conn(), ...){
  if(!is.neuron(new.neuron)){
    new.neuron = catmaid::read.neuron.catmaid(new.neuron, pid = pid, conn = conn, ...)
  }
  if(!is.neuron(old.neuron)){
    old.neuron = catmaid::read.neuron.catmaid(new.neuron, pid = pid2, conn = conn2, ...)
  }
  if(!is.null(old.neuron$tags)&length(old.neuron$tags)){
    labels = unlist(sapply(1:length(old.neuron$tags),function(t) rep(names(old.neuron$tags[t]),length(old.neuron$tags[[t]]))))
    tagged.points = unlist(old.neuron$tags)
    tagged.points = tagged.points[tagged.points%in%old.neuron$d$PointNo]
    tagged.points = old.neuron$d[match(tagged.points,old.neuron$d$PointNo),]
    tagged.points$tag = labels
    near = nabor::knn(query = tagged.points[,c('X','Y', 'Z')], data =nat::xyzmatrix(new.neuron), k=1, radius=ifelse(offset,search.range.nm,0))$nn.idx
    near = near[near!=0]
    labels = tagged.points$tag[near!=0]
    nodes = new.neuron$d$PointNo[near]
    if(length(nodes)!=length(labels)){
      warning("error with assigning tag information")
    }
    if(length(nodes)>0){
      for(i in 1:length(nodes)){
        catmaid::catmaid_set_labels(node = nodes[i] , labels = labels[i], type= "treenode", pid = pid, conn = conn, ...)
      }
    }
  }
}
#' @export
#' @rdname fafbseg_transfer
fafbseg_transfer_connectors<- function(new.neuron, old.neuron,
                                       offset = FALSE, search.range.nm = 1000, links = TRUE,
                                       pid = 1, conn = NULL, pid2 = 1, conn2 = fafb_seg_conn(), ...){
  if(!is.neuron(new.neuron)){
    new.neuron = catmaid::read.neuron.catmaid(new.neuron, pid = pid, conn = conn, ...)
  }
  if(!is.neuron(old.neuron)){
    old.neuron = catmaid::read.neuron.catmaid(new.neuron, pid = pid2, conn = conn2, ...)
  }
  if(!is.null(old.neuron$connector)&length(old.neuron$connector)){
    c.df = old.neuron$connector
    near = nabor::knn(query = nat::xyzmatrix(c.df), data = nat::xyzmatrix(new.neuron), k=1, radius=ifelse(offset,search.range.nm,0))$nn.idx
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
          link_type = ifelse(new.c.df[i,"prepost"]==0,"presynaptic","postsynaptic")
          catmaid_link_connectors(treenode_id = new.c.df[i,"near_id"],
                                  connector_id = res$connector_id,
                                  link_type = link_type,
                                  pid = pid, conn = conn,...)
        }
      }
    }
  }
}

#' Get the UTC creation / edit time for a CATMAID node
#'
#' @description Get the UTC creation / edit time for a CATMAID treenode or connector. Useful for makign 'state' arguments to be passed to other functions that edit data on a CATMAID server.
#' @param id a treenode or connector ID
#' @param connector_id the connector ID to/from which a link is to be made
#' @param time whether to return the creation_time or edition_time
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_set_labels
#' @export
#' @rdname catmaid_node_time
catmaid_node_time <- function(id, time = c("creation_time", "edition_time"), pid = 1, conn = NULL, ...){
  time = match.arg(time)
  id = as.numeric(id)
  post_data = list()
  post_data["node_ids"] = id
  path = sprintf("/%d/node/user-info", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                simplifyVector = T, conn = conn, ...)
  if(!is.null(res$error)){
    stop(res$error)
  }else{
    res[[1]][[time]]
  }
}

#' Add synaptic links between connectors and treenodes
#'
#' @description Add links between connectors and treenodes to designate pre- and postsynaptic connections in CATMAID
#' @param treenode_id the treenode ID to/from which a link is to be made
#' @param connector_id the connector ID to/from which a link is to be made
#' @param from_type whether the from_id is a connector or a treenode. to_type will be set to the other type from from_id
#' @param link_type whether the link is presynaptic or postsynaptic
#' @param verbose whether or not to report if a link is successfully made
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_link_connectors
catmaid_link_connectors <- function(treenode_id, connector_id,
                                    link_type = c("presynaptic","postsynaptic"),
                                    verbose = TRUE,
                                    pid = 1, conn = NULL, ...){
  link_type = match.arg(link_type)
  link_type = paste0(link_type,"_to")
  ### Get states ###
  from_time = catmaid_node_time(id=treenode_id, time = "edition_time", pid=pid, conn=conn, ...)
  to_time = catmaid_node_time(id=connector_id, time = "edition_time", pid=pid, conn=conn, ...)
  ### Make link ###
  post_data = list()
  post_data["from_id"] = treenode_id
  post_data["to_id"] = connector_id
  post_data["link_type"] = link_type
  post_data["state"] = sprintf('[[%s, "%s"], [%s, "%s"]]',treenode_id, from_time, connector_id, to_time)
  path = sprintf("/%d/link/create", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                simplifyVector = T, conn = conn, ...)
  invisible(catmaid:::catmaid_error_check(res2))
  if(verbose){
    if(!is.null(res$error)){
      stop(res$error)
    }else{
      message(res$message, " linking treenode ", treenode_id," and connector ", connector_id)
    }
  }
}


#' Find the location of specified tags for a CATMAID neuron
#'
#' @description  Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
#' @param x a neuron or neuronlist object
#' @param tag a single character specifying which tag to look for. Defaults to TODO.
#' @param url if TRUE (defualt) a list of URLs pertaining to specified tag locations are returned. If FALSE, a data.frame subsetted from x$d is returned, reporting treenode ID and X,Y,Z positions for specified tags
#' @export
#' @rdname catmaid_get_tag
catmaid_get_tag<-function(x, tag = "TODO", url = TRUE) UseMethod("catmaid_get_tag")
catmaid_get_tag.neuron <- function(x, tag = "TODO", url = TRUE){
  TODO = unique(unlist(x$tags[[tag]]))
  if(is.null(TODO)){
    NULL
  }else{
    df = subset(x$d,PointNo%in%TODO)
    if(url){
      catmaid_url = paste0(catmaid_get_server(conn), "?pid=",pid)
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
catmaid_get_tag.neuronlist <- function(x, tag = "TODO", url = TRUE){
  if(url){
    unlist(lapply(x,catmaid_get_tag.neuron, url=url, tag= tag))
  }else{
    do.call(rbind,lapply(x,catmaid_get_tag.neuron, url=url, tag = tag))
  }
}

#' Find the TODO tagged merge sites between a given neuron and the CATMAID database at large
#'
#' @description  Find the location of tags in a CATMAID neuron, either as URLs to the location of a TODO tag in CATMAID or as a data.frame reporting the location and skeleton treenode locations of specified tags.
#' @param TODO a data.frame, as returned by catmaid_get_tag with tag = "TODO" (or a different tag indicative of a merge point) and url = FALSE
#' @param skid the skeleton ID of a neuron for which merges are to be found
#' @param fafbseg whether or not to use fafbseg::read_brainmaps_meshes on the TODO tag locations and restrict possible merges to neurons within the search range and within the cognate auto-segmented volume
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_find_likely_merge
catmaid_find_likely_merge <- function(TODO, fafbseg = FALSE, min_nodes = 2, search.range.nm = 1000, pid=1, conn = conn, ...){
  possible.merges = data.frame()
  for(i in 1:nrow(TODO)){
    todo = TODO[i,]
    skid = rownames(todo)
    bbx = rbind(todo[,c('X','Y', 'Z')]-search.range.nm,todo[,c('X','Y', 'Z')]+search.range.nm)
    skids.bbx = catmaid_skeletons_in_bbx(bbx=bbx, min_nodes=min_nodes, pid=pid, conn = conn, ...)
    neuron.bbx = catmaid::read.neurons.catmaid(skids.bbx, OmitFailures = TRUE, pid=pid, conn = conn, ...)
    neuron.bbx = neuron.bbx[setdiff(names(neuron.bbx),skid)]
      if(fafbseg){
        require(fafbseg)
        seg = tryCatch(fafbseg::brainmaps_xyz2id(todo[,c('X','Y', 'Z')]),error=function(e)NULL)
        s = fafbseg::find_merged_segments(seg)
        vol = tryCatch(fafbseg::read_brainmaps_meshes(s), error = function(e) "read error")
        in.vol = sapply(neuron.bbx,function(n) sum(nat::pointsinside(nat::xyzmatrix(n),surf=vol))>0)
        neuron.bbx = neuron.bbx[in.vol]
        message(length(neuron.bbx), " fragments found in fafb auto-traced segment corresponding with merge tag")
      }else{
        message(length(neuron.bbx), " fragments found within ", search.range.nm," nm of merge tag")
      }
      neuron.bbx.d = do.call(rbind,lapply(1:length(neuron.bbx),function(n) cbind(neuron.bbx[[n]]$d,skid=names(neuron.bbx[n]))))
      near = nabor::knn(query = todo[,c('X','Y', 'Z')], data =nat::xyzmatrix(neuron.bbx.d),k=10, radius=search.range.nm)$nn.idx
      near = near[near!=0]
      near = neuron.bbx.d[near,]
      near$TODO = todo$PointNo
      near$merger = skid
      possible.merges = rbind(possible.merges,near)
  }
  possible.merges = possible.merges[,setdiff(colnames(possible.merges),"Parent")]
  colnames(possible.merges) = c("upstream.node","label","X","Y","Z","W","upstream.skid","downstream.node","downstream.skid")
  possible.merges
}

#' Interactively choose to join neurons in CATMAID, visualising with rgl in R
#'
#' @description Interactively choose to join neurons in CATMAID via rgl in R. Use with caution in live tracing environments.
#' @param possible.merges a data frame of possible mergers between tree nodes for neurons in a CATMAID database, as returned via catmaid_find_likely_merge
#' @param downstream.neurons neurons (typically smaller, new tracings) to be joined into other (typically larger) neurons, which will be pulled from CATMAID. If NULL, the downstream neurons are also pulled from CATMAID.
#' @param brain the brain to plot while visualising potential mergers using rgl. Defaults to NULL, no brain plotted.
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_interactive_join
catmaid_interactive_join <- function(possible.merges, downstream.neurons = NULL, brain = NULL, pid = 1, conn = conn, ...){
  if(!is.null(downstream.neurons)&!nat::is.neuronlist(downstream.neurons)){
    stop("downstream.neurons must either be a neuronlist,
         or left NULL if you want to fetch them from CATMAID within this function")
  }
  neurons = catmaid::read.neurons.catmaid(skids=unique(possible.merges$upstream.skid),pid=pid,conn=conn,...)
  for(todo in 1:length(unique(possible.merges$downstream.node))){
    TODO.possible = subset(possible.merges,downstream.node==unique(possible.merges$downstream.node)[todo])
    skid = unique(TODO.possible[,"downstream.skid"])[todo]
    todo = unique(possible.merges$downstream.node)[todo]
    if(!is.null(downstream.neurons)){
      neuron = downstream.neurons[as.character(skid)]
    }else{
      neuron = catmaid::read.neurons.catmaid(skid, pid = pid, conn = conn, ...)
    }
    downstream.node = catmaid::catmaid_get_treenodes_detail(tnids = todo, pid = pid, conn = conn, ...)
    continue = FALSE
    i = 1
    while(!continue){
      message("Option ",i," of ", nrow(TODO.possible), " for ", skid)
      rgl::clear3d()
      if(!is.null(brain)){
        rgl::plot3d(brain,alpha=0.1,col="grey")
      }
      merger = catmaid::catmaid_get_treenodes_detail(tnids = TODO.possible[i,"upstream.node"], pid = pid, conn = conn, ...)
      merger.neuron = neurons[as.character(merger$skid)]
      rgl::plot3d(neuron,col="red",lwd=2,soma=T)
      rgl::plot3d(merger.neuron,col="black",lwd=2,soma=T)
      rgl::spheres3d(nat::xyzmatrix(downstream.node), col="orange", radius = 10)
      rgl::spheres3d(nat::xyzmatrix(TODO.possible[i,]), col="grey", radius = 10)
      message(paste(merger.neuron[,],collapse = " "))
      catmaid_url = paste0(catmaid_get_server(conn), "?pid=",pid)
      catmaid_url = paste0(catmaid_url, "&zp=", TODO.possible[i,"Z"])
      catmaid_url = paste0(catmaid_url, "&yp=", TODO.possible[i,"Y"])
      catmaid_url = paste0(catmaid_url, "&xp=", TODO.possible[i,"X"])
      catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
      catmaid_url = paste0(catmaid_url, "&active_skeleton_id=", skid)
      catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
      message("See merge site in CATMAID: ", catmaid_url)
      progress = readline("Make merge? y = yes, n = no, c = cycle : ")
      if(progress=="y"){
        nat::npop3d()
        rgl::plot3d(merger.neuron,col="green",lwd=2,soma=T)
        sure = readline("Sure? y = yes, n = no : ")
        if(sure=="y"){
          catmaid_join_skeletons(from_treenode_id = possible.merges[i,"downstream.node"], to_treenode_id = possible.merges[i,"upstream.node"], pid = pid, conn = conn, ...)
          continue=TRUE
        }else{
          i = ifelse(i==nrow(TODO.possible),i,i+1)
        }
      }else if(progress=="n"){
        TODO.possible = TODO.possible[-i,]
        if(nrow(TODO.possible)==0){
          continue = TRUE
        }else if (i > nrow(TODO.possible)){
          i = nrow(TODO.possible)
        }else{
          i = ifelse(i==nrow(TODO.possible),1,i+1)
        }
      }else{
        i = ifelse(i==nrow(TODO.possible),1,i+1)
      }
    }
  }
}

#' Interactively choose to join neurons in CATMAID, visualising with rgl in R
#'
#' @description Interactively choose to join neurons in CATMAID via rgl in R. Use with caution in live tracing environments.
#' @param connector_id a data frame of possible mergers between tree nodes for neurons in a CATMAID database, as returned via catmaid_find_likely_merge
#' @param node the brain to plot while visualising potential mergers using rgl. Defaults to NULL, no brain plotted.
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch and catmaid::catmaid_get_treenode_detail
#' @export
#' @rdname catmaid_connector_nodes
catmaid_connector_nodes <- function(connector_id, node = c("presynaptic","postsynaptic"),
                                            pid=1, conn = conn, ...){
  node = match.arg(node)
  connector_id = as.numeric(connector_id)
  if(length(connector_id)!=1){
    stop("connector_id must be a single connector_id")
  }
  post_data = list()
  post_data["connector_ids[0]"] = connector_id
  path = sprintf("/%d/connector/skeletons", pid)
  res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                               simplifyVector = T, conn = conn,...)
  tnids = res[[1]][[2]][[paste0(node,"_to_node")]]
  detail = catmaid::catmaid_get_treenodes_detail (tnids=tnids, pid = pid, conn = conn, ...)
  detail$connector_id = connector_id
  detail
}


fafbseg_join_connectors_in_brainmaps_volumes <- function(neurons,
                                                                       putatively.connected.skids,
                                                                       connector.range.nm = 1000,
                                                                       pid=1,conn = NULL, ...){
  require(fafbseg)
  if(length(putatively.connected.skids)>10){
    stop("A maximum 10 potentially connected skeleton IDs can be given at any one time")
  }
  if(!nat::is.neuronlist(neurons)){
    stop("Neurons must be a neuronlist object")
  }
  ### Get ngl IDs corresponding to potential presynapses ###
  connectors.pre = catmaid::catmaid_get_connector_table(skids=putatively.connected.skids, direction = "incoming",
                                                       get_partner_nodes = TRUE,
                                                       pid = pid, conn = conn, ...)
  connectors.pre = connectors.pre[!duplicated(connectors.pre$connector_id),]
  connectors.pre[is.na(connectors.pre$partner_skid)&is.na(connectors.pre$partner_nodes),"partner_nodes"] = 0
  connectors.pre = subset(connectors.pre,partner_nodes<10)
  connector.pre.segs = fafbseg::brainmaps_xyz2id(nat::xyzmatrix(connectors.pre))
  df.pre.connectors = data.frame(connector_id = connectors.pre$connector_id, connector.skid = connectors.pre$skid,
                             x = connectors.pre$x, y = connectors.pre$y, z = connectors.pre$z,
                             ngl_id = connector.pre.segs, link_type = "presynaptic", partner_nodes = connectors.pre$partner_nodes)
  ### Get ngl IDs corresponding to potential postsynapses ###
  connectors.post = catmaid::catmaid_get_connector_table(skids=putatively.connected.skids, direction = "outgoing",
                                                       get_partner_nodes = TRUE,
                                                       pid = pid, conn = conn, ...)
  connectors.post = connectors.post[!duplicated(connectors.post$connector_id),]
  connectors.post[is.na(connectors.post$partner_skid)&is.na(connectors.post$partner_nodes),"partner_nodes"] = 0
  connectors.post = subset(connectors.post,partner_nodes<10&partner_nodes>0)
  for(i in 1:nrow(connectors.post)){
    catmaid_connector_nodes(connector_id = df[i,"connector_id"],node = df[i,"link_type"],pid=pid,conn=conn,...)
  }
  connector.post.segs = fafbseg::brainmaps_xyz2id(nat::xyzmatrix(connectors.post))
  df.post.connectors = data.frame(connector_id = connectors.post$connector_id, connector.skid = connectors.post$skid,
                             x = connectors.post$x, y = connectors.post$y, z = connectors.post$z,
                             ngl_id = connector.post.segs, link_type = "postsynaptic", partner_nodes = connectors.post$partner_nodes)
  df.connectors = rbind(df.pre.connectors,df.post.connectors)
  df.connectors = subset(df.connectors,ngl_id!=0)
  ### Get ngl IDs corresponding to the skeletons ###
  segs = fafbseg::brainmaps_xyz2id(nat::xyzmatrix(neurons))
  skids = unlist(sapply(1:length(neurons),function(n) rep(names(neurons)[n],nrow(nat::xyzmatrix(neurons[[n]])))))
  df = data.frame(skid=skids,ngl_id=segs)
  df = subset(df,ngl_id!=0)
  df = merge(df,df.connectors,all.x = FALSE, all.y = FALSE)
  df = dplyr::distinct(df)
  ### For connectors and neurons in the same seg, find nearest treenode to connectors ###
  for(i in 1:nrow(df)){
    neuron = neurons[df[i,"skid"]][[1]]
    near = nabor::knn(data = nat::xyzmatrix(neuron), query = nat::xyzmatrix(df[i,]),k=1,radius = connector.range.nm)$nn.idx
    tnid = neuron$d[near,]
    ### Make links ###
    if(nrow(tnid)==1){
      if(df[i,"partner_nodes"]==0){
        catmaid_link_connectors(treenode_id = tnid, connector_id = df[i,"connector_id"],
                                link_type = df[i,"link_type"],
                                pid = pid, conn = conn, verbose =FALSE, ...)
      }else{
        catmaid_connector_nodes(connector_id = df[i,"connector_id"],node = df[i,"link_type"],
                                pid=pid,conn=conn, ...)

      }

    }
  }

}

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
#' @param skid the skeleton ID that corresponds to neuron
#' @param tolerance how many potentially duplictated nodes you will tolerate. If NULL, skeleton IDs for potentially duplictated neurons are returned.
#' @param search.range.nm determines the size of the bounding box around each node in neuron to search for a duplicated. Defaults to 1 nm
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @return If tolerance is NULL, then a list, with entires duplicated.nodes (a TRUE or FALSE for ever node in neuron, TRUE is there is another skeleton within duplication.range.nm)
#' and overlapping.skids, a list of skeleton IDs for potentially overlapping neurons. If tolerance is a numeric value between 0 and 1, TRUE or FALSe is returned.
#' @export
#' @rdname catmaid_skeletons_in_bbx
catmaid_duplicated <- function(neuron, skid = 0, tolerance = NULL, duplication.range.nm = 10, pid = 1, conn = NULL, ...){
  duplicated = c()
  skids = list()
  pb <- utils::txtProgressBar(min = 0, max = nrow(neuron$d), style = 3)
  for(i in 1:nrow(neuron$d)){
    bbx = rbind(nat::xyzmatrix(neuron$d[i,])-duplication.range.nm, nat::xyzmatrix(neuron$d[i,])+duplication.range.nm)
    skids.bbx = catmaid_skeletons_in_bbx(bbx=bbx, min_nodes = 2, pid = pid, conn = conn, ...)
    skids.bbx = setdiff(skids.bbx,as.integer(skid))
    duplicated = c(duplicated,length(skids.bbx)>0)
    skids[[i]] = unique(skids.bbx)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  if(is.null(tolerance)){
    list(duplicated.nodes = duplicated, overlapping.skids = skids)
  }else{
    calc = (sum(duplicated)/nrow(duplicated))
    ifelse(length(calc)>0,calc,0) > tolerance
  }
}

#' Function used to delete single neurons from CATMAID based on a given skeleton ID
#'
#' @description  Delete a single neuron from a given CATMAID instance.
#' Deletes a neuron if and only if two things are the case: 1. The user
#' owns all treenodes of the skeleton modeling the neuron in question and
#' 2. The neuron is not annotated by other users.
#' Use with extreme caution as you may be significantly affecting others' work.
#' You will first be shown the neuron you want to delete and who has worked on it, and then asked whether or not you want to continue.
#' @param skid the skeleton ID for a single neuron you want to delete
#' @param delete_ids CATMAID connector IDs to be deleted, as long as they do not connect to any skeletons
#' @param delete_connectors if TRUE, connectors attached to the deleted skeleton, and no other skeleton (even if it is a single node) are also deleted
#' @param plot whether or not the plot the neuron you are considering deleting, before deciding to end it
#' @param brain the brain you want to plot alongside the neuron you are considering deleting.
#' @param max.nodes the maximum number of nodes a neuron can have, and still be successfully deleted. Helps prevent the accidental deletion of large neurons
#' control if TRUE (default) you will be asked at each stage if you wish to continue. Helps prevent accidental deletions. If FALSE, then will delete without asking.
#' Really should only be set to FALSE if you are deleting from your own CATMAID instance
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_delete
catmaid_delete_neuron <- function(skid,
                                  delete_connectors = TRUE,
                                  plot = TRUE, brain = NULL,
                                  max.nodes = 100, control = TRUE,
                                  pid = 1, conn = NULL, ...){
  if(length(skid)>1){
    stop("length(skid)>1 - You may only attempt to delete one neuron at a time.")
  }
  message("You are about to delete a NEURON from a CATMAID instance.
          This may affect people's hard work.
          You will first be shown the neuron you want to delete and who has worked on it,
          and then asked whether or not you want to continue.")
  if(control){
    progress = readline(paste0("Do you want to continue with: ", skid," in CATMAID instance ", server, " y=yes, n=no : "))
  }else{
    progress = "y"
  }
  if(progress=="y"){
    nid = catmaid_get_neuronid(skids = skid, pid = pid, conn = conn, ...)
    neuron = tryCatch(catmaid::read.neurons.catmaid(skid, pid = pid, conn = conn, ...), error = function(e) "ERROR")
    if(neuron=="ERROR"){
     warning("Neuron with skeleton ID ", skid, " could not be read with pid ", pid, " in ", catmaid_get_server(conn))
    }else{
      connector_ids = unique(catmaid::connectors(neuron)$connector_id)
      xyz = nat::xyzmatrix(neuron)
      meta = summary(neuron)
      if(meta$nodes>max.nodes){
        stop("The neuron flagged for deletion has ", meta$nodes,
             " which is greater than the given argument max.nodes, ", max.nodes, ". Aborting.")
      }
      catmaid_url = paste0(catmaid_get_server(conn), "?pid=",pid)
      catmaid_url = paste0(catmaid_url, "&zp=", xyz[1,"Z"])
      catmaid_url = paste0(catmaid_url, "&yp=", xyz[1,"Y"])
      catmaid_url = paste0(catmaid_url, "&xp=", xyz[1,"X"])
      catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
      catmaid_url = paste0(catmaid_url, "&active_skeleton_id=", skid)
      catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
      ### Plot ###
      if(plot){
        rgl::clear3d()
        if(!is.null(brain)){
          rgl::plot3d(brain,alpha=0.1,col="grey")
        }
        rgl::plot3d(neuron,lwd=3,soma=T,WithConnectors = TRUE, col ="black")
      }
      ### Decide ###
      print(cbind(neuron[,], meta))
      print(paste0("See neuron in CATMAID at: ", catmaid_url))
      decide = "n"
      if(control){
        decide = readline(paste0("Are you SURE you want to delete: ", skid," in CATMAID instance ", server, "? y=yes, n=no : "))
      }else{
        decide = "y"
      }
      if(decide=="y"){
        path = sprintf("/%d/neuron/%d/delete", pid, nid)
        res = catmaid_fetch(path, body = NULL, include_headers = F,
                            conn = conn, ...)
        if(!is.null(res$success)){
          message(res$success)
          if(delete_connectors&!is.null(connector_ids)){
            message("Now deleting free connectors")
            catmaid_delete_connectors(connector_ids = connector_ids, pid = pid, conn = conn, ...)
          }
        }else{
          message(res$error)
        }
        message("This activity has been logged.")
      }
    }
  }
}
#' @export
#' @rdname catmaid_delete
catmaid_delete_connectors <- function(connector_ids, pid = pid, conn = conn, ...){
  delete = c()
  for(i in 1:length(connector_ids)){
    post_data = list()
    post_data["connector_ids[0]"] = as.numeric(connector_ids[i])
    path = sprintf("/%d/connector/skeletons", pid)
    res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                 simplifyVector = T, conn = conn,...)
    invisible(catmaid:::catmaid_error_check(res2))
    del = ifelse(length(res),FALSE,TRUE)
    delete = c(delete, del)
  }
  connector_ids_delete = connector_ids[delete]
  if(length(connector_ids_delete)){
    for(connector_id_delete in connector_ids_delete){
      ctime = catmaid_node_time(id=connector_id_delete, time = "edition_time", pid=pid, conn=conn, ...)
      post_data = list()
      post_data["connector_id"] = connector_id_delete
      post_data["state"] =  sprintf('{"edition_time":"%s","c_links": []}',ctime)
      path = sprintf("/%d/connector/delete", pid)
      res = catmaid::catmaid_fetch(path, body = post_data, include_headers = F,
                                 simplifyVector = T, conn = conn,...)
      invisible(catmaid:::catmaid_error_check(res2))
      if(is.null(res$error)){
        message(res$message)
      }else{
        warning(res$error)
      }
    }
  }
}
#' @export
#' @rdname catmaid_delete
catmaid_delete_neurons <- function(skid,
                                  delete_connectors = TRUE,
                                  plot = TRUE, brain = NULL,
                                  max.nodes = 100,
                                  pid = 1, conn = NULL, ...){
  deletions = sapply(skid,catmaid_delete_neuron,delete_connectors=delete_connectors,plot=plot,brain=brain,max.nodes=max.nodes,pid=pid,conn=conn,...)
}

#' Get the CATMAID neuron ID that corresponds to the skeleton ID
#'
#' @description Retrieve the neuron IDs for given skeleton IDs. This is typically the skeleton ID + 1, and is often, but not always accurately, kept by CATMAID tracers in the name of a neuron.
#' @param skids a vector of skeleton IDs or argument applicable to catmaid::catmaid_get_neuronid
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_get_neuronid
catmaid_get_neuronid <- function(skids, pid = 1, conn = NULL, ...){
  skids = catmaid_skids(skids, conn = conn, pid = pid, ...)
  if (any(duplicated(skids))) {
    uskids = unique(skids)
    unids = catmaid_get_neuronid(uskids, pid = pid, conn = conn,...)
    res = unids[match(skids, uskids)]
    return(res)
  }
  skids[is.na(skids)] = -1L
  res = lapply(skids,function(skid)
    catmaid::catmaid_fetch(sprintf("/%d/skeleton/%s/neuronname", pid, skid), body = NULL, include_headers = F,
                           conn = conn, ...)$neuronid)
  res = sapply(res,function(r) ifelse(is.null(r),NA,r))
  names(res) = skids
  res
}

#' Lock or unlock a CATMAID neuron reconstruction
#'
#' @description  Lock or unlock a CATMAID neuron reconstruction by adding or removing a 'locked' annotation to a set of skeleton IDs (skids). A locked neuron cannot be edited until it is unlocked.
#' @param skids the skeleton IDs neurons you wish to lock / unlock
#' @param server the CATMAID server from which to delete the neuron identified by skid
#' @param pid project id. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param ... methods passed to catmaid::catmaid_fetch
#' @export
#' @rdname catmaid_lock_neurons
catmaid_lock_neurons <- function(skids, pid = 1, conn = NULL, ...){
  skids = catmaid::catmaid_skids(skid, pid=pid,conn=conn,...)
  catmaid::catmaid_set_annotations_for_skeletons(skids, annotations = "locked", pid = pid,
                                        conn = conn, ...)
}
#' @export
#' @rdname catmaid_lock_neurons
catmaid_unlock_neurons <- function(skids, pid = 1, conn = NULL, ...){
  skids = catmaid::catmaid_skids(skid, pid=pid,conn=conn,...)
  catmaid::catmaid_remove_annotations_for_skeletons(skids, annotations = "locked", pid = pid,
                                                 conn = conn, ...)
}


# A helper function, not exported
catmaid_convert_time <- function(utc){
  t = format(as.POSIXlt(utc,tz="GMT",origin="1970-01-01"), "%Y-%m-%d %H:%M:%OS3")
  s = unlist(strsplit(t," "))
  t = paste0(s[1],"T",s[2],"Z")
}

#' Upload neuron(s) to CATMAID
#'
#' @description  Uploads neurons to CATMAID, names them and annotates them, from the environment specified with conn2 to that specified by conn2.
#' Please use with caution, as you could be heavily adding to a live tracing environment.
#' @param x either skeleton IDs in the environment specified by conn2 (by default, the v14-seg CATMAID instance), or a neuronlist object to upload to the CATMAID instance specified when you use catmaid::catmaid_login(), if conn is NULL, else specified by conn$server
#' @param tolerance how many potentially duplictated nodes you will tolerate
#' @param name whatever you want to name your uploaded neurons. If a single character, then it will be added to all uploaded neurons. Else, can be a character vector the same length as swc.
#' @param annotations a character vector of annotations, to be added to all of the uploaded neurons
#' @param include.tags whether of not to transfer each neuron's tags to its newly uploaded cognate
#' @param include.connectors whether of not to transfer each neuron's connectors to its newly uploaded cognate
#' @param join whether or not to attempt to join each uploaded neuron to neurons within the new CATMAID instance, based on the placement of join.tags specified using the next argument
#' @param join.tag a single character specifying a tag that has been used to signify a potential merge point
#' @param min_nodes the minimum number of nodes a potential merger skeleton needs to have
#' @param search.range.nm the maimum distance from which to search from the TODO point tag to find potential mergers
#' @param fafbseg whether or not to use fafbseg::read_brainmaps_meshes on the TODO tag locations and restrict possible merges to neurons within the search range and within the cognate auto-segmented volume
#' @param brain the brain to plot while visualising potential mergers using rgl. Defaults to NULL, no brain plotted.
#' @param pid project id for conn. Defaults to 1
#' @param pid project id for conn2. Defaults to 1
#' @param conn CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param conn2 CATMAID connection object, see ?catmaid::catmaid_login for details
#' @param max.upload the maximum number of files that the function will allow you to upload at once
#' @param ... methods passed to catmaid::catmaid_fetch #' @examples
#' @examples
#' # This function, first, gets the neurons we want from the v14-seg CATMAID instance
#' # Then it checks that we have not already uploaded to v14-seg using the annotations you specify with the avoid argument
#' # Then it will seek to upload them in a controlled way, that gives us options to have their tags, connectors, and make joins, as well as trying to check for possible duplication.
#' # Be aware that the interactive join will often suggest the same neuron multiple times at different join sites.
#' # The join sites are indicated by spheres - you may need to zoom in in order to see them.
#' uploaded = catmaid_controlled_upload(x = "name:ASB Tester", join = TRUE, name = "ASB Tester from v14-seg",
#' search.range.nm = 1000, annotations = "ASB Test v14-seg Upload", brain = elmr::FAFB14)
#' # Note that CATMAID links will also be supplied, so you can inspect a merge site in CATMAID. If you join in CATMAID, do not join in the interactive window, ad this will throw an errror, just keep hitting 'n' for no, until all options are exhausted.
#' # let's lock the neurons we just uploaded, to lessen the chance someone else will connect stuff to them and want to upload them AGAIN later...
#' catmaid_lock_neurons(skids = uploaded$downloaded.skids, conn = fafb_seg_conn())
#' # Oops, did you make a mistake in uploading this neuron?
#' delete.skids = catmaid::catmaid_skids("annotation:ASB Tester from v14-seg")
#' catmaid_delete_neurons(delete.skids)
#' # Be careful wen deleting, especially if you have merged your fragment during the interactive join.
#' # Phew.
#' @export
#' @rdname catmaid_controlled_upload
catmaid_controlled_upload <- function(x, tolerance = 0.5, name = "v14-seg neuron upload",
                                      annotations = "v14-seg upload", avoid = "v14",
                                      include.tags = TRUE, include.connectors = TRUE,
                                      search.range.nm = 1000, join = FALSE, join.tag = "TODO",
                                      fafbseg = FALSE, min_nodes = 2,
                                      brain = NULL, return.uploaded.skids = TRUE,
                                      pid = 1, conn = NULL, pid2 = 1, conn2 = fafb_seg_conn(),  ...){
  if(!is.neuronlist(x)){
    neurons = catmaid::read.neurons.catmaid(x,pid=pid2,conn=conn2, ...)
    anns = catmaid::catmaid_get_annotations_for_skeletons(x,pid=pid2,conn=conn2,...)
    already.there = subset(anns,annotation%in%avoid)$skid
    neurons = neurons[setdiff(names(neurons),already.there)]
    if(length(already.there)){
      warning(length(already.there), " neurons have an annotation that indicates that they should not be uploaded: ", avoid)
    }
  }else{
    neurons = x
  }
  if(length(name)!=1&length(name)!=length(neurons)){
    stop("The name argument must be a chracter vector the same length as neurons,
         or if all uploaded neurons are to be named the same, a character vector of length one")
  }
  uploaded.new = c()
  uploaded.old = c()
  nat::nopen3d()
  for(i in 1:length(neurons)){
    neuron = neurons[[i]]
    old.skid = neurons[i,"skid"]
    if(length(name)==1){
      nam = names
    }else{
      nam = names[i]
    }
    rgl::clear3d()
    if(!is.null(brain)){
      rgl::plot3d(brain,alpha=0.1,col="grey")
    }
    rgl::plot3d(neuron,WithConnectors = TRUE)
    message("Assessing whether upload of neuron ",i, " may cause a duplication:" )
    dupe = catmaid_duplicated(neuron, tolerance = NULL, pid = pid, conn = conn, ...)
    calc = (sum(dupe$duplicated.nodes)/nrow(dupe$duplicated.nodes))
    calc = ifelse(length(calc)>0,calc,0)
    dupe = calc > tolerance
    message("The neuron flagged for upload seems to have ", calc*100, "% of its nodes already in the targetted CATMAID instance")
    progress = "n"
    progress = readline("Do you want to upload this neuron? y=yes, n=no: ")
    if(dupe){
      warning("Neuron ", i, " with skid ", old.skid, " appears to already exist in the CATMAID instance to which you are seeking to upload.
              Upload for neuron ", i, " aborted.")
      next
    }else if(progress=="y"){
      message("Uploading neuron ", i)
      new.skid = catmaid_upload_neurons(swc=neuron,name=name,annotations= c(annotations, paste0("v14-seg: ",old.skid), avoid),
                                        include.tags=include.tags,include.connectors=include.connectors,
                                        search.range.nm=10, return.new.skids = TRUE,
                                        conn = conn, pid = pid, max.upload = 1, ...)
      new.neuron = read.neurons.catmaid(new.skid, conn=conn, pid=pid, ...)
      message("Upload successful, neuron ", new.skid, " created, named: ", new.neuron[,"name"])
      uploaded.new = c(uploaded.new,new.skid)
      uploaded.old = c(uploaded.old,old.skid)
      if(join){
        message("Finding tagged potential join points ...")
        TODO = catmaid_get_tag(x = new.neuron, tag = join.tag, url = FALSE)
        if(nrow(TODO)>0){
          possible.merges = catmaid_find_likely_merge(TODO = TODO, pid=pid, fafbseg = fafbseg, conn = conn,
                                                       min_nodes = min_nodes, search.range.nm = search.range.nm, ...)
            if(nrow(possible.merges)>1){
              message("Choose join sites interactively: ")
              catmaid_interactive_join(possible.merges=possible.merges, downstream.neurons = new.neuron, brain = brain, pid = pid, conn = conn, ...)
            }else{
              message("No potential join sites found for new neuron ", new.skid)
              next
            }
        }else{
          message("No potential join sites found for new neuron ", new.skid)
          next
        }
      }
    }
  }
  catmaid::catmaid_set_annotations_for_skeletons(skids = uploaded.old, annotations = avoid, pid = pid2, conn = conn2, ...)
  if(return.uploaded.skids){
    list(uploaded.skids = uploaded.new, downloaded.skids = uploaded.old)
  }
}

# Hidden function to connect to a local CATMAID server
local_conn <- function(){
  catmaid::catmaid_login(server = "http://localhost:8000/",
                token = "5c93cd0d5a75427aac7d1e39f3deb0cd59de19e7",
                Cache = FALSE)
}

# Get CATMAID server
catmaid_get_server<-function(conn=NULL,...){
  if(is.null(conn)){
    conn = catmaid::catmaid_login()
  }
  conn$server
}


# Test these functions
#frags = read.neurons.fafbseg("annotation:ASB downseg")
uploaded = catmaid_upload_neurons(swc = frags[398:length(frags)], name ="neuron SWC upload", annotations = NULL,
                                  include.tags = TRUE,include.connectors = TRUE,
                                  search.range.nm = 100, return.new.skids = TRUE,
                                  pid = 4, conn = local_conn(), max.upload = 1000)
dels1 = catmaid_skids(x = "annotation:SWC upload",conn=local_conn(),pid=4)
dels2 = catmaid_skids(x = "annotation:neuron upload",conn=local_conn(),pid=4)
dels = c(dels1,dels2)
catmaid_delete_neurons(skid=dels,conn=local_conn(),pid=4, max.nodes = 5000)



