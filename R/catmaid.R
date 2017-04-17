#' Resamle a CATMAID neuron
#'
#' @description Resample a catmaid neuron so that connector information is retained
#'
#' @param someneuronlist a neuronlist or neuron object
#' @param neuron A neuron object
#' @param stepsize The new spacing along the tracing
#' @param ... additional arguments passed to methods.
#' @export
#' @alias resample
#' @importFrom nat resample
resample.catmaidneuron<-function(neuron,stepsize=1){
  r = NextMethod()
  c = connectors(neuron)
  c$treenode_id = nabor::knn(
    data = nat::xyzmatrix(r),
    query = nat::xyzmatrix(c),
    k = 1
  )$nn.idx
  r$connectors = c
  r
}

#' @export
#' @alias resample
#' @importFrom nat resample
resample.catmaidneuronlist<-function(someneuronlist,stepsize=1){
  nl = nlapply(someneuronlist,resample.catmaid.neuron,stepsize=stepsize)
  nl
}

# Some extra catmaid related functions
#open_catmaid<-function(x, s=select3d(), mirror=FALSE, sample=FCWB, scale=1) {
#  xyz=xyzmatrix(x)
  # calculate centroid
#  cent=colMeans(xyz[s(xyz),, drop=F])
#  cent=matrix(cent, ncol=3)*scale

 # xyzi=as.integer(cent)
#  url=sprintf("https://neurocean.janelia.org/catmaidL1/?pid=1&zp=%d&yp=%d&xp=%d&tool=tracingtool&sid0=1&s0=1",
 #             xyzi[3], xyzi[2], xyzi[1])
  #browseURL(url)
#}


#unique.connections <- function(someneuronlist, anotherneuronlist, direction, omit = NULL, minimum_synapses = 2, min_nodes = 1000){
#  results = list()
#  reverse = 2
#  if (direction == 2) { reverse = 1} # 1 is outgoing, 2 is incoming.
#  for (n in 1:length(someneuronlist)){
#    print(n)
#    sneuron = someneuronlist[n]
#    hits = c(names(sneuron))
#    tryCatch({cn.sneuron = catmaid_query_connected(names(sneuron), minimum_synapses = minimum_synapses, boolean_op = "OR")}, error=function(e){cn.sneuron = NULL})
#    sneuron.inputs = (unique(subset(cn.sneuron[[direction]]$partner, !cn.sneuron[[direction]]$partner%in%names(omit) & cn.sneuron[[direction]]$num_nodes>min_nodes)))
#    if (is.null(cn.sneuron) == F){
#      for (neuron in sneuron.inputs){
#        tryCatch({cn.targeter = catmaid_query_connected(neuron, minimum_synapses = minimum_synapses, boolean_op = "OR")}, error=function(e){cn.targeter = NULL})
#        if (is.null(cn.targeter) == F){
#          if (all(!cn.targeter[[reverse]]$partner%in%names(anotherneuronlist[-n])) == T){
#            hits = c(hits, neuron)
#          }
#        }
#      }
#    }
#    if (length(hits) > 0){results[[length(results)+1]] <- hits}
#  }
#  return(results)
#}

#connectivity <- function(skids, min_nodes = 1000, min_synapses = 1, direction = 1, names = F){
#  skids=catmaid_skids(skids, conn = conn)
#  cn = catmaid_query_connected(skids)
#  cn = subset(cn[[direction]], cn[[direction]]$partner%in%skids & cn[[direction]]$num_nodes >= min_nodes & cn[[direction]]$syn.count >= min_synapses)
#  m = matrix(0,nrow = length(skids), ncol = length(skids))
#  rownames(m) <- colnames(m) <- skids
#  for (skid in as.character(skids)){
#    for (partner in as.character(unique(cn$partner))){
#      c = subset(cn$syn.count, cn$partner == partner & cn$skid == skid)
#      m[skid,partner] <- ifelse(length(c) > 0, c, 0)
#    }
#  }
#  if (names == T){rownames(m) <- colnames(m) <- catmaid_get_neuronnames(skids)}
#  blacklist = c(names(which(rowSums(m) == 0)), names(which(colSums(m) == 0)))
#  blacklist = blacklist[duplicated(blacklist)]
#  m = m[!rownames(m)%in%blacklist,!colnames(m)%in%blacklist]
#  return(m)
#}


#' Meta-annotate CATMAID annotations
#'
#' @description Meta-annotate a grops of CATMAID annotations
#'
#' @param annotations annotations to meta-annotate
#' @param meta_annotations meta-annotation to add
#' @param conn a catmaid_connection objection returned by catmaid_login. If NULL (the default) a new connection object will be generated using the values of the catmaid.* package options as described in the help for catmaid_login
#' @param pid project id (default 1)
#' @param ... additional arguments passed to methods.
#'
#' @export
#' @rdname catmaid_set_meta_annotations
catmaid_set_meta_annotations<-function(meta_annotations,annotations,pid=1,conn=NULL,...){
  post_data = list()
  post_data[sprintf("annotates[%d]", seq_along(annotations))] = as.list(annotations)
  path = sprintf("/%d/annotations/add", pid)
  post_data[sprintf("annotations[%d]", seq_along(meta_annotations))] = as.list(meta_annotations)
  res = catmaid_fetch(path, body = post_data, include_headers = F,
                      simplifyVector = T, conn = conn, ...)
}

#' Annotate CATMAID partners
#'
#' @description Intended to use to annotate CATMAID left-right cognates, and fetch them
#'
#' @param partners a vector of two left-right cognates
#' @param conn a catmaid_connection objection returned by catmaid_login. If NULL (the default) a new connection object will be generated using the values of the catmaid.* package options as described in the help for catmaid_login
#' @param pid project id (default 1)
#' @param skids CATMAID skeleton IDs
#' @param name a vector of neuron names
#' @param ... additional arguments passed to methods.
#'
#' @export
#' @rdname catmaid_annotate_partners
catmaid_annotate_partners<-function(partners,pid=1,conn=NULL,...){
  if(is.vector(partners)){
    if (length(partners)!=2){
      message("Too many skids supplied")
      break
    }else{
      an = paste0("paired with #",partners[1])
      catmaid_set_annotations_for_skeletons(partners[2],annotations = an,pid=pid,conn=conn)
      an = paste0("paired with #",partners[2])
      catmaid_set_annotations_for_skeletons(partners[1],annotations = an,pid=pid,conn=conn)
    }
  }
}

#' @rdname catmaid_annotate_partners
catmaid_get_annotated_partner<-function(skids,pid=1,conn=NULL,...){
  an = sapply(skids,function(x) paste0("paired with #",x))
  ids = catmaid_query_by_annotation(an)$skid
  read.neurons.catmaid(unique(ids))
}

#' @export
#' @rdname catmaid_annotate_partners
reverse.name.side<-function(names){
  for (n in names){
    if (grepl("left|Left|_l|L$",n)){
      left.sub = gsub("left","right",n)
      left.sub = gsub("Left","Right",left.sub)
      left.sub = gsub("_l","_r",left.sub)
      return(gsub("L$","R",left.sub))
    }else if (grepl("right|Right|_r|R$",n)){
      right.sub = gsub("right","left",n)
      right.sub = gsub("Right","Left",right.sub)
      right.sub = gsub("_r","_l",right.sub)
      return(gsub("R$","L",right.sub))
    }
  }
}


#' Update a local neuronlist with new CATMAID data
#'
#' @description Use to update a large neuronlist quickly by pulling just certain neurons from CATMAID
#'
#' @param skids sekeleton IDs
#' @param someneuronlist a neuronlist object
#' @param ... additional arguments passed to read.neurons.catmaid
#'
#' @export
#' @rdname update.neuronlist
update.neuronlist<-function(someneuronlist,skids,...){
  someneuronlist[skids] = read.neurons.catmaid(skids)
  someneuronlist
}

#' @export
#' @rdname update.neuronlist
update.names.neuronlist<-function(someneuronlist,skids = names(someneuronlist),...){
  someneuronlist[,"name"] = catmaid_get_neuronnames(skids)
  someneuronlist
}
