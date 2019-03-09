#' Get data and plots describing tracer activity statistics
#'
#' @description Some functions that quickly visualise tracer effort, or retrieve data to do with tracer effort.
#'
#' @param skids The skids of the neurons to consider
#' @param names Full names of tracers
#' @param calc.date The date from which to retrieve data for the graphs
#' @param from.date The date from which to plot the graphs
#' @param to.date The date to which to plot the graphs
#' @param cumulative Whether to generate a cumulative plot. Defaults to true.
#' @param value Whether to look at tree node, pre synapse connector or post synaptic node data.
#' @param ... additional arguments passed to methods.
#'
#' @return tracer.plot() returns a plot for cable length traced and connectors placed by named tracers. tracer.treenodes.plot() visualises data on node generation.
#' @export
#' @rdname tracer.plot
#' @importFrom dplyr filter
tracer.plot <- function(names = c("Alex Bates","Ruairi Roberts", "Greg Jefferis", "Clement Hallou", "Philipp Ranft", "Philipp Schlegel", "Fiona Love", "Amelia Edmondson-Stait"), calc.date = "2016-04-01", from.date = "2016-04-01", to.date = Sys.Date(), cumulative = T, ...){
  if(!requireNamespace('ggplot2', quietly = TRUE))
    stop("You must install the suggested ggplot2 package to use this function!")
  if(!requireNamespace('easyGgplot2', quietly = TRUE))
    stop("You must install the suggested easyGgplot2 package to use this function!")

  dataf = catmaid::catmaid_user_history(from = calc.date)
  dataf[is.na(dataf)] = 0
  dataf <- dataf %>%
    dplyr::group_by(full_name) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(cable.cum.sum = cumsum(new_cable), connectors.cum.sum = cumsum(new_connectors))
  if (!cumulative){
    p1 <- ggplot2::qplot(date, new_cable, col=full_name,
              data=filter(dataf, full_name%in%names),  ylim=c(0, max(filter(dataf, full_name%in%names)$new_cable)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      ggplot2::geom_point()+
      ggplot2::geom_smooth()+
      ggplot2::theme(legend.position="none")
    p2 <- ggplot2::qplot(date, new_connectors, col=full_name,
              data=filter(dataf, full_name%in%names),  ylim=c(0, max(filter(dataf, full_name%in%names)$new_connectors)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      ggplot2::geom_point()+
      ggplot2::geom_smooth()+
      ggplot2::theme(legend.position="top")
  }else{
    p1 <- ggplot2::qplot(date, cable.cum.sum, col=full_name,
                data=filter(dataf, full_name%in%names),  ylim=c(0, max(filter(dataf, full_name%in%names)$cable.cum.sum)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      ggplot2::geom_point()+
      ggplot2::geom_path()+
      ggplot2::theme(legend.position="none")
    p2 <- ggplot2::qplot(date, connectors.cum.sum, col=full_name,
                data=filter(dataf, full_name%in%names),  ylim=c(0, max(filter(dataf, full_name%in%names)$connectors.cum.sum)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      ggplot2::geom_point()+
      ggplot2::geom_path()+
      ggplot2::theme(legend.position="top")
  }
  easyGgplot2::ggplot2.multiplot(p2,p1,cols=1)
}

#' @export
#' @rdname tracer.plot
tracer.treenodes.plot <- function(skids = NULL, names = c("Alex Bates","Ruairi Roberts", "Greg Jefferis", "Clement Hallou", "Philipp Ranft", "Philipp Schlegel", "Fiona Love", "Amelia Edmondson-Stait"), from.date = "2016-04-01", to.date = Sys.Date(), ...){
  if (is.null(skids)){
    skids=as.integer(catmaid::catmaid_fetch(paste("/1/skeletons/?nodecount_gt=1")))
  }
  treenodes=nat::nlapply(skids, catmaid_get_treenode_table, OmitFailures = T)
  ul=catmaid::catmaid_get_user_list()
  #as.Date(as.POSIXct(z, origin = "1970-01-01"))
  if (is.null(names)){
    names = ul$full_name
  }
  user_ids = subset(ul, full_name%in%names)$id
  treenodes %>%
    dplyr::bind_rows() %>%
    dplyr::filter(user_id%in%user_ids) %>%
    dplyr::arrange(last_modified) %>%
    dplyr::mutate(last_modified = as.Date(as.POSIXct(last_modified, origin = "1970-01-01"))) %>%
    dplyr::select(-one_of(c("id","parent_id","confidence","x","y","z","r", "reviewer_id")))%>%
    dplyr::mutate(node_count = 1) %>%
    dplyr::group_by(user_id, last_modified) %>%
    dplyr::summarize(new_treenodes = n()) %>%
    dplyr::group_by(user_id) %>%
    dplyr::mutate(full_name = unique(subset(ul$full_name, ul$id==user_id))) %>%
    dplyr::mutate(treenode.cum.sum = cumsum(new_treenodes))->
    treenodes.names
  treenodes.names$full_name = ul$full_name[match(treenodes.names$user_id, ul$id)]
  p1 <- ggplot2::qplot(last_modified, new_treenodes, col=full_name,
                       data=treenodes.names,  ylim=c(0, max(treenodes.names$new_treenodes)), xlim = c(as.Date(from.date), as.Date(to.date)))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::theme(legend.position="none")
  p2 <- ggplot2::qplot(last_modified, treenode.cum.sum, col=full_name,
                data=treenodes.names,  ylim=c(0, max(treenodes.names$treenode.cum.sum)), xlim = c(as.Date(from.date), as.Date(to.date)))+
    ggplot2::geom_point()+
    ggplot2::geom_path()+
    ggplot2::theme(legend.position="top")
  easyGgplot2::ggplot2.multiplot(p2,p1,cols=1)
}

#' @export
#' @rdname tracer.plot
tracer.neuron.stats <- function(skids, value = c("nodes","pre","post"), ...){
  neuronds.stats=catmaid::catmaid_get_contributor_stats(skids)
  ul=catmaid::catmaid_get_user_list()
  uls=ul[,c('full_name','id')]
  neuronds.stats$node_contributors
  if (value == "nodes"){
    dplyr::inner_join(neuronds.stats$node_contributors, uls, by='id') %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(pct=n/sum(n)*100, cpct=cumsum(pct))
  }else if (value == "pre"){
    dplyr::inner_join(neuronds.stats$pre_contributors, uls, by='id') %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(pct=n/sum(n)*100)
  }else if (value == "post"){
    dplyr::inner_join(neuronds.stats$post_contributors, uls, by='id') %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(pct=n/sum(n)*100)
  }
}

#applyTransform.neuron <- function(neuron, trafo, inverse = F){
#  xyzmatrix(neuron$d)<-Morpho::applyTransform(xyzmatrix(neuron$d), trafo = trafo, inverse = inverse)
#  neuron
#}

#applyTransform.neuronlist <- function(someneuronlist, trafo, inverse = F){
#  nlapply(someneuronlist, applyTransform.neuron, trafo, inverse)
#}


#' Get data tracers' contributions to skeletons
#'
#' @description Get data tracers' contributions to ndoes of skeletons of CATMAID neurons, where the skeletons have been edited e.g. branches removed or added
#'
#' @param x a neuronlist where the name of each neuron is its skeleton ID
#' @param exclude.authors a vector of full author names to exclude from result. The result will then print the contributions of non-authors that meet the threshold specified by arguments to suggest_authorship
#' @param direction for synapse.contribution.to.skeletons, whether to look at upstream or downstream connections
#' @param ... argument supplied to suggest_authorship
#' @export
#' @rdname contribution.to.skeletons
#' @importFrom dplyr mutate_ arrange_
node.contribution.to.skeletons <-function(x, exclude.authors = NULL, ...){
  if(!requireNamespace('elmr', quietly = TRUE))
    stop("Please install suggested package elmr to use this function!")

  skids = names(x)
  treenodes=nat::nlapply(skids, catmaid_get_treenode_table, OmitFailures = T)
  treenodes = do.call(rbind,treenodes)
  ul=catmaid::catmaid_get_user_list()
  uls = ul[, c("full_name", "id")]
  nodes = unlist(sapply(x,function(x) x$d$PointNo))
  trimmed.treenodes = subset(treenodes, id%in%nodes)
  t = sort(table(trimmed.treenodes$user_id),decreasing = TRUE)
  fns = uls[match(as.numeric(names(t)),uls$id),]$full_name
  df = data.frame(id = names(t), n = c(t), full_name = fns)
  df <- df %>% dplyr::arrange_(~desc(n)) %>%
    dplyr::mutate_(pct = ~n/sum(n) * 100, cpct = ~cumsum(pct))
  attr(df,"type") = "nodes"
  if(!is.null(exclude.authors)){
    s = df[!df$full_name%in%exclude.authors,]
    s = elmr::suggest_authorship(s, ...)
    s$action = "ack"
    elmr::write_ack(s)
  }
  df
}

#' @export
#' @rdname contribution.to.skeletons
synapse.contribution.to.skeletons <-function(x, exclude.authors = NULL, direction = NULL){
  skids = names(x)
  c = do.call(rbind,lapply(x,function(x) x$connectors))
  connectors = do.call(rbind,lapply(skids,catmaid_get_connector_table))
  if(!is.null(direction)){
    connectors = subset(connectors,direction==direction)
  }
  ul=catmaid::catmaid_get_user_list()
  uls = ul[, c("full_name", "id")]
  cids = unlist(sapply(x,function(x) x$connectors$connector_id))
  trimmed.connectors = subset(connectors, connector_id%in%cids)
  # Remove duplicated connectors
  trimmed.connectors = trimmed.connectors[!duplicated(trimmed.connectors$connector_id),]
  t = sort(table(trimmed.connectors$user_id),decreasing = TRUE)
  fns = uls[match(as.numeric(names(t)),uls$id),]$full_name
  df = data.frame(id = names(t), n = c(t), full_name = fns)
  df <- df %>% dplyr::arrange_(~desc(n)) %>%
    dplyr::mutate_(pct = ~n/sum(n) * 100, cpct = ~cumsum(pct))
  attr(df,"type") = "synapses"
  if(!is.null(exclude.authors)){
    s = df[!df$full_name%in%exclude.authors,]
    s = elmr::suggest_authorship(s,..)
    s$action = "ack"
    write_ack(s)
  }
  df
}







