#' Get data and plots describing tracer activity statistics
#'
#' @description Some functions that quickly visualise tracer effort, or retrieve data to do with tracer effort.
#'
#' @param skids The skids of the neurons to consider
#' @param names Full names of tracers
#' @param clac.date The date from which to retrieve data for the graphs
#' @param from.date The date from which to plot the graphs
#' @param to.date The date to which to plot the graphs
#' @param Cumulative Whether to generate a cumulative plot. Defaults to true.
#' @param value Whether to look at tree node, pre synapse connector or post synaptic node data.
#' @param ... additional arguments passed to methods.
#'
#' @return tracer.plot() returns a plot for cable length traced and connectors placed by named tracers. tracer.treenodes.plot() visualises data on node generation.
#' @export
#' @rdname tracer.plot
tracer.plot <- function(names = c("Alex Bates","Ruairi Roberts","istvan taisz", "Greg Jefferis", "Adam Heath", "Clement Hallou", "Philipp Ranft", "Philipp Schlegel"), calc.date = "2013-04-01", from.date = "2016-04-01", to.date = Sys.Date(), cumulative = T, ...){
  require(ggplot2)
  require(easyGgplot2)
  require(lubridate)
  require(dplyr)
  df = catmaid::catmaid_user_history(from = calc.date)
  df[is.na(df)] = 0
  df <- df %>%
    group_by(full_name) %>%
    arrange(date) %>%
    mutate(cable.cum.sum = cumsum(new_cable), connectors.cum.sum = cumsum(new_connectors))
  if (!cumulative){
    p1 <- ggplot2::qplot(date, new_cable, col=full_name,
              data=filter(df, full_name%in%names),  ylim=c(0, max(filter(df, full_name%in%names)$new_cable)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      geom_point()+
      geom_smooth()+
      theme(legend.position="none")
    p2 <- ggplot2::qplot(date, new_connectors, col=full_name,
              data=filter(df, full_name%in%names),  ylim=c(0, max(filter(df, full_name%in%names)$new_connectors)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      geom_point()+
      geom_smooth()+
      theme(legend.position="top")
  }else{
    p1 <- ggplot2::qplot(date, cable.cum.sum, col=full_name,
                data=filter(df, full_name%in%names),  ylim=c(0, max(filter(df, full_name%in%names)$cable.cum.sum)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      geom_point()+
      geom_path()+
      theme(legend.position="none")
    p2 <- ggplot2::qplot(date, connectors.cum.sum, col=full_name,
                data=filter(df, full_name%in%names),  ylim=c(0, max(filter(df, full_name%in%names)$connectors.cum.sum)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      geom_point()+
      geom_path()+
      theme(legend.position="top")
  }
  easyGgplot2::ggplot2.multiplot(p2,p1,cols=1)
}


#' @export
#' @rdname tracer.plot
summarise_contribution <- function(skids, auth=5.0, ack=3000, ...) {
  ul=catmaid_get_user_list()
  uls=ul[,c('full_name','id')]
  stats=catmaid_get_contributor_stats(skids, ...)
  stats.summ <- inner_join(stats$node_contributors, uls, by='id') %>%
    arrange(desc(n)) %>%
    mutate(pct=n/sum(n)*100, cpct=cumsum(pct)) %>%
    filter(n>=ack) %>%
    mutate(action=ifelse(pct>=auth, "auth", "ack"))
  stats.summ
}

#' @export
#' @rdname tracer.plot
write_ack <- function(x, ...) {
  with(subset(x, action=='ack'),
       cat("We thank", paste(full_name, collapse=", "),
           "for contributing", round(sum(pct), digits = 1),
           "% of reconstructed arbour cable."))
}

#' @export
#' @rdname tracer.plot
tracer.treenodes.plot <- function(skids = NULL, names = c("Alex Bates","Ruairi Roberts","istvan taisz", "Greg Jefferis", "Adam Heath", "Clement Hallou", "Philipp Ranft", "Philipp Schlegel", "Melissa Ryan"), from.date = "2016-04-01", to.date = Sys.Date(), ...){
  require(ggplot2)
  require(easyGgplot2)
  require(lubridate)
  require(dplyr)
  if (is.null(skids)){
    skids=as.integer(catmaid::catmaid_fetch(paste("/1/skeletons/?nodecount_gt=1")))
  }
  treenodes=nlapply(skids, catmaid_get_treenode_table, OmitFailures = T)
  ul=catmaid_get_user_list()
  as.Date(as.POSIXct(z, origin = "1970-01-01"))
  if (is.null(names)){
    names = ul$full_name
  }
  user_ids = subset(ul, full_name%in%names)$id
  treenodes %>%
    bind_rows %>%
    filter(user_id%in%user_ids) %>%
    arrange(last_modified) %>%
    mutate(last_modified = as.Date(as.POSIXct(last_modified, origin = "1970-01-01"))) %>%
    select(-one_of(c("id","parent_id","confidence","x","y","z","r", "reviewer_id")))%>%
    mutate(node_count = 1) %>%
    group_by(user_id, last_modified) %>%
    summarize(new_treenodes = n()) %>%
    group_by(user_id) %>%
    mutate(full_name = unique(subset(ul$full_name, ul$id==user_id))) %>%
    mutate(treenode.cum.sum = cumsum(new_treenodes))->
    treenodes.names
  treenodes.names$full_name = ul$full_name[match(treenodes.names$user_id, ul$id)]
  p1 <- ggplot2::qplot(last_modified, new_treenodes, col=full_name,
                       data=treenodes.names,  ylim=c(0, max(treenodes.names$new_treenodes)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      geom_point()+
      geom_smooth()+
      theme(legend.position="none")
  p2 <- ggplot2::qplot(last_modified, treenode.cum.sum, col=full_name,
                data=treenodes.names,  ylim=c(0, max(treenodes.names$treenode.cum.sum)), xlim = c(as.Date(from.date), as.Date(to.date)))+
      geom_point()+
      geom_path()+
      theme(legend.position="top")
  easyGgplot2::ggplot2.multiplot(p2,p1,cols=1)
}

#' @export
#' @rdname tracer.plot
tracer.neuron.stats <- function(skids, value = c("nodes","pre","post"), ...){
  neuronds.stats=catmaid_get_contributor_stats(skids)
  ul=catmaid_get_user_list()
  uls=ul[,c('full_name','id')]
  neuronds.stats$node_contributors
  if (value == "nodes"){
    inner_join(neuronds.stats$node_contributors, uls, by='id') %>%
    arrange(desc(n)) %>%
    mutate(pct=n/sum(n)*100, cpct=cumsum(pct))
  }else if (value == "pre"){
    inner_join(neuronds.stats$pre_contributors, uls, by='id') %>%
    arrange(desc(n)) %>%
    mutate(pct=n/sum(n)*100)
  }else if (value == "post"){
    inner_join(neuronds.stats$post_contributors, uls, by='id') %>%
    arrange(desc(n)) %>%
    mutate(pct=n/sum(n)*100)
  }
}



