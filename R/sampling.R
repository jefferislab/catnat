

randomly_sample_connections <- function(someneuronlist, direction = 0, all.connections = T, experiments = 10, clustering = c("None","flow centrality","Synapse Clusters"), soma = T, min.nodes = 1000, min.strength = 1){
  c = catmaid_get_connectors(subset(connectors(neuron)$connector_id,connectors(neuron)$prepost==direction))
  all.post = read.neurons.catmaid(unique(c$post), OmitFailures = T)
  post = subset(all.post, sapply(all.post, function(x) x$NumPoints>min.nodes))
  post = su
  if(soma){post=post[sapply(post, function(x) !is.null(x$tags$soma))]}
  frags = c$post[c$post%in%post]
  chance.its.new
}

#partner_progression <- function(someneuronlist, direction = 0, all.connections = T, soma = T, min_nodes = 1000, min_strength = 1){
#  connected = get_connected_skeletons(someneuronlist, soma=soma,prepost=direction,min_nodes=min_nodes,min_synapse = min_strength)
#  treenodes=lapply(names(connected), catmaid_get_treenode_table, OmitFailures = T)
#  c = catmaid_get_connector_table(names(someneuronlist), direction = "outgoing")
#  c = c[c$connector_id%in%connectors(someneuronlist)$connector_id,]
#  ids = sapply(1:nrow(c),function(x) connected[as.character(c[x,10])][[1]]$d[connected[as.character(c[x,10])][[1]]$d$PointNo==c[x,8],]$Parent)


#  get_2nd_point_from_synapse <- function()

#  sapply(treenodes, function(x))


#  treenodes %>%
#    dplyr::arrange(last_modified) %>%
#    dplyr::mutate(last_modified = as.Date(as.POSIXct(last_modified, origin = "1970-01-01")))->
#    treenodes.data

#  treenodes %>%
#    dplyr::bind_rows() %>%
#    dplyr::arrange(last_modified) %>%
#    dplyr::mutate(last_modified = as.Date(as.POSIXct(last_modified, origin = "1970-01-01"))) %>%
#    dplyr::select(-one_of(c("r", "reviewer_id")))%>%
#    dplyr::mutate(node_count = 1) %>%
#    dplyr::group_by(last_modified) %>%
#    dplyr::summarize(new_treenodes = n()) %>%
#    dplyr::mutate(treenode.cum.sum = cumsum(new_treenodes))->
#    treenodes.data
#}
