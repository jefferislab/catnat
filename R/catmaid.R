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


