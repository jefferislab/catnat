#' Generate URLs to CATMAID connector location
#'
#' @description   Generate URLs to CATMAID connector location
#' @param df data frame that must contain a connector_id column, and x, y, z columns
#' @param server CATMAID instance
#' @export
#' @rdname connector_URL
#' @export
connector_URL <- function(df, server = "https://neuropil.janelia.org/tracing/fafb/v14/"){
  if(!is.data.frame(df))
    stop("Please give me a data frame!")
  catmaid_url = paste0(server, "?pid=1")
  catmaid_url = paste0(catmaid_url, "&zp=", df[["z"]])
  catmaid_url = paste0(catmaid_url, "&yp=", df[["y"]])
  catmaid_url = paste0(catmaid_url, "&xp=", df[["x"]])
  catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
  if(is.null(df["connector_id"])) {
    id = df[["partner_skid"]]
    catmaid_url = paste0(catmaid_url, "&active_skeleton_id=", id)
  } else {
    id = df[["connector_id"]]
    catmaid_url = paste0(catmaid_url, "&active_node_id=", id)
  }
  catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
  invisible(catmaid_url)
}

# Hidden
update_tracing_sheet <- function(df, prepost = 1, polypre = FALSE){
  if(prepost==1){
    df$URL = connector_URL(df)
    message("Reading information on presynaptic partners...")
    df$pre_name = catmaid::catmaid_get_neuronnames(as.integer(df$partner_skid))
    if(sum(is.na(df$pre_name))>0){
      df[which(is.na(df$pre_name)),]$partner_skid = sapply(df[which(is.na(df$pre_name)),]$connector_id, function(x) catmaid::catmaid_get_connectors(x)$pre[1])
      df$URL = connector_URL(df)
    }
    df$pre_nodes = 1
    df$partner_skid = as.character(df$partner_skid)
    partner_skids = unique(df$partner_skid)[unique(df$partner_skid)!="NULL"]
    df$pre_soma = FALSE
    if(length(partner_skids)>0){
      pre = catmaid::read.neurons.catmaid(partner_skids,OmitFailures= T)
      result = summary(pre)
      df$pre_nodes[df$partner_skid%in%names(pre)] = result[as.character(df$partner_skid[df$partner_skid%in%names(pre)]),]$nodes
      df$pre_soma[df$partner_skid%in%names(pre)] = result[as.character(df$partner_skid[df$partner_skid%in%names(pre)]),]$nsoma
    }
    df = unique.neurons.trace(df,prepost=prepost,polypre = polypre)
  }else if (polypre==FALSE){
    message("Reading information on postsynaptic partners...")
    df = df[,1:6]
    df.c = lapply(df$connector_id, catmaid::catmaid_get_connectors)
    df$connections.laid = lapply(df.c, nrow)
    max.post = max(unlist(lapply(df.c, nrow)))
    df = cbind(df, matrix(NA,ncol=max.post,nrow=nrow(df)))
    df$nsoma = 0
    skids = sapply(df.c,function(x) x$post)
    c = 0
    for(n in 1:length(df.c)){
      message(paste0("Updating information on connector ",n," of ",length(df.c)))
      x = df.c[[n]]$post
      if(is.null(x)) {
        df[n,]$nsoma = 0
      }else{
        p = read.neurons.catmaid(unique(x),OmitFailures= FALSE,.progress="none")
        xn = catmaid::catmaid_get_neuronnames(as.integer(x))
        max.post = max(unlist(lapply(df.c, nrow)))
        length(xn) <- max.post
        df[n,8:(max.post+7)] = xn
        df[n,]$nsoma = tryCatch(sum(summary(p)$nsoma>0),error=function(e) 0)
      }
    }
    df = unique.neurons.trace(df,prepost=0,polypre = FALSE)
  }else if (polypre==TRUE){
    message("Reading information on postsynaptic partners...")
    df$connections.laid = sapply(df$connector_id, function(x) sum(df$connector_id==x) )
    df$partner_skid = as.character(df$partner_skid)
    df$post_name = catmaid::catmaid_get_neuronnames(as.integer(df$partner_skid)) #sapply(df$partner_skid, function(x) tryCatch(catmaid::catmaid_get_neuronnames(x),error = function(e) "neuron"))
    if(sum(is.na(df$post_name))>0){
      df[which(is.na(df$post_name)),]$partner_skid = as.character(sapply(df[which(is.na(df$post_name)),]$connector_id, function(x) catmaid::catmaid_get_connectors(x)$post[1]))
      df$URL = connector_URL(df)
    }
    df$post_nodes = 1
    partner_skids = unique(df$partner_skid)[unique(df$partner_skid)!="NULL"]
    df$post_soma = FALSE
    if(length(partner_skids)>0){
      post = catmaid::read.neurons.catmaid(partner_skids,OmitFailures= T)
      result = summary(post)
      df$post_nodes[df$partner_skid%in%names(post)] = result[as.character(df$partner_skid[df$partner_skid%in%names(post)]),]$nodes
      df$post_soma = FALSE
      df$post_soma[df$partner_skid%in%names(post)] = result[as.character(df$partner_skid[df$partner_skid%in%names(post)]),]$nsoma
    }
    df = unique.neurons.trace(df,prepost=prepost,polypre = polypre)
  }
  df
}

# Hidden
unique.neurons.trace <- function(df, prepost = 1, polypre = FALSE){
  if(nrow(df)>0){
    count = c()
    c = 0
    if(prepost==0&!polypre){
      df.c = lapply(df$connector_id, catmaid::catmaid_get_connectors)
      skids = sapply(df.c,function(x) x$post)
      for(n in 1:length(df.c)){
        x = df.c[[n]]$post
        if(is.null(x)) {
          count = c(count,c)
        }else{
          c= c + ifelse(n==1,1,sum(!x%in%unlist(skids[1:(n-1)])))
          count = c(count,c)
        }
      }
    }else{
      for(n in 1:length(df$partner_skid)){
        c = c + ifelse(n==1,1,!df$partner_skid[n]%in%df$partner_skid[1:(n-1)])
        count = c(count,c)
      }
    }
    df$unique.neurons = count
  }
  df
}

#' Create or update a googlesheet of randomised tracing targets for a given neuron
#'
#' @description Creates or updates a googlesheet in your drive for a neuron. Worksheets in the googlesheet containing randomised lists of synapses and CATMAID URLs to those synapses,
#' as well as other meta information to inform tracing. Different worksheets contain separate lists for pre and post synapses, in the whole neuron and in labelled axonic or
#' dendritic compartments. Randomised pre synapse lists for both presynaptic connector objects for the given neuron and presynaptic connections can be given.
#'
#' @param neuron a neuron object. If neuron = NULL is given to the update function, data will be updated in the google sheet but any addition or removal of nodes/synapses to the CATMAID neuron will be ignored. If axon.dendrite.split = TRUE, neuron should be passed through catnat::flow.centrality first in order to identify its axon and dendrite
#' @param sheet_title title of the googlesheet to be created or updated
#' @param folder location to save tracing sheets. Defaults to 'googlesheet' which creates a tracing sheet in your home folder on googledrive. It can take a while to write a googlesheet. Alternatively, CSVs can be created in a local folder.
#' @param skid CATMAID skeleton ID for the neuron
#' @param polypre if TRUE, then randomised worksheets for connections laid at the neuron's presynapses are generated
#' @param axon.dendrite.split if TRUE and the neuron is labelled in its neuron$d data frame in SWC fashion, then separate worksheets for randomised synapse lists in axon and dendrite are generated
#' @param randomise whether to randomise the synapse sampling list (recommended). If false, the sheet is organised by strongest known synaptic partners
#' @param ws the individual google worksheet or path to .csv file to update
#'
#' @examples
#' \dontrun{
#' # Load the neurons for which we want to create sampling sheets
#' wedpns.chosen = read.neurons.catmaid("annotation:^WED-PN Complete PDP$")
#' # Assign their cell types to these neurons
#' # This should be changed from being manual to depending on an external gogolesheet or the name in CATMAID at some point
#' wedpns.chosen[,"cell.type"]  = c("WED-PN3","WED-PN4","WED-PN1","WED-PN5","WED-PN6","WED-PN2")
#' # Run the flow-centrality axon-dendrite split algorithm
#' wedpns.chosen.flow = catnat::flow.centrality(wedpns.chosen, polypre= FALSE, mode = "centrifugal")
#' # Create a new folder for seach of these cell types locally
#' dir.create("Data/sampling")
#' for(ct in wedpns.chosen.flow[,"cell.type"]){
#'   message(ct)
#'   if(dir.exists(paste0("Data/sampling/",ct))){
#'     message("Sampling sheets for this cell type ought already to be present!")
#'   }else{
#'     dir.create(paste0("Data/sampling/",ct))
#'     wedpn = subset(wedpns.chosen.flow,cell.type==ct)[[1]]
#'     catnat::create_tracing_samplesheet(neuron=wedpn,sheet_title = ct, folder = paste0("Data/sampling/",ct,"/"), polypre = TRUE, axon.dendrite.split = TRUE, randomise = TRUE)
#'   }
#' }
#' #Update the sampling sheets as connections have been sampled / the subject neuron has been further modified
#' for(ct in wedpns.chosen.flow[,"cell.type"]){
#'   message("Working on ",ct)
#'   wedpn = subset(wedpns.chosen.flow,cell.type==ct)[[1]]
#'   catnat::update_tracing_samplesheets(neuron=wedpn,sheet_title = ct, folder = paste0("Data/sampling/",ct,"/"), polypre = TRUE)
#' }
#' }
#' @export
#' @rdname create_tracing_samplesheet
create_tracing_samplesheet <-function(neuron, sheet_title = "mystery_neuron", folder = "googlesheet", skid = neuron$skid, polypre = TRUE, axon.dendrite.split = FALSE, randomise = TRUE){
  if(is.null(skid)){
    stop("Please provide a CATMAID skeleton ID")
  }
  message("Reminder: Ideally, a CATMAID neuron should be as complete as possible before one starts to sample its partners")
  if(folder=="googlesheet"){
    message("NOTE: You are directly writing to a googlesheet on your drive. This can take a while to complete. Enjoy a coffee in the meantime.
            It is quicker to write a .csv directly to a local folder on your machine by changing the argument 'folder' to the path to some local folder.")
    googlesheets::gs_auth(verbose=TRUE)
    # Create the googlesheet object
    gs = googlesheets::gs_new(title = sheet_title, ws_title = "whole neuron input")
    if(axon.dendrite.split){
        googlesheets::gs_ws_new(gs, ws_title = "dendritic input", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "axonal input", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "other input", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "dendritic output connectors", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "axonal output connectors", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "other output connectors", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "dendritic output connections", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "axonal output connections", verbose = TRUE)
        googlesheets::gs_ws_new(gs, ws_title = "other output connections", verbose = TRUE)
    }
    googlesheets::gs_ws_new(gs, ws_title = "whole neuron output connectors", verbose = TRUE)
    if(polypre)googlesheets::gs_ws_new(gs, ws_title = "whole neuron output connections", verbose = TRUE)
    gs = googlesheets::gs_title(sheet_title)
  }
  # Get neuron synaptic data
  df = catmaid::catmaid_get_connector_table(skid)
  df = subset(df,partner_treenode_id%in%neuron$d$PointNo) # remove points not in neuron skeleton
  df$treenode_id = lapply(df$connector_id, function(y) ifelse(y%in%neuron$connectors$connector_id,neuron$connectors$treenode_id[neuron$connectors$connector_id==y][1],0) )
  df = subset(df, treenode_id%in%neuron$d$PointNo)
  df.post =  df[df$direction == "incoming",]
  dend = subset(neuron$d, Label==3)$PointNo
  axon = subset(neuron$d, Label==2)$PointNo
  if(nrow(df.post)>2){
    message("Adding input synapse list...")
    df.post = update_tracing_sheet(df=df.post, prepost = 1) # Add more information to this sheet
    if(randomise){
      df.post = df.post[sample(nrow(df.post)),] # Randomise synapses
    }else{
      df.post = df.post[order(df.post$pre_nodes,decreasing=TRUE),] # Order by strongest inputs
    }
    df.post[] = lapply(df.post, as.character)
    if(axon.dendrite.split){
      message("Adding input synapse list in separate dendrite and axon worksheets...")
      df.dend = subset(df.post,treenode_id%in%dend)
      if(nrow(df.dend)>0){df.dend$running.completion = (1:nrow(df.dend))/nrow(df.dend)}
      df.axon = subset(df.post,treenode_id%in%axon)
      if(nrow(df.axon)>0){df.axon$running.completion = (1:nrow(df.axon))/nrow(df.axon)}
      df.other = subset(df.post,!treenode_id%in%c(dend,axon))
      df.dend[] = lapply(df.dend, as.character)
      df.axon[] = lapply(df.axon, as.character)
      if(nrow(df.other)>0){
        df.other[] = lapply(df.other, as.character)
        df.other = subset(df.post,!treenode_id%in%c(dend,axon))
        df.other$running.completion = (1:nrow(df.other))/nrow(df.other)
        if(nrow(df.other)>0){
          if(folder=="googlesheet"){
            googlesheets::gs_edit_cells(gs, ws = "other input", input = unique.neurons.trace(df.other, prepost=1), col_names = TRUE)
          }else{
            utils::write.csv(x = unique.neurons.trace(df.other, prepost = 1),file = paste0(folder,sheet_title,"_inputs_other",".csv"),row.names=FALSE)
          }
        }
      }
      if(folder=="googlesheet"){
        googlesheets::gs_edit_cells(gs, ws = "dendritic input", input = unique.neurons.trace(df.dend, prepost = 1), col_names = TRUE)
        googlesheets::gs_edit_cells(gs, ws = "axonal input", input = unique.neurons.trace(df.axon, prepost = 1), col_names = TRUE)
      }else{
        utils::write.csv(x = unique.neurons.trace(df.axon, prepost = 1),file = paste0(folder,sheet_title,"_inputs_axon",".csv"),row.names=FALSE)
        utils::write.csv(x = unique.neurons.trace(df.dend, prepost = 1),file = paste0(folder,sheet_title,"_inputs_dendrite",".csv"),row.names=FALSE)
      }
    }
    df.post$running.completion = (1:nrow(df.post))/nrow(df.post)
    if(folder=="googlesheet"){
      googlesheets::gs_edit_cells(gs, ws = "whole neuron input", input = unique.neurons.trace(df.post), col_names = TRUE)
    }else{
      utils::write.csv(x = unique.neurons.trace(df.post, prepost = 1),file = paste0(folder,sheet_title,"_inputs_whole",".csv"),row.names=FALSE)
    }
  }
  df.pre =  subset(neuron$connectors, prepost==0)
  if(nrow(df.pre)>2){
    message("Adding output connector list...")
    df.pre = update_tracing_sheet(df=df.pre, prepost = 0, polypre = FALSE)
    if(randomise){
      df.pre = df.pre[sample(nrow(df.pre)),] # Randomise synapses
    }else{
      df.pre[is.na(df.pre)] = 0
      df.pre = df.pre[order(df.pre$nsoma,decreasing=TRUE),]
    }
    df.pre[] = lapply(df.pre, as.character)
    if(axon.dendrite.split){
      message("Adding output connector list in separate dendrite and axon worksheets...")
      df.pre.dend = subset(df.pre,treenode_id%in%dend)
      if(nrow(df.pre.dend)>0){df.pre.dend$running.completion = (1:nrow(df.pre.dend))/nrow(df.pre.dend)}
      df.pre.axon = subset(df.pre,treenode_id%in%axon)
      if(nrow(df.pre.axon)>0){df.pre.axon$running.completion = (1:nrow(df.pre.axon))/nrow(df.pre.axon)}
      df.pre.other = subset(df.pre,!treenode_id%in%c(dend,axon))
      df.pre.dend[] = lapply(df.pre.dend, as.character)
      df.pre.axon[] = lapply(df.pre.axon, as.character)
      if(nrow(df.pre.other)>0){
        df.pre.other[] = lapply(df.pre.other, as.character)
        if(nrow(df.pre.other)>0){
          df.pre.other$running.completion = (1:nrow(df.pre.other))/nrow(df.pre.other)
          if(folder=="googlesheet"){
            googlesheets::gs_edit_cells(gs, ws = "other output connectors", input = unique.neurons.trace(df.pre.other,prepost = 0, polypre = FALSE), col_names = TRUE)
          }else{
            utils::write.csv(x = unique.neurons.trace(df.pre.other, prepost = 0, polypre = FALSE),file = paste0(folder,sheet_title,"_ouput_connectors_other",".csv"),row.names=FALSE)
          }
        }
      }
      if(folder=="googlesheet"){
        googlesheets::gs_edit_cells(gs, ws = "dendritic output connectors", input = unique.neurons.trace(df.pre.dend,prepost = 0, polypre = FALSE), col_names = TRUE)
        googlesheets::gs_edit_cells(gs, ws = "axonal output connectors", input = unique.neurons.trace(df.pre.axon,prepost = 0, polypre = FALSE), col_names = TRUE)
      }else{
        utils::write.csv(x = unique.neurons.trace(df.pre.axon,prepost = 0, polypre = FALSE),file = paste0(folder,sheet_title,"_output_connectors_axon",".csv"),row.names=FALSE)
        utils::write.csv(x = unique.neurons.trace(df.pre.dend,prepost = 0, polypre = FALSE),file = paste0(folder,sheet_title,"_ouput_connectors_dendrite",".csv"),row.names=FALSE)
      }
    }
    df.pre$running.completion = (1:nrow(df.pre))/nrow(df.pre)
    if(folder=="googlesheet"){
      googlesheets::gs_edit_cells(gs, ws = "whole neuron output connectors", input = unique.neurons.trace(df.pre,prepost=0,polypre = FALSE), col_names = TRUE)
    }else{
      utils::write.csv(x = unique.neurons.trace(df.pre,prepost = 0, polypre = FALSE),file = paste0(folder,sheet_title,"_ouput_connectors_whole",".csv"),row.names=FALSE)
    }
    if(polypre){
      message("Adding output connections list...")
      df.polypre =  df[df$direction == "outgoing",]
      df.polypre = update_tracing_sheet(df=df.polypre, prepost = 0, polypre = TRUE)
      if(randomise){
        df.polypre = df.polypre[sample(nrow(df.polypre)),] # Shuffle
      }else{
        df.polypre = df.polypre[order(df.polypre$post_nodes,decreasing=TRUE),] # Order by size of already-traced fragments
      }
      df.polypre[] = lapply(df.polypre, as.character)
      if(axon.dendrite.split){
        message("Adding randomised output connections list in separate dendrite and axon worksheets...")
        df.dend = subset(df.polypre,treenode_id%in%dend)
        if(nrow(df.dend)>0){df.dend$running.completion = (1:nrow(df.dend))/nrow(df.dend)}
        df.axon = subset(df.post,treenode_id%in%axon)
        if(nrow(df.axon)>0){df.axon$running.completion = (1:nrow(df.axon))/nrow(df.axon)}
        df.other = subset(df.polypre,!treenode_id%in%c(dend,axon))
        df.dend[] = lapply(df.dend, as.character)
        df.axon[] = lapply(df.axon, as.character)
        if(nrow(df.other)>0){
          df.other[] = lapply(df.other, as.character)
          df.other = subset(df.post,!treenode_id%in%c(dend,axon))
          if(nrow(df.other)>0){
            df.other$running.completion = (1:nrow(df.other))/nrow(df.other)
            if(folder=="googlesheet"){
              googlesheets::gs_edit_cells(gs, ws = "other output connections", input = unique.neurons.trace(df.other, prepost=0, polypre = TRUE), col_names = TRUE)
            }else{
              utils::write.csv(x = unique.neurons.trace(df.other, prepost=0, polypre = TRUE),file = paste0(folder,sheet_title,"_ouput_connections_other",".csv"),row.names=FALSE)
            }
          }
        }
        if(folder=="googlesheet"){
          googlesheets::gs_edit_cells(gs, ws = "dendritic output connections", input = unique.neurons.trace(df.dend, prepost = 0, polypre = TRUE), col_names = TRUE)
          googlesheets::gs_edit_cells(gs, ws = "axonal output connections", input = unique.neurons.trace(df.axon, prepost = 0, polypre = TRUE), col_names = TRUE)
        }else{
          utils::write.csv(x = unique.neurons.trace(df.axon, prepost=0, polypre = TRUE),file = paste0(folder,sheet_title,"_output_connections_axon",".csv"),row.names=FALSE)
          utils::write.csv(x = unique.neurons.trace(df.dend, prepost=0, polypre = TRUE),file = paste0(folder,sheet_title,"_ouput_connections_dendrite",".csv"),row.names=FALSE)
          utils::write.csv(x = unique.neurons.trace(df.other, prepost=0, polypre = TRUE),file = paste0(folder,sheet_title,"_ouput_connections_other",".csv"),row.names=FALSE)
        }
      }
      df.polypre$running.completion = (1:nrow(df.polypre))/nrow(df.polypre)
      if(folder=="googlesheet"){
        googlesheets::gs_edit_cells(gs, ws = "whole neuron output connections", input = unique.neurons.trace(df.polypre, prepost = 0, polypre = TRUE), col_names = TRUE)
      }else{
        utils::write.csv(x = unique.neurons.trace(df.polypre, prepost=0, polypre = TRUE),file = paste0(folder,sheet_title,"_output_connections_whole",".csv"),row.names=FALSE)
      }
    }
  }
  if(folder=="googlesheet"){
    message("Tracing sheet created in your home google drive directory. Please move to desired location on google drive")
  }else{
    message("Tracing sheet CSVs created in your chosen local folder")
  }
}

#' @export
#' @rdname create_tracing_samplesheet
#' @importFrom utils read.csv
update_single_samplesheet <- function(sheet_title = "mystery_neuron", folder = "googlesheet", neuron = NULL, skid = neuron$skid, ws="whole neuron input"){
  if(folder=="googlesheet"){
    if(!is.character(ws)){
      stop("Worksheet must be named")
    }
    googlesheets::gs_auth(verbose=TRUE)
    gs = googlesheets::gs_title(sheet_title)
    gss = googlesheets::gs_read(gs, ws = ws, range = NULL, literal = TRUE, verbose = TRUE, col_names = TRUE)
  }else{
    file = paste0(folder,ws)
    gss = read.csv(file=file)
    gss[] = lapply(gss, as.character)
  }
  if(grepl("dend",ws)&!is.null(neuron)){
    neuron = dendritic_cable(neuron)
  }
  if(grepl("axon",ws)&!is.null(neuron)){
    neuron = axonic_cable(neuron)
  }
  if(grepl("input",ws)){
    prepost = 1
  }else{
    prepost = 0
  }
  if(grepl("connection",ws)){
    polypre = TRUE
  }else{
    polypre = FALSE
  }
  if(is.null(neuron)){
    gss.final = update_tracing_sheet(gss, polypre = polypre, prepost = prepost)
  }else{
    message("Updating sampling sheet based on supplied neuron")
    if(polypre==T|prepost==1){
      if(is.null(skid)){
        stop("Please provide a CATMAID skeleton ID")
      }
      new.df = catmaid::catmaid_get_connector_table(skid)
      new.df$treenode_id = lapply(new.df$connector_id, function(y) ifelse(y%in%neuron$connectors$connector_id,neuron$connectors$treenode_id[neuron$connectors$connector_id==y][1],0) )
      new.df = subset(new.df, treenode_id%in%neuron$d$PointNo)
      p = ifelse(prepost==1,"incoming","outgoing")
      new.df = subset(new.df, direction==p)
      new = apply(new.df[,c("connector_id","partner_treenode_id","treenode_id")],1,paste,collapse="")
      old = apply(gss[,c("connector_id","partner_treenode_id","treenode_id")],1,paste,collapse="")
      add = which(!new%in%old)
      add.df = new.df[add,]
      gss = gss[which(old%in%new),] # remove connectors that no longer exist
    }else{
      p = prepost
      new = subset(neuron$connectors, prepost==p)$connector_id
      old = gss$connector_id
      add = new[!new%in%old]
      add.df = neuron$connectors[neuron$connectors$connector_id%in%add,]
      gss = gss[old%in%new,] # remove connectors that no longer exist
    }
    if(nrow(gss)>0){
      gss.updated = update_tracing_sheet(gss, polypre = polypre, prepost = prepost)     # Update old
      if(nrow(add.df)>0){
        add.updated =  update_tracing_sheet(add.df,polypre = polypre, prepost = prepost)
        random.rows = base::sample(x = 1:(nrow(add.updated)+nrow(gss.updated)),size = nrow(add.updated))
        gss.final = as.data.frame(matrix(0,ncol = ncol(gss.updated), nrow = (nrow(add.updated)+nrow(gss.updated))))
        colnames(gss.final) = colnames(gss.updated)
        for(r in 1:(nrow(add.updated)+nrow(gss.updated))){
          if(r%in%random.rows){
            gss.final[r,] = add.updated[1,]
            add.updated = add.updated[-1,]
          }else{
            gss.final[r,] = gss.updated[1,]
            gss.updated = gss.updated[-1,]
          }
        }
        gss.final = unique.neurons.trace(gss.final,prepost = prepost, polypre = polypre)
        gss.final$running.completion = (1:nrow(gss.final))/nrow(gss.final)
      }else{
        gss.final = gss.updated
        gss.final = unique.neurons.trace(gss.final,prepost = prepost, polypre = polypre)
        gss.final$running.completion = (1:nrow(gss.final))/nrow(gss.final)
      }
    }
  }
  if(folder=="googlesheet"){
    googlesheets::gs_ws_rename(gs, from = ws, to = paste0("old ",ws), verbose = TRUE)
    gs = googlesheets::gs_title(sheet_title)
    googlesheets::gs_ws_new(gs, verbose = TRUE, ws_title = ws)
    gs = googlesheets::gs_title(sheet_title)
    message("Writing googlesheet...")
    googlesheets::gs_edit_cells(gs, ws = ws, input = gss.final, col_names = TRUE)
    googlesheets::gs_ws_delete(gs, ws = paste0("old ",ws), verbose = TRUE)
    message(paste0("Neuron tracing worksheet ",ws," fully updated in ",sheet_title))
    message("Updated google worksheet")
  }else{
    gss.final[] = lapply(gss.final, as.character)
    utils::write.csv(x=gss.final,file=file,row.names=FALSE)
    message("Updated .csv file")
  }
}

#' @export
#' @rdname create_tracing_samplesheet
update_tracing_samplesheets <-function(sheet_title = "mystery_neuron", folder = "googlesheet", neuron = NULL, skid = neuron$skid, polypre = TRUE){
  if(folder=="googlesheet"){
    gs = googlesheets::gs_title(sheet_title)
    googlesheets::gs_copy(from=gs, to = paste0(sheet_title,"_copy"), verbose = TRUE)
    message("A copy of the old spreadsheet has been made and left in your home directory on google drive, just in case...")
    worksheets = googlesheets::gs_ws_ls(gs)
    possible = c("whole neuron input","dendritic input","axonal input", "other input", "whole neuron output connectors","dendritic output connectors",
                 "axonal output connectors", "other output connectors",  "whole neuron output connections",  "dendritic output connections",
                 "axonal output connections", "other output connections")
    if(!polypre){
      possible = possible[!grepl("connections", possible)]
    }
    worksheets = worksheets[worksheets%in%possible]
    for(ws in worksheets){
      update_single_samplesheet(sheet_title = sheet_title, folder = folder, neuron = neuron, skid = skid, ws = ws)
    }
  }else{
    files = list.files(folder)
    files = files[grep(sheet_title,files)]
    for(file in files){
      message("updating: ",file)
      update_single_samplesheet(sheet_title = sheet_title, folder = folder, neuron = neuron, skid = skid, ws = file)
    }
  }
}

