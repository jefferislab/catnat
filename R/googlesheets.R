# Hidden
connector_URL <- function(df){
  if(!is.data.frame(df))
    stop("Please give me a data frame!")
  base = "https://neuropil.janelia.org/tracing/fafb/v14"
  catmaid_url = paste0(base, "?pid=1")
  catmaid_url = paste0(catmaid_url, "&zp=", df[["z"]])
  catmaid_url = paste0(catmaid_url, "&yp=", df[["y"]])
  catmaid_url = paste0(catmaid_url, "&xp=", df[["x"]])
  catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
  id = if(is.null(df["partner_skid"])) {
    df[["connector_id"]]
  } else {
    df[["partner_skid"]]
  }
  catmaid_url = paste0(catmaid_url, "&active_skeleton_id=", id)
  catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
  invisible(catmaid_url)
}

# Hidden
catmaid_urls <- function (df) {
  if (!is.data.frame(df)){
    stop("Please give me a data frame!")
  }
  base = "https://neuropil.janelia.org/tracing/fafb/v14"
  catmaid_url = paste0(base, "?pid=1")
  catmaid_url = paste0(catmaid_url, "&zp=", df[["z"]])
  catmaid_url = paste0(catmaid_url, "&yp=", df[["y"]])
  catmaid_url = paste0(catmaid_url, "&xp=", df[["x"]])
  catmaid_url = paste0(catmaid_url, "&tool=tracingtool")
  id = if (is.null(df$partner_skid)) {
    df$connector_id
  }else {
    df$partner_skid
  }
  catmaid_url = paste0(catmaid_url, "&active_skeleton_id=",
                       id)
  catmaid_url = paste0(catmaid_url, "&sid0=5&s0=0")
  catmaid_url
}

# Hidden
update_tracing_sheet <- function(df, prepost = 1, polypre = FALSE){
  df$URL = catmaid_urls(df)
  if(prepost==1){
    message("Reading information on presynaptic partners...")
    df$pre_name = catmaid::catmaid_get_neuronnames(as.integer(df$partner_skid))
    if(sum(is.na(df$pre_name))>0){
      df[which(is.na(df$pre_name)),]$partner_skid = sapply(df[which(is.na(df$pre_name)),]$connector_id, function(x) catmaid::catmaid_get_connectors(x)$pre[1])
      #if(sum(is.na(df$pre_name))>0){
      #  newskids = df[which(is.na(df$pre_name)),]$partner_skid
      #  df[which(is.na(df$pre_name))[!sapply(newskids, is.null)],]$pre_name = catmaid::catmaid_get_neuronnames(unlist(newskids))
      #}
      df$URL = connector_URL(df)
    }
    df$pre_nodes = 1
    pre = catmaid::read.neurons.catmaid(unique(unlist(df$partner_skid)),OmitFailures= T)
    result = summary(pre)
    df$pre_nodes[df$partner_skid%in%names(pre)] = result[as.character(df$partner_skid[df$partner_skid%in%names(pre)]),]$nodes
    df$pre_soma = FALSE
    df$pre_soma[df$partner_skid%in%names(pre)] = result[as.character(df$partner_skid[df$partner_skid%in%names(pre)]),]$nsoma
    df = unique.neurons.trace(df,prepost=prepost,polypre = polypre)
  }else if (!polypre){
    message("Reading information on postsynaptic partners...")
    df.c = lapply(df$connector_id, catmaid::catmaid_get_connectors)
    df$connections.laid = lapply(df.c, nrow)
    max.post = max(unlist(lapply(df.c, nrow)))
    df = cbind(df, matrix(NA,ncol=max.post,nrow=nrow(df)))
    df$nsoma = 0
    skids = sapply(df.c,function(x) x$post)
    count = c()
    c = 0
    for(n in 1:length(df.c)){
      message(paste0("Updating information on connector ",n," of ",length(df.c)))
      x = df.c[[n]]$post
      if(is.null(x)) {
        df[n,]$nsoma = 0
        count = c(count,c)
      }else{
        p = read.neurons.catmaid(unique(x),OmitFailures= TRUE,.progress="none")
        xn = catmaid::catmaid_get_neuronnames(as.integer(x))
        max.post = max(unlist(lapply(df.c, nrow)))
        length(xn) <- max.post
        df[n,9:(ncol(df)-1)] = xn
        df[n,]$nsoma = sum(summary(p)$nsoma>0)
        c= c + ifelse(n==1,1,sum(!x%in%unlist(skids[1:(n-1)])))
        count = c(count,c)
      }
    }
    df$unique.neurons = count
  }else if (polypre&prepost==0){
    message("Reading information on presynaptic partners...")
    df$connections.laid = sapply(df$connector_id, function(x) sum(df$connector_id==x) )
    df$post_name = catmaid::catmaid_get_neuronnames(as.integer(df$partner_skid)) #sapply(df$partner_skid, function(x) tryCatch(catmaid::catmaid_get_neuronnames(x),error = function(e) "neuron"))
    if(sum(is.na(df$post_name))>0){
      df[which(is.na(df$post_name)),]$partner_skid = sapply(df[which(is.na(df$post_name)),]$connector_id, function(x) catmaid::catmaid_get_connectors(x)$post[1])
      #if(sum(is.na(df$post_name))>0){
      #  newskids = df[which(is.na(df$post_name)),]$partner_skid
      #  df[which(is.na(df$post_name))[!sapply(newskids, is.null)],]$post_name = catmaid::catmaid_get_neuronnames(unlist(newskids))
      #}
      df$URL = connector_URL(df)
    }
    df$post_nodes = 1
    post = catmaid::read.neurons.catmaid(unique(df$partner_skid),OmitFailures= T)
    result = summary(post)
    df$post_nodes[df$partner_skid%in%names(post)] = result[as.character(df$partner_skid[df$partner_skid%in%names(post)]),]$nodes
    df$post_soma = FALSE
    df$post_soma[df$partner_skid%in%names(post)] = result[as.character(df$partner_skid[df$partner_skid%in%names(post)]),]$nsoma
    df = unique.neurons.trace(df,prepost=prepost,polypre = polypre)
  }
  df
}

# Hidden
unique.neurons.trace <- function(df, prepost = 1, polypre = FALSE){
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
  df
}

#' Create or update a googlesheet of randomised tracing targets for a given neuron
#'
#' @description Creates or updates a googlesheet in your drive for a neuron. Worksheets in the googlesheet containin randomised lists of synapses and CATMAID URLs to those synapses,
#' as well as other meta information to inform tracing. Different worksheets contain separate lists for pre and post synapses, in the whole neuron and in labelled axonic or
#' dendritic compartments. Randomised pre synapse lists for both presynaptic connector objects for the given neuron and presynaptic connections can be given.
#'
#' @param neuron a neuron object. If neuron = NULL is given to the update function, data will be updated in the google sheet but any addition or removal of nodes/synapses to the CATMAID neuron will be ignored
#' @param sheet_title title of the googlesheet to be created or updated
#' @param polypre if TRUE, then randomised worksheets for connections laid at the neuron's presynapses are generated
#' @param axon.dendrite.split if TRUE and the neuron is labelled in its neuron$d data frame in SWC fashion, then separate worksheets for randomised synapse lists in axon and dendrite are generated
#' @param ws name of worksheet in tracing googlesheet to specifically update
#' @param skid CATMAID skeleton ID for the neuron
#' @export
#' @rdname create_tracing_googlesheet
create_tracing_googlesheet <-function(sheet_title, neuron, skid = neuron$skid, polypre = FALSE, axon.dendrite.split = FALSE){
  if(is.null(skid)){
    stop("Please provide a CATMAID skeleton ID")
  }
  message("Reminder: Ideally, a CATMAID neuron should be as complete as possible before one starts to sample its partners")
  googlesheets::gs_auth(verbose=TRUE)
  message("This might take a small while, get a coffee or look at Facebook or something")
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
  googlesheets::gs_ws_new(gs, ws_title = "whole neuron output connections", verbose = TRUE)
  gs = googlesheets::gs_title(sheet_title)
  # Get neuron synaptic data
  df = catmaid_get_connector_table(skid)
  df = subset(df,partner_treenode_id%in%neuron$d$PointNo) # remove points not in neuron skeleton
  df$treenode_id = lapply(df$connector_id, function(y) ifelse(y%in%neuron$connectors$connector_id,neuron$connectors$treenode_id[neuron$connectors$connector_id==y][1],0) )
  df = subset(df, treenode_id%in%neuron$d$PointNo)
  df.post =  df[df$direction == "incoming",]
  if(nrow(df.post)>2){
    message("Adding randomised input synapse list...")
    df.post = update_tracing_sheet(df=df.post, prepost = 1)
    df.post = df.post[sample(nrow(df.post)),] # Randomise synapses
    if(axon.dendrite.split){
      message("Adding randomised input synapse list in separate dendrite and axon worksheets...")
      dend = subset(neuron$d, Label==3)$PointNo
      df.dend = subset(df.post,treenode_id%in%dend)
      df.dend$running.completion = (1:nrow(df.dend))/nrow(df.dend)
      axon = subset(neuron$d, Label==2)$PointNo
      df.axon = subset(df.post,treenode_id%in%axon)
      df.axon$running.completion = (1:nrow(df.axon))/nrow(df.axon)
      df.other = subset(df.post,!treenode_id%in%c(dend,axon))
      if(nrow(df.other)>0){
        if(nrow(df.other)>0){
          df.other$running.completion = (1:nrow(df.other))/nrow(df.other)
          googlesheets::gs_edit_cells(gs, ws = "other input", input = unique.neurons.trace(df.other), col_names = TRUE)
        }
      }
      googlesheets::gs_edit_cells(gs, ws = "dendritic input", input = unique.neurons.trace(df.dend), col_names = TRUE)
      googlesheets::gs_edit_cells(gs, ws = "axonal input", input = unique.neurons.trace(df.axon), col_names = TRUE)
    }
    df.post$running.completion = (1:nrow(df.post))/nrow(df.post)
    googlesheets::gs_edit_cells(gs, ws = "whole neuron input", input = unique.neurons.trace(df.post), col_names = TRUE)
  }
  df.pre =  subset(neuron$connectors, prepost==0)
  if(nrow(df.pre)>2){
    message("Adding randomised output connector list...")
    df.pre = update_tracing_sheet(df=df.pre, prepost = 0, polypre = FALSE)
    df.pre= df.pre[sample(nrow(df.pre)),] # Randomise synapses
    df.pre$status = "TODO"
    if(axon.dendrite.split){
      message("Adding randomised output connector list in separate dendrite and axon worksheets...")
      df.dend = subset(df.pre,treenode_id%in%dend)
      df.dend$running.completion = (1:nrow(df.dend))/nrow(df.dend)
      df.axon = subset(df.pre,treenode_id%in%axon)
      df.axon$running.completion = (1:nrow(df.axon))/nrow(df.axon)
      df.other = subset(df.pre,!treenode_id%in%c(dend,axon))
      if(nrow(df.other)>0){
        df.other = subset(df.post,!treenode_id%in%c(dend,axon))
        if(nrow(df.other)>0){
          df.other$running.completion = (1:nrow(df.other))/nrow(df.other)
          googlesheets::gs_edit_cells(gs, ws = "other input", input = unique.neurons.trace(df.other, prepost=0), col_names = TRUE)
        }
      }
      googlesheets::gs_edit_cells(gs, ws = "dendritic output connectors", input = unique.neurons.trace(df.dend,prepost=0), col_names = TRUE)
      googlesheets::gs_edit_cells(gs, ws = "axonal output connectors", input = unique.neurons.trace(df.axon,prepost=0), col_names = TRUE)
    }
    df.pre$running.completion = (1:nrow(df.pre))/nrow(df.pre)
    googlesheets::gs_edit_cells(gs, ws = "whole neuron output connectors", input = unique.neurons.trace(df.pre,prepost=0), col_names = TRUE)
    if(polypre){
      message("Adding randomised output connections list...")
      df.polypre =  df[df$direction == "outgoing",]
      df.polypre = update_tracing_sheet(df=df.polypre, prepost = 0, polypre = TRUE)
      df.polypre = df.polypre[sample(nrow(df.polypre)),] # Shuffle
      if(axon.dendrite.split){
        message("Adding randomised output connections list in separate dendrite and axon worksheets...")
        df.dend = subset(df.polypre,treenode_id%in%dend)
        df.dend$running.completion = (1:nrow(df.dend))/nrow(df.dend)
        df.axon = subset(df.polypre,treenode_id%in%axon)
        df.axon$running.completion = (1:nrow(df.axon))/nrow(df.axon)
        df.other = subset(df.polypre,!treenode_id%in%c(dend,axon))
        if(nrow(df.other)>0){
          df.other = subset(df.post,!treenode_id%in%c(dend,axon))
          if(nrow(df.other)>0){
            df.other$running.completion = (1:nrow(df.other))/nrow(df.other)
            googlesheets::gs_edit_cells(gs, ws = "other input", input = unique.neurons.trace(df.other, prepost = 0, polypre = TRUE), col_names = TRUE)
          }
        }
        googlesheets::gs_edit_cells(gs, ws = "dendritic output connections", input = unique.neurons.trace(df.dend, prepost = 0, polypre = TRUE), col_names = TRUE)
        googlesheets::gs_edit_cells(gs, ws = "axonal output connections", input = unique.neurons.trace(df.axon, prepost = 0, polypre = TRUE), col_names = TRUE)
      }
      df.polypre$running.completion = (1:nrow(df.polypre))/nrow(df.polypre)
      googlesheets::gs_edit_cells(gs, ws = "whole neuron output connections", input = unique.neurons.trace(df.polypre, prepost = 0, polypre = TRUE), col_names = TRUE)
    }
  }
  message("Tracing sheet created in your home google drive directory. Please move to desired location on google drive")
}

#' @export
#' @rdname create_tracing_googlesheet
update_tracing_worksheet <- function(sheet_title, neuron = NULL, skid = neuron$skid, ws = "whole neuron input"){
  googlesheets::gs_auth(verbose=TRUE)
  if(!is.character(ws)){
    stop("Worksheet must be named")
  }
  if(grepl("dend",ws)&!is.null(neuron)){
    dend = rownames(subset(neuron$d, Label==3))
    neuron = nat::prune_vertices(neuron, dend, invert = TRUE)
  }
  if(grepl("axon",ws)&!is.null(neuron)){
    axon = rownames(subset(neuron$d, Label==2))
    neuron = nat::prune_vertices(neuron, axon, invert = TRUE)
  }
  if(grepl("input",ws)){
    prepost = 1
  }else{
    prepost = 0
  }
  if(grepl("connections",ws)){
    polypre = TRUE
  }else{
    polypre = FALSE
  }
  gs = googlesheets::gs_title(sheet_title)
  gss = googlesheets::gs_read(gs, ws = ws, range = NULL, literal = TRUE, verbose = TRUE, col_names = TRUE)
  if(is.null(neuron)){
    gss.final = update_tracing_sheet(gss, polypre = polypre, prepost = prepost)
  }else{
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
    # Update old
    gss.updated = update_tracing_sheet(gss, polypre = polypre, prepost = prepost)
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
      gss.updated = gss.final
      gss.final = unique.neurons.trace(gss.final,prepost = prepost, polypre = polypre)
      gss.final$running.completion = (1:nrow(gss.final))/nrow(gss.final)
    }
  }
  googlesheets::gs_ws_rename(gs, from = ws, to = paste0("old ",ws), verbose = TRUE)
  gs = googlesheets::gs_title(sheet_title)
  googlesheets::gs_ws_new(gs, verbose = TRUE, ws_title = ws)
  gs = googlesheets::gs_title(sheet_title)
  message("Writing googlesheet...")
  googlesheets::gs_edit_cells(gs, ws = ws, input = gss.final, col_names = TRUE)
  googlesheets::gs_ws_delete(gs, ws = paste0("old ",ws), verbose = TRUE)
  message(paste0("Neuron tracing worksheet ",ws," fully updated in ",sheet_title))
}

#' @export
#' @rdname create_tracing_googlesheet
update_tracing_googlesheet <-function(sheet_title, neuron = NULL, polypre = FALSE){
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
    update_tracing_worksheet(sheet_title = sheet_title, neuron = neuron, ws = ws)
  }
}








