## ----setup, include=FALSE, eval=F, echo=T--------------------------------
#  knitr::opts_chunk$set(echo = TRUE)

## ----load.catnat, eval=F, echo=T-----------------------------------------
#  # installation
#  if (!require("devtools")) install.packages("devtools")
#  if (!require("catnat")) devtools::install_github("alexanderbates/catnat")
#  if (!require("fafbseg")) install.packages("jefferis/fafbseg")
#  
#  # Load catnat
#  library(catnat)
#  library(fafbseg)

## ----get.neurons, eval=F, echo=T-----------------------------------------
#  # Login to CATMAID
#  # see ?catmaid::catmaid_login() for details
#  
#  # Read some interesting neurons from CATMAID
#  a2sc = read.neurons.catmaid("name:MBON a2sc")
#  
#  # Let's see what we hsve
#  print(a2sc[,])

## ----tracing.sheets, eval=F, echo=T--------------------------------------
#  # Generate tracing sheet
#  tl.incoming = fafb_seg_tracing_list(skids = names(a2sc),connector_ids = NULL, direction = "incoming", unique = FALSE)
#  tl.outgoing= fafb_seg_tracing_list(skids = names(a2sc),connector_ids = NULL, direction = "outgoing", unique = FALSE)

## ----subset.for.mbona2sc.r, eval=F, echo=T-------------------------------
#  # Generate tracing sheet
#  tl.incoming = subset(tl.incoming,skid==names(a2sc)[1])
#  tl.outgoing = subset(tl.outgoing,skid==names(a2sc)[1])

## ----subset.ngl_id, eval=F, echo=T---------------------------------------
#  # Generate tracing sheet
#  tl.incoming = tl.incoming[!duplicated(tl.incoming$ngl_id),]
#  tl.outgoing = tl.outgoing[!duplicated(tl.outgoing$ngl_id),]

## ----assign.neuropil, eval=F, echo=T-------------------------------------
#  # Assign points to neuropils
#  pin.in = points_in_neuropil(x=tl.incoming[,c("x","y","z")],brain = elmr::FAFB14NP.surf, alpha = 30000)
#  pin.out = points_in_neuropil(x=tl.outgoing[,c("x","y","z")],brain = elmr::FAFB14NP.surf, alpha = 30000)
#  
#  # Add to data frame
#  tl.incoming$neuropil = pin.in$neuropil
#  tl.outgoing$neuropil = pin.out$neuropil
#  
#  # And now if you want, you can subset
#  tl.incoming = subset(tl.incoming, neuropil=="LH_R")
#  tl.outgoing = subset(tl.outgoing, neuropil=="LH_R")

## ----save.csv, include = FALSE, eval=F, echo=T---------------------------
#  write.csv2(x = tl.incoming,file = "/Users/abates/projects/centrifugal/data/tracing/MBONa2scRight_In.csv")
#  write.csv2(x = tl.outgoing,file = "/Users/abates/projects/centrifugal/data/tracing/MBONa2scRight_Out.csv")

## ----neuronvolume, eval=F, echo=T----------------------------------------
#  neuron = read.neurons.catmaid("1299700")
#  neuronvolume = fafb_segs_stitch_volumes(neuron = neuron, map = TRUE)

## ----neuronvolume3d, eval=F, echo=T--------------------------------------
#  nopen3d()
#  neuronvolume3d(neuronvolume)

## ----update.radii, eval=F, echo=T----------------------------------------
#  # This updates CATMAID. In this example, it updates the fafbseg `v14-seg` envrionment neuron.
#  fafbseg_update_node_radii(x = neuron, conn = fafb_seg_conn(), method = "nearest.mesh.point")
#  # This is rather slow, but you could run it overnight for a bunch of neurons of interest if you wanted
#  
#  # If you made a mistake, you can return all radii to -1 (default ) by
#  catmaid_update_radius(tnids = neuron[[1]]$d$PointNo, radii = -1)

## ----controlled.upload.to.CATMAID, eval=F, echo=T------------------------
#  uploaded = catmaid_controlled_upload(x = "name:ASB Tester",name = "ASB Tester from v14-seg",
#                                       search.range.nm = 1000, annotations = "ASB Test v14-seg Upload",
#                                       fafbseg = TRUE,  join = TRUE, join,tag = "TODO",
#                                       brain = elmr::FAFB14, lock = TRUE)

## ----auto.attach.connectors, eval=F, echo=T------------------------------
#  fafbseg_join_connectors_in_ngl_volumes("annotation:ASB Test v14-seg Upload",
#                                         putatively.connected.skids="annotation:WTPN2017_AL_PN")

## ----batch.delete, eval=F, echo=T----------------------------------------
#  catmaid_delete_neurons("annotation:ASB Test v14-seg Upload")

