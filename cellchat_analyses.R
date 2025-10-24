getCellChat = function(seurat_object, the.ccDB, the_label_column = "main.short"){
    future::plan("multisession", workers = 8) # do parallel
    data.input = as(seurat_object@assays$RNA$data, "CsparseMatrix")
    meta = seurat_object@meta.data 
    meta$labels = meta[,the_column]
    cell.use = rownames(meta)
    data.input = data.input[, cell.use]
    meta = meta[cell.use, ]
    #unique(meta$labels) # check the cell labels
    the.cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    #
    the.cellChat <- addMeta(the.cellChat, meta = meta)
    the.cellChat <- setIdent(the.cellChat, ident.use = "labels") # set "labels" as default cell identity
    #levels(the.cellChat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(the.cellChat@idents)) # number of cells in each cell group
    the.cellChat@DB <- the.ccDB
    #
    options(stringsAsFactors = FALSE)
    plan("multisession", workers = 20)
    # subset the expression data of signaling genes for saving computation cost
    the.cellChat <- subsetData(the.cellChat) # This step is necessary even if using the whole database
    the.cellChat <- identifyOverExpressedGenes(the.cellChat)
    the.cellChat <- identifyOverExpressedInteractions(the.cellChat)
    the.cellChat <- computeCommunProb(the.cellChat, type = "triMean")
    the.cellChat <- filterCommunication(the.cellChat, min.cells = 10)
    the.cellChat <- computeCommunProbPathway(the.cellChat)
    #
    the.cellChat <- aggregateNet(the.cellChat)
    the.cellChat <- netAnalysis_computeCentrality(the.cellChat, slot.name = "netP")
    return(the.cellChat)
}


# example usage
# start from a Seurat Object integrated already.
# always save object at each step! In case of error in pipeline, you can load the objects again to save time!
# that is there object list 1-4 and final!

the.seu.objs = lapply(c("gControl","gSema4ab"), function(x){
    my.obj = subset(seurat_object, sample == x)
    return(my.obj)
})
names(the.seu.objs) = c("gControl","gSema4ab")

# I suggest to run this part as a job as it may take 24 hours!
# ccaWprefix; is a column in meta.data by adding prefixes to cca_clusters.


library(NMF)
library(ggalluvial)
object.list = foreach(trt = c("gControl","gSema4ab")) %dopar%{
    if(!file.exists(paste0("rdsFiles/",trt,"_the.cellchat.objs.rds"))){
        the.cellchat.obj = getCellChat(the.seu.objs[[trt]], CellChatDB.use, the_label_column = "ccaWprefix")
        saveRDS(the.cellchat.obj, file = paste0("rdsFiles/",trt,"_the.cellchat.objs.rds"))
    }else{
        the.cellchat.obj = readRDS(paste0("rdsFiles/",trt,"_the.cellchat.objs.rds"))
    }
    return(the.cellchat.obj)
}
names(object.list) = c("gControl","gSema4ab")


if(!file.exists("object.list2.rds")){
    setCPU(20)
    the_names = names(object.list)
    object.list = foreach(trt = names(object.list)) %dopar% {
        cellchat = object.list[[trt]]
        pdf(paste0("SelectK",trt,"_Outgoing.pdf"), width = 12, height = 5)
        print(selectK(cellchat, pattern = "outgoing"))
        dev.off()
        #object.list[[trt]] = cellchat
        return(cellchat)
    }
    names(object.list) = the_names
    
    the_names = names(object.list)
    object.list = foreach(trt = names(object.list)) %dopar%{
        cellchat = object.list[[trt]]
        pdf(paste0("SelectK",trt,"_Incoming.pdf"), width = 12, height = 5)
        print(selectK(cellchat, pattern = "incoming"))
        dev.off()
        #object.list[[trt]] = cellchat
        return(cellchat)
    }
    names(object.list) = the_names
    saveRDS(object.list, "object.list2.rds")
}else{
    object.list = readRDS("object.list2.rds")
}

# run this after analyses of incoming and outgoing patterns calculated above (SelectK)!
iPatterns = c(4,3)
oPatterns = c(5,3)

ii = 1
for(trt in names(object.list)){
    pdf(paste0("OutGoing_CommPropPattern_",trt,".pdf"), width = 6, height = 24)
    object.list[[trt]] <- identifyCommunicationPatterns(object.list[[trt]], pattern = "outgoing", k = oPatterns[ii], width = 6,height = 28)   # Compute the network centrality scores
    dev.off()
    pdf(paste0("InComing_CommPropPattern_",trt,".pdf"), width = 6, height = 24)
    object.list[[trt]] <- identifyCommunicationPatterns(object.list[[trt]], pattern = "incoming", k = iPatterns[ii], width = 6,height = 28)   # Compute the network centrality scores
    dev.off()
    ii = ii + 1
}
saveRDS(object.list, "object.list2.rds")

future::plan("multisession", workers = 4) 

the_names = names(object.list)
object.list = foreach(trt = names(object.list)) %dopar% {
    cellchat = object.list[[trt]]
    cellchat <- computeNetSimilarity(cellchat, type = "functional")
    cellchat <- netEmbedding(cellchat, type = "functional")
    #> Manifold learning of the signaling networks for a single dataset
    cellchat <- netClustering(cellchat, type = "functional")
    #> Classification learning of the signaling networks for a single dataset
    # Visualization in 2D-space
    netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
    #object.list[[trt]] = cellchat
    return(cellchat)
}
names(object.list) = the_names
saveRDS(object.list, "object.list3.rds")

the_names = names(object.list)
object.list = foreach(trt = names(object.list)) %dopar% {
    cellchat = object.list[[trt]]
    cellchat <- computeNetSimilarity(cellchat, type = "structural")
    cellchat <- netEmbedding(cellchat, type = "structural")
    #> Manifold learning of the signaling networks for a single dataset
    cellchat <- netClustering(cellchat, type = "structural")
    #> Classification learning of the signaling networks for a single dataset
    # Visualization in 2D-space
    netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
    #object.list[[trt]] = cellchat
    return(cellchat)
}
names(object.list) = the_names
saveRDS(object.list, "object.list4.rds")

#group.new = unique(c(levels(object.list[["gControl"]]@idents),
#table(levels(object.list[["gControl"]]@idents) %in% group.new)
#table(levels(object.list[["gSema4abR1"]]@idents) %in% group.new)
#table(levels(object.list[["gSema4abR2"]]@idents) %in% group.new)
group.new = levels(obj1.3$ccaWprefix)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
#weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets


the_names = names(object.list)
object.list = foreach(trt = names(object.list)) %dopar%{
    the.cellChat = object.list[[trt]]
    the.cellChat = aggregateNet(the.cellChat)
    #object.list[[i]] = the.cellChat
    return(the.cellChat)
}
names(object.list) = the_names

the_group1 = levels(obj1.3$ccaWprefix)
the_group = gsub("#[0-9]+","",the_group1)
the_group[the_group == "c"] = the_group1[the_group == "c"]
the_group
names(the_group) = the_group1
the_group

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], label.size = 6,dot.size = 10, show.legend = T, #group = the_group,
                title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,90)+ylim(0,70)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
oo(20,10)
patchwork::wrap_plots(plots = gg)


print(".. saving ...")
saveRDS(object.list, "object.list.final.rds")
print(".. done ..")

#group.cellType <- unique(seurat_object$main.subs)
group.cellType <- factor(group.new, levels = group.new)
group.cellType
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

# perform down stream analyses and plotting using CelllChat tutorial!
