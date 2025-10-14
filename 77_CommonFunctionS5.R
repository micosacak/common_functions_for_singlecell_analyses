###### COMMON FUNCTIONS

crtDfull <- function(the_fullpath, current_path = "") {
    path_parts <- unlist(strsplit(the_fullpath, "/"))
    for(part in path_parts){
      current_path <- if (current_path == "") part else paste0(current_path, "/", part)
      if (!dir.exists(current_path)) {
        dir.create(current_path)
      }
    }
    return(the_fullpath)
}
tff = function(input_file, width = 5, height = 5, res = 300, units = "in", out = NULL){
    if(!is.null(out)) input_file = paste0(out, "/", input_file)
    tiff(input_file, width = width, height = height, res = res, units = units)
}
crtD = function(theDir = NULL, setDir = F){
    if(!is.null(theDir)){
        dir.create(theDir, showWarnings = F)
        if(setDir){
            setwd(theDir)
            dir.create("rdsFiles", showWarnings = F)
        }
    }else{
        theDir = gsub("[- :\\.]","",as.character(Sys.time()))
        crtD(theDir)
        if(setDir){
            setwd(theDir)
            dir.create("rdsFiles", showWarnings = F)
        }
    }
}
getQuantiles = function(the.values, the.groups, prb.values = 0.98){
    the.name = paste0(prb.values*100,"%")
    quantiles_list = list()
    unique_groups = unique(the.groups)
    indexes_above = c() 
    #
    for (group in unique_groups) {
        group_values = the.values[the.groups == group]
        quantiles = quantile(group_values, probs = prb.values)
        quantiles_list[[group]] = quantiles
        above_quant = group_values > quantiles[the.name]
        indexes_above = c(indexes_above, above_quant)
    }
    return(list(the.quants = quantiles_list, the.index = indexes_above))
}
getVln <- function(data, x_var = NULL, y_var = NULL, ...) {
  if (is.null(x_var) & is.null(y_var)) {
    x_var <- names(data)[1]  # Default to the first column for x variable
    y_var <- names(data)[2]  # Default to the second column for y variable
  } else if (is.null(x_var)) {
    x_var <- names(data)[1]  # Default to the first column for x variable
  } else if (is.null(y_var)) {
    y_var <- names(data)[2]  # Default to the second column for y variable
  }
  
  # Create violin plot
  ggplot(data, aes_string(x = x_var, y = y_var, fill = x_var)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), color = "black") +
    geom_jitter(...) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}
ldFnx = function(nCPUs = 120, mRAM = 8000) { ##
    set.seed(123)
    source("/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/77_CommonFunctionS5.R")
    setwd("/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/002_AnalysesOutputs")
    tmp_folder <<- "/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/000_jobs_commonfiles"
    dir.create(tmp_folder, showWarnings = F)
  
    nCPUs <<- nCPUs
    mRAM <<- mRAM
    library("tidyverse")
    library("foreach")
    library("doMC")  
    library("future")
    library("Seurat")
    library("Signac")
    library("SeuratWrappers")
    library("DoubletFinder")
    library("Matrix")
    library("SingleCellExperiment")
    library("scDblFinder")
    library("monocle3")
    library("CellChat")
    library("patchwork")
    library("BiocParallel")
    library("scCustomize")
    library("biomaRt")
    library("dplyr")
    library("ggforce")
    library("ggplot2")
    library("GenomeInfoDb")
    library("BSgenome.Drerio.UCSC.danRer11")
    library("JASPAR2020")
    library("TFBSTools")
    library("FigR")
    library("ComplexHeatmap")
    library("grid")
    library("motifmatchr")
    library("ggseqlogo")
    library("AnnotationHub")
    library("FNN")
    library("pheatmap")
    library("RColorBrewer")  # For default color palette
    library("igraph")
    library("ggraph")
    library("scales")
    library("ggrepel")  # For better label repulsion
    library("rtracklayer")
    library("GenomicRanges")
    library("Biostrings")
    library("Nebulosa")
    library("cowplot")
    library("viridis")
    library("data.table")
    library("pagoda2")
    library("dplyr")
    library("tidyr")
    library("ggplot2")
    library("dplyr")
    library("reshape2")
    library("tibble")
    library("org.Dr.eg.db")
    library("doParallel")
    library("cicero")
    registerDoMC(cores = nCPUs)
    plan("multicore", workers = nCPUs)
    options(stringsAsFactors = FALSE)
    #options(future.globals.maxSize = 8000*1024^2)
    #
    ldData(tmp_folder)
} 
ldData = function(tmp_folder){

  if(!file.exists(paste0(tmp_folder,"/hsa_ensDb_genes.rds"))){
    hsa_ensDb_genes = parseGTF(the_gtf_path = paste0(tmp_folder,"/Homo_sapiens.GRCh38.109.gtf.gz"), orgId = "hsa")
    saveRDS(hsa_ensDb_genes, file = paste0(tmp_folder,"/hsa_ensDb_genes.rds"))
  }else{
    print("Loading ...")
    hsa_ensDb_genes = readRDS(paste0(tmp_folder,"/hsa_ensDb_genes.rds"))
  }
  
  hsa_ensDb_genes <<- hsa_ensDb_genes

  if(!file.exists(paste0(tmp_folder,"/dre_ensDb_genes.rds"))){
    dre_ensDb_genes = parseGTF(the_gtf_path = paste0(tmp_folder,"/Danio_rerio.GRCz11.109.gtf.gz"), orgId = "dre")
    saveRDS(dre_ensDb_genes, file = paste0(tmp_folder,"/dre_ensDb_genes.rds"))
  }else{
    print("Loading ...")
    dre_ensDb_genes = readRDS(paste0(tmp_folder,"/dre_ensDb_genes.rds"))
  }
  
  dre_ensDb_genes <<- dre_ensDb_genes
  
  if(!file.exists(paste0(tmp_folder,"/hsa2mmu2dre.orthos.rds"))){
    hsa2dre = getOrthos("hsa","dre")
    hsa2mmu = getOrthos("hsa","mmu")
    mmu2dre = getOrthos("mmu","dre")
    library(dplyr)

    # Example structure assuming getOrthos returns data frames with columns like hsaENS, hsaSYM, dreENS, dreSYM
    # Ensure the column names in your actual data frames align with these expected names or adjust the names accordingly
    
    # Merge human and mouse data on human Ensembl IDs
    human_mouse <- merge(hsa2dre[, c("hsaENS", "hsaSYM", "dreENS", "dreSYM")], 
                         hsa2mmu[, c("hsaENS", "hsaSYM", "mmuENS", "mmuSYM")], 
                         by = c("hsaENS", "hsaSYM"), 
                         all = TRUE)
    
    # Rename for clarity if needed
    colnames(human_mouse)[colnames(human_mouse) == "dreENS.x"] <- "dreENS"
    colnames(human_mouse)[colnames(human_mouse) == "dreSYM.x"] <- "dreSYM"
    
    # Merge with zebrafish data on mouse Ensembl IDs
    common_data_frame <- merge(human_mouse, 
                               mmu2dre[, c("mmuENS", "mmuSYM", "dreENS", "dreSYM")], 
                               by.x = c("mmuENS", "mmuSYM"), 
                               by.y = c("mmuENS", "mmuSYM"), 
                               all = TRUE)
    
    # Clean up column names if there are duplicates or conflicts
    # This may involve selecting which columns to keep or renaming them appropriately
    common_data_frame <- common_data_frame[, !duplicated(colnames(common_data_frame))]
    
    # View the resulting data frame
    head(common_data_frame)
    hsa2mmu2dre.orthos = common_data_frame
    hsa2mmu2dre.orthos = hsa2mmu2dre.orthos[, c(1,2,3,4,5,6)]
    colnames(hsa2mmu2dre.orthos)[5:6] = c("dreENS","dreSYM")
    idx = paste(hsa2mmu2dre.orthos$mmuENS, hsa2mmu2dre.orthos$hsaENS,hsa2mmu2dre.orthos$dreENS, sep = "@")
    hsa2mmu2dre.orthos = hsa2mmu2dre.orthos[!duplicated(idx),]
    dim(hsa2mmu2dre.orthos)
    saveRDS(hsa2mmu2dre.orthos, file = paste0(tmp_folder,"/hsa2mmu2dre.orthos.rds"))
  }else{
    print("Loading ...")
    hsa2mmu2dre.orthos = readRDS(paste0(tmp_folder,"/hsa2mmu2dre.orthos.rds"))
  }

  hsa2mmu2dre.orthos <<- hsa2mmu2dre.orthos

  if(!file.exists(paste0(tmp_folder,"/the_core_pfm.rds"))){
    # Get a list of motif position frequency matrices from the JASPAR database
    the_core_pfm <- getMatrixSet(
      x = JASPAR2020,
      opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
    )
    saveRDS(the_core_pfm, file = paste0(tmp_folder,"/the_core_pfm.rds"))
  }else{
    print("Loading ...")
    the_core_pfm = readRDS(paste0(tmp_folder,"/the_core_pfm.rds"))
  }

  the_core_pfm <<- the_core_pfm

  if(!file.exists(paste0(tmp_folder,"/the_hsa_pfm.rds"))){
    the_hsa_pfm <- getMatrixSet(
      x = JASPAR2020,
      opts = list(species = 9606, all_versions = FALSE)
    )
    saveRDS(the_hsa_pfm, file = paste0(tmp_folder,"/the_hsa_pfm.rds"))
  }else{
    print("Loading ...")
    the_hsa_pfm = readRDS(file = paste0(tmp_folder,"/the_hsa_pfm.rds"))
  }

  the_hsa_pfm <<- the_hsa_pfm

  if(!file.exists(paste0(tmp_folder,"/the_mmu_pfm.rds"))){
    the_mmu_pfm <- getMatrixSet(
      x = JASPAR2020,
      opts = list(species = 10090, all_versions = FALSE)
    )
    saveRDS(the_mmu_pfm, file = paste0(tmp_folder,"/the_mmu_pfm.rds"))
  }else{
    print("Loading ...")
    the_mmu_pfm = readRDS(file = paste0(tmp_folder,"/the_mmu_pfm.rds"))
  }

  the_mmu_pfm <<- the_mmu_pfm

  if(file.exists(paste0(tmp_folder,"/annotation.dre.rds"))){
    annotation.dre  = readRDS(paste0(tmp_folder,"/annotation.dre.rds"))
  }else{
      ## Load the annotation resource.
      ah <- AnnotationHub()
      ## Query for all available EnsDb databases
      #query(ah, "EnsDb")
      #ahDb <- query(ah, pattern = c("Danio rerio", "EnsDb", 109))
      ## What have we got
      #ahDb
      annotation.dre <- ah[["AH109573"]]
      saveRDS(annotation.dre, file = paste0(tmp_folder,"/annotation.dre.rds"))
  }

  annotation.dre <<- annotation.dre

  if(file.exists(paste0(tmp_folder,"/annotation.zfx.rds"))){
      annotation.zfx = readRDS(paste0(tmp_folder,"/annotation.zfx.rds"))
  }else{
      # get gene annotations for hg38
      annotation.zfx <- GetGRangesFromEnsDb(ensdb = annotation.dre)
      saveRDS(annotation.zfx, file = paste0(tmp_folder,"/annotation.zfx.rds"))
  }

  annotation.zfx <<- annotation.zfx

  if(file.exists(paste0(tmp_folder,"/bsgenome.zfx.rds"))){
    bsgenome.zfx = readRDS(paste0(tmp_folder,"/bsgenome.zfx.rds"))
    print(bsgenome.zfx)
  }else{
      print(annotation.zfx)
      #
      nw.seq = seqlevels(annotation.zfx)
      nw.seq = gsub("^chr|^chrUn","",nw.seq)
      nw.seq = gsub("_alt$","",nw.seq)
      nw.seq = gsub("v1$",".1",nw.seq)
      nw.seq = gsub("v2$",".2",nw.seq)
      nw.seq = gsub("^_K","K",nw.seq)
      nw.seq = gsub("^M$","MT",nw.seq)
      #
      seqlevels(annotation.zfx) = nw.seq
      bsgenome.zfx = BSgenome.Drerio.UCSC.danRer11
      nw.seq = seqlevels(bsgenome.zfx)
      nw.seq = gsub("^chr|^chrUn","",nw.seq)
      nw.seq = gsub("_alt$","",nw.seq)
      nw.seq = gsub("v1$",".1",nw.seq)
      nw.seq = gsub("v2$",".2",nw.seq)
      nw.seq = gsub("^_K","K",nw.seq)
      nw.seq = gsub("^M$","MT",nw.seq)
      seqlevels(bsgenome.zfx) = nw.seq
      print(bsgenome.zfx)
      saveRDS(bsgenome.zfx,"bsgenome.zfx.rds")
  }

  bsgenome.zfx <<- bsgenome.zfx

  all_names = c(names(the_core_pfm), names(the_hsa_pfm),names(the_mmu_pfm))
  all_pfm = c(the_core_pfm, the_hsa_pfm, the_mmu_pfm)
  all_pfm = all_pfm[!duplicated(all_names)]
  length(all_pfm)
  the_pfm = getPFM(tmp_folder)
  head(the_pfm,1)
  head(all_pfm,1)
  the_pfm = all_pfm
  the_names <- unlist(lapply(the_pfm, function(x) {
    if (inherits(x, "PFMatrixList")) {
      return(x@listData$remap_tf_name)
    } else {
      return(x@tags$remap_tf_name)
    }
  }))

  the_symbols = names(the_names)
  the_orthos = hsa2mmu2dre.orthos
  the_orthos$rwns = NULL
  #dim(the_orthos)
  aa = the_names[!the_names %in% the_orthos$hsaSYM]
  aa = sort(aa)
  the_orthos = the_orthos[the_orthos$hsaSYM %in% the_names,]
  idx = match(the_orthos$hsaSYM, the_names)
  the_orthos$tfSYM = names(the_names)[idx] 

  new_tfxs = 
  "ENSDARG00000044301,atf1,ENSG00000123268,ATF1,MA0604.1
  ENSDARG00000036073,cebpg,ENSG00000153879,CEBPG,MA0838.1
  ENSDARG00000036073,cebpg,ENSG00000153879,CEBPG,MA1636.1
  ENSDARG00000037421,egr1,ENSG00000120738,EGR1,MA0162.4
  ENSDARG00000020759,elf1,ENSG00000120690,ELF1,MA0473.3
  ENSDARG00000069763,etv5a,ENSG00000244405,ETV5,MA0765.2
  ENSDARG00000044511,etv5b,ENSG00000244405,ETV5,MA0765.2
  ENSDARG00000012788,foxa3,ENSG00000170608,FOXA3,MA1683.1
  ENSDARG00000057680,foxj2,ENSG00000170608,FOXA3,MA0614.1
  ENSDARG00000004843,foxp1a,ENSG00000114861,FOXP1,MA0481.3
  ENSDARG00000014181,foxp1b,ENSG00000114861,FOXP1,MA0481.3
  ENSDARG00000013477,gata1a,ENSG00000102145,GATA1,MA0035.4
  ENSDARG00000008818,hsf1,ENSG00000185122,HSF1,MA0486.2
  ENSDARG00000017400,klf1,ENSG00000105610,KLF1,MA0493.1
  ENSDARG00000061368,klf13,ENSG00000169926,KLF13,MA0657.1
  ENSDARG00000062420,nfia,ENSG00000147862,NFIB,MA1643.1
  ENSDARG00000038687,nfkb2,ENSG00000077150,NFKB2,MA0778.1
  ENSDARG00000078280,nkx3-1,ENSG00000167034,NKX3_1,MA0124.2
  ENSDARG00000019717,pbx2,ENSG00000204304,PBX2,MA1113.2
  ENSDARG00000002779,pdx1,ENSG00000139515,PDX1,MA0132.2
  ENSDARG00000056783,raraa,ENSG00000131759,RARA,MA0159.1
  ENSDARG00000034893,rarab,ENSG00000131759,RARA,MA0159.1
  ENSDARG00000056783,raraa,ENSG00000131759,RARA,MA0730.1
  ENSDARG00000034893,rarab,ENSG00000131759,RARA,MA0730.1
  ENSDARG00000056175,scrt2,ENSG00000215397,SCRT2,MA0744.2
  ENSDARG00000040046,snai2,ENSG00000019549,SNAI2,MA0745.2
  ENSDARG00000101576,tbxta,ENSG00000164458,TBXT,MA0009.2
  ENSDARG00000039806,tbxtb,ENSG00000164458,TBXT,MA0009.2
  ENSDARG00000028159,tead1a,ENSG00000074219,TEAD2,MA1121.1
  ENSDARG00000059483,tead1b,ENSG00000074219,TEAD2,MA1121.1
  ENSDARG00000059020,thap1,ENSG00000131931,THAP1,MA0597.1
  ENSDARG00000017953,tp73,ENSG00000078900,TP73,MA0861.1
  ENSDARG00000030402,twist1a,ENSG00000122691,TWIST1,MA1123.2
  ENSDARG00000039899,zbtb7a,ENSG00000178951,ZBTB7A,MA0750.2
  ENSDARG00000074453,zfx,ENSG00000005889,ZFX,MA0146.2
  ENSDARG00000007823,atf3,ENSG00000162772,ATF3,MA0605.1
  ENSDARG00000007823,atf3,ENSG00000162772,ATF3,MA0605.2
  ENSDARG00000102899,crema,ENSG00000095794,CREM,MA0609.1
  ENSDARG00000102899,cremb,ENSG00000095794,CREM,MA0609.1
  ENSDARG00000102899,crema,ENSG00000095794,CREM,MA0609.2
  ENSDARG00000102899,cremb,ENSG00000095794,CREM,MA0609.2"
  new_tfxs = read.table(textConnection(new_tfxs), sep = ",", header = FALSE, stringsAsFactors = FALSE)
  colnames(new_tfxs) = c("dreENS","dreSYM","hsaENS","hsaSYM","tfSYM")

  is_nas = subset(the_orthos, is.na(dreENS))
  head(is_nas)

  qq = merge(is_nas, new_tfxs, by = c("hsaENS"), all = TRUE)
  qq$tfSYM.x[is.na(qq$tfSYM.x)] = qq$tfSYM.y[is.na(qq$tfSYM.x)]
  qq$tfSYM.y[is.na(qq$tfSYM.y)] = qq$tfSYM.x[is.na(qq$tfSYM.y)]
  qq$dreENS.x = NULL
  qq$dreSYM.x = NULL
  colnames(qq) = gsub("[.x|.y]","",colnames(qq))
  qq = qq[,colnames(the_orthos)]
  the_orthos = the_orthos[!is.na(the_orthos$dreENS),]
  the_orthos = rbind(the_orthos, qq)

  dre.orgDb.genes = AnnotationDbi::select(
      org.Dr.eg.db,
      keys = keys(org.Dr.eg.db),
      columns = c("ENSEMBL","SYMBOL","GENENAME","GO"))

  dre.all.trxs = subset(dre.orgDb.genes, GO == "GO:0003700")

  table(duplicated(dre.all.trxs$ENSEMBL))


  not_found = dre.all.trxs$ENSEMBL[!dre.all.trxs$ENSEMBL %in% dre_ensDb_genes$dreENS]
  not_found = not_found[!is.na(not_found)]

  #dre.all.trxs$ENSEMBL[!dre.all.trxs$ENSEMBL %in% the_orthos$dreENS]

  new_pfm = the_pfm
  new_pfm = the_pfm[names(the_pfm) %in% the_orthos$tfSYM[!is.na(the_orthos$tfSYM) & !is.na(the_orthos$dreSYM)]]
  rem_pfm = the_pfm[!names(the_pfm) %in% the_orthos$tfSYM[!is.na(the_orthos$tfSYM) & !is.na(the_orthos$dreSYM)]]
  idx = match(the_orthos$tfSYM, names(new_pfm))
  new_pfm = new_pfm[idx[!is.na(idx)]]
  names(new_pfm) = the_orthos$dreSYM[!is.na(idx)]
  table(is.na(names(new_pfm)))

  the_names <- unlist(lapply(rem_pfm, function(x) {
    if (inherits(x, "PFMatrixList")) {
      return(x@listData$remap_tf_name)
    } else {
      return(x@tags$remap_tf_name)
    }
  }))
  length(the_names)

  length(new_pfm)

  ## add a motif for sall1b

  my_mt = data.frame(
      "A" = c(0,0,0,120,0,0,0,0,0,0,0,0),
      "C" = c(0,0,0,0,0,0,0,0,0,0,0,0),
      "G" = c(0,0,0,0,0,0,0,0,0,120,0,0),
      "T" = c(0,0,120,0,120,120,120,120,120,0,0,0)
  )
  my_mt = t(my_mt)
    my_tfs_matrix = as.matrix(my_mt)
  # Create PFMatrix object
  sall1b_motif <- PFMatrix(
    ID = "MA9999.0",              # Custom MA ID
    name = "sall1b",              # Custom name
    matrixClass = "Unknown",      # Placeholder (customize if known, e.g., "Zinc finger")
    strand = "+",                 # Strand orientation
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),  # Equal background frequencies
    tags = list(
      alias = "-",                # Placeholder
      centrality_logp = "-",      # Placeholder
      description = "sall1b transcription factor",  # Custom description
      family = "Unknown",         # Placeholder (customize if known)
      medline = "-",              # Placeholder
      remap_tf_name = "sall1b",   # Same as name
      source = "Custom",          # Indicate custom source
      symbol = "sall1b",          # Same as name
      tax_group = "vertebrates",  # Placeholder (customize if known)
      tfbs_shape_id = "-",        # Placeholder
      type = "Unknown",           # Placeholder (e.g., "ChIP-seq" if applicable)
      unibind = "-",              # Placeholder
      collection = "CORE",        # JASPAR collection
      species = c("9606" = "Danio rerio"),  # Placeholder (customize if known)
      acc = "-"                   # Placeholder for accession
    ),
    profileMatrix = my_tfs_matrix  # Your position frequency matrix
  )

  #print(pfm)

  new_pfm$sall1b = sall1b_motif
  the_pfm$MA9999.1 = sall1b_motif

  new_pfm <<- new_pfm
  the_pfm <<- the_pfm
  sall1b_motif <<- sall1b_motif
  gtf_data <<- rtracklayer::import("/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/001_Alignments/Genomes/dre/22_Multiome/dre_GRCz11_atacARC_custom/genes/genes.gtf.gz")
}
setDir = function(theFolder = "LAB5061", theOutput = "20240208"){
    theFolder <<- theFolder 
    dir.create(theFolder, showWarnings = F)
    setwd(theFolder)
    theOutput <<- theOutput
    dir.create(theOutput, showWarnings = F)
    setwd(theOutput)
    dir.create("rdsFiles", showWarnings = F)
}
setRAM = function(mRAM) options(future.globals.maxSize = mRAM*1024^3)
pL = function() print("... loading ...")
setCPU = function(x) plan("multicore", workers = as.integer(x))
ldm = function(x) lapply(x, function(x) dim(x))
NL = Seurat::NoLegend()
RA = Seurat::RotatedAxis()

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
get_vCols = function(x = "D") return(scale_color_viridis_c(option = x))
oo = function(the.wd = 5, the.ht = 5){
    options(repr.plot.width=the.wd, repr.plot.height=the.ht)
}
getQuantiles = function(the.values, the.groups, prb.values = 0.98){
    the.name = paste0(prb.values*100,"%")
    quantiles_list = list()
    unique_groups = unique(the.groups)
    indexes_above = c() 
    #
    for (group in unique_groups) {
        group_values = the.values[the.groups == group]
        quantiles = quantile(group_values, probs = prb.values)
        quantiles_list[[group]] = quantiles
        above_quant = group_values > quantiles[the.name]
        indexes_above = c(indexes_above, above_quant)
    }
    return(list(the.quants = quantiles_list, the.index = indexes_above))
}
getVln <- function(data, x_var = NULL, y_var = NULL, ...) {
  if (is.null(x_var) & is.null(y_var)) {
    x_var <- names(data)[1]  # Default to the first column for x variable
    y_var <- names(data)[2]  # Default to the second column for y variable
  } else if (is.null(x_var)) {
    x_var <- names(data)[1]  # Default to the first column for x variable
  } else if (is.null(y_var)) {
    y_var <- names(data)[2]  # Default to the second column for y variable
  }
  
  # Create violin plot
  ggplot(data, aes_string(x = x_var, y = y_var, fill = x_var)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), color = "black") +
    geom_jitter(...) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}
getOrthos = function(org_id1 = "hsa", org_id2 = "dre", ensVersion = 105){
    org.info = list("dre" = "drerio", "nfu" = "nfurzeri", "mmu" = "mmusculus", "hsa" = "hsapiens")              
    org1_ensembl <- useEnsembl(biomart = 'genes', dataset = paste0(org.info[[org_id1]],'_gene_ensembl'), version = ensVersion)
    org2_ensembl <- useEnsembl(biomart = 'genes', dataset = paste0(org.info[[org_id2]],'_gene_ensembl'), version = ensVersion)
    org1_genes = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = org1_ensembl)
    org2_genes = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = org2_ensembl)
    ensList = list("org1" = org1_genes, "org2" = org2_genes)
    
    org1TOorg2.orthos = getLDS(attributes=c("ensembl_gene_id","external_gene_name"),mart=org1_ensembl,attributesL=c("ensembl_gene_id","external_gene_name"), martL=org2_ensembl)
    colnames(org1TOorg2.orthos) = c(paste0(org_id1, "ENS"), paste0(org_id1, "SYM"), paste0(org_id2, "ENS"), paste0(org_id2, "SYM"))
    idx1 = org1TOorg2.orthos[,2] == ""
    org1TOorg2.orthos[idx1, 2] = org1TOorg2.orthos[idx1, 1]
    idx2 = org1TOorg2.orthos[,4] == ""
    org1TOorg2.orthos[idx2, 4] = org1TOorg2.orthos[idx2, 3]
    org1TOorg2.orthos$rwns = paste(org1TOorg2.orthos[,1], org1TOorg2.orthos[,2], org1TOorg2.orthos[,3], org1TOorg2.orthos[,4], sep = "qq")
    return(org1TOorg2.orthos)
}

find_match_improved <- function(the_gene, the_row_names, max.distance = 0.01, topN = 10) {
    the_gene_lower <- tolower(the_gene)
    row_names_lower <- tolower(the_row_names)
    
    exact_match <- which(the_gene_lower == row_names_lower)
    if (length(exact_match) > 0) {
        return(the_row_names[exact_match])
    }
    
    prefix_matches <- grep(paste0("^", the_gene_lower), row_names_lower, ignore.case = TRUE)
    if (length(prefix_matches) > 0) {
        pos_genes <- the_row_names[prefix_matches]
        return(list("pos_genes" = pos_genes, "first_gene" = pos_genes[1]))
    }

    approx_matches <- grep(paste0("^", the_gene_lower), row_names_lower, ignore.case = TRUE)
    if (length(approx_matches) == 0) {
        matching_indices <- agrepl(paste0("^", the_gene_lower), row_names_lower, 
        ignore.case = TRUE, max.distance = max.distance)
        matching_indices <- which(matching_indices)
    } else {
        matching_indices <- approx_matches
    }
    
    if (length(matching_indices) > 0) {
        pos_genes <- the_row_names[matching_indices]
        return(list("pos_genes" = pos_genes, "first_gene" = pos_genes[1]))
    }
  
    gene_length <- nchar(the_gene_lower)
    while (gene_length >= 3) {
        shortened_gene <- substr(the_gene_lower, 1, gene_length - 1)
        matching_indices <- grep(paste0("^", shortened_gene), row_names_lower, ignore.case = TRUE)
        if (length(matching_indices) > 0) {
            pos_genes <- the_row_names[matching_indices]
            return(list("pos_genes" = pos_genes, "first_gene" = pos_genes[1]))
        }
        gene_length <- gene_length - 1
    }

    pos_genes <- grep(paste0("^", substr(the_gene_lower, 1, 3)), row_names_lower, ignore.case = TRUE)
    if (length(pos_genes) > 0) {
        pos_genes <- the_row_names[pos_genes]
        pos_genes <- unique(pos_genes)[order(nchar(pos_genes))]  # Sort by length
        if (length(pos_genes) > topN) {
            pos_genes <- head(pos_genes, topN)  # Limit to topN default 10
        }
        return(list("pos_genes" = pos_genes, "first_gene" = NULL))
    }
    
    return(NULL)
}

find_match =  function(the_gene, the_row_names, max.distance = 0.01) {
  # TO DO: improve to find the closest gene
  exact_match <- which(tolower(the_gene) == tolower(the_row_names))
  if (!identical(exact_match, integer(0))) {
    return(the_row_names[exact_match])
  } else {
    case_match <- which(tolower(the_gene) == tolower(the_row_names))
    #
    if (length(case_match) > 0) {
      return(the_row_names[case_match[1]])
    } else {
      #max.distance <- 0.01  # Adjust this value based on your desired threshold for approximate matching
      matching_indices <- agrepl(the_gene, the_row_names, ignore.case = TRUE, max.distance = max.distance)
      matching_indices <- which(matching_indices)
      #
      if (length(matching_indices) == 0) {
        gene_length <- nchar(the_gene)
        while (gene_length >= 3) {
          the_gene <- substr(the_gene, 1, gene_length - 1)
          matching_indices <- agrepl(the_gene, the_row_names, ignore.case = TRUE, max.distance = max.distance)
          matching_indices <- which(matching_indices)
          if (length(matching_indices) > 0) {
            return(list("pos_genes" = the_row_names[matching_indices], "first_gene" = the_row_names[matching_indices[1]]))
          }
          gene_length <- nchar(the_gene)
        }
        #
        pos_genes <- agrepl(paste0("^", the_gene), the_row_names, ignore.case = TRUE, max.distance = max.distance)
        pos_genes <- which(pos_genes)
        pos_genes <- unique(pos_genes)[order(nchar(the_row_names[pos_genes]))]  # Sort by gene name length
        if (length(pos_genes) > 20) {
          pos_genes <- head(pos_genes, 20)  # Limit to 20 possible gene names
        }
        return(list("pos_genes" = the_row_names[pos_genes], "first_gene" = NULL))
      } else {
        return(list("pos_genes" = the_row_names[matching_indices], "first_gene" = the_row_names[matching_indices[1]]))
      }
    }
  }
}

###### FUNCTIONS TO KEEP
get_factors = function(factor_list, level_order, sep = "") {
    lengths <- sapply(factor_list, length)
    if (length(unique(lengths)) != 1) stop("All vectors must have the same length")
    if (!all(level_order %in% names(factor_list))) stop("level_order must contain only names from factor_list")
    if (length(level_order) != length(factor_list)) stop("level_order must include all factor names")
    levels_list <- lapply(factor_list, levels)
    level_combinations <- do.call(expand.grid, levels_list[level_order])
    ordered_levels <- do.call(paste, c(level_combinations, sep = sep))
    concatenated <- do.call(paste, c(factor_list[level_order], sep = sep))
    result <- factor(concatenated, levels = ordered_levels)
    return(result)
}
get_network_data_old = function (
    figR.d, 
    score.cut = 1, 
    TFs = NULL, 
    DORCs = NULL, 
    weight.edges = FALSE, 
    TFnodecol = "Tomato", 
    DORCnodecol = "Sky Blue", 
    posEdgecol = "Forest Green", 
    negEdgecol = "Purple", 
    labelSize = 13, 
    TFnodeSize = 5, 
    DORCnodeSize = 5, 
    edgeWidthScaleMin = 0.0,
    edgeWidthScaleMax = 1,
    edgeAlpha = 0.6,  # Transparency for edge clarity
    verbose = FALSE,   # Option for debugging output
    the_layout = "tree", # options; sugiyama, linear, tree, bipartite
    the_circular = TRUE) {
     # Load required packages
    require(dplyr)
    require(igraph)
    require(ggraph)
    require(ggplot2)
    require(scales)
    require(ggrepel)  # For better label repulsion
    
    # Validate inputs
    if (is.null(TFs)){ #stop("TFs list (my_TFs) must be provided.")
        TFs = figR.d$Motif      
    }   
    if (is.null(DORCs)){ #stop("DORCs list (my_dorcs) must be provided.")
        DORCs = figR.d$DORC
    }
    
    # Check for missing values in critical columns
    if (any(is.na(figR.d$Score)) || any(is.na(figR.d$Corr))) {
        stop("Missing values detected in Score or Corr columns.")
    }
    
    # Validate color inputs
    validate_color <- function(col, name) {
        tryCatch(
            gplots::col2hex(col),
            error = function(e) stop("Invalid color for ", name, ": ", col)
        )
    }
    validate_color(TFnodecol, "TFnodecol")
    validate_color(DORCnodecol, "DORCnodecol")
    validate_color(posEdgecol, "posEdgecol")
    validate_color(negEdgecol, "negEdgecol")
    
    # Dynamic score cutoff if score.cut is NULL
    if (is.null(score.cut)) {
        score.cut <- quantile(abs(figR.d$Score), 0.9, na.rm = TRUE)
        message("Using score.cut = ", round(score.cut, 3), " (90th percentile of abs(Score)).")
    }
    
    # Filter data by score threshold and provided TFs/DORCs
    net.dat <- figR.d %>% 
        dplyr::filter(
            abs(Score) >= score.cut,
            Motif %in% TFs,  # Restrict Motif to TFs only
            DORC %in% DORCs  # DORC must be in provided DORCs
        )
    
    # Warn if no interactions remain after filtering
    if (nrow(net.dat) == 0) {
        warning("No interactions meet score.cut = ", score.cut, ". Consider lowering the threshold.")
        stop("No valid interactions found after filtering. Check TFs, DORCs, or score.cut.")
    }
    
    # Format node and edge names (add "." for consistency with original)
    net.dat$Motif <- paste0(net.dat$Motif, ".")
    net.dat$DORC <- paste0(net.dat$DORC, ".")
    
    # Define nodes: distinguish TFs and DORCs
    all_nodes <- unique(c(net.dat$Motif, net.dat$DORC))
    tf_names <- paste0(TFs, ".")
    nodes <- data.frame(
        name = all_nodes,
        group = ifelse(all_nodes %in% tf_names, "TF", "DORC"),
        size = ifelse(all_nodes %in% tf_names, TFnodeSize, DORCnodeSize)
    )
    
    # Debugging output for nodes
    if (verbose) {
        #message("Node classifications:")
        #print(nodes)
    }
    
    # Define edges
    edges <- as.data.frame(net.dat)
    links <- data.frame(
        from = edges$Motif,
        to = edges$DORC,
        corr = edges$Corr,
        enrichment = edges$Enrichment.P,
        weight = if(weight.edges) scales::rescale(log1p(abs(edges$Score))) * edgeWidthScaleMax else 1
    )
    
    # Debugging output for edges
    if (verbose) {
        #message("Edge correlations:")
        #print(links[, c("from", "to", "corr")])
    }
    
    # Create igraph object
    g <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)
    return(list("g" = g, "edges" = edges, "links" = links, "nodes" = nodes))
}
DoHeatmap_custom <- function(the.obj, 
                        features, 
                        group.by = "seurat_clusters", 
                        color = colorRampPalette(c("blue", "white", "red"))(100), 
                        scale = "row", 
                        title = "Heatmap with Gene Groups",
                        color_limits = NULL,
                        cluster_colors = NULL,
                        row_names_bold = TRUE,
                        use_raster = FALSE) {
  # Input validation
  if (!inherits(the.obj, "Seurat")) stop("the.obj must be a Seurat object")
  if (!is.list(features)) stop("features must be a named list of gene vectors")
  if (!all(sapply(features, is.character))) stop("Each element in features must be a character vector of genes")
  if (!group.by %in% colnames(the.obj@meta.data)) stop("group.by not found in Seurat object metadata")
  
  # Combine all genes into a single vector
  all_genes <- unlist(features)
  fontface_row <- if (row_names_bold) "bold" else "plain"
  
  # Check for missing genes
  scaled_data <- GetAssayData(the.obj, slot = "scale.data")
  missing_genes <- all_genes[!all_genes %in% rownames(scaled_data)]
  if (length(missing_genes) > 0) {
    warning("The following genes are missing in the Seurat object: ", paste(missing_genes, collapse = ", "))
    all_genes <- all_genes[all_genes %in% rownames(scaled_data)]
  }
  
  # Exit if no valid genes remain
  if (length(all_genes) == 0) stop("No valid genes found in the Seurat object")
  
  # Extract scaled data for the selected genes
  heatmap_data <- scaled_data[all_genes, , drop = FALSE]
  
  # Define gaps between gene groups (cumulative number of genes)
  gaps_row <- cumsum(lengths(features)[-length(features)])
  
  # Extract cell annotations and ensure group.by is treated as factor
  cell_metadata <- as.factor(the.obj@meta.data[[group.by]])
  names(cell_metadata) <- colnames(the.obj)
  
  # Order cells by factor levels of group.by
  unique_levels <- levels(cell_metadata)  # Get factor levels
  cell_order <- order(factor(cell_metadata, levels = unique_levels))
  heatmap_data <- heatmap_data[, cell_order, drop = FALSE]  # Reorder columns by factor levels
  
  # Create column annotation for cell groups
  cell_annotation <- data.frame(
    Cluster = cell_metadata[cell_order],
    row.names = colnames(heatmap_data)
  )
  
  # Define gaps between cell groups (based on factor levels)
  group_counts <- table(factor(cell_metadata, levels = unique_levels))[unique_levels]
  gaps_col <- cumsum(as.numeric(group_counts)[-length(group_counts)])
  
  # Create column labels: one label per cluster, centered over the group
  col_labels <- rep("", ncol(heatmap_data))  # Initialize empty labels
  cum_counts <- cumsum(as.numeric(group_counts))
  start_counts <- c(1, cum_counts[-length(cum_counts)] + 1)  # Start index of each group
  label_positions <- floor((start_counts + cum_counts) / 2)  # Center of each group
  col_labels[label_positions] <- unique_levels  # Place label at center of each group
  
  # Apply custom color limits if provided
  breaks <- NA
  if (!is.null(color_limits)) {
    if (!is.numeric(color_limits) || length(color_limits) != 2) {
      stop("color_limits must be a numeric vector of length 2 (e.g., c(-1.5, 1))")
    }
    if (scale == "none") {
      # Use breaks directly when scale = "none"
      breaks <- seq(color_limits[1], color_limits[2], length.out = length(color) + 1)
    } else {
      # For scaled data (row or column), adjust breaks to match scaled range
      scaled_range <- range(heatmap_data, na.rm = TRUE)
      breaks <- seq(color_limits[1], color_limits[2], length.out = length(color) + 1)
      warning("Using color_limits with scale = '", scale, "'. Ensure color_limits match the scaled data range (", 
              scaled_range[1], " to ", scaled_range[2], ").")
    }
  }
  
  # Define cluster colors for column annotations
  annotation_colors <- NULL
  if (!is.null(cluster_colors)) {
    # Validate cluster_colors
    if (!is.vector(cluster_colors) || !all(unique_levels %in% names(cluster_colors))) {
      stop("cluster_colors must be a named vector with colors for each level of group.by: ", paste(unique_levels, collapse = ", "))
    }
    annotation_colors <- list(Cluster = cluster_colors[unique_levels])
  } else {
    # Use default color palette (RColorBrewer "Set2" or similar)
    n_clusters <- length(unique_levels)
    default_colors <- brewer.pal(min(n_clusters, 8), "Set2")
    if (n_clusters > 8) {
      default_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_clusters)
    }
    annotation_colors <- list(Cluster = setNames(default_colors, unique_levels))
  }
  
  # Plot heatmap using pheatmap
  pheatmap(
    heatmap_data,
    cluster_rows = FALSE,  # Preserve gene order
    cluster_cols = FALSE,  # Disable column clustering
    annotation_col = cell_annotation,  # Cell group annotation colors
    annotation_colors = annotation_colors,  # Custom or default cluster colors
    gaps_row = gaps_row,  # White space between gene groups
    gaps_col = gaps_col,  # White space between cell groups
    scale = scale,  # Scale data (row, column, or none)
    color = color,  # Color palette for heatmap
    breaks = if (all(is.na(breaks))) NULL else breaks,  # Use breaks only if defined
    show_rownames = TRUE,
    show_colnames = TRUE,  # Enable column names
    labels_col = col_labels,  # One label per cluster, centered
    angle_col = "45",  # Rotate column labels by 45 degrees
    row_names_side = "left",  # Gene names on the left side
    main = title,
    fontsize_row = 8,
    fontsize_col = 8,
    fontface_row = fontface_row,
    use_raster = use_raster,
    annotation_names_col = FALSE  # Hide annotation name
  )
}
getCellMatrixH5 = function(the_counts){
    out <- emptyDrops(the_counts)
    is.cell <- out$FDR <= 0.001 & !is.na(out$FDR)
    sum(is.cell, na.rm=TRUE)
    the_counts = the_counts[,is.cell[!is.na(is.cell)]]
    return(list("the_counts" = the_counts, "out" = out))
}

###### scCUSTOMIZE Plotting customized

get_DP = function(the_gene){
    
    # Function to extract the max density value from a density plot
    extract_max_density <- function(plot) {
      ggplot_build(plot)$data[[1]]$density %>% max()
    }
    # Generate density plots for each sample
    p0 <- Plot_Density_Custom(seurat_object = subset(obj1.3, orig.ident == "ControlR4"), features =  the_gene, aspect_ratio = 1)
    p1 <- Plot_Density_Custom(seurat_object = subset(obj1.3, orig.ident == "Sema4abR3"), features = the_gene, aspect_ratio = 1)
    p2 <- Plot_Density_Custom(seurat_object = subset(obj1.3, orig.ident == "Sema4abR4"), features = the_gene, aspect_ratio = 1)
    
    # Extract maximum density values
    max_density <- max(extract_max_density(p0), extract_max_density(p1), extract_max_density(p2))
    
    # Update plots to have a consistent color scale limit
    p0 <- p0 + scale_fill_viridis_c(limits = c(0, max_density))
    p1 <- p1 + scale_fill_viridis_c(limits = c(0, max_density))
    p2 <- p2 + scale_fill_viridis_c(limits = c(0, max_density))
    
    oo(21,7)
    # Combine the plots using cowplot
    combined_plot <- plot_grid(p0, p1, p2, ncol = 3)
    return(combined_plot)
}
get_DP1 <- function(the_gene, grouping_var) {
    # Function to extract the max density value from a density plot
    extract_max_density <- function(plot) {
      ggplot_build(plot)$data[[1]]$density %>% max()
    }
  # Get unique levels of the grouping variable
  the_levels <- levels(obj1.3[[grouping_var]])
  
  # Generate density plots for each level
  plots <- lapply(the_levels, function(level) {
    subset_obj <- subset(obj1.3, subset = eval(as.name(grouping_var)) == level)
    Plot_Density_Custom(seurat_object = subset_obj, features = the_gene, aspect_ratio = 1)
  })
  
  # Extract maximum density values
  max_density <- max(sapply(plots, extract_max_density))
  
  # Update plots to have a consistent color scale limit
  plots <- lapply(plots, function(plot) {
    plot + scale_fill_viridis_c(limits = c(0, max_density))
  })
  
  # Combine the plots using cowplot
  combined_plot <- plot_grid(plotlist = plots, ncol = length(levels))
  return(combined_plot)
}
getDensPlot = function(seurat_obj){
    # Extract UMAP embeddings and metadata
    umap_embeddings <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
    umap_embeddings$sample <- seurat_obj$sample  # Replace 'sample' with the relevant column name
    # Get unique sample identifiers
    sample_ids <- unique(umap_embeddings$sample)

    # Function to create density plot for each sample
    create_density_plot <- function(sample_id) {
      ggplot(umap_embeddings[umap_embeddings$sample == sample_id, ], aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(alpha = 0.3) +
        geom_density_2d() +
        ggtitle(paste("Density plot for sample:", sample_id)) +
        theme_minimal()
    }
    # Create list of density plots
    density_plots <- lapply(sample_ids, create_density_plot)
    
    # Combine all density plots into one figure
    combined_plot <- plot_grid(plotlist = density_plots, ncol = 2)
    
    # Display the combined plot
    print(combined_plot)
    return(combined_plot)
}
getDensPlot1 = function(seurat_obj){
    # Extract UMAP embeddings and metadata
    umap_embeddings <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
    umap_embeddings$sample <- seurat_obj$sample  # Replace 'sample' with the relevant column name
    # Get unique sample identifiers
    sample_ids <- unique(umap_embeddings$sample)

    # Function to create density plot for each sample
    create_density_plot <- function(sample_id) {
      ggplot(umap_embeddings[umap_embeddings$sample == sample_id, ], aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(alpha = 0.3) +
        geom_density_2d() +
        ggtitle(paste("Density plot for sample:", sample_id)) +
        theme_minimal()
    }
    # Create list of density plots
    density_plots <- lapply(sample_ids, create_density_plot)
    
    # Combine all density plots into one figure
    combined_plot <- plot_grid(plotlist = density_plots, ncol = 2)
    
    # Display the combined plot
    print(combined_plot)
    return(combined_plot)
}
getDensPlot2 <- function(seurat_obj, group_by = "sample", color_option = "C", ncol = NULL) {
  # Extract UMAP embeddings and metadata
  umap_embeddings <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_embeddings$group <- seurat_obj[[group_by, drop = TRUE]]  # Use the specified group by column
  
  # Ensure the group_by column is a factor
  umap_embeddings$group <- as.factor(umap_embeddings$group)
  
  # Get unique group identifiers
  group_ids <- levels(umap_embeddings$group)
  
  # Function to create density plot for each group with heatmap-like effect
  create_density_plot <- function(group_id) {
    ggplot(umap_embeddings[umap_embeddings$group == group_id, ], aes(x = UMAP_1, y = UMAP_2)) +
      stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE, alpha = 0.8) +
      scale_fill_viridis_c(option = color_option) +
      ggtitle(paste("Density plot for group:", group_id)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
      )
  }
  
  # Remove rows with non-finite values in UMAP embeddings
  umap_embeddings <- umap_embeddings[is.finite(rowSums(umap_embeddings[, c("UMAP_1", "UMAP_2")])), ]
  
  # Create list of density plots
  density_plots <- lapply(group_ids, create_density_plot)
  
  # Determine the number of columns for the grid if not specified
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(length(group_ids)))
  }
  
  # Combine all density plots into one figure
  combined_plot <- plot_grid(plotlist = density_plots, ncol = ncol)
  
  # Display the combined plot
  print(combined_plot)
  
  return(combined_plot)
}

####### GENE REGULATORY NETWORKS (FigR)
runGenePeakcorr_dre = function (ATAC.se, RNAmat, genome, 
        geneList = NULL, 
        windowPadSize = 50000, 
        normalizeATACmat = TRUE, 
        nCores = 4, 
        keepPosCorOnly = TRUE, 
        keepMultiMappingPeaks = FALSE, 
        n_bg = 100, 
        p.cut = NULL,
        pairsPerChunk = 10000,
        largeChunkSize = 100000
    ) {
        stopifnot(inherits(ATAC.se, "RangedSummarizedExperiment"))
        stopifnot(inherits(RNAmat, c("Matrix", "matrix")))
        if (!all.equal(ncol(ATAC.se), ncol(RNAmat))) stop("Input ATAC and RNA objects must have same number of cells")
        message("Assuming paired scATAC/scRNA-seq data ..")
        peakRanges.OG <- granges(ATAC.se)
        rownames(ATAC.se) <- paste0("Peak", 1:nrow(ATAC.se))
        ATACmat <- assay(ATAC.se)
        if (normalizeATACmat) ATACmat <- centerCounts(ATACmat)
        if (is.null(rownames(RNAmat))) stop("RNA matrix must have gene names as rownames")
        if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
            message("Peaks with 0 accessibility across cells exist ..")
            message("Removing these peaks prior to running correlations ..")
            peaksToKeep <- Matrix::rowSums(assay(ATAC.se)) != 0
            ATAC.se <- ATAC.se[peaksToKeep, ]
            ATACmat <- ATACmat[peaksToKeep, ]
            message("Important: peak indices in returned gene-peak maps are relative to original input SE")
        }
        peakRanges <- granges(ATAC.se)
        if (any(Matrix::rowSums(RNAmat) == 0)) {
            message("Genes with 0 expression across cells exist ..")
            message("Removing these genes prior to running correlations ..")
            genesToKeep <- Matrix::rowSums(RNAmat) != 0
            RNAmat <- RNAmat[genesToKeep, ]
        }
        cat("Number of peaks in ATAC data:", nrow(ATACmat), "\n")
        cat("Number of genes in RNA data:", nrow(RNAmat), "\n")
        if (!genome %in% c("hg19", "hg38", "mm10", "GRCz11")) stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
        switch(genome, 
               hg19 = {TSSg <- FigR::hg19TSSRanges}, 
               hg38 = {TSSg <- FigR::hg38TSSRanges}, 
               mm10 = {TSSg <- FigR::mm10TSSRanges}, 
               GRCz11 = {TSSg = grcZ11TSSRanges}
              )
        names(TSSg) <- as.character(TSSg$gene_name)
        if (!is.null(geneList)) {
            if (length(geneList) == 1) 
                stop("Please specify more than 1 valid gene symbol")
            if (any(!geneList %in% names(TSSg))) {
                cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
                cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
                cat("\n")
                stop()
            }
            TSSg <- TSSg[geneList]
        }
        genesToKeep <- intersect(names(TSSg), rownames(RNAmat))
        cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ", length(genesToKeep), "\n")
        RNAmat <- RNAmat[genesToKeep, ]
        TSSg <- TSSg[genesToKeep]
        TSSflank <- GenomicRanges::flank(TSSg, width = windowPadSize, both = TRUE)
        cat("\nTaking peak summits from peak windows ..\n")
        peakSummits <- resize(peakRanges, width = 1, fix = "center")
        cat("Finding overlapping peak-gene pairs ..\n")
        genePeakOv <- findOverlaps(query = TSSflank, subject = peakSummits)
        numPairs <- length(genePeakOv)
        cat("Found ", numPairs, "total gene-peak pairs for given TSS window ..\n")
        cat("Number of peak summits that overlap any gene TSS window: ", length(unique(subjectHits(genePeakOv))), "\n")
        cat("Number of gene TSS windows that overlap any peak summit: ", length(unique(queryHits(genePeakOv))), "\n\n")
        set.seed(123)
        cat("Determining background peaks ..\n")
        if (is.null(rowData(ATAC.se)$bias)) {
            if (genome %in% "hg19") myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
            if (genome %in% "mm10") myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
            if (genome %in% "hg38") myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
            if (genome %in% "GRCz11") myGenome <- bsgenome.zfx
            ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
        }
        cat("Using ", n_bg, " iterations ..\n\n")
        set.seed(123)
        bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
        cat("Computing gene-peak correlations ..\n")
        startingPoint <- 1
        chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
        chunkEnds <- chunkStarts + largeChunkSize - 1
        chunkEnds[length(chunkEnds)] <- numPairs

        dorcList <- list()
        for (i in 1:length(chunkStarts)) {
            cat("Running pairs: ", chunkStarts[i], "to", chunkEnds[i], "\n")
            if(genome != "GRCz11"){
              ObsCor <- FigR::PeakGeneCor(ATAC = ATACmat, RNA = RNAmat, OV = genePeakOv[chunkStarts[i]:chunkEnds[i]], chunkSize = pairsPerChunk, ncores = nCores, bg = bg)
            }else{
              message("Runing PeakGeneCor_dre ...")
              ObsCor <- PeakGeneCor_dre(ATAC = ATACmat, RNA = RNAmat, OV = genePeakOv[chunkStarts[i]:chunkEnds[i]], chunkSize = pairsPerChunk, ncores = nCores, bg = bg)
            }
            gc()
            dorcList[[i]] <- ObsCor
        }
        cat("\nMerging results ..\n")
        dorcTab <- bind_rows(dorcList)
        cat("Performing Z-test for correlation significance ..\n")
        permCols <- 4:(ncol(bg) + 3)
        if (keepPosCorOnly) {
            cat("Only considering positive correlations ..\n")
            dorcTab <- dorcTab %>% dplyr::filter(rObs > 0)
        }
        if (!keepMultiMappingPeaks) {
            cat("Keeping max correlation for multi-mapping peaks ..\n")
            dorcTab <- dorcTab %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs == max(rObs))
        }
        dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]
        dorcTab$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTab$Peak], "Peak", 2))
        dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[, permCols]))
        dorcTab$rBgMean <- rowMeans(dorcTab[, permCols])
        dorcTab$pvalZ <- 1 - stats::pnorm(q = dorcTab$rObs, mean = dorcTab$rBgMean, sd = dorcTab$rBgSD)
        cat("\nFinished!\n")
        if (!is.null(p.cut)) {
            cat("Using significance cut-off of ", p.cut, " to subset to resulting associations\n")
            dorcTab <- dorcTab[dorcTab$pvalZ <= p.cut, ]
        }
        dorcTab$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[dorcTab$Peak])), 
                                    paste(start(peakRanges.OG[dorcTab$Peak]), end(peakRanges.OG[dorcTab$Peak]), sep = "-"), sep = ":")
        return(as.data.frame(dorcTab[, c("Peak", "PeakRanges", "Gene", "rObs", "pvalZ")], stringsAsFactors = FALSE))
}
run_FigGRN_dre_final = function(ATAC.se = NULL, dorcMat = NULL,  rnaMat = NULL, dorcTab = NULL,  
                          DORC.knn = NULL, dorcGenes = NULL, dorcK = 30, n_bg = 100, genome = "GRCz11", nCores = 4){
    stopifnot(all.equal(ncol(dorcMat), ncol(rnaMat)))
    if (!all(c("Peak", "Gene") %in% colnames(dorcTab))) stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")
    #
    if (all(grepl("chr", dorcTab$Peak, ignore.case = TRUE))) {
        usePeakNames <- TRUE
        message("Detected peak region names in Peak field")
        if (!(all(grepl("chr", rownames(ATAC.se), ignore.case = TRUE)))) 
        stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")
        if (!all(dorcTab$Peak %in% rownames(ATAC.se))) 
        stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
    }else{
        usePeakNames <- FALSE
        message("Assuming peak indices in Peak field")
        if (max(dorcTab$Peak) > nrow(ATAC.se)) 
        stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
    }
    #
    zero_var_cols <- apply(Matrix::t(dorcMat), 2, var) == 0
    if (any(zero_var_cols)) {
      dorcMat <- dorcMat[!zero_var_cols,]
    }
    #
    if(is.null(dorcGenes)){
        dorcGenes <- rownames(dorcMat)
    }else{
        cat("Using specified list of dorc genes ..\n")
        if (!(all(dorcGenes %in% rownames(dorcMat)))) {
        cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
        cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], 
        sep = ", ")
        cat("\n")
        stop()
        }
    }
    # edited on 09.04.2025.
    # DORC.knn can be given based on the Seurat object.
    if(is.null(DORC.knn)){
        DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))), k = dorcK)$nn.index
        rownames(DORC.knn) <- rownames(dorcMat)
    }
    #
    if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
        if (genome %in% "hg19") myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        if (genome %in% "mm10") myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
        if (genome %in% "hg38") myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        if (genome %in% "GRCz11"){ 
        myGenome <- bsgenome.zfx
        print("using zebrafish genome ...")    
        }
        ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
    }
    packagePath <- find.package("FigR", lib.loc = NULL, quiet = TRUE)
    if(grepl("hg", genome)){
        print("Using human motifs ...")
        pwm <- readRDS(paste0(packagePath, "/data/cisBP_human_pfms_2021.rds"))
    }else if(grepl("mm", genome)){
        print("Using mouse motifs ...")
        pwm <- readRDS(paste0(packagePath, "/data/cisBP_mouse_pfms_2021.rds"))
    }else if(grepl("GRCz11", genome)) {
        print("Using zebrafish motifs ...")
        pwm <- readRDS(paste0(packagePath, "/data/cisBP_zebrafish_pfms_2021.rds"))
    }else{
        stop("Could not find the pwm data ...")
    }
    all_open = showConnections(all = TRUE)
    #closeAllConnections()
    if(nrow(all_open) > 4) closeAllConnections()
    #
    gc()
    pwm = new_pfm
    nCores = nCPUs
    if (all(grepl("_", names(pwm), fixed = TRUE))) names(pwm) <- FigR::extractTFNames(names(pwm))
    message("Removing genes with 0 expression across cells ..\n")
    rnaMat <- rnaMat[Matrix::rowSums(rnaMat) != 0, ]
    myGeneNames <- gsub(x = rownames(rnaMat), pattern = "-", replacement = "")
    rownames(rnaMat) <- myGeneNames
    motifsToKeep <- intersect(names(pwm), myGeneNames)
    cat("Getting peak x motif matches ..\n")
    motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se, pwms = pwm[motifsToKeep], genome = bsgenome.zfx)
    motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix)) != 0]
    cat("Determining background peaks ..\n")
    cat("Using ", n_bg, " iterations ..\n\n")
    if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
        ATAC.mat <- assay(ATAC.se)
        ATAC.mat <- cbind(ATAC.mat, 1)
        ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = ATAC.mat), rowRanges = granges(ATAC.se))
        set.seed(123)
        bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
    }else {
        set.seed(123)
        bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
    }
    # added on 09.04.2025
    if(grepl("GRCz11", genome)){
        genome(bsgenome.zfx) <- "GRCz11"
        genome(atac_se) = "GRCz11"
    }

    all_open = showConnections(all = TRUE)
    #closeAllConnections()
    if(nrow(all_open) > 4) closeAllConnections()
    #
    gc()
    nCores = nCPUs
    setCPU(nCPUs)
    gc()
    cat("Testing ", length(motifsToKeep), " TFs\n")
    cat("Testing ", nrow(dorcMat), " DORCs\n")
    
    if (nCores > 1) message("Running FigR using ", nCores, " cores ..\n")
    opts <- list()
    pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
    cl <- parallel::makeCluster(nCores)
    clusterEvalQ(cl, .libPaths())
    doSNOW::registerDoSNOW(cl)
    mZtest.list <- foreach(g = dorcGenes, .options.snow = opts, .packages = c("FigR", "dplyr", "Matrix", "Rmpfr")) %dopar% {
        DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% c(g, rownames(dorcMat)[DORC.knn[g, ]])])
        if (usePeakNames) DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks)
        mZ <- FigR::motifPeakZtest(peakSet = DORCNNpeaks, bgPeaks = bg, tfMat = assay(motif_ix))
        mZ <- mZ[, c("gene", "z_test")]
        colnames(mZ)[1] <- "Motif"
        colnames(mZ)[2] <- "Enrichment.Z"
        mZ$Enrichment.P <- 2 * pnorm(abs(mZ$Enrichment.Z), lower.tail = FALSE)
        mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
        mZ <- cbind(DORC = g, mZ)
        corr.r <- cor(dorcMat[g, ], t(as.matrix(rnaMat[mZ$Motif, ])), method = "spearman")
        stopifnot(all.equal(colnames(corr.r), mZ$Motif))
        mZ$Corr <- corr.r[1, ]
        mZ$Corr.Z <- scale(mZ$Corr, center = TRUE, scale = TRUE)[, 1]
        mZ$Corr.P <- 2 * pnorm(abs(mZ$Corr.Z), lower.tail = FALSE)
        mZ$Corr.log10P <- sign(mZ$Corr.Z) * -log10(mZ$Corr.P)
        return(mZ)
    }
    cat("Finished!\n")
    cat("Merging results ..\n")
    TFenrich.d <- do.call("rbind", mZtest.list)
    dim(TFenrich.d)
    rownames(TFenrich.d) <- NULL
    TFenrich.d <- TFenrich.d %>% dplyr::mutate(Score = sign(Corr) * as.numeric(-log10(1 - (1 - Rmpfr::mpfr(Enrichment.P, 100)) * (1 - Rmpfr::mpfr(Corr.P, 100)))))
    TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
    return(TFenrich.d)
}

PeakGeneCor_dre = function (ATAC, RNA, OV, ncores = 4, chunkSize = 200, metric = "spearman", 
    bg = NULL) 
{
    cat(as.character(Sys.time()))
    cat("\n")
    stopifnot(ncol(ATAC) == ncol(RNA))
    if (chunkSize > 10000) 
        stop("Do not specify very large chunk sizes. Please use chunkSize < 1000")
    n <- length(OV)
    starts <- seq(1, n, chunkSize)
    ends <- starts + chunkSize - 1
    ends[length(ends)] <- n
    OVd <- OV %>% as.data.frame() %>% dplyr::rename(Gene = "queryHits", 
        Peak = "subjectHits")
    chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
    time_elapsed <- Sys.time()
    cat("Running in parallel using ", ncores, "cores ..\n")
    cat("Computing observed correlations ..\n")
    corList <- pbmcapply::pbmclapply(X = chunkList, FUN = function(x) {
        FigR::chunkCore(chunk = x, A = ATAC, R = RNA, O = OVd, 
            met = metric)
    }, mc.cores = ncores)
    if (any(unlist(sapply(corList, is.null)))) {
        message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
        message("Please check to see you have enough cores/memory allocated")
        message("Also make sure you have filtered down to non-zero peaks/genes")
    }
    OVd$rObs <- unlist(corList)
    cat("Finished!\n")
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)), 
        "\n\n")
    if (!is.null(bg)) {
        n_iter <- ncol(bg)
        cat("Computing background correlations ..\n")
        time_elapsed <- Sys.time()
        bgCor <- foreach(i = 1:n_iter, .combine = "cbind", .export = c("chunkCore", 
            "t"), .packages = c("pbmcapply", "FigR", "Matrix")) %do% 
            {
                OVdBg <- OVd[, 1:2]
                OVdBg$Peak <- bg[OVdBg$Peak, i]
                bgCorList <- pbmcapply::pbmclapply(X = chunkList, 
                  FUN = function(x) {
                    chunkCore(chunk = x, A = ATAC, R = RNA, O = OVdBg, 
                      met = metric)
                  }, mc.cores = ncores)
                unlist(bgCorList)
            }
        if (sum(is.null(bgCor)) != 0 | sum(is.na(bgCor)) != 0) 
            stop("One or more of the chunk processes failed unexpectedly (returned NULL) .. Please check to see you have enough cores/m\n           emory allocated")
        time_elapsed <- Sys.time() - time_elapsed
        cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)), 
            "\n\n")
        colnames(bgCor) <- paste0("rBg", 1:ncol(bgCor))
        OVd <- cbind(OVd, bgCor)
    }
    cat(as.character(Sys.time()))
    cat("\n")
    return(OVd)
}

##### start of figR functions

complete_df <- function(df, all_cmn_genes) {
    complete <- data.frame(cmn_genes = all_cmn_genes)
    merged <- left_join(complete, df, by = "cmn_genes")
    merged$DORC <- sapply(strsplit(merged$cmn_genes, "@"), `[`, 1)
    merged$Motif <- sapply(strsplit(merged$cmn_genes, "@"), `[`, 2)
    merged$Score <- ifelse(is.na(merged$Score), runif(nrow(merged), -0.1, 0.1), merged$Score)
    merged <- merged[, c("DORC", "Motif", "Score")]
    return(merged)
}

get_tf_dorc = function(figR.d, score.cut = 1.5){
    DORCsToKeep <- figR.d %>% 
        dplyr::filter(abs(Score) >= score.cut) %>% 
        pull(DORC) %>% 
        unique()
    #      
    TFsToKeep <- figR.d %>% 
        dplyr::filter(abs(Score) >= score.cut) %>% 
        pull(Motif) %>% 
        unique()
    figR.d = figR.d[figR.d$DORC %in% DORCsToKeep & figR.d$Motif %in% TFsToKeep,]
    return(unique(paste(figR.d$DORC,figR.d$Motif, sep = "@")))
}

process_df <- function(df, all_cmn_genes) {
    # Subset the data frame to include only DORC@Motif pairs in all_cmn_genes
    df$cmn_genes <- paste(df$DORC, df$Motif, sep = "@")
    df_subset <- df[df$cmn_genes %in% all_cmn_genes, ]
    
    # Use complete_df to ensure all_cmn_genes are present and fill missing scores
    complete <- complete_df(df_subset, all_cmn_genes)
    return(complete)
}

df_to_matrix <- function(dff) {
  net <- dff %>%
    reshape2::dcast(DORC ~ Motif) %>%
    tibble::column_to_rownames("DORC") %>%
    as.matrix()
  return(net)
}

plotMarker2D_custom = function (df, markers, markerMat, pointSize = 0.5, plotClean = TRUE, 
    rasteRize = TRUE, splitBy = NULL, minCutoff = NULL, maxCutoff = NULL, 
    colorPalette = "solar_extra", legend.position = "bottom", 
    showDataRange = TRUE, combine = TRUE, colorRange = "0-1", ...) 
{
    stopifnot(class(df) == "data.frame")
    if (!is.null(splitBy)) 
        if (!all(splitBy %in% colnames(df))) 
            stop("One or more of the splitBy variables is not a column in the provided dataframe\n")
    dimnames <- colnames(df)[1:2]
    if (is.null(markers)) {
        stop("Must provide at least 1 marker to plot")
    }
    else {
        if (!(all(markers %in% rownames(markerMat)))) {
            message("One or more of the features supplied is not present in the matrix provided: \n")
            message(markers[!markers %in% rownames(markerMat)], 
                sep = ", ")
            cat("\n")
            stop()
        }
    }
    markerMat <- markerMat[markers, , drop = FALSE]
    if (length(markers) > 20 & combine) 
        stop("Plotting this many features at once can lead to squished graphics ..\n")
    
    # Validate colorRange
    if (!colorRange %in% c("0-1", "1-100")) 
        stop("colorRange must be either '0-1' or '1-100'\n")
    
    transformScores <- function(markerScore, minCutoff, maxCutoff, colorRange) {
        # Apply cutoffs
        if (!is.null(minCutoff)) {
            if (grepl(pattern = "^q.*", minCutoff)) {
                quantile.cutMin <- as.numeric(strsplit(minCutoff, 
                  "q", fixed = TRUE)[[1]][2])
                if (!quantile.cutMin > 0 & quantile.cutMin < 1) 
                  stop("Must provide a fractional value to use for quantile cutoff. See help page for details\n")
                minCutoff <- quantile(markerScore, quantile.cutMin)
            }
            if (grepl(pattern = "^sd.*", minCutoff)) {
                sd.cutMin <- as.numeric(strsplit(minCutoff, "sd", 
                  fixed = TRUE)[[1]][2])
                if (sd.cutMin > 0 & sd.cutMin < 1) 
                  stop("SD cut-off must be an integer number representing the number of standard deviations away from the mean to use. See help page for details\n")
                minCutoff <- mean(markerScore) - (sd.cutMin * 
                  sd(markerScore))
            }
        }
        else {
            minCutoff <- min(markerScore)
        }
        if (!is.null(maxCutoff)) {
            if (grepl(pattern = "^q.*", maxCutoff)) {
                quantile.cutMax <- as.numeric(strsplit(maxCutoff, 
                  "q", fixed = TRUE)[[1]][2])
                if (!quantile.cutMax > 0 & quantile.cutMax < 1) 
                  stop("Must provide a fractional value to use for quantile cutoff. See help page for details\n")
                maxCutoff <- quantile(markerScore, quantile.cutMax)
            }
            if (grepl(pattern = "^sd.*", maxCutoff)) {
                sd.cutMax <- as.numeric(strsplit(maxCutoff, "sd", 
                  fixed = TRUE)[[1]][2])
                if (sd.cutMax > 0 & sd.cutMax < 1) 
                  stop("SD cut-off must be an integer number representing the number of standard deviations away from the mean to use. See help page for details\n")
                maxCutoff <- mean(markerScore) + (sd.cutMax * 
                  sd(markerScore))
            }
        }
        else {
            maxCutoff <- max(markerScore)
        }
        
        # Clip values
        markerScore[markerScore <= minCutoff] <- minCutoff
        markerScore[markerScore >= maxCutoff] <- maxCutoff
        
        # Normalize to 0-1
        if (maxCutoff != minCutoff) {
            markerScore <- (markerScore - minCutoff) / (maxCutoff - minCutoff)
        } else {
            markerScore[] <- 0 # If min=max, set all to 0 to avoid division by zero
        }
        
        # Scale to 1-100 if specified
        if (colorRange == "1-100") {
            markerScore <- markerScore * 99 + 1 # Maps 0-1 to 1-100
        }
        
        markerScore
    }
    
    gglist <- lapply(markers, function(i) {
        cat("Plotting ", i, "\n")
        mScore <- transformScores(markerScore = markerMat[i, ], 
                                 minCutoff = minCutoff, 
                                 maxCutoff = maxCutoff, 
                                 colorRange = colorRange)
        if (sum(mScore) == 0) 
            warning("No counts detected for marker: ", i, "\n")
        if (anyNA(mScore)) 
            warning("NAs detected for marker: ", i, "\n")
        i <- gsub(pattern = "-", replacement = "", x = i)
        df[, i] <- mScore
        if (legend.position != "bottom") {
            if (!legend.position %in% c("top", "right", "left", "none")) 
                stop("Must specify a valid legend position specification .. See function options for details\n")
        }
        if (rasteRize) {
            require(ggrastr)
            pointfun <- geom_point_rast
        }
        else {
            pointfun <- geom_point
        }
        
        # Set color scale with fixed range
        if (colorRange == "0-1") {
            col_limits <- c(0, 1)
            col_breaks <- if (showDataRange) scales::breaks_pretty(n = 3)(c(0, 1)) else c(0, 1)
            col_labels <- if (showDataRange) NULL else c("min", "max")
        } else {
            col_limits <- c(1, 100)
            col_breaks <- if (showDataRange) scales::breaks_pretty(n = 3)(c(1, 100)) else c(1, 100)
            col_labels <- if (showDataRange) NULL else c("min", "max")
        }
        
        myColfun <- scale_color_gradientn(
            colours = BuenColors::jdb_palette(colorPalette), 
            na.value = "orange", 
            limits = col_limits, 
            breaks = col_breaks, 
            labels = col_labels
        )
        
        g <- ggplot(df, aes_string(dimnames[1], dimnames[2], color = i)) + 
            pointfun(size = pointSize) + 
            theme_classic() + 
            labs(title = i) + 
            theme(plot.title = element_text(hjust = 0.5), 
                  legend.position = legend.position, 
                  legend.key.size = unit(0.3, "cm"), 
                  legend.text = element_text(size = 5), 
                  legend.title = element_blank()) + 
            myColfun
        if (!is.null(splitBy)) 
            g <- g + facet_wrap(as.formula(paste("~", splitBy))) + 
                theme(strip.background = element_blank())
        if (plotClean) 
            g <- g + theme(axis.line = element_blank(), 
                           axis.ticks = element_blank(), 
                           axis.title = element_blank(), 
                           axis.text = element_blank())
        g
    })
    
    if (length(gglist) == 1) {
        return(gglist[[1]])
    }
    else if (combine) {
        cat("Merging plots ..\n")
        return(cowplot::plot_grid(plotlist = gglist, align = "hv", ...))
    }
    else {
        cat("Returning list of plots\n")
        names(gglist) <- markers
        return(gglist)
    }
}

get_figR_rds_files <- function(the_dir) {
    if (!dir.exists(the_dir)) {
        warning("Directory does not exist: ", the_dir)
        return(NULL)
    }
    the_files <- dir(the_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(the_files) == 0) {
    return(NULL)
    }
    the_rds <- lapply(the_files, readRDS)
    names(the_rds) <- basename(the_files)
    return(the_rds)
}

get_figR_folders <- function(the_folder, thePattern = NULL) {
    if (!dir.exists(the_folder)) {
        warning("Directory does not exist: ", the_folder)
        return(NULL)
    }
    # Get all entries in the folder
    the_entries <- dir(the_folder, full.names = TRUE)
    if (length(the_entries) == 0) {
        return(NULL)
    }
    # Process each entry
    result <- lapply(the_entries, 
        function(entry){
            # Check if the entry is a directory
            if (file.info(entry)$isdir){
                # Try to read .rds files directly
                rds_files <- get_figR_rds_files(entry)
                if(!is.null(thePattern)) rds_files <- rds_files[grepl(thePattern, names(rds_files))]
                if (!is.null(rds_files)) {
                    return(rds_files) # Return .rds files if found
                } else {
                # Recursively process sub-folders
                    return(get_figR_folders(entry))
                }
            } else {
                return(NULL) # Skip non-directories
            }
        })
    # Name the results by entry names
    names(result) <- basename(the_entries)
    # Filter out NULL entries (non-folders or empty folders)
    result <- result[!sapply(result, is.null)]
    # Return NULL if no valid results
    if (length(result) == 0) {
        return(NULL)
    }
    return(result)
}

# these three functions, will extract DORC regions, and plot motifs on the genomic regions
getDorcPeaks <- function(the_tfenrich = NULL, the_footprints = NULL, the_genome = NULL, dorc_name = "lin28a", score.cut = 1) {
  if (is.null(the_tfenrich) || is.null(the_footprints) || is.null(the_genome)) {
    stop("All input parameters (the_tfenrich, the_footprints, the_genome) must be provided.")
  }
  
  dorc_peaks <- the_footprints[the_footprints$Gene == dorc_name, ]$PeakRanges
  if (length(dorc_peaks) == 0) {
    stop(paste("No peaks found in footprints for DORC:", dorc_name))
  }    
  message("Found ", length(dorc_peaks), " peaks for DORC: ", dorc_name)
  message("First few peak ranges:")
  message(head(dorc_peaks))
  
  tf_hits <- the_tfenrich[the_tfenrich$DORC == dorc_name & the_tfenrich$Enrichment.P < 0.05 & abs(the_tfenrich$Score) >= score.cut, ]
  the_neg <- tf_hits$Motif[tf_hits$Score < 0]
  valid_tfs <- unique(tf_hits$Motif)
  if (length(valid_tfs) == 0) {
    stop(paste("No significantly enriched TFs found for DORC:", dorc_name))
  }
  message("Enriched TFs: ", paste(valid_tfs, collapse = ", "))
  
  peak_granges <- Signac::StringToGRanges(dorc_peaks, sep = c(":", "-"))
  message("Peak GRanges: ", paste(Signac::GRangesToString(peak_granges), collapse = ", "))
  
  peak_sequences <- Biostrings::getSeq(the_genome, peak_granges)
  names(peak_sequences) <- Signac::GRangesToString(peak_granges)
  if (length(peak_sequences) == 0) {
    stop("Failed to retrieve sequence for peak.")
  }
  message("Retrieved sequence for peak, length: ", length(peak_sequences[[1]]), " bp")
  
  return(list("valid_tfs" = valid_tfs, "peak_sequences" = peak_sequences, "peak_granges" = peak_granges, "the_neg" = the_neg))
}
# Modified plot_motif_hits function
plot_motif_hits <- function(hits, peak_granges, the_neg = NULL, dorc_name = "lin28a") {
  parsePeaks <- function(x) {
    aa <- unlist(strsplit(x, split = "-"))
    aa[1] <- as.integer(aa[1])
    aa[2] <- as.integer(aa[2])
    aa[3] <- as.integer(aa[3])
    aa
  }
  
  peak <- peak_granges[1]
  peak_str <- Signac::GRangesToString(peak)
  peak_start <- start(peak)
  peak_end <- end(peak)
  
  # Subset hits to this peak
  hits_in_peak <- subsetByOverlaps(hits, peak)
  if (length(hits_in_peak) == 0) {
    message("No motif hits for peak: ", peak_str, ". Skipping plot.")
    return(NULL)
  }
  
  # Prepare data for plotting
  df <- as.data.frame(hits_in_peak)
  df$start_rel <- df$start - peak_start
  df$end_rel <- df$end - peak_start
  df$midpoint <- df$start_rel + (df$end_rel - df$start_rel) / 2
  
  # Peak line data
  peak_df <- data.frame(
    start_rel = 0,
    end_rel = peak_end - peak_start,
    y = "Peak"
  )
  
  # Grid data for ticks on peak line
  grid_positions <- seq(0, peak_end - peak_start, by = 100)
  grid_df <- data.frame(
    x = grid_positions,
    y = "Peak",
    label = as.character(grid_positions)
  )
  
  # Define TFs to color
  if (is.null(the_neg)) {
    colored_tfs <- c("fosl2", "fosab", "fosl1a", "fosl1b","jun","junba","junbb","jund", "tead1a","tead1b")
    tf_colors <- setNames(rep("grey10", length(unique(df$TF))), unique(df$TF))
    tf_colors[colored_tfs[colored_tfs %in% unique(df$TF)]] <- c("red", "red","red","red", "cyan3","cyan3","cyan3","purple", "purple")[1:sum(colored_tfs %in% unique(df$TF))]
  } else {
    colored_tfs <- c("fosl2", "fosab", "fosl1a", "fosl1b","jun","junba","junbb","jund", "tead1a","tead1b")
    tf_colors <- setNames(rep("grey10", length(unique(df$TF))), unique(df$TF))
    tf_colors[the_neg[the_neg %in% unique(df$TF)]] <- "gray80"
    tf_colors[colored_tfs[colored_tfs %in% unique(df$TF)]] <- c("red", "red","red","red", "cyan3","cyan3","cyan3","purple", "purple")[1:sum(colored_tfs %in% unique(df$TF))]
  }


  # Define TFs and their fixed color mapping
  colored_tfs <- c("fosl2", "fosab", "fosl1a", "fosl1b", "jun", "junba", "junbb", "jund", "tead1a", "tead1b")
  # Fixed color assignments for each TF
  tf_color_map <- c(
    fosl2 = "red", fosab = "red", fosl1a = "red", fosl1b = "red",
    jun = "cyan3", junba = "cyan3", junbb = "cyan3",
    jund = "purple", tead1a = "purple", tead1b = "purple"
  )

  # Initialize tf_colors with default color (grey10) for all unique TFs in df$TF
  tf_colors <- setNames(rep("grey10", length(unique(df$TF))), unique(df$TF))

  # If the_neg is not NULL, assign gray80 to specified TFs
  if (!is.null(the_neg)) {
    tf_colors[the_neg[the_neg %in% unique(df$TF)]] <- "gray80"
  }

  # Assign fixed colors to colored_tfs present in df$TF
  present_tfs <- colored_tfs[colored_tfs %in% unique(df$TF)]
  tf_colors[present_tfs] <- tf_color_map[present_tfs]



  
  # Plot
  p <- ggplot() +
    geom_segment(data = peak_df, aes(x = start_rel, xend = end_rel, y = y, yend = y),
                 color = "black", size = 5) +
    geom_segment(data = grid_df, aes(x = x, xend = x, y = y, yend = y),
                 color = "black", size = 0.5, linetype = "solid") +
    geom_text(data = grid_df, aes(x = x, y = y, label = label),
              size = 8, angle = 90, vjust = 0.5, hjust = 1.5, color = "black") +
    geom_segment(data = df, aes(x = midpoint, xend = midpoint, y = "Peak", yend = TF, color = TF),
                 size = 0.5, linetype = "dashed") +
    ggrepel::geom_text_repel(data = df, aes(x = midpoint, y = TF, label = TF, color = TF),
                             size = 8, angle = 90, hjust = -0.5, vjust = -0.5,
                             segment.size = 0, min.segment.length = 0,
                             box.padding = 0.5, point.padding = 0.5,
                             max.overlaps = Inf, direction = "both",
                             position = position_dodge(width = 10)) +
    scale_color_manual(values = tf_colors) +
    theme_minimal() +
    labs(title = paste("TF Motifs for Peak:", peak_str), subtitle = paste("DORC:", dorc_name)) +
    scale_x_continuous(limits = c(0, peak_end - peak_start)) +
    scale_y_discrete(limits = c("Peak", unique(df$TF)), expand = c(0, 0.1)) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt"))
  
  return(p)
}
getUpstreamDownstreamPeaks = function(the_gene, gtf_data, the_peaks, 
                                     direction = "upstream", 
                                     upstream_distance = 5000, 
                                     downstream_distance = 5000) {
  # Validate direction parameter
  if (!direction %in% c("upstream", "downstream", "both")) {
    stop("Direction must be 'upstream', 'downstream', or 'both'.")
  }
  
  # Filter for the gene, handling NAs
  gene_gtf <- gtf_data[!is.na(gtf_data$gene_name) & gtf_data$gene_name == the_gene & gtf_data$type == "gene"]
  
  # Check if gene is found
  if (length(gene_gtf) == 0) {
    stop("Gene ", the_gene, " not found in GTF file. Check gene name or GTF content.")
  }
  
  # Extract TSS based on strand
  tss <- ifelse(strand(gene_gtf) == "+", start(gene_gtf), end(gene_gtf))
  
  # Initialize result list
  result_peaks <- list(upstream = character(0), downstream = character(0))
  
  # Define upstream range if requested
  if (direction %in% c("upstream", "both")) {
    upstream_range <- GRanges(
      seqnames = seqnames(gene_gtf),
      ranges = IRanges(
        start = ifelse(strand(gene_gtf) == "+", tss - upstream_distance, tss),
        end = ifelse(strand(gene_gtf) == "+", tss, tss + upstream_distance)
      ),
      strand = strand(gene_gtf)
    )
  }
  
  # Define downstream range if requested
  if (direction %in% c("downstream", "both")) {
    downstream_range <- GRanges(
      seqnames = seqnames(gene_gtf),
      ranges = IRanges(
        start = ifelse(strand(gene_gtf) == "+", tss, tss - downstream_distance),
        end = ifelse(strand(gene_gtf) == "+", tss + downstream_distance, tss)
      ),
      strand = strand(gene_gtf)
    )
  }
  
  # Parse peaks into GRanges object
  peak_data <- strsplit(the_peaks, "-")
  peak_gr <- GRanges(
    seqnames = sapply(peak_data, `[`, 1),
    ranges = IRanges(
      start = as.numeric(sapply(peak_data, `[`, 2)),
      end = as.numeric(sapply(peak_data, `[`, 3))
    )
  )
  
  # Find overlapping peaks
  if (direction %in% c("upstream", "both")) {
    upstream_overlaps <- findOverlaps(peak_gr, upstream_range)
    result_peaks$upstream <- the_peaks[queryHits(upstream_overlaps)]
  }
  
  if (direction %in% c("downstream", "both")) {
    downstream_overlaps <- findOverlaps(peak_gr, downstream_range)
    result_peaks$downstream <- the_peaks[queryHits(downstream_overlaps)]
  }
  
  # Return results
  return(result_peaks)
}
# a custom heatmap function to plot TFs and DORCs from figR TFenrich data frame.
plotfigRHeatmap_C1 <- function(figR.d, score.cut = 1, score.pos = NULL, my_cols = NULL, DORCs = NULL, TFs = NULL, heatmap_name = "Score", rText = 10, cText = 10, ...) {
  if(is.null(myCols)) myCols <- circlize::colorRamp2(seq(-2, 2, length.out = 9), colors = BuenColors::jdb_palette("solar_flare"))
  
  message("Using absolute score cut-off of: ", score.cut, " ..\n")
  if(is.null(score.pos)){
    DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% 
      pull(DORC) %>% unique()
    TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% 
      pull(Motif) %>% unique()
  
  if (!is.null(DORCs)) {
    if (!all(DORCs %in% figR.d$DORC)) 
      stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    DORCsToKeep <- intersect(DORCsToKeep, DORCs)
    TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & 
                                          DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }
  
  if (!is.null(TFs)) {
    if (!all(TFs %in% figR.d$Motif)) 
      stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    TFsToKeep <- intersect(TFsToKeep, TFs)
    DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & 
                                            Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }
  }else{
    DORCsToKeep <- figR.d %>% dplyr::filter(Score >= score.pos) %>% 
      pull(DORC) %>% unique()
    TFsToKeep <- figR.d %>% dplyr::filter(Score >= score.pos) %>% 
      pull(Motif) %>% unique()
  
  if (!is.null(DORCs)) {
    if (!all(DORCs %in% figR.d$DORC)) 
      stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    DORCsToKeep <- intersect(DORCsToKeep, DORCs)
    TFsToKeep <- figR.d %>% dplyr::filter(Score >= score.pos & 
                                          DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }
  
  if (!is.null(TFs)) {
    if (!all(TFs %in% figR.d$Motif)) 
      stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    TFsToKeep <- intersect(TFsToKeep, TFs)
    DORCsToKeep <- figR.d %>% dplyr::filter(Score >= score.pos & 
                                            Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }
  }
  net.d <- figR.d %>% dplyr::filter(DORC %in% DORCsToKeep & 
                                    Motif %in% TFsToKeep) %>% 
    reshape2::dcast(DORC ~ Motif) %>% 
    tibble::column_to_rownames("DORC") %>% as.matrix()
  
  message("Plotting ", nrow(net.d), " DORCs x ", ncol(net.d), " TFs\n")

  myHeat <- ComplexHeatmap::Heatmap(net.d, 
                                   col = myCols, 
                                   clustering_distance_rows = "pearson", 
                                   clustering_distance_columns = "pearson", 
                                   name = heatmap_name,  # Use custom heatmap name
                                   border = TRUE, 
                                   row_names_gp = gpar(fontsize = rText, fontface = "italic"),
                                   column_names_gp = gpar(fontsize = rText, fontface = "italic"),
                                   ...)
  myHeat
}
plotfigRHeatmap_C3 <- function(figR.d, score.cut = 1, DORCs = NULL, TFs = NULL, myCols = NULL, heatmap_name = "Score", rText = 10, cText = 10, the_col_max = 1, ...) {
  if(is.null(myCols)) myCols <- circlize::colorRamp2(seq(-2, 2, length.out = 9), colors = BuenColors::jdb_palette("solar_flare"))
  
  message("Using absolute score cut-off of: ", score.cut, " ..\n")
  
  # Use provided DORCs and TFs, or fall back to all unique DORCs and TFs in the data
  DORCsToKeep <- if (!is.null(DORCs)) {
    if (!all(DORCs %in% figR.d$DORC)) {
      warning("Some DORCs specified are not in the data frame. Adding them with Score = 0.001.")
    }
    DORCs
  } else {
    unique(figR.d$DORC)
  }
  
  TFsToKeep <- if (!is.null(TFs)) {
    if (!all(TFs %in% figR.d$Motif)) {
      warning("Some TFs specified are not in the data frame. Adding them with Score = 0.001.")
    }
    TFs
  } else {
    unique(figR.d$Motif)
  }
  
  # Create a complete data frame with all DORC-TF pairs
  complete_pairs <- expand.grid(DORC = DORCsToKeep, Motif = TFsToKeep, stringsAsFactors = FALSE)
  
  # Merge with original data, filling missing pairs with Score = 0.001
  # Add a flag to track padded values
  figR.d_padded <- complete_pairs %>%
    left_join(figR.d, by = c("DORC", "Motif")) %>%
    mutate(Score = ifelse(is.na(Score), 0.001, Score),
           is_padded = is.na(Score))  # Flag padded values
  
  # Filter based on score cut-off, but always include padded values (Score = 0.001)
  net.d <- figR.d_padded %>%
    dplyr::filter(abs(Score) >= score.cut | is_padded) %>%
    reshape2::dcast(DORC ~ Motif, value.var = "Score") %>%
    tibble::column_to_rownames("DORC") %>%
    as.matrix()
  
  # Check if the resulting matrix is empty
  if (nrow(net.d) == 0 || ncol(net.d) == 0) {
    stop("No data remains after filtering. Check score.cut or input data.")
  }
  
  message("Plotting ", nrow(net.d), " DORCs x ", ncol(net.d), " TFs\n")
  
  # Define color scale
  myCols <- circlize::colorRamp2(seq(-the_col_max, the_col_max, length.out = 9), 
                                 colors = BuenColors::jdb_palette("solar_flare"))
  
  # Create heatmap
  myHeat <- ComplexHeatmap::Heatmap(net.d, 
                                   col = myCols, 
                                   clustering_distance_rows = "pearson", 
                                   clustering_distance_columns = "pearson", 
                                   name = heatmap_name, 
                                   border = TRUE, 
                                   row_names_gp = grid::gpar(fontsize = rText, fontface = "italic"),
                                   column_names_gp = grid::gpar(fontsize = cText, fontface = "italic"),
                                   ...)
  myHeat
}
# these functions will generate Gene Regulatory Networks
get_network_data_old1 = function (
    figR.d, 
    score.cut = 1, 
    TFs = NULL, 
    DORCs = NULL, 
    weight.edges = FALSE, 
    TFnodecol = "Tomato", 
    DORCnodecol = "Sky Blue", 
    posEdgecol = "Forest Green", 
    negEdgecol = "Purple", 
    labelSize = 13, 
    TFnodeSize = 5, 
    DORCnodeSize = 5, 
    edgeWidthScaleMin = 0.0,
    edgeWidthScaleMax = 1,
    edgeAlpha = 0.6,  # Transparency for edge clarity
    verbose = FALSE,   # Option for debugging output
    the_layout = "tree", # options; sugiyama, linear, tree, bipartite
    the_circular = TRUE) {
     # Load required packages
    require(dplyr)
    require(igraph)
    require(ggraph)
    require(ggplot2)
    require(scales)
    require(ggrepel)  # For better label repulsion
    
    # Validate inputs
    if (is.null(TFs)){ #stop("TFs list (my_TFs) must be provided.")
        TFs = figR.d$Motif      
    }   
    if (is.null(DORCs)){ #stop("DORCs list (my_dorcs) must be provided.")
        DORCs = figR.d$DORC
    }
    
    # Check for missing values in critical columns
    if (any(is.na(figR.d$Score)) || any(is.na(figR.d$Corr))) {
        stop("Missing values detected in Score or Corr columns.")
    }
    
    # Validate color inputs
    validate_color <- function(col, name) {
        tryCatch(
            gplots::col2hex(col),
            error = function(e) stop("Invalid color for ", name, ": ", col)
        )
    }
    validate_color(TFnodecol, "TFnodecol")
    validate_color(DORCnodecol, "DORCnodecol")
    validate_color(posEdgecol, "posEdgecol")
    validate_color(negEdgecol, "negEdgecol")
    
    # Dynamic score cutoff if score.cut is NULL
    if (is.null(score.cut)) {
        score.cut <- quantile(abs(figR.d$Score), 0.9, na.rm = TRUE)
        message("Using score.cut = ", round(score.cut, 3), " (90th percentile of abs(Score)).")
    }
    
    # Filter data by score threshold and provided TFs/DORCs
    net.dat <- figR.d %>% 
        dplyr::filter(
            abs(Score) >= score.cut,
            Motif %in% TFs,  # Restrict Motif to TFs only
            DORC %in% DORCs  # DORC must be in provided DORCs
        )
    
    # Warn if no interactions remain after filtering
    if (nrow(net.dat) == 0) {
        warning("No interactions meet score.cut = ", score.cut, ". Consider lowering the threshold.")
        stop("No valid interactions found after filtering. Check TFs, DORCs, or score.cut.")
    }
    
    # Format node and edge names (add "." for consistency with original)
    net.dat$Motif <- paste0(net.dat$Motif, ".")
    net.dat$DORC <- paste0(net.dat$DORC, ".")
    
    # Define nodes: distinguish TFs and DORCs
    all_nodes <- unique(c(net.dat$Motif, net.dat$DORC))
    tf_names <- paste0(TFs, ".")
    nodes <- data.frame(
        name = all_nodes,
        group = ifelse(all_nodes %in% tf_names, "TF", "DORC"),
        size = ifelse(all_nodes %in% tf_names, TFnodeSize, DORCnodeSize)
    )
    
    # Debugging output for nodes
    if (verbose) {
        #message("Node classifications:")
        #print(nodes)
    }
    
    # Define edges
    edges <- as.data.frame(net.dat)
    links <- data.frame(
        from = edges$Motif,
        to = edges$DORC,
        corr = edges$Corr,
        enrichment = edges$Enrichment.P,
        weight = if(weight.edges) scales::rescale(log1p(abs(edges$Score))) * edgeWidthScaleMax else 1
    )
    
    # Debugging output for edges
    if (verbose) {
        #message("Edge correlations:")
        #print(links[, c("from", "to", "corr")])
    }
    
    # Create igraph object
    g <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)
    return(list("g" = g, "edges" = edges, "links" = links, "nodes" = nodes))
}


create_network_plots_old1 <- function(
    the_list,           # List containing TFenrich.d.rds for each condition
    conditions,        # Named vector of conditions
    my_TFs = NULL,            # List of transcription factors
    my_dorcs = NULL,          # List of DORCs
    score.cut = 0.1,   # Score cutoff for filtering
    weight.edges = TRUE,
    verbose = TRUE,
    the_layout = "sugiyama",
    the_circular = FALSE,
    labelSize = 13,
    edgeWidthScaleMin = 0.0,
    edgeWidthScaleMax = 1,
    posEdgecol = "Forest Green",
    negEdgecol = "Red",
    missingEdgecol = "Gray",
    TFnodecol = "Tomato", 
    DORCnodecol = "Sky Blue", 
    edgeAlpha = 0.6,
    dashedEdgeWidth = 0.2,
    arrLength = 0.1,
    arrRation = 2,
    TFnodeSize = 5, 
    DORCnodeSize = 5
) {
    # Load required packages
    require(dplyr)
    require(igraph)
    require(ggraph)
    require(ggplot2)
    require(scales)
    require(ggrepel)
    require(patchwork)
    require(gplots)
    
    # Step 1: Generate combined graph data
    combined_figR.d <- do.call(rbind, lapply(conditions, function(x) the_list[[x]]$TFenrich.d.rds))
    
    # Validate combined data
    if (!is.data.frame(combined_figR.d)) stop("Combined figR.d is not a data frame.")
    if (!all(c("Motif", "DORC", "Score", "Corr", "Enrichment.P") %in% colnames(combined_figR.d))) {
        stop("Required columns missing in combined figR.d.")
    }
    
    combined_data <- get_network_data_old1(
        figR.d = combined_figR.d,
        score.cut = score.cut,
        TFs = my_TFs,
        DORCs = my_dorcs,
        weight.edges = weight.edges,
        verbose = verbose,
        the_layout = the_layout,
        the_circular = the_circular,
        TFnodeSize = TFnodeSize, 
        DORCnodeSize = DORCnodeSize
    )
    
    # Extract combined graph components
    combined_nodes <- combined_data$nodes
    combined_links <- combined_data$links
    combined_g <- combined_data$g
    
    # Set edge attributes for combined graph
    combined_links$status <- "present"
    combined_links$edge_linetype <- "solid"
    combined_links$edge_color_category <- factor(
        ifelse(combined_links$corr > 0, "Positive", "Negative"),
        levels = c("Positive", "Negative", "Missing")
    )
    
    # Step 2: Generate condition-specific data with missing edges
    condition_data_list <- lapply(names(conditions), function(cond) {
        figR.d <- the_list[[conditions[cond]]]$TFenrich.d.rds
        if (!is.data.frame(figR.d)) stop("figR.d is not a data frame for condition: ", cond)
        if (!all(c("Motif", "DORC", "Score", "Corr", "Enrichment.P") %in% colnames(figR.d))) {
            stop("Required columns missing in figR.d for condition: ", cond)
        }
        
        cond_data <- get_network_data_old1(
            figR.d = figR.d,
            score.cut = score.cut,
            TFs = my_TFs,
            DORCs = my_dorcs,
            weight.edges = weight.edges,
            verbose = verbose,
            the_layout = the_layout,
            the_circular = the_circular,
            TFnodeSize = TFnodeSize, 
            DORCnodeSize = DORCnodeSize
        )
        
        cond_links <- cond_data$links
        all_links <- combined_links %>%
            dplyr::mutate(
                status = ifelse(paste0(from, "-", to) %in% paste0(cond_links$from, "-", cond_links$to), 
                               "present", "missing"),
                edge_linetype = ifelse(status == "present", "solid", "dashed"),
                weight = ifelse(status == "present", 
                               cond_links$weight[match(paste0(from, "-", to), 
                                                      paste0(cond_links$from, "-", cond_links$to))], 
                               dashedEdgeWidth),
                corr = ifelse(status == "present", 
                             cond_links$corr[match(paste0(from, "-", to), 
                                                  paste0(cond_links$from, "-", cond_links$to))], 
                             0),
                enrichment = ifelse(status == "present", 
                                   cond_links$enrichment[match(paste0(from, "-", to), 
                                                              paste0(cond_links$from, "-", cond_links$to))], 
                                   NA),
                edge_color_category = factor(
                    case_when(
                        status == "missing" ~ "Missing",
                        status == "present" & corr > 0 ~ "Positive",
                        status == "present" & corr <= 0 ~ "Negative"
                    ),
                    levels = c("Positive", "Negative", "Missing")
                )
            )
        
        nodes <- combined_nodes
        g <- igraph::graph_from_data_frame(d = all_links, vertices = nodes, directed = TRUE)
        
        if (verbose) {
            #message("Edge linetypes for ", cond, ":")
            #print(table(all_links$edge_linetype))
            #message("Edge status for ", cond, ":")
            #print(table(all_links$status))
            #message("Edge color categories for ", cond, ":")
            #print(table(all_links$edge_color_category))
        }
        
        return(list(g = g, links = all_links, nodes = nodes))
    }) %>% setNames(names(conditions))
    
    # Check if any condition has edges
    if (all(sapply(condition_data_list, function(x) nrow(x$links) == 0))) {
        stop("No valid interactions found across all conditions. Check TFs, DORCs, or score.cut.")
    }
    
    # Step 3: Create condition-specific plots
    plots <- lapply(names(conditions), function(cond) {
        #message("Generating plot for condition: ", cond)
        
        links <- condition_data_list[[cond]]$links
        
        ggraph(combined_g, layout = the_layout, circular = the_circular) +
            geom_edge_diagonal(
                aes(
                    edge_width = links$weight,
                    edge_colour = links$edge_color_category,
                    edge_linetype = links$edge_linetype,
                    edge_alpha = edgeAlpha
                ),
                arrow = arrow(length = unit(arrLength, "inches"), type = "closed"),
                end_cap = circle(max(combined_nodes$size) / arrRation, "points")
            ) +
            geom_node_point(
                aes(size = size, color = group)
            ) +
            geom_text_repel(
                aes(x = x, y = y, label = name),
                size = labelSize / 3,
                min.segment.length = 0,
                box.padding = 0.5,
                max.overlaps = 20,
                nudge_y = 0.2,
                nudge_x = 0.1,
                hjust = 0.5,
                vjust = 0.5
            ) +
            scale_edge_width(range = c(edgeWidthScaleMin, edgeWidthScaleMax)) +
            scale_edge_colour_manual(
                values = c(
                    "Positive" = gplots::col2hex(posEdgecol),
                    "Negative" = gplots::col2hex(negEdgecol),
                    "Missing" = gplots::col2hex(missingEdgecol)
                ),
                labels = c("Positive" = "Positive", "Negative" = "Negative", "Missing" = "Missing"),
                breaks = c("Positive", "Negative", "Missing"),
                drop = FALSE
            ) +
            scale_edge_linetype_manual(
                values = c("solid" = "solid", "dashed" = "dashed"),
                labels = c("solid" = "Present", "dashed" = "Missing")
            ) +
            scale_color_manual(
                values = c("TF" = gplots::col2hex(TFnodecol), "DORC" = gplots::col2hex(DORCnodecol)),
                labels = c("DORC", "TF")
            ) +
            scale_edge_alpha_identity() +
            scale_size_identity() +
            theme_void() +
            theme(
                legend.position = "right",
                legend.title = element_text(size = 10),
                plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
                legend.text = element_text(size = 8)
            ) +
            guides(
                edge_colour = guide_legend(
                    title = "Edge Correlation/Status",
                    override.aes = list(
                        edge_colour = c(gplots::col2hex(posEdgecol), 
                                        gplots::col2hex(negEdgecol), 
                                        gplots::col2hex(missingEdgecol)),
                        linetype = c("solid", "solid", "dashed")
                    )
                ),
                color = guide_legend(title = "Node Type"),
                size = "none",
                edge_width = guide_legend(title = "Edge Weight"),
                edge_alpha = "none",
                edge_linetype = guide_legend(title = "Edge Status")
            ) +
            ggtitle(cond)
    }) %>% setNames(names(conditions))
    
    # Step 4: Plot combined graph
    p_combined <- ggraph(combined_g, layout = the_layout, circular = the_circular) +
        geom_edge_diagonal(
            aes(
                edge_width = combined_links$weight,
                edge_colour = combined_links$edge_color_category,
                edge_linetype = combined_links$edge_linetype,
                edge_alpha = edgeAlpha
            ),
            arrow = arrow(length = unit(arrLength, "inches"), type = "closed"),
            end_cap = circle(max(combined_nodes$size) / arrRation, "points")
        ) +
        geom_node_point(
            aes(size = size, color = group)
        ) +
        geom_text_repel(
            aes(x = x, y = y, label = name),
            size = labelSize / 3,
            min.segment.length = 0,
            box.padding = 0.5,
            max.overlaps = 20,
            nudge_y = 0.2,
            nudge_x = 0.1,
            hjust = 0.5,
            vjust = 0.5
        ) +
        scale_edge_width(range = c(edgeWidthScaleMin, edgeWidthScaleMax)) +
        scale_edge_colour_manual(
            values = c(
                "Positive" = gplots::col2hex(posEdgecol),
                "Negative" = gplots::col2hex(negEdgecol),
                "Missing" = gplots::col2hex(missingEdgecol)
            ),
            labels = c("Positive" = "Positive", "Negative" = "Negative", "Missing" = "Missing"),
            breaks = c("Positive", "Negative", "Missing"),
            drop = FALSE
        ) +
        scale_edge_linetype_manual(
            values = c("solid" = "solid", "dashed" = "dashed"),
            labels = c("solid" = "Present", "dashed" = "Missing")
        ) +
        scale_color_manual(
                values = c("TF" = gplots::col2hex(TFnodecol), "DORC" = gplots::col2hex(DORCnodecol)),
                labels = c("DORC", "TF")
            ) +
        scale_edge_alpha_identity() +
        scale_size_identity() +
        theme_void() +
        theme(
            legend.position = "right",
            legend.title = element_text(size = 10),
            plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
            legend.text = element_text(size = 8)
        ) +
        guides(
            edge_colour = guide_legend(
                title = "Corr",
                override.aes = list(
                    edge_colour = c(gplots::col2hex(posEdgecol), 
                                    gplots::col2hex(negEdgecol), 
                                    gplots::col2hex(missingEdgecol)),
                    linetype = c("solid", "solid", "dashed")
                )
            ),
            color = guide_legend(title = "Node Type"),
            size = "none",
            edge_width = guide_legend(title = "Edge Weight"),
            edge_alpha = "none",
            edge_linetype = guide_legend(title = "Edge Status")
        ) +
        ggtitle("Combined")
    
    # Return list of plots
    return(c(plots, list(Combined = p_combined)))
}
get_network_data = function (
    figR.d, 
    score.cut = 1, 
    TFs = NULL, 
    DORCs = NULL, 
    score.min = 1,
    weight.edges = FALSE, 
    TFnodecol = "Tomato", 
    DORCnodecol = "Sky Blue", 
    posEdgecol = "Forest Green", 
    negEdgecol = "Purple", 
    labelSize = 13, 
    TFnodeSize = 5, 
    DORCnodeSize = 5, 
    edgeWidthScaleMin = 0.0,
    edgeWidthScaleMax = 1,
    edgeAlpha = 0.6,  # Transparency for edge clarity
    verbose = FALSE,   # Option for debugging output
    the_layout = "tree", # options; sugiyama, linear, tree, bipartite
    the_circular = TRUE,
    all_trx = NULL,
    keepDorc = TRUE
    ) {
  
    # Validate inputs
    if (is.null(TFs) & !keepDorc){ #stop("TFs list (my_TFs) must be provided.")
        TFs = figR.d$Motif      
    }   
    if (is.null(DORCs)){ #stop("DORCs list (my_dorcs) must be provided.")
        DORCs = figR.d$DORC
    }
    
    # Check for missing values in critical columns
    if (any(is.na(figR.d$Score)) || any(is.na(figR.d$Corr))) {
        stop("Missing values detected in Score or Corr columns.")
    }
    
    # Validate color inputs
    validate_color <- function(col, name) {
        tryCatch(
            gplots::col2hex(col),
            error = function(e) stop("Invalid color for ", name, ": ", col)
        )
    }
    validate_color(TFnodecol, "TFnodecol")
    validate_color(DORCnodecol, "DORCnodecol")
    validate_color(posEdgecol, "posEdgecol")
    validate_color(negEdgecol, "negEdgecol")
    
    # Dynamic score cutoff if score.cut is NULL
    if (is.null(score.cut)) {
        score.cut <- quantile(abs(figR.d$Score), 0.9, na.rm = TRUE)
        message("Using score.cut = ", round(score.cut, 3), " (90th percentile of abs(Score)).")
    }
    
    # Filter data by score threshold and provided TFs/DORCs
    if(keepDorc){
        if(is.null(TFs)){
            tmp.dat <- figR.d %>% 
                dplyr::filter(
                abs(Score) >= score.cut,
                DORC %in% DORCs  # DORC must be in provided DORCs
            )
            TFs = tmp.dat$Motif
        }else{
            tmp.dat <- figR.d %>% 
            dplyr::filter(
                abs(Score) >= score.cut,
                Motif %in% TFs,
                DORC %in% DORCs  # DORC must be in provided DORCs
            )
            TFs = c(TFs,tmp.dat$Motif)
        }
        DORCs = c(DORCs,unique(TFs))
        net.dat <- figR.d[figR.d$Motif %in% DORCs & figR.d$DORC %in% DORCs & abs(figR.d$Score) >= score.min,]
    }else{
        # Filter data by score threshold and provided TFs/DORCs
        net.dat <- figR.d %>% 
            dplyr::filter(
                abs(Score) >= score.cut,
                Motif %in% TFs,  # Restrict Motif to TFs only
                DORC %in% DORCs  # DORC must be in provided DORCs
            )
    }
    # Warn if no interactions remain after filtering
    if (nrow(net.dat) == 0) {
        warning("No interactions meet score.cut = ", score.cut, ". Consider lowering the threshold.")
        stop("No valid interactions found after filtering. Check TFs, DORCs, or score.cut.")
    }
    
    # Format node and edge names (add "." for consistency with original)
    net.dat$Motif <- paste0(net.dat$Motif, ".")
    net.dat$DORC <- paste0(net.dat$DORC, ".")
    
    # Define nodes: distinguish TFs and DORCs
    all_nodes <- unique(c(net.dat$Motif, net.dat$DORC))
    tf_names <- paste0(TFs, ".")
    nodes <- data.frame(
        name = all_nodes,
        group = ifelse(all_nodes %in% tf_names, "TF", "DORC"),
        size = ifelse(all_nodes %in% tf_names, TFnodeSize, DORCnodeSize)
    )
    
    # Debugging output for nodes
    if (verbose) {
        #message("Node classifications:")
        #print(nodes)
    }
    
    # Define edges
    edges <- as.data.frame(net.dat)
    links <- data.frame(
        from = edges$Motif,
        to = edges$DORC,
        corr = edges$Corr,
        enrichment = edges$Enrichment.P,
        weight = if(weight.edges) scales::rescale(log1p(abs(edges$Score))) * edgeWidthScaleMax else 1
    )
    
    # Debugging output for edges
    if (verbose) {
        #message("Edge correlations:")
        #print(links[, c("from", "to", "corr")])
    }
    
    # Create igraph object
    g <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)
    return(list("g" = g, "edges" = edges, "links" = links, "nodes" = nodes))
}
create_network_plots <- function(
    the_list,                 # List containing TFenrich.d.rds for each condition
    conditions = NULL,        # Named vector of conditions
    my_TFs = NULL,            # List of transcription factors
    my_dorcs = NULL,          # List of DORCs
    score.cut = 0.1,          # Score cutoff for filtering
    score.min = 0.1,
    weight.edges = TRUE,
    verbose = TRUE,
    the_layout = "sugiyama", #eigen, tree, kk, linear, partition, hive, circlepack, treemap, cactustree, unrooted
    the_circular = FALSE,
    labelSize = 8,
    edgeWidthScaleMin = 0.0,
    edgeWidthScaleMax = 1,
    TFnodecol = "purple", 
    DORCnodecol = "orange", 
    posEdgecol = "cyan3",
    negEdgecol = "red",
    missingEdgecol = "gray80",
    showDORClabels = TRUE,    # Parameter to control DORC labels
    edgeAlpha = 0.6,
    dashedEdgeWidth = 0.2,
    arrLength = 0.12,
    arrRation = 0.7,
    TFnodeSize = 4, 
    DORCnodeSize = 4,
    the_TFcolumn = "TFenrich.d.rds",
    all_trx = NULL,
    keepDorc = TRUE
) {
    
    if(is.null(conditions)) conditions = names(the_list)

    # Step 1: Generate combined graph data
    combined_figR.d <- do.call(rbind, lapply(conditions, function(x) the_list[[x]][[the_TFcolumn]]))
    if(is.null(all_trx)) all_trx = unique(combined_figR.d$Motif)
    
    # Validate combined data
    if (!is.data.frame(combined_figR.d)) stop("Combined figR.d is not a data frame.")
    if (!all(c("Motif", "DORC", "Score", "Corr", "Enrichment.P") %in% colnames(combined_figR.d))) {
        stop("Required columns missing in combined figR.d.")
    }
    
    combined_data <- get_network_data(
        figR.d = combined_figR.d,
        score.cut = score.cut,
        TFs = my_TFs,
        DORCs = my_dorcs,
        score.min = score.min,
        weight.edges = weight.edges,
        verbose = verbose,
        the_layout = the_layout,
        the_circular = the_circular,
        TFnodeSize = TFnodeSize, 
        DORCnodeSize = DORCnodeSize,
        all_trx = all_trx,
        keepDorc = keepDorc
    )
    
    # Extract combined graph components
    combined_nodes <- combined_data$nodes
    combined_links <- combined_data$links
    combined_g <- combined_data$g
    
    # Add label column to combined_nodes based on showDORClabels
    combined_nodes$label <- ifelse(showDORClabels | combined_nodes$group == "TF", 
                                  combined_nodes$name, 
                                  "")
    
    # Set edge attributes for combined graph
    combined_links$status <- "present"
    combined_links$edge_linetype <- "solid"
    combined_links$edge_color_category <- factor(
        ifelse(combined_links$corr > 0, "Positive", "Negative"),
        levels = c("Positive", "Negative", "Missing")
    )
    
    # Step 2: Generate condition-specific data with missing edges
    condition_data_list <- lapply(names(conditions), function(cond) {
        figR.d <- the_list[[conditions[cond]]][[the_TFcolumn]]
        if (!is.data.frame(figR.d)) stop("figR.d is not a data frame for condition: ", cond)
        if (!all(c("Motif", "DORC", "Score", "Corr", "Enrichment.P") %in% colnames(figR.d))) {
            stop("Required columns missing in figR.d for condition: ", cond)
        }
        
        cond_data <- get_network_data(
            figR.d = figR.d,
            score.cut = score.cut,
            TFs = my_TFs,
            DORCs = my_dorcs,
            score.min = score.min,
            weight.edges = weight.edges,
            verbose = verbose,
            the_layout = the_layout,
            the_circular = the_circular,
            TFnodeSize = TFnodeSize, 
            DORCnodeSize = DORCnodeSize,
        all_trx = all_trx,
        keepDorc = keepDorc
        )
        
        cond_links <- cond_data$links
        all_links <- combined_links %>%
            dplyr::mutate(
                status = ifelse(paste0(from, "-", to) %in% paste0(cond_links$from, "-", cond_links$to), 
                               "present", "missing"),
                edge_linetype = ifelse(status == "present", "solid", "dashed"),
                weight = ifelse(status == "present", 
                               cond_links$weight[match(paste0(from, "-", to), 
                                                      paste0(cond_links$from, "-", cond_links$to))], 
                               dashedEdgeWidth),
                corr = ifelse(status == "present", 
                             cond_links$corr[match(paste0(from, "-", to), 
                                                  paste0(cond_links$from, "-", cond_links$to))], 
                             0),
                enrichment = ifelse(status == "present", 
                                   cond_links$enrichment[match(paste0(from, "-", to), 
                                                              paste0(cond_links$from, "-", cond_links$to))], 
                                   NA),
                edge_color_category = factor(
                    case_when(
                        status == "missing" ~ "Missing",
                        status == "present" & corr > 0 ~ "Positive",
                        status == "present" & corr <= 0 ~ "Negative"
                    ),
                    levels = c("Positive", "Negative", "Missing")
                )
            )
        
        nodes <- combined_nodes
        g <- igraph::graph_from_data_frame(d = all_links, vertices = nodes, directed = TRUE)
        
        if (verbose) {
            #message("Edge linetypes for ", cond, ":")
            #print(table(all_links$edge_linetype))
            #message("Edge status for ", cond, ":")
            #print(table(all_links$status))
            #message("Edge color categories for ", cond, ":")
            #print(table(all_links$edge_color_category))
        }
        
        return(list(g = g, links = all_links, nodes = nodes))
    }) %>% setNames(names(conditions))
    
    # Check if any condition has edges
    if (all(sapply(condition_data_list, function(x) nrow(x$links) == 0))) {
        stop("No valid interactions found across all conditions. Check TFs, DORCs, or score.cut.")
    }
    
    # Step 3: Create condition-specific plots
    plots <- lapply(names(conditions), function(cond) {
        #message("Generating plot for condition: ", cond)
        
        links <- condition_data_list[[cond]]$links
        
        # Create the ggraph object
        g <- ggraph(combined_g, layout = the_layout, circular = the_circular)
        
        # Get node data with layout coordinates from ggraph
        node_data <- g$data
        # Merge with label column from combined_nodes
        node_data <- node_data %>% 
            dplyr::left_join(dplyr::select(combined_nodes, name, label), by = "name")
        
        g +
            geom_edge_diagonal(
                aes(
                    edge_width = links$weight,
                    edge_colour = links$edge_color_category,
                    edge_linetype = links$edge_linetype,
                    edge_alpha = edgeAlpha
                ),
                arrow = arrow(length = unit(arrLength, "inches"), type = "closed"),
                end_cap = circle(max(combined_nodes$size) / arrRation, "points")
            ) +
            geom_node_point(
                aes(size = size, color = group)
            ) +
            geom_text_repel(
                data = node_data,  # Use node_data with layout coordinates and labels
                aes(x = x, y = y, label = label),
                size = labelSize / 3,
                min.segment.length = 0,
                box.padding = 0.5,
                max.overlaps = 20,
                nudge_y = 0.2,
                nudge_x = 0.1,
                hjust = 0.5,
                vjust = 0.5
            ) +
            scale_edge_width(range = c(edgeWidthScaleMin, edgeWidthScaleMax)) +
            scale_edge_colour_manual(
                values = c(
                    "Positive" = gplots::col2hex(posEdgecol),
                    "Negative" = gplots::col2hex(negEdgecol),
                    "Missing" = gplots::col2hex(missingEdgecol)
                ),
                labels = c("Positive" = "Positive", "Negative" = "Negative", "Missing" = "Missing"),
                breaks = c("Positive", "Negative", "Missing"),
                drop = FALSE
            ) +
            scale_edge_linetype_manual(
                values = c("solid" = "solid", "dashed" = "dashed"),
                labels = c("solid" = "Present", "dashed" = "Missing")
            ) +
            scale_color_manual(
                values = c("TF" = gplots::col2hex(TFnodecol), "DORC" = gplots::col2hex(DORCnodecol)),
                labels = c("DORC", "TF")
            ) +
            scale_edge_alpha_identity() +
            scale_size_identity() +
            theme_void() +
            theme(
                legend.position = "right",
                legend.title = element_text(size = 10),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                legend.text = element_text(size = 8)
            ) +
            guides(
                edge_colour = guide_legend(
                    title = "Edge Correlation/Status",
                    override.aes = list(
                        edge_colour = c(gplots::col2hex(posEdgecol), 
                                        gplots::col2hex(negEdgecol), 
                                        gplots::col2hex(missingEdgecol)),
                        linetype = c("solid", "solid", "dashed")
                    )
                ),
                color = guide_legend(title = "Node Type"),
                size = "none",
                edge_width = guide_legend(title = "Edge Weight"),
                edge_alpha = "none",
                edge_linetype = guide_legend(title = "Edge Status")
            ) +
            ggtitle(cond)
    }) %>% setNames(names(conditions))
    
    # Step 4: Plot combined graph
    # Create the ggraph object
    p_combined <- ggraph(combined_g, layout = the_layout, circular = the_circular)
    
    # Get node data with layout coordinates from ggraph
    node_data <- p_combined$data
    # Merge with label column from combined_nodes
    node_data <- node_data %>% 
        dplyr::left_join(dplyr::select(combined_nodes, name, label), by = "name")
    
    p_combined <- p_combined +
        geom_edge_diagonal(
            aes(
                edge_width = combined_links$weight,
                edge_colour = combined_links$edge_color_category,
                edge_linetype = combined_links$edge_linetype,
                edge_alpha = edgeAlpha
            ),
            arrow = arrow(length = unit(arrLength, "inches"), type = "closed"),
            end_cap = circle(max(combined_nodes$size) / arrRation, "points")
        ) +
        geom_node_point(
            aes(size = size, color = group)
        ) +
        geom_text_repel(
            data = node_data,  # Use node_data with layout coordinates and labels
            aes(x = x, y = y, label = label),
            size = labelSize / 3,
            min.segment.length = 0,
            box.padding = 0.5,
            max.overlaps = 20,
            nudge_y = 0.2,
            nudge_x = 0.1,
            hjust = 0.5,
            vjust = 0.5
        ) +
        scale_edge_width(range = c(edgeWidthScaleMin, edgeWidthScaleMax)) +
        scale_edge_colour_manual(
            values = c(
                "Positive" = gplots::col2hex(posEdgecol),
                "Negative" = gplots::col2hex(negEdgecol),
                "Missing" = gplots::col2hex(missingEdgecol)
            ),
            labels = c("Positive" = "Positive", "Negative" = "Negative", "Missing" = "Missing"),
            breaks = c("Positive", "Negative", "Missing"),
            drop = FALSE
        ) +
        scale_edge_linetype_manual(
            values = c("solid" = "solid", "dashed" = "dashed"),
            labels = c("solid" = "Present", "dashed" = "Missing")
        ) +
        scale_color_manual(
            values = c("TF" = gplots::col2hex(TFnodecol), "DORC" = gplots::col2hex(DORCnodecol)),
            labels = c("DORC", "TF")
        ) +
        scale_edge_alpha_identity() +
        scale_size_identity() +
        theme_void() +
        theme(
            legend.position = "right",
            legend.title = element_text(size = 10),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.text = element_text(size = 8)
        ) +
        guides(
            edge_colour = guide_legend(
                title = "Corr",
                override.aes = list(
                    edge_colour = c(gplots::col2hex(posEdgecol), 
                                    gplots::col2hex(negEdgecol), 
                                    gplots::col2hex(missingEdgecol)),
                    linetype = c("solid", "solid", "dashed")
                )
            ),
            color = guide_legend(title = "Node Type"),
            size = "none",
            edge_width = guide_legend(title = "Edge Weight"),
            edge_alpha = "none",
            edge_linetype = guide_legend(title = "Edge Status")
        ) +
        ggtitle("Combined")
    
    # Return list of plots
    return(c(plots, list(Combined = p_combined)))
}

####### MARKER GENES

the_neuron_markers = list("V1b" = c("en1b"),
    "V2b/s" = c( "gata3", "irx1a","nkx1.2lb"), #"tal1",
    "V2a" = c("vsx2","shox2"),
    "V0" = c("evx1", "evx2"),
    "V1a" = c("sp9"), #, "foxd3"),
    "dl2" = c("slc17a6b"), #"foxd3", 
    "DoLa" = c("pnoca"),
    "CoLo" = c("chga"),
    "dl1" = c("barhl2"),
    "KA" = c("pkd1l2a", "tal1","tal2"),
    "OPC" = c("aplnra", "aplnrb"),
    "OD"  = c("cldnk", "mpz"),
    "FB" = c("pdgfrb", "col12a1a"),
    "ENS" = c("phox2bb"),
    "Schw" = c("slc35c2","cldnh","epcam"),
    "MN-sec" =c("nr2f1a"),
    "ISN" = c("tph2","ddc"),
    #"MN-pri" = c("chga"),
    "MN-V3-SRN" = c("slc18a3a"),
    "V3" = c("sim1a"),
    "dl3" = c("otpa", "otpb", "isl1"),
    "dl6" = c("wt1a"),"immN2"= c("nhlh2"),"immN1"=c("pcna","fabp7a")
    )
the_main_markers =list(
    "NT" = c("chata","gad1b","gad2","slc6a5","slc17a6a","slc17a6b"),
    "RP" = c("lmx1a","lmx1bb","lmx1ba"),
    "FP" = c("shha","shhb","foxa2","arxa"),
    "LFP" = c("nkx2.2a","nkx2.2b","wnt4b","nkx2.9"),
    "ERG" = c("gfap","fabp7a","sox2"),
    "p2" = c("vsx1","foxn4"),
    "immN" = c("nhlh2","nfia", "nfic","neurod1","neurog1"),
    "dl1" = c("barhl2", "barhl1a", "lhx2b", "lhx9"),
    "dl2" = c("foxd3"),
    "dl3" = c("tlx1","tlx2","tlx3b"), # https://neuraldevelopment.biomedcentral.com/articles/10.1186/s13064-016-0059-9; 
    "dl4" = c("otpa", "otpb"),
    "dl5" = c("lbx1a","lbx1b"),
    "dl6" = c("wt1a","dmrt3a", "bhlhe22"),
    "V0c" = c("evx1","evx2","pitx2"), # https://neuraldevelopment.biomedcentral.com/articles/10.1186/s13064-023-00176-w
    "V0v" = c("en1b"), # https://neuraldevelopment.biomedcentral.com/articles/10.1186/s13064-023-00176-w
    "V1" = c("sp9"),
    "V2a" = c("vsx2","shox2"),
    "V2b" = c("gata2a","gata3","tal1"),
    "V2s" = c("sox1a","sox1b","foxp2"),
    "MN" = c("mnx1","isl1"),
    "DoLa" = c("pnoca"),
    "V3" = c("sim1a"),
    "ISN" = c("tph2","fev","ddc"),
    "CSF-cNs" = c("sst1.1","urp1","pkd1l2a"),
    "OPC" = c("aplnra","aplnrb"),
    "OD" = c("mbpb","olig2","mpz"),
    "Schw" = c("slc35c2","erbb3b"), # https://www.cell.com/current-biology/fulltext/S0960-9822(05)00171-5
    "NC" = c("cmn", "tgfb3"),
    "RB" = c("isl2b","hpca","p2rx3a"),
    "ENS" = c("phox2bb","phox2a"),
    "FB" = c("col12a1a", "pdgfrb", "pfn1"),
    "MC" = c("mylz3","meox1"),
    "rps" = c("rpl32","rps27a"),
    "IC" = c("lcp1","mpx"),
    "RBC" = c("hbbe1.3","hbbe2")
    )


cell_cycle_genes = list(
g1_s_dr = c("mcm5","pcna","tyms","fen1","mcm2","mcm4","rrm1","unga","gins2","mcm6","cdca7a","dtl","prim1","uhrf1",
  "si:dkey-185e18.7","hells","rfc2","rpa2","nasp","rad51ap1","gmnn","wdr76","slbp","ccne2","ubr7","pold3",
  "msh2","atad2","rad51","rrm2","cdc45","cdc6","exo1","tipin","dscc1","blm","casp8ap2","usp1","pola1","chaf1b",
  "brip1","e2f8"),
g2_m_dr = c("hmgb2a","cdk1","nusap1","ube2c","birc5a","tpx2","top2a","ndc80","cks2","nuf2","cks1b","mki67","tmpoa","cenpf",
  "tacc3","smc4","ccnb2","ckap2l","aurkb","bub1","kif11","anp32e","tubb4b","gtse1","kif20ba","si:ch211-69g19.2",
  "jpt1a","cdc20","ttk","kif2c","rangap1a","ncapd2","dlgap5","si:ch211-244o22.2","cdca8","ect2","kif23","hmmr",
  "aurka","anln","lbr","ckap5","cenpe","ctcf","nek2","g2e3","gas2l3","cbx5","selenoh")
)

####### PLOTTING FOR DEG ANALYSES

get_fold_RectPlot <- function(
        the_df, 
        x_column = "Percentage",
        y_column = "Clusters", 
        f_column = "Samples", 
        the_colors = NULL, 
        x_limits = NULL,
        x_breaks = NULL, 
        fold_threshold = 1.25,
        the_angle = 0,
        vjs = NULL,
        hjs = NULL,
        the_line1_shape = "dashed",
        the_line1_color = "black",
        the_line1_size = 1,
        the_line2_shape = "dashed",
        the_line2_color = "red",
        the_line2_size = 1,
        showLegend = TRUE,           # Keep this for backward compatibility
        legend_position = "top",     # New parameter for legend position
        x_title = NULL,
        y_title = NULL,
        legend_title_size = 16,
        x_title_size = 16,
        y_title_size = 16,
        x_text_size = 16,
        y_text_size = 16,
        legend_text_size = 16,
        x_breaks_by = 1,
        y_breaks_by = 1
        ){
    
    # example
    #p0_8 = get_fold_RectPlot(the_counts, x_column = "Fold", y_column = "Clusters", f_column = "Folds", 
    #        the_colors = the_main_colors, the_angle = 90, vjs = 0.5, hjs = 1, x_title = "Fold Change",
    #        showLegend = T,legend_position = "right", the_line1_shape = "solid", the_line1_size = 0.5,
    #        the_line2_shape = "dashed", the_line2_size = 0.5,
    #        x_title_size = 24, x_text_size = 22,y_title_size = 24, y_text_size = 22, x_breaks_by = 1)

    if(is.null(x_title)) x_title = x_column
    if(is.null(y_title)) y_title = y_column
    
    # rename the columns for internal use
    colnames(the_df)[colnames(the_df) == x_column] = "x_column"
    colnames(the_df)[colnames(the_df) == y_column] = "y_column"
    colnames(the_df)[colnames(the_df) == f_column] = "f_column"
    
    #set limits and breaks
    if(is.null(x_limits)){
        x_min = abs(min(the_df$x_column, na.rm = TRUE))+0.1
        x_max = max(the_df$x_column, na.rm = TRUE)+0.1
        x_limits = c(-x_min, x_max)
    }
                           
    if(is.null(x_breaks)){
        xx_left = abs(x_limits[1]-0.1)
        yy_right = abs(x_limits[2]+0.1)
        dff = xx_left-floor(xx_left)
        if(dff > 0 & dff < 0.5){
            dff = 0.5
        }else{
            dff = 1
        }
        xx_left = floor(xx_left)+dff
        #
        dff = yy_right-floor(yy_right)
        if(dff > 0 & dff < 0.5){
            dff = 0.5
        }else{
            dff = 1
        }
        yy_right = floor(yy_right)+dff
        
        x_breaks = c(-xx_left, seq(ceiling(x_limits[1]), -1, by = y_breaks_by), seq(1, floor(x_limits[2]), by = x_breaks_by),yy_right)
        print(x_breaks)
    }

    the_df$Start <- ifelse(the_df$x_column >= 1, 1, -1)  # Bar starts at 1 or -1
    the_df$End <- the_df$x_column  # Bar ends at the actual value
    if(sum(the_df$x_column == Inf) != 0){
        idx = which(the_df$x_column == Inf)
        print(idx)
        the_df$x_column[idx] = 1
        the_df$Start[idx] = 1
        the_df$End[idx] = 1
        the_df$Deviation[idx] = 0
    }

    if(is.null(the_colors)){
        if(class(the_df$f_column) == "factor"){
            the_levels = levels(the_df$f_column)
        }else{
            the_levels = unique(the_df$f_column)
        }
        the_levels = unique(the_df$f_column)
        print(the_levels)
        the_colors = scPalette(8*length(the_levels))
        print(the_colors)
        the_colors = the_colors[seq(1, length(the_colors),8)]
        the_colors = setNames(the_colors, the_levels)
    }
    
    if(is.null(vjs) | is.null(hjs)){
            if(the_angle != 0){
                vjs = 1
                hjs = 0.5
            }else{
                vjs = 1
                hjs = 1
            }
    }
    
    # Create the plot with geom_rect
    pp <- ggplot(the_df, aes(y = y_column, fill = f_column)) +
        geom_rect(aes(xmin = Start, xmax = End, 
                  ymin = as.numeric(factor(y_column)) - 0.44, 
                  ymax = as.numeric(factor(y_column)) + 0.44), 
              position = position_dodge(width = 0.9)) +
        labs(x = x_title, y = y_title) +
        theme_classic() +
        scale_fill_manual(values = the_colors) +
        scale_x_continuous(breaks = x_breaks, limits = x_limits) +
        theme(
            axis.ticks = element_line(colour = "black"), 
            legend.title = element_text(face = "bold", colour = "black", size = legend_title_size),
            axis.text.x = element_text(angle = the_angle, vjust = vjs, hjust = hjs, face = "bold", size = x_text_size, color = "black"),
            axis.text.y = element_text(face = "bold", colour = "black", size = y_text_size),
            axis.title.x = element_text(face = "bold", size = x_title_size),
            axis.title.y = element_text(face = "bold", size = y_title_size),
            panel.background = element_rect(fill = "white", colour = "white"), 
            plot.background = element_rect(fill = "white", colour = "white"), 
            legend.text = element_text(face = "bold", size = legend_text_size)
        ) +
        geom_vline(xintercept = c(-1, 1), linetype = the_line1_shape, color = the_line1_color, size = the_line1_size) +
        geom_vline(xintercept = c(-fold_threshold, fold_threshold), linetype = the_line2_shape, color = the_line2_color, size = the_line2_size) +
        guides(fill = guide_legend(title = NULL))

    # Handle legend positioning and visibility
    if(!showLegend || legend_position == "none") {
        pp <- pp + theme(legend.position = "none")
    } else {
        # Validate legend_position
        valid_positions <- c("top", "right", "bottom", "left")
        if(!legend_position %in% valid_positions) {
            warning("Invalid legend_position. Using 'top' instead.")
            legend_position <- "top"
        }
        
        pp <- pp + theme(
            legend.position = legend_position,
            legend.justification = "center"
        )
    }
    
    print(x_limits)
    return(pp)
}
get_fold_RectPlot2 <- function(
  the_df, 
  x_column = "Percentage",
  y_column = "Clusters", 
  f_column = "Samples",  # Used for patterns
  the_colors = NULL,     # Now applies to y_column (Clusters)
  the_patterns = NULL,   # Applies to f_column
  x_limits = NULL,
  x_breaks = NULL, 
  fold_threshold = 1.25,
  the_angle = 0,
  vjs = NULL,
  hjs = NULL,
  the_line1_shape = "dashed",
  the_line1_color = "black",
  the_line1_size = 1,
  the_line2_shape = "dashed",
  the_line2_color = "red",
  the_line2_size = 1,
  showLegend = TRUE,
  legend_position = "top",
  x_title = NULL,
  y_title = NULL,
  legend_title_size = 16,
  x_title_size = 16,
  y_title_size = 16,
  x_text_size = 16,
  y_text_size = 16,
  legend_text_size = 16,
  x_breaks_by = 1,
  y_breaks_by = 1,
  bar_border_color = "black",  # New parameter for bar border color
  bar_border_size = 0.5        # New parameter for bar border thickness
) {

    # example usage
    #p0_8 <- get_fold_RectPlot2(
    #the_df = the_counts,
    #x_column = "Fold",
    #y_column = "Clusters",
    #f_column = "Folds",
    #the_colors = the_main_colors,  # Colors for Clusters
    #the_patterns = c("gControl_Les" = "none", "gSall1b_Les" = "circle", "gSall1b_gControl" = "crosshatch"),
    #the_angle = 90,
    #vjs = 0.5,
    #hjs = 1,
    #x_title = "Fold Change",
    #showLegend = TRUE,
    #legend_position = "right",
    #the_line1_shape = "solid",
    #the_line1_size = 0.5,
    #the_line2_shape = "dashed",
    #the_line2_size = 0.5,
    #x_title_size = 24,
    #x_text_size = 22,
    #y_title_size = 24,
    #y_text_size = 22,
    #x_breaks_by = 1
    #)
  
  if(is.null(x_title)) x_title = x_column
  if(is.null(y_title)) y_title = y_column
  
  # Rename columns for internal use
  colnames(the_df)[colnames(the_df) == x_column] = "x_column"
  colnames(the_df)[colnames(the_df) == y_column] = "y_column"
  
  # Handle multiple f_column inputs for patterns
  if(length(f_column) > 1) {
    the_df$f_column <- interaction(the_df[, f_column], sep = "_")
  } else {
    colnames(the_df)[colnames(the_df) == f_column] = "f_column"
  }
  
  # Set limits and breaks
  if(is.null(x_limits)){
    x_min = abs(min(the_df$x_column, na.rm = TRUE)) + 0.1
    x_max = max(the_df$x_column, na.rm = TRUE) + 0.1
    x_limits = c(-x_min, x_max)
  }
  
  if(is.null(x_breaks)){
    xx_left = abs(x_limits[1] - 0.1)
    yy_right = abs(x_limits[2] + 0.1)
    dff = xx_left - floor(xx_left)
    if(dff > 0 & dff < 0.5){
      dff = 0.5
    } else {
      dff = 1
    }
    xx_left = floor(xx_left) + dff
    dff = yy_right - floor(yy_right)
    if(dff > 0 & dff < 0.5){
      dff = 0.5
    } else {
      dff = 1
    }
    yy_right = floor(yy_right) + dff
    x_breaks = c(-xx_left, seq(ceiling(x_limits[1]), -1, by = y_breaks_by), 
                 seq(1, floor(x_limits[2]), by = x_breaks_by), yy_right)
    print(x_breaks)
  }
  
  # Define bar start and end points
  the_df$Start <- ifelse(the_df$x_column >= 1, 1, -1)
  the_df$End <- the_df$x_column
  if(sum(the_df$x_column == Inf) != 0){
    idx = which(the_df$x_column == Inf)
    print(idx)
    the_df$x_column[idx] = 1
    the_df$Start[idx] = 1
    the_df$End[idx] = 1
  }
  
  # Set default colors for y_column (Clusters) if not provided
  if(is.null(the_colors)){
    the_levels = unique(the_df$y_column)
    print(the_levels)
    the_colors = scPalette(8 * length(the_levels))  # Assumes scPalette exists
    the_colors = the_colors[seq(1, length(the_colors), 8)]
    the_colors = setNames(the_colors, the_levels)
  }
  
  # Set default patterns for f_column if not provided
  if(is.null(the_patterns)){
    the_levels = unique(the_df$f_column)
    the_patterns = rep("none", length(the_levels))
    the_patterns = setNames(the_patterns, the_levels)
  }
  
  # Adjust text justification based on angle
  if(is.null(vjs) | is.null(hjs)){
    if(the_angle != 0){
      vjs = 1
      hjs = 0.5
    } else {
      vjs = 1
      hjs = 1
    }
  }
  
  # Create the plot with geom_rect_pattern
  pp <- ggplot(the_df, aes(y = y_column, fill = y_column)) +  # Fill by y_column (Clusters)
    geom_rect_pattern(
      aes(xmin = Start, xmax = End, 
          ymin = as.numeric(factor(y_column)) - 0.44, 
          ymax = as.numeric(factor(y_column)) + 0.44,
          pattern = f_column),  # Pattern by f_column
      position = position_dodge(width = 0.9),
      pattern_fill = "black",
      pattern_angle = 45,
      colour = bar_border_color,  # Add border color to each bar
      size = bar_border_size      # Add border thickness to each bar
    ) +
    labs(x = x_title, y = y_title) +
    theme_classic() +
    scale_fill_manual(values = the_colors) +      # Colors for Clusters
    scale_pattern_manual(values = the_patterns) + # Patterns for Folds
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    theme(
      axis.ticks = element_line(colour = "black"), 
      legend.title = element_text(face = "bold", colour = "black", size = legend_title_size),
      axis.text.x = element_text(angle = the_angle, vjust = vjs, hjust = hjs, face = "bold", size = x_text_size, color = "black"),
      axis.text.y = element_text(face = "bold", colour = "black", size = y_text_size),
      axis.title.x = element_text(face = "bold", size = x_title_size),
      axis.title.y = element_text(face = "bold", size = y_title_size),
      panel.background = element_rect(fill = "white", colour = "white"), 
      plot.background = element_rect(fill = "white", colour = "white"), 
      legend.text = element_text(face = "bold", size = legend_text_size)
      #panel.border = element_rect(colour = "black", fill = NA, size = 1  # Add square border)
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = the_line1_shape, color = the_line1_color, size = the_line1_size) +
    geom_vline(xintercept = c(-fold_threshold, fold_threshold), linetype = the_line2_shape, color = the_line2_color, size = the_line2_size) +
    guides(fill = guide_legend(title = "Clusters"), pattern = guide_legend(title = "Folds"))
  
  # Handle legend positioning and visibility
  if(!showLegend || legend_position == "none") {
    pp <- pp + theme(legend.position = "none")
  } else {
    valid_positions <- c("top", "right", "bottom", "left")
    if(!legend_position %in% valid_positions) {
      warning("Invalid legend_position. Using 'top' instead.")
      legend_position <- "top"
    }
    pp <- pp + theme(
      legend.position = legend_position,
      legend.justification = "center"
    )
  }
  
  print(x_limits)
  return(pp)
}

rename_idents <- function(current_idents, current_idents_vector, new_idents_vector) {
  # Check if the lengths of the vectors match
  if (length(current_idents_vector) != length(new_idents_vector)) {
    stop("The lengths of the current and new identities vectors must match.")
  }
  rename_map <- setNames(new_idents_vector, current_idents_vector)
  current_idents <- current_idents
  new_idents <- rename_map[as.character(current_idents)]
  return(new_idents )
}

process_results <- function(results_list, control_condition, p_val_adj_cutoff = 0.1, log2fc_cutoff = 0.25, colors = NULL, the_size = 4){
  filtered_results <- lapply(results_list, function(df) {
    if (!is.null(df)) {
      subset(df, p_val_adj < p_val_adj_cutoff & abs(avg_log2FC) >= 
               log2fc_cutoff)
    }
    else {
      NULL
    }
  })
  process_comparison <- function(comp_name, df) {
    y <- unlist(strsplit(comp_name, split = "@"))
    the_comp <- y[1]
    y <- unlist(strsplit(comp_name, split = "__vs__"))
    the_cluster <- gsub(paste0(control_condition, "@"), "", y[2])
    if (!is.null(df)) {
      ups <- sum(df$avg_log2FC > 0)
      dws <- sum(df$avg_log2FC < 0)
    }
    else {
      ups <- 0
      dws <- 0
    }
    data.frame(upS = ups, dwS = dws, Cluster = the_cluster, Comps = the_comp)
  }
  degs_list <- mapply(process_comparison, names(filtered_results), filtered_results, SIMPLIFY = FALSE)
  degs_df <- do.call(rbind, degs_list)
  upS <- degs_df[, c("upS", "Cluster", "Comps")]
  dwS <- degs_df[, c("dwS", "Cluster", "Comps")]
  upS$upDws <- "UP"
  dwS$upDws <- "DW"
  colnames(upS) <- c("Counts", "Cluster", "Comps", "upDws")
  colnames(dwS) <- c("Counts", "Cluster", "Comps", "upDws")
  plot_data <- rbind(upS, dwS)
  gg <- ggplot(plot_data, aes(x = Cluster, y = Counts, shape = upDws)) + 
    geom_point(aes(color = Comps), size = the_size) + theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust = 1), plot.background = element_rect(fill = "white", 
    colour = "white"))
  if (!is.null(colors)) {
    gg <- gg + scale_color_manual(values = colors)
  }
  return(list(plot = gg, data = degs_df))
}

getDEGs <- function(the_obj, the_treatment, the_control, the_cluster_column, the_sample = "orig.ident", min_pct = 0.01) {
  Idents(the_obj) <- paste(the_obj@meta.data[[the_sample]], the_obj@meta.data[[the_cluster_column]], sep = "@")  
  comps <- c(
    paste(
        paste(the_treatment, unique(the_obj@meta.data[[the_cluster_column]]), sep = "@"), 
        paste(the_control, unique(the_obj@meta.data[[the_cluster_column]]), sep = "@"), 
        sep = "__vs__")
  )
  the.counts <- table(Idents(the_obj))
  degs <- foreach(the_comp = comps) %dopar% {
    the_idents <- unlist(strsplit(the_comp, split = "__vs__"))
    if (sum(the_idents %in% names(the.counts)) == 2) {
      if (the.counts[names(the.counts) == the_idents[1]] >= 3 & the.counts[names(the.counts) == the_idents[2]] >= 3) {
        dgs <- FindMarkers(the_obj, ident.1 = the_idents[1], ident.2 = the_idents[2], min.pct = min_pct)
        return(dgs)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }
  names(degs) <- comps
  return(degs)
}

remGenes = function(the_dfr, pct_1 = 0.1, pct_2 = 0.1){
    the_dfr = the_dfr[the_dfr$pct.1 >= pct_1 | the_dfr$pct.2 >= pct_2,]
    return(the_dfr)
}

generate_deg_plot = function (
    results_list, 
    control_condition, 
    p_val_adj_cutoff = 0.1, p_val_cutoff = NULL,
    log2fc_cutoff = 0.25, 
    colors = c("gray60", "gray80"), 
    the_size = 4, 
    cluster_order = NULL, 
    axis_text_size = 12, 
    axis_title_size = 14, 
    keep_all = NULL, 
    lg_size = NULL, 
    lg_symbol_size = 4) 
{
    filtered_results <- lapply(results_list, function(df) {
        if (!is.null(df)) {
            if(is.null(p_val_cutoff)){
                subset(df, p_val_adj < p_val_adj_cutoff & abs(avg_log2FC) >= log2fc_cutoff)
            }else{
                subset(df, p_val < p_val_cutoff & abs(avg_log2FC) >= log2fc_cutoff)
            }
        } else {
            NULL
        }
    })
    
    process_comparison <- function(comp_name, df) {
        y <- unlist(strsplit(comp_name, split = "@"))
        the_comp <- y[1]
        y <- unlist(strsplit(comp_name, split = "__vs__"))
        the_cluster <- gsub(paste0(control_condition, "@"), "", y[2])
        if (!is.null(df)) {
            ups <- sum(df$avg_log2FC > 0)
            dws <- sum(df$avg_log2FC < 0)
        } else {
            ups <- 0
            dws <- 0
        }
        data.frame(upS = ups, dwS = dws, Cluster = the_cluster, Comps = the_comp)
    }
    
    degs_list <- mapply(process_comparison, names(filtered_results), filtered_results, SIMPLIFY = FALSE)
    degs_df <- do.call(rbind, degs_list)
    
    upS <- degs_df[, c("upS", "Cluster", "Comps")]
    dwS <- degs_df[, c("dwS", "Cluster", "Comps")]
    upS$UpDown <- "Up"
    dwS$UpDown <- "Down"
    colnames(upS) <- c("Counts", "Cluster", "Comparisons", "UpDown")
    colnames(dwS) <- c("Counts", "Cluster", "Comparisons", "UpDown")

    if(is.null(lg_size)) lg_size = axis_text_size
    if (!is.null(keep_all)) {
        keep_cells <- (upS$Counts >= keep_all) | (dwS$Counts >= keep_all)
        upS <- upS[keep_cells, ]
        dwS <- dwS[keep_cells, ]
    }

    plot_data <- rbind(upS, dwS)
    cluster_order = cluster_order[cluster_order %in% plot_data$Cluster]
    if (!is.null(cluster_order)) {
        plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)
    }

    gg <- ggplot(plot_data, aes(x = Cluster, y = Counts, shape = UpDown)) + 
        geom_point(aes(color = Comparisons, fill = Comparisons), size = the_size) + 
        scale_shape_manual(values = c(Up = 24, Down = 25)) + 
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, face = "bold", colour = "black", size = axis_text_size), 
              axis.text.y = element_text(angle = 90, vjust = 0, hjust = 1, face = "bold", colour = "black", size = axis_text_size), 
              axis.ticks = element_line(colour = "black"), 
              panel.background = element_rect(fill = "white", colour = "white"), 
              plot.background = element_rect(fill = "white", colour = "white"), 
              panel.grid.major = element_line(colour = "gray90"), 
              panel.grid.minor = element_line(colour = "gray90"), 
              legend.text = element_text(face = "bold", colour = "black", size = lg_size), 
              legend.title = element_text(face = "bold", colour = "black"), 
              axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)), 
              axis.title.y = element_text(size = axis_title_size, face = "bold")) + 
        ylim(0, NA) + 
        guides(shape = guide_legend(override.aes = list(fill = NA, size = lg_symbol_size)))
    
    if (!is.null(colors)) {
        gg <- gg + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
    }
    
    return(list(gg,plot_data))
}

get_DEG_plot_data <- function (
    results_list, 
    control_condition, 
    p_val_adj_cutoff = 0.1, 
    p_val_cutoff = NULL,
    min_percent_cutoff = NULL,
    log2fc_cutoff = 0.25, 
    colors = c("gray60", "gray80"), 
    the_size = 4, 
    cluster_order = NULL, 
    axis_text_size = 12, 
    axis_title_size = 14, 
    keep_all = NULL, 
    lg_size = NULL, 
    lg_symbol_size = 4, 
    h_line = NULL) 
{
    filtered_results <- lapply(results_list, function(df) {
        if (!is.null(df)) {
            if(is.null(p_val_cutoff)){
                the_res = subset(df, p_val_adj < p_val_adj_cutoff & abs(avg_log2FC) >= log2fc_cutoff)
            } else {
                the_res = subset(df, p_val < p_val_cutoff & abs(avg_log2FC) >= log2fc_cutoff)
            }
        } else {
            the_res = NULL
        }
        if(!is.null(the_res) & !is.null(min_percent_cutoff)){
            the_res = subset(the_res, pct.1 >= min_percent_cutoff |pct.2 >= min_percent_cutoff) 
        }
        return(the_res)
    })
    #
    process_comparison <- function(comp_name, df) {
        y <- unlist(strsplit(comp_name, split = "@"))
        the_comp <- y[1]
        y <- unlist(strsplit(comp_name, split = "__vs__"))
        the_cluster <- gsub(paste0(control_condition, "@"), "", y[2])
        if (!is.null(df)) {
            ups <- sum(df$avg_log2FC > 0)
            dws <- sum(df$avg_log2FC < 0)
        } else {
            ups <- 0
            dws <- 0
        }
        data.frame(upS = ups, dwS = dws, Cluster = the_cluster, Comps = the_comp)
    }
    #    
    degs_list <- mapply(process_comparison, names(filtered_results), filtered_results, SIMPLIFY = FALSE)
    degs_df <- do.call(rbind, degs_list)
    #    
    upS <- degs_df[, c("upS", "Cluster", "Comps")]
    dwS <- degs_df[, c("dwS", "Cluster", "Comps")]
    upS$UpDown <- "Up"
    dwS$UpDown <- "Down"
    colnames(upS) <- c("Counts", "Cluster", "Comparisons", "UpDown")
    colnames(dwS) <- c("Counts", "Cluster", "Comparisons", "UpDown")
    #
    if(is.null(lg_size)) lg_size = axis_text_size
    if (!is.null(keep_all)) {
        keep_cells <- (upS$Counts >= keep_all) | (dwS$Counts >= keep_all)
        upS <- upS[keep_cells, ]
        dwS <- dwS[keep_cells, ]
    }
    #
    plot_data <- rbind(upS, dwS)
    cluster_order = cluster_order[cluster_order %in% plot_data$Cluster]
    if (!is.null(cluster_order)) {
        plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)
    }
    #
    # Calculate the sum of Up and Down for each Cluster and Comparison
    sum_data <- aggregate(Counts ~ Cluster + Comparisons, data = plot_data, sum)
    sum_data$UpDown <- "Total"
    plot_data <- rbind(plot_data, sum_data)
    #
    return(plot_data)
}

performDEG = function(the_obj, the_comps, the_columns, n_cpus = 2, the_sample = "orig.ident"){
    # Outer foreach loop
    the_degs <- foreach(comp = names(the_comps), .combine = 'c') %dopar% { 
    # .combine = "c" will return all results in one list
    setCPU(n_cpus)
    deg <- getDEGs(the_obj, the_comps[[comp]][1], the_comps[[comp]][2], the_columns, the_sample = the_sample)
    return(deg)
    }
    return(the_degs)
}

getGeneDfr <- function(results_list, gene_name) {
    gene_rows <- list()
    for (comp_name in names(results_list)){
        comp_df <- results_list[[comp_name]]
        if (is.null(comp_df)) {
            next
        }
    if (gene_name %in% rownames(comp_df)) {
        gene_row <- comp_df[gene_name, , drop = FALSE]
    } else {
        gene_row <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(comp_df)))
        colnames(gene_row) <- colnames(comp_df)
        for (col in colnames(comp_df)) {
            if (!is.numeric(comp_df[[col]])){
                gene_row[[col]] <- NA
            }
        }
    }
    gene_row$the_comp <- comp_name
    gene_row$gene_name <- gene_name
    gene_rows[[comp_name]] <- gene_row
  }
  gene_dfr <- do.call(rbind, gene_rows)
  return(gene_dfr)
}
generate_deg_plot =function (results_list, control_condition, p_val_adj_cutoff = 0.1, 
    log2fc_cutoff = 0.25, colors = c("gray60", "gray80"), the_size = 4, 
    cluster_order = NULL, axis_text_size = 12, axis_title_size = 14) 
{
    filtered_results <- lapply(results_list, function(df) {
        if (!is.null(df)) {
            subset(df, p_val_adj < p_val_adj_cutoff & abs(avg_log2FC) >= 
                log2fc_cutoff)
        }
        else {
            NULL
        }
    })
    process_comparison <- function(comp_name, df) {
        y <- unlist(strsplit(comp_name, split = "@"))
        the_comp <- y[1]
        y <- unlist(strsplit(comp_name, split = "__vs__"))
        the_cluster <- gsub(paste0(control_condition, "@"), "", 
            y[2])
        if (!is.null(df)) {
            ups <- sum(df$avg_log2FC > 0)
            dws <- sum(df$avg_log2FC < 0)
        }
        else {
            ups <- 0
            dws <- 0
        }
        data.frame(upS = ups, dwS = dws, Cluster = the_cluster, 
            Comps = the_comp)
    }
    degs_list <- mapply(process_comparison, names(filtered_results), 
        filtered_results, SIMPLIFY = FALSE)
    degs_df <- do.call(rbind, degs_list)
    upS <- degs_df[, c("upS", "Cluster", "Comps")]
    dwS <- degs_df[, c("dwS", "Cluster", "Comps")]
    upS$UpDown <- "Up"
    dwS$UpDown <- "Down"
    colnames(upS) <- c("Counts", "Cluster", "Comparisons", "UpDown")
    colnames(dwS) <- c("Counts", "Cluster", "Comparisons", "UpDown")
    plot_data <- rbind(upS, dwS)
    if (!is.null(cluster_order)) {
        plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)
    }
    gg <- ggplot(plot_data, aes(x = Cluster, y = Counts, shape = UpDown)) + 
        geom_point(aes(color = Comparisons, fill = Comparisons), 
            size = the_size) + scale_shape_manual(values = c(Up = 24, 
        Down = 25)) + theme(axis.text.x = element_text(angle = 60, 
        vjust = 1, hjust = 1, face = "bold", colour = "black", size = axis_text_size), 
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 1, 
            face = "bold", colour = "black", size = axis_text_size), 
        axis.ticks = element_line(colour = "black"), 
        panel.background = element_rect(fill = "white", colour = "white"), 
        plot.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_line(colour = "gray90"), 
        panel.grid.minor = element_line(colour = "gray90"), 
        legend.text = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(face = "bold", colour = "black"), 
        axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)), 
        axis.title.y = element_text(size = axis_title_size, face = "bold")) + 
        ylim(0, NA) + guides(shape = guide_legend(override.aes = list(fill = NA)))
    if (!is.null(colors)) {
        gg <- gg + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
    }
    return(gg)
}

###### GO-KEGG ANALYSES

# topGO functions
topDiffGenes <- function(x) return(x > 0)
topGOanalysis <- function(uniGenes = NULL, output_folder = NULL, minNodeSize = 5, annotFunction = annFUN.org, ids = "ensembl", degDfr = NULL, firstSigNodes = 10, useInfo = "all", algoComp = NULL){
    for(goTerm in c("BP","CC","MF")){   
        print(goTerm)
        sampleGOdata <- new(
        "topGOdata", 
        description = output_folder,
        ontology = goTerm,
        allGenes = uniGenes,
        geneSel = topDiffGenes,
        nodeSize = minNodeSize,
        annot = annotFunction,
        mapping = orgInfo$orgDb,
        ID = ids)
        print(c("Now analysing",goTerm," by topGO ..."))
        runTopGO(sampleGOdata, goTerm, output_folder, degDfr = degDfr, firstSigNodes = firstSigNodes, useInfo = "all", algoComp = algoComp)
    }
}
runTopGO <- function(sampleGOdata, goTerm, output_folder, degDfr = NULL, firstSigNodes = firstSigNodes, useInfo = "all", algoComp = NULL){
    dir.create(output_folder ,showWarnings = FALSE)
    # algorithms and statistics that can be used together!
    #algoComp <- rbind(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    #               c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    #               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #               c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    #               c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    #               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
    # algorithms and statistics that have been tested!
    #algoComp <- 
    #   rbind(
    #   c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
    #   c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
    #   c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #   c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
    #   c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
    #   c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
    #   )
    # algorithms and statistics chosen to be used!
    if(is.null(algoComp)){
        algoComp <- 
        rbind(
        c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
        )
        rownames(algoComp) <- c("classic", "elim", "weight", "weight01", "lea", "parentchild")
        colnames(algoComp) <- c("fisher", "z", "ks", "t", "globaltest", "category", "sum", "ks.ties")
    }
    #algoComp
    if(sum(algoComp) == 0){
        stop("ERROR: no algorithms and test found to use!")
    }else{
        results <- vector("list", sum(algoComp))
    }
    h = 0
    for(i in 1:nrow(algoComp)){
        for(j in 1:ncol(algoComp)){
            if(algoComp[i,j]){
                h = h + 1
                names(results)[h] <- paste(rownames(algoComp)[i], colnames(algoComp)[j], sep = "")
                print(names(results)[h])
                results[h] <- runTest(sampleGOdata, algorithm = rownames(algoComp)[i], statistic = colnames(algoComp)[j], cutOff = 0.01)
            }
        }
    }
    allGO = usedGO(object = sampleGOdata) 
    # There are 48 possible comparisons with different statics and algorithms! User must define the allRes length using the following template!!!
    #allRes <- GenTable(sampleGOdata, 
    #               results[[1]], results[[2]], results[[3]], results[[4]], results[[5]],
    #               results[[6]], results[[7]], results[[8]], results[[9]], results[[10]],
    #               results[[11]], results[[12]], results[[13]], results[[14]], results[[15]],
    #               results[[16]], results[[17]], results[[18]],results[[19]], results[[20]], 
    #               results[[21]], results[[22]], results[[23]],results[[24]], results[[25]],
    #               results[[26]], results[[27]], results[[28]],results[[29]], results[[30]],
    #               results[[31]], results[[32]], results[[33]],results[[34]], results[[35]],
    #               results[[36]], results[[37]], results[[38]],results[[39]], results[[40]],
    #               results[[41]], results[[42]], results[[43]],results[[44]], results[[45]],
    #               results[[46]], results[[47]], results[[48]]
    #                   orderBy = 1, ranksOf = 1, topNodes = length(allGO), numChar = 100)
    allRes <- GenTable(sampleGOdata, results[[1]], results[[2]], orderBy = 1, ranksOf = 1, topNodes = length(allGO), numChar = 100)
    colnames(allRes)[6:ncol(allRes)] <- names(results)
    for(h in 6:ncol(allRes)) allRes[allRes[,h] == "< 1e-30",h] = 1e-30 
    for(h in 6:ncol(allRes)) class(allRes[,h]) = "numeric" 
    for(i in 1:sum(algoComp)){
        goAnalysis <- names(results[i])
        outputPie <- paste("topGO_", goTerm,"_", goAnalysis, "_pieChart.pdf", sep = "")
        outputNodes <- paste("topGO_", goTerm,"_", goAnalysis, "_nodes.pdf", sep = "")
        goDfr <- allRes
        colnames(goDfr)[1] <- "GO"
        goDfr <- goDfr[goDfr[,which(colnames(allRes) == goAnalysis)] < 0.05,]
        writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder,  resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("topGO_",goTerm,"_", goAnalysis))        
        drawPie(allRes, outputName = outputPie, output_folder = output_folder, 
                pvalCol = which(colnames(allRes) == goAnalysis), 
                labCol = 2, countCol = 4, totalCol = 3)
        
        pdf(paste(output_folder,outputNodes, sep = "/"), width = 21, height = 21, pointsize = 5)
        showSigOfNodes(sampleGOdata, score(results[[i]]), firstSigNodes = firstSigNodes, useInfo = useInfo)
        dev.off()
    }
    output = paste(output_folder,"/topGO_",goTerm,"_allResults.tab",sep = "")
    write.table(allRes, output,sep = "\t")
    #return(allRes)
}
writeGOs <- function(goDfr = NULL, orgDb = NULL, output_folder = NULL, resultDfr = NULL, geneIDtype = "ENSEMBL", goType = NULL){
    orgDb <- eval(parse(text = orgDb))
    keyid2GO <- select(orgDb, keys = keys(orgDb,keytype = geneIDtype), column = c("GO"), keytype = geneIDtype)
    keyid2GO <- keyid2GO[,-c(3:4)]
    head(keyid2GO)
    colnames(keyid2GO) <- c("geneIds","GO")
    tmpVector = as.character(unname(unlist(goDfr[1,])))
    colnames(goDfr)[which(startsWith(tmpVector,"GO:"))] <- "GO"
    colnames(goDfr)[which(tolower(colnames(goDfr)) == "term")] <- "Term"
    if(nrow(goDfr) != 0){
        lt <- vector("list",length = length(goDfr$GO))
        names(lt) <- gsub("[/:]","",goDfr$GO)
        for(i in 1:length(names(lt))){
            lt[[i]] <- data.frame()
        }
        for(j in 1:nrow(goDfr)){
            goId <- goDfr$GO[j]
            geneIds <- keyid2GO$geneIds[keyid2GO$GO == goId]        
            tmpRes <- resultDfr[rownames(resultDfr) %in% geneIds,]
            tmpGoId <- gsub(":","",goId)
        if(nrow(tmpRes) > 0){
            prefix <- goDfr$Term[goDfr$GO == goId]
            outFile <- paste0(tmpGoId,"_",goDfr$Term[j])
            outFile  <- gsub("[/ :]","", outFile) #replace any special character with ""
            outFile <- strtrim(outFile,30)        #trimming might be better to prevent of writing into files!   
            tmpRes <- merge(keyid2GO, tmpRes, all.x = FALSE, all.y = TRUE, by.x = "geneIds", by.y = "row.names")
            tmpRes <- tmpRes[tmpRes$GO == goId,]
            tmpRes[,"Term"] <- rep(goDfr$Term[j], each = nrow(tmpRes))
            lt[tmpGoId][[1]] <- tmpRes
            names(lt)[which(names(lt) == tmpGoId)] <- outFile
        }else{
            lt[which(names(lt) == tmpGoId)] <- NULL
        }
    }
    lt <- lt[nchar(names(lt)) != 9]
    if(length(lt) >= 1) write.xlsx(lt, file = paste0(output_folder,"/",paste0(goType,"_genes.xlsx")))
    }
}

### Plotting GO and KEGG outputs
addBold = function() return(theme(text = element_text(face = "bold")))
### GOKEGG Plots
get_GOKEGG_dfr = function(the_out_df, row_cluster = T, col_cluster = T){
    get_top10 <- function(data) {
        top10 <- data %>%
        arrange(desc(pval)) %>%
        slice_head(n = 50)
        return(top10)
    }
    reshaped_df <- the_out_df %>%
    pivot_longer(
        cols = starts_with("Pvalue") | starts_with("Count"), 
        names_to = c(".value", "CellID"), 
        names_pattern = "(Pvalue|Count)_(.*)"
    )
    colnames(reshaped_df) = c("Term","GOKEGGID","cluster","pval","ratio")
    reshaped_df = reshaped_df[,-1]
    #reshaped_df$pval[reshaped_df$pval >= 0.05] = 10
    reshaped_df$pval = -log10(reshaped_df$pval)
    # Check the reshaped data frame
    #head(reshaped_df)
    # 
    # Apply the function to each cluster
    df_top10 <- reshaped_df %>%
      group_by(cluster) %>%
      group_modify(~ get_top10(.x)) %>%
      ungroup()
    #
    # Pivot the data to wide format for clustering
    df_wide <- df_top10 %>%
      dplyr::select(GOKEGGID, cluster, pval) %>%
      pivot_wider(names_from = cluster, values_from = pval)
    #
    # Replace NAs with a large value (e.g., the max value in the data) for clustering purposes
    max_value <- max(df_wide %>% dplyr::select(-GOKEGGID), na.rm = TRUE)
    df_wide[is.na(df_wide)] <- max_value
    #
    # Perform hierarchical clustering only on numeric columns (excluding GOKEGGID)
    numeric_data <- df_wide %>% dplyr::select(-GOKEGGID)
    #
    #
    if(row_cluster){
        row_order <- hclust(dist(numeric_data))$order
        # Reorder GOKEGGID based on clustering
        df_wide <- df_wide[row_order,]

    }
    if(col_cluster){
        ordered_GOKEGGID <- sort(df_wide$GOKEGGID)
        # Reorder df_long_unique based on clustered GOKEGGID
        df_long_clustered <- df_top10 %>%
            mutate(GOKEGGID = factor(GOKEGGID, levels = ordered_GOKEGGID))
    }
    df_long_clustered$ratio[df_long_clustered$ratio >= 20] = 21

    #
    # Check the updated data frame
    #head(df_long_clustered)
    #
    # Assuming df_long_clustered is already available
    # Filter out rows where ratio is 0 and adjust the color based on -log10(p-value)
    df_long_clustered <- df_long_clustered %>%
      filter(ratio != 0) %>%
      mutate(
        Color = case_when(
          pval < 1.30102999566398 ~ "black",                 # p-value >= 0.05
          pval >= 1.30102999566398 & pval < 2 ~ "red",       # 0.05 > p-value >= 0.01
          pval >= 2 ~ "blue"                      # p-value < 0.01
        ),
        pval_range = case_when(
          pval < 1.30102999566398 ~ "p_val >= 0.05",
          pval >= 1.30102999566398 & pval < 2 ~ "0.05 > p_val >= 0.01",
          pval >= 2 ~ "p_val < 0.01"
        )
      )
    return(list("df_long_clustered" = df_long_clustered, "df_top10" = df_top10))
}
getFiles = function(base_folder,
                    the_subs = c("allDEGs", "upDEGs", "downDEGs"),
                    file_patterns = c("GOstats_BP_Up.tab", 
                   "GOstats_CC_Up.tab", 
                   "GOstats_MF_Up.tab", 
                   "GOstats_kegg_over.tab",
                   "GOstats_BP_Down.tab", 
                   "GOstats_CC_Down.tab", 
                   "GOstats_MF_Down.tab", 
                   "GOstats_kegg_under.tab")){
    parse_folder <- function(folder) {
        treatment_control <- strsplit(folder, "__vs__")[[1]]
        treatment_cluster <- strsplit(treatment_control[1], "@")[[1]]
        control_cluster <- strsplit(treatment_control[2], "@")[[1]]
        
        data.frame(
            Folder = folder,
            Treatment = treatment_cluster[1],
            Control = control_cluster[1],
            Cluster = treatment_cluster[2],
            stringsAsFactors = FALSE
        )
    }
    combine_folders_files <- function(folder_df, combined_paths) {
        combined_df <- do.call(rbind, lapply(seq_len(nrow(folder_df)), function(i) {
            folder_info <- folder_df[i, ]
            paths <- data.frame(do.call(rbind, lapply(combined_paths, function(path) {
                path_split <- strsplit(path, "/")[[1]]
                data.frame(
                    Folder = folder_info$Folder,
                    Treatment = folder_info$Treatment,
                    Control = folder_info$Control,
                    Cluster = folder_info$Cluster,
                    DEGs = path_split[1],
                    File = path_split[2],
                    stringsAsFactors = FALSE
                )
            })))
            paths
        }))
        return(combined_df)
    }
    # Function to extract category from File
    extract_category <- function(file) {
      if (grepl("BP_Up", file)) return("BP_Up")
      if (grepl("BP_Down", file)) return("BP_Down")
      if (grepl("CC_Up", file)) return("CC_Up")
      if (grepl("CC_Down", file)) return("CC_Down")
      if (grepl("MF_Up", file)) return("MF_Up")
      if (grepl("MF_Down", file)) return("MF_Down")
      if (grepl("kegg_over", file)) return("KEGG_Up")
      if (grepl("kegg_under", file)) return("KEGG_Down")
      return(NA)
    }
    the_folders <- dir(base_folder)
    combinations <- expand.grid(the_subs, file_patterns, stringsAsFactors = FALSE)
    combined_paths <- apply(combinations, 1, function(x) paste(x, collapse = "/"))
    folder_df <- do.call(rbind, lapply(the_folders, parse_folder))
    combinations_df <- combine_folders_files(folder_df, combined_paths)
    combinations_df <- combinations_df %>% mutate(the_full_path = paste(Folder, DEGs, File, sep = "/"))
    combinations_df$the_full_path = paste(base_folder,combinations_df$the_full_path, sep = "/")
    # Extract category and create a column
    combinations_df <- combinations_df %>%
      mutate(Category = sapply(File, extract_category))
    return(combinations_df)
}                         
# Function to create data frames for each combination
create_combined_df <- function(combinations_df, deg_type, category) {
    # Function to read tab file and return a data frame
    read_tab_file <- function(file_path) {
        if (file.exists(file_path)) {
            dff <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
            dff = dff[,c(1,2,5)]
        } else {
            dff <- data.frame(GO = character(), Pvalue = numeric(), Count = integer(), stringsAsFactors = FALSE)
        }
        return(dff)
    }
    subset_df <- combinations_df %>% filter(DEGs == deg_type, Category == category)
    #subset_df
    combined_data <- list()
    for (i in seq_len(nrow(subset_df))) {
        file_path <- subset_df$the_full_path[i]
        cluster <- subset_df$Cluster[i]
        data <- read_tab_file(file_path) %>% mutate(Cluster = cluster)
        combined_data[[i]] <- data
    }
    # Filter out empty data frames
    combined_data <- combined_data[sapply(combined_data, nrow) > 0]
    #
    if (length(combined_data) == 0) {
        return(NULL)
    }
    combined_df <- bind_rows(combined_data)
    if(colnames(combined_df)[1] %in% c("KEGGID")){
        combined_df <- combined_df %>% complete(KEGGID, Cluster, fill = list(Pvalue = 1, Count = 0))
        combined_df <- combined_df %>% 
                dplyr::select(KEGGID, Pvalue, Count, Cluster) %>% 
                pivot_wider(names_from = Cluster, values_from = c(Pvalue, Count), values_fill = list(Pvalue = 2, Count = -1))
    }else{
        combined_df <- combined_df %>% complete(GO, Cluster, fill = list(Pvalue = 1, Count = 0))
        combined_df <- combined_df %>% 
                dplyr::select(GO, Pvalue, Count, Cluster) %>% 
                pivot_wider(names_from = Cluster, values_from = c(Pvalue, Count), values_fill = list(Pvalue = 2, Count = -1))
    }
  return(combined_df)
}
upDw_plot = function(the_dfr, the.levels = NULL, y.max = 300, the.cols = NULL){
    if(is.null(the.cols)) the.cols = c("Up" = "gray70", "Down" = "gray50","Total"="gray30")
    custom_axis_text <- function(breaks, the_size = 14, the_levels = NULL) {
        element_text(
            angle = 45, vjust = 1, hjust = 1, face = "bold", size = the_size,
            colour = ifelse(breaks %in% c("SRE", "SRN"), 
            "red",
        ifelse(breaks %in% c("NN#2"), 
            "red", "black"))
        )
    }
    #
    if(!is.null(the.levels)){
        the_dfr$UpDown = factor(the_dfr$UpDown, levels = the.levels)
    }else{
        if(sum(c("Total","Down","Up") %in% the_dfr$UpDown) == 3){
            the_dfr$UpDown = factor(the_dfr$UpDown, levels = c("Total","Down","Up"))
        }else if(sum(c("Down","Up") %in% the_dfr$UpDown) == 2){
            the_dfr$UpDown = factor(the_dfr$UpDown, levels = c("Down","Up"))
        }else{
            stop(paste("Unknown levelss:", unique(the_dfr$UpDown)))
        }
    }
    # Create a bar plot with UpDown as part of the fill aesthetics
    ggplot(the_dfr, aes(x = Cluster, y = Counts, fill = UpDown)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.88) +
        #geom_bar(stat = "identity", position = "dodge",width = 1) +
        labs(y = "Number of DEGs", x = "Clusters") +
        theme_classic() + 
        scale_fill_manual(values = the.cols) +
        theme(
            axis.ticks = element_line(colour = "black"), 
            legend.title = element_text(face = "bold", colour = "black", size = 16),
            axis.text.x =  custom_axis_text(levels(the_res$Cluster), the_size = 16),
            axis.text.y = element_text(face = "bold",colour = "black",size = 16),
            axis.title.x = element_text(face = "bold", size = 20),
            axis.title.y = element_text(face = "bold", size = 20),
            panel.background = element_rect(fill = "white", colour = "white"), 
            plot.background = element_rect(fill = "white", colour = "white"), 
            #panel.grid.major.y = element_line(colour = "gray90"), 
            #panel.grid.minor = element_line(colour = "gray90"),
            legend.position = "top",
            legend.justification = "center",
            legend.text = element_text(face = "bold", size = 14)) + #ylim(0, 100) + 
            guides(shape = guide_legend(override.aes = list(fill = NA, size = 12))) + guides(fill = guide_legend(title = NULL)
            ) + geom_hline(yintercept = y.max, linetype = "dashed", color = "red", size = 1)
} 

###### VELOCYTO functions

### Velocyto
getVeloData = function(the_obj, the_loom_data){
    for(trt in names(the_loom_data)){
        x = the_loom_data[[trt]]
        x = x[, colnames(x) %in% colnames(the_obj)]
        rownames(x) = getUniqs(rownames(x))
        idx = match(colnames(the_obj), colnames(x))
        ids = rownames(x) %in% rownames(the_obj)
        the_loom_data[[trt]] = x[ids,idx]

    }
    return(the_loom_data)
}
runVelocyto = function(the_obj, loom_data, nCPUs = 1, the_out_folder = NULL, custom_colors = NULL, cluster2use = "main.cells", useMe = FALSE){

    spliced_matrix <- loom_data$spliced
    emb = the_obj@reductions$umap@cell.embeddings
    tff("histo.tiff", width = 5, height = 5, out = the_out_folder)
    hist(log10(colSums(spliced_matrix)), col='wheat', xlab='cell size')
    dev.off()
    vc.obj <- Pagoda2$new(spliced_matrix, modelType='plain',trim=10,log.scale=T)
    tff("adjustedVariance.tiff", width = 8, height = 4, out = tod)
    vc.obj$adjustVariance(plot=T, do.par=T, gam.k=10)
    dev.off()
    #
    vc.obj$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
    vc.obj$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
    #
    vc.obj$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
    #vc.obj$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T, n.cores = nCPUs)  
    #
    vc.obj$getEmbedding(type='PCA',embeddingType='UMAP',perplexity=50,verbose=T, n.cores = nCPUs)  
    #
    vc.obj$embeddings$PCA$umap = the_obj@reductions$umap@cell.embeddings
    #
    vc_clusters = vc.obj$clusters$PCA$multilevel
    sc_clusters = the_obj@meta.data[,cluster2use]
    names(sc_clusters) = colnames(the_obj)
    print(table(sc_clusters))
    if(useMe){
        the_clusters = sc_clusters
        color2use = custom_colors
        sffx = "sc"
        vc.obj$clusters$PCA$multilevel = sc_clusters
    }else{
        the_clusters = vc_clusters
        color2use = sccore::fac2col(the_clusters)
        sffx = "vc"
    }
    #
    tff(paste0("pca_plots_cellclusters_UMAP_",sffx,".tiff"), width = 5, height = 5, out = the_out_folder)
    vc.obj$plotEmbedding(type='PCA',embeddingType='UMAP',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
    dev.off()
    
    tff(paste0("pca_plots_depth_UMAP_",sffx,".tiff"), width = 5, height = 5, out = the_out_folder)
    vc.obj$plotEmbedding(type='PCA',embeddingType='UMAP',colors=vc.obj$depth, main='depth', n.cores = nCPUs)
    dev.off()

    tff(paste0("pca_plots_cellclusters_umaps_",sffx,".tiff"), width = 5, height = 5, out = the_out_folder)
    vc.obj$plotEmbedding(type='PCA',embeddingType='umap',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
    dev.off()
    
    tff(paste0("pca_plots_depth_umaps_",sffx,".tiff"), width = 5, height = 5, out = the_out_folder)
    vc.obj$plotEmbedding(type='PCA',embeddingType='umap',colors=vc.obj$depth, main='depth', n.cores = nCPUs)
    dev.off()
    #
    spliced_matrix <- loom_data$spliced
    unspliced_matrix <- loom_data$unspliced
    spliced_matrix <- spliced_matrix[,rownames(vc.obj$counts)] 
    unspliced_matrix <- unspliced_matrix[,rownames(vc.obj$counts)]
    #
    emb_umap <- vc.obj$embeddings$PCA$umap
    cell.dist <- as.dist(1-armaCor(t(vc.obj$reductions$PCA)))
    #
    spliced_matrix <- filter.genes.by.cluster.expression(spliced_matrix,the_clusters,min.max.cluster.average = 0.2)
    unspliced_matrix <- filter.genes.by.cluster.expression(unspliced_matrix,the_clusters,min.max.cluster.average = 0.05)
    #
    fit.quantile <- 0.02
    vc.vel <- gene.relative.velocity.estimates(spliced_matrix,unspliced_matrix,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = nCPUs, , do.par=T)
    #
    tff(paste0("velo_embded_umap_",sffx,".tiff"), width = 5, height = 5, out = tod)
    show.velocity.on.embedding.cor(vc.obj$embeddings$PCA$UMAP, vc.vel, return.details = F,xlab = "umap_1", ylab = "umap_2", 
                               n=200,scale='sqrt',cell.colors=color2use,cex=0.8,arrow.scale=3,
                               show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,
                               cell.border.alpha = 0.1, n.cores = nCPUs)
    dev.off()
    #
    tff(paste0("velo_embded_umaps_",sffx,".tiff"), width = 5, height = 5, out = tod)
    show.velocity.on.embedding.cor(vc.obj$embeddings$PCA$umap, vc.vel, return.details = F,xlab = "umap_1", ylab = "umap_2", 
                               n=200,scale='sqrt',cell.colors=color2use,cex=0.8,arrow.scale=3,
                               show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,
                               cell.border.alpha = 0.1, n.cores = nCPUs)
    dev.off()
    return(list("vc_obj" = vc.obj, "vc_vel" = vc.vel))
}
getLoomPaths <- function(prfx = "/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/rawData_featureCounts/", sffx = NULL) {
    if(is.null(sffx)) stop("please provide a suffix!")
    
    # Example usage
    #root_paths <- c(
    #  "Cavone_etal_2021_E-MTAB/",
    #  "LAB4933/",
    #  "LAB4996/"
    #)

    #loom_files <- getLoomPaths(sffx = root_paths)
    #loom_files

    root_paths = paste0(prfx, sffx)
    # Ensure root_paths is a character vector
  if (is.character(root_paths)) {
    root_paths <- as.list(root_paths) # Convert to a list if single string
  }
  
  # Initialize an empty list to store velocyto paths
  all_velocyto_paths <- list()
  
  # Loop over each root path
  for (the_root_path in root_paths) {
    # Get all directories recursively
    velocyto_paths <- list.dirs(path = the_root_path, recursive = TRUE, full.names = TRUE)
    # Filter only directories containing "velocyto"
    velocyto_paths <- grep("velocyto", velocyto_paths, value = TRUE)
    # Append to the collection
    all_velocyto_paths <- c(all_velocyto_paths, velocyto_paths)
  }
  
  # Extract loom file paths
  loom_files <- unlist(lapply(all_velocyto_paths, function(path) {
    dir_contents <- dir(path)
    loom_paths <- file.path(path, dir_contents)
    return(loom_paths)
  }))
  
  return(loom_files)
}


###### MONOCLE3
### start of functions to add pseudotime to a seurat object.
                            
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
#
get_earliest_principal_node <- function(cds, time_bin = "ERG", the_col_data = "main.cells"){
  cell_ids <- which(colData(cds)[, the_col_data] == time_bin)
  #cell_ids = which(colData(cds)[, "the_cells"] %in% time_bin)
  #cell_ids <- which(colData(cds)$main.cells == "ERG")
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}
#
getColors = function(the.colors, the.vector){
  idx = match(the.vector, names(the.colors))
  return(the.colors[idx])
}

runMonocle = function(my.seurat.obj, the_root_cell = "10", the_col_data = "cca_clusters", add_more_var = NULL){
    the_var_genes = rownames(my.seurat.obj@reductions[["pca"]]@feature.loadings)
    if(!is.null(add_more_var)){ 
        my.seurat.obj = FindVariableFeatures(my.seurat.obj, nfeatures = add_more_var) #,layer = "data")
        the_var_genes = union(the_var_genes, VariableFeatures(my.seurat.obj))
    }
    #gene_annotation <- as.data.frame(rownames(my.seurat.obj@reductions[["pca"]]@feature.loadings), row.names = rownames(my.seurat.obj@reductions[["pca"]]@feature.loadings))
    gene_annotation = as.data.frame(the_var_genes, row.names = the_var_genes)
    colnames(gene_annotation) <- "gene_short_name"
    cell_metadata <- my.seurat.obj@meta.data
    New_matrix <- my.seurat.obj[["RNA"]]$counts
    New_matrix <- New_matrix[the_var_genes, ]
    expression_matrix <- New_matrix
    cds_from_seurat <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
    recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
    names(recreate.partition) <- cds_from_seurat@colData@rownames
    recreate.partition <- as.factor(recreate.partition)
    cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
    list_cluster <- my.seurat.obj@active.ident
    names(list_cluster) <- my.seurat.obj[["RNA"]]$data@Dimnames[[2]]
    cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
    cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- my.seurat.obj@reductions[["umap.cca"]]@cell.embeddings
    cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)
    cds_from_seurat <- order_cells(cds_from_seurat, 
                root_pr_nodes = get_earliest_principal_node(cds_from_seurat,time_bin = the_root_cell, the_col_data = the_col_data))
    return(cds_from_seurat)
}
### end
getTheCells = function(seurat_obj, gene, min.exp = 0){
    if(length(gene) == 1){
        tt = seurat_obj@assays$RNA$counts[rownames(seurat_obj) == gene,] > min.exp
    }else{
        tt = colSums(seurat_obj@assays$RNA$counts[rownames(seurat_obj) %in% gene,]) > min.exp
    }
    return(colnames(seurat_obj)[tt])
}
find_common_and_unique <- function(...){
  # Capture the input vectors as a list
  vec_list <- list(...)
  # Find elements that are common in at least two vectors
  all_elements <- unlist(vec_list)
  element_counts <- table(all_elements)
  common_elements <- names(element_counts[element_counts >= 2]) 
  # Find elements that are unique in only one vector
  unique_elements <- names(element_counts[element_counts == 1])
  # Return results as a list
  return(common_elements)
}
# General function to get cells and find common and unique cells
get_cells_by_genes <- function(the_obj, gene_list) {
    cell_groups <- list()
    cell_counter <- list()
    for (genes in gene_list) {
        min_exp <- tail(genes, n = 1)
        gene_names <- head(genes, n = -1)
        current_cells <- getTheCells(the_obj, gene_names, min.exp = min_exp)
        group_name <- paste(gene_names, collapse = "_")
        cell_groups[[group_name]] <- current_cells
        for (cell in current_cells) {
            if (is.null(cell_counter[[cell]])) {
                cell_counter[[cell]] <- 1
            }else{
                cell_counter[[cell]] <- cell_counter[[cell]] + 1
            }
        }
    }
  common_cells <- names(which(sapply(cell_counter, function(x) x > 1)))
  for (group in names(cell_groups)) {
    cell_groups[[group]] <- setdiff(cell_groups[[group]], common_cells)
  }
  cell_groups[["common"]] <- common_cells
  return(cell_groups)
}
# Function to get evenly spaced palette indices based on the length of the cells
get_spaced_indices <- function(palette_length, num_cells) {
  # Generate equally spaced indices
  indices <- seq(1, palette_length, length.out = num_cells)
  # Round indices to the nearest integer
  indices <- round(indices)
  return(indices)
}
get_distinct_colors <- function(num_cells) {
  # Use RColorBrewer's "Set1" palette for distinct colors, or "Paired"
  max_colors <- min(num_cells, brewer.pal.info["Set1", "maxcolors"]) # Set a limit to the max available colors
  
  # Generate the palette with the required number of distinct colors
  colors <- brewer.pal(max_colors, "Set1")
  
  # If more colors are needed than available in the palette, recycle or generate more
  if (num_cells > max_colors) {
    colors <- colorRampPalette(colors)(num_cells)  # Interpolate more colors if needed
  }
  
  return(colors)
}
get_earliest_principal_node <- function(cds, time_bin = "ERG", the_col_data = "main.cells"){
  cell_ids <- which(colData(cds)[, the_col_data] == time_bin)
  #cell_ids = which(colData(cds)[, "the_cells"] %in% time_bin)
  #cell_ids <- which(colData(cds)$main.cells == "ERG")
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}
getColors = function(the.colors, the.vector){
  idx = match(the.vector, names(the.colors))
  return(the.colors[idx])
}
get_monocle = function(the_seu_object, the_dims = 15, the_group = "orig.ident", 
                       n_max_comp = 4, red_method = "UMAP", the_bin = "10", the_col_data = "cca_clusters"){
    gene_annotation = as.data.frame(rownames(the_seu_object))
    colnames(gene_annotation) = "gene_short_name"
    rownames(gene_annotation) = gene_annotation$gene_short_name
    #
    cds <- new_cell_data_set(
      the_seu_object$RNA$counts,
      cell_metadata = the_seu_object@meta.data,
      gene_metadata = gene_annotation
      )
    #
    cds <- preprocess_cds(cds, num_dim = the_dims)
    cds <- align_cds(cds, alignment_group = the_group)
    cds <- reduce_dimension(cds, max_components = n_max_comp)
    cds <- cluster_cells(cds,reduction_method = red_method)
    cds <- learn_graph(cds)
    #
    cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, time_bin = the_bin, the_col_data = the_col_data))
    #cell_ids <- which(colData(cds)[, the_col_data] == the_bin)
    return(cds)
}
get_monocle <- function(the_seu_object, the_dims = 15, the_group = "orig.ident", 
                       n_max_comp = 4, red_method = "UMAP", the_bin = "10", the_col_data = "cca_clusters") {
    # Extract gene annotations
    gene_annotation <- as.data.frame(rownames(the_seu_object), stringsAsFactors = FALSE)
    colnames(gene_annotation) <- "gene_short_name"
    rownames(gene_annotation) <- gene_annotation$gene_short_name
    
    # Create a CellDataSet with the scaled data
    cds <- new_cell_data_set(
      expression_data = the_seu_object@assays$RNA$scale.data,
      cell_metadata = the_seu_object@meta.data,
      gene_metadata = gene_annotation
    )
    
    # Optionally, preprocess the data without normalization (skipping scaling)
    # Note: preprocess_cds does normalization and scaling, but we skip this
    # If you want to include PCA without scaling:
    cds <- preprocess_cds(cds, num_dim = the_dims, norm_method = "none", use_genes = VariableFeatures(the_seu_object))
    
    # Align and reduce dimensions
    cds <- align_cds(cds, alignment_group = the_group)
    cds <- reduce_dimension(cds, max_components = n_max_comp, reduction_method = red_method)
    
    # Cluster cells
    cds <- cluster_cells(cds, reduction_method = red_method)
    
    # Learn the trajectory graph
    cds <- learn_graph(cds)
    
    # Order cells along the trajectory
    cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, time_bin = the_bin, the_col_data = the_col_data))
    
    return(cds)
}

###### CELLCHAT

convert2csv = function(the_database){
  #the_database <- the_database.zebrafish # set the_database <- the_database.human if working on the human dataset
  interaction_input <- the_database$interaction
  complex_input <- the_database$complex
  cofactor_input <- the_database$cofactor
  geneInfo <- the_database$geneInfo
  write.csv(interaction_input, file = "interaction_input_the_database.csv")
  write.csv(complex_input, file = "complex_input_the_database.csv")
  write.csv(cofactor_input, file = "cofactor_input_the_database.csv")
  write.csv(geneInfo, file = "geneInfo_input_the_database.csv")
}
convert2db = function(){
  options(stringsAsFactors = FALSE)
  interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
  complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
  cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
  geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
  new_database <- list()
  new_database$interaction <- interaction_input
  new_database$complex <- complex_input
  new_database$cofactor <- cofactor_input
  new_database$geneInfo <- geneInfo
  return(new_database)
}

###### SIGNAC

parseGTF = function(the_gtf_path = NULL, orgId = "hsa"){
    if(is.null(the_gtf_path) | is.null(orgId)){
        stop("Please provide both; the GTF path and org ID")
    }
    # Load necessary libraries

    # Read the GTF file (assuming it's gzipped)
    gtf_data <- fread(the_gtf_path, header = FALSE, col.names = c("seqname", "source", 
                                                                  "feature", "start", "end", "score", "strand", "frame", "attribute"))
    # Filter for rows corresponding to genes
    genes <- gtf_data[feature == "gene"]
    # Extract Ensembl IDs and Gene Names
    ensembl_ids <- sapply(strsplit(genes$attribute, ";"), function(x) {
      ensembl <- grep("gene_id", x, value = TRUE)
      ensembl_id <- gsub("gene_id \"|\"", "", ensembl)
      return(ensembl_id)
    })
    #
    gene_names <- sapply(strsplit(genes$attribute, ";"), function(x) {
      gene <- grep("gene_name", x, value = TRUE)
      gene_name <- gsub("gene_name \"|\"", "", gene)
        if(length(gene_name) == 0) gene_name = "@"
      return(gene_name)
    })
    the_genes = data.frame("ENS" = ensembl_ids, "SYM" = gene_names)
    idx = grep("@", the_genes$SYM)
    the_genes$SYM[idx] = the_genes$ENS[idx]
    colnames(the_genes) = paste0(orgId, colnames(the_genes))
    return(the_genes)
}
# these 2 functions; getPFM, getTrx; are for multiome data
getPFM = function(tmp_path){
    
    # Get a list of motif position frequency matrices from the JASPAR database
    if(file.exists(paste0(tmp_path,"/the_pfm_all.rds"))){
        the_pfm = readRDS(paste0(tmp_path,"/the_pfm_all.rds"))
    }else{
        the_pfm <- getMatrixSet(
          x = JASPAR2020,
          opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
        )
        saveRDS(the_pfm, file = paste0(tmp_path,"/the_pfm_all.rds"))
    }
    return(the_pfm)
}
getTrx = function(tmp_path){
    
    # Get a list of motif position frequency matrices from the JASPAR database
    if(file.exists(paste0(tmp_path,"/the_pfm_all.rds"))){
        the_pfm = readRDS(paste0(tmp_path,"/the_pfm_all.rds"))
    }else{
        the_pfm = getPFM()
        saveRDS(the_pfm, file = paste0(tmp_path,"/the_pfm_all.rds"))
    }
    # Extract names for each motif
    the_names <- unlist(lapply(the_pfm, function(x) {
      if (inherits(x, "PFMatrixList")) {
        return(x@listData$remap_tf_name)
      } else {
        return(x@tags$remap_tf_name)
      }
    }))
    the_symbols = names(the_names)
    the_orthos = readRDS(paste0(tmp_path,"/dre2hsa.orthos.rds"))
    the_orthos$rwns = NULL
    #dim(the_orthos)
    the_orthos = the_orthos[the_orthos$hsaSYM %in% the_names,]
    idx = match(the_orthos$hsaSYM, the_names)
    the_orthos$tfSYM = names(the_names)[idx] 
    return(list("the_orthos" = the_orthos, "the_pfm" = the_pfm))
}
getPeaks <- function(coordinates_vector, query_coordinate, uPs = 5000, dWs = 5000) {
  # Split the query coordinate
  query_parts <- unlist(strsplit(query_coordinate, "-"))
  query_chromosome <- query_parts[1]
  query_start <- as.numeric(query_parts[2])
  query_end <- as.numeric(query_parts[3])
  
  # Calculate range
  range_start <- query_start - uPs
  range_end <- query_end + dWs
  
  # Split coordinates vector
  split_coordinates <- strsplit(coordinates_vector, "-")
  
  # Initialize empty vector to store matching coordinates
  matching_coordinates <- c()
  
  # Iterate through coordinates vector
  for (i in 1:length(split_coordinates)) {
    # Extract chromosome, start, and end
    chrom <- split_coordinates[[i]][1]
    start <- as.numeric(split_coordinates[[i]][2])
    end <- as.numeric(split_coordinates[[i]][3])
    
    # Check if within range and same chromosome
    if (chrom == query_chromosome && start >= range_start && end <= range_end) {
      matching_coordinates <- c(matching_coordinates, coordinates_vector[i])
    }
  }
  
  return(matching_coordinates)
}
getRange <- function(coordinates) {
  # Extract start and end values from each coordinate
  start_values <- sapply(strsplit(coordinates, "-"), function(x) as.numeric(x[2]))
  end_values <- sapply(strsplit(coordinates, "-"), function(x) as.numeric(x[3]))
  
  # Find the smallest start value and the largest end value
  range_start <- min(start_values)
  range_end <- max(end_values)
  
  # Extract chromosome number from the first coordinate
  chromosome_number <- unlist(strsplit(coordinates[1], "-"))[1]
  
  # Create the range
  coordinate_range <- paste(chromosome_number, range_start, range_end, sep = "-")
  
  return(coordinate_range)
}

###### SEURAT

con2CsMatrix = function(the.data){
  if(class(the.data) != "list"){
      the.data = as(the.data, "CsparseMatrix")
    }else{
      for(trt in names(the.data)){
        tmp = the.data[[trt]]
        tmp = as(tmp, "CsparseMatrix")
        the.data[[trt]] = tmp
      }
    return(the.data)
    }
  }
getMarkers = function(the.obj, the.idents, onlyPos = T, the.levels = NULL){
  if(is.null(the.levels)){ 
    the.levels = names(table(the.idents)) 
  }
  the.mrkrs = foreach(the.ident = the.levels) %dopar% {
    print(the.ident)
    nw.idents = the.idents
    nw.idents[nw.idents != the.ident] = "rem.cells"
    Idents(the.obj) = nw.idents
    mrkrs = FindMarkers(the.obj, ident.1 = the.ident, ident.2 = "rem.cells", only.pos = onlyPos)
    mrkrs$cluster = the.ident
    mrkrs$gene = rownames(mrkrs)
    return(mrkrs)
  }
  names(the.mrkrs) = names(table(the.idents))
  #
  tmp.mrkrs = the.mrkrs[[1]]
  tmp.mrkrs = tmp.mrkrs[order(tmp.mrkrs$avg_log2FC, decreasing = T),]
  for(i in seq(2, length(the.mrkrs))){
    tmp2 = the.mrkrs[[i]]
    tmp2 = tmp2[order(tmp2$avg_log2FC, decreasing = T),]
    tmp.mrkrs = rbind(tmp.mrkrs, tmp2)
  }
  #tmp.mrkrs = tmp.mrkrs[order(tmp.mrkrs$avg_log2FC, decreasing = T),]
  return(tmp.mrkrs)
}
find_doublets_by_scDblFinder = function(the.seu.obj, the.clusters = "seurat_clusters", the.dbr = 0.1, nreps = 12, min.reps = 3, n.iter = 6, n.dims = NULL, the.cpus = 4){
  # the.clusters; seurat_clusters (default), or use another type of cluster (e.g., main.cells)
  the.seeds = list(121,122,123,124,125,126,127,128,129,131,132,133)
  if(is.null(n.dims)) n.dims = the.seu.obj@commands$RunPCA.RNA$npcs
  the.counts = as(the.seu.obj[["RNA"]]$counts, "sparseMatrix")
  sce = SingleCellExperiment(list(counts = the.counts), colData = the.seu.obj@meta.data) #, rowData=DataFrame() ' one can add ensebembl and gene names'
  the.dfs = vector("list", nreps)
  names(the.dfs) = paste("rep", seq(1, nreps))
  ii = 1
  for(rplcs in names(the.dfs)){
    set.seed(the.seeds[[ii]])
    ii = ii + 1
    the.doublets = scDblFinder(sce, samples = "orig.ident",dbr = the.dbr, clusters = the.clusters, iter = n.iter, dims = n.dims,BPPARAM=MulticoreParam(the.cpus))
    the.doublets = the.doublets$scDblFinder.class
    the.dfs[[rplcs]] = the.doublets
  }
  the.counts = NULL
  sce = NULL
  gc()
  the.dfs = as.data.frame(the.dfs)
  the.dfs = rowSums(the.dfs == "singlet")
  the.dfs = the.dfs >= nreps-min.reps
  the.dfs[the.dfs == TRUE] = "Singlet"
  the.dfs[the.dfs != "Singlet"] = "Doublet"
  return(the.dfs)
}
getObjs = function(rawData, use_all_genes = TRUE, trt = NULL, n.dims = 20, resL = 0.5, vars2regress = c("nCount_RNA"), nVars = 2000, mito.prfx = "^mt-"){
    obj.list <- CreateSeuratObject(counts = rawData, project = trt, min.cells = 3, min.features = 200)
    obj.list[["percent.mt"]] <- PercentageFeatureSet(obj.list, pattern = mito.prfx)
    obj.list <- subset(obj.list, subset = nFeature_RNA > 200 & percent.mt < 20)
    print("Normalizing ...")
    obj.list <- NormalizeData(obj.list)
    print("Scaling ...") 
    obj.list <- FindVariableFeatures(obj.list, selection.method = "vst", nfeatures = nVars)
    if(use_all_genes){
      all.genes <- rownames(obj.list)
      obj.list <- ScaleData(obj.list, features = all.genes, vars.to.regress = vars2regress)
    }else{
      obj.list <- ScaleData(obj.list, vars.to.regress = vars2regress)
    }
    print("PCA ...")
    obj.list <- RunPCA(obj.list, features = VariableFeatures(object = obj.list))
    obj.list <- FindNeighbors(obj.list, dims = 1:n.dims)
    obj.list <- FindClusters(obj.list, resolution = resL)
    print("Umap ...")
    obj.list <- RunUMAP(obj.list, dims = 1:n.dims)
    return(obj.list)
}
getIntegrated = function(the.raw.data, use_all_genes = TRUE, minCells = 5, minFtr = 200, npcas = 50, ndims = 30, varGenes = 2000, vars2regress = c("nCount_RNA"), resLs = c(0.5,1,1.5,2), mito.prefix = "^mt-", ribo.prefix = NULL){
    the.obj <<- CreateSeuratObject(the.raw.data, min.cells = minCells, min.features = minFtr)
    if(!is.null(mito.prefix)) the.obj[["percent.mt"]] <<- PercentageFeatureSet(the.obj, pattern = mito.prefix)
    if(!is.null(ribo.prefix)) the.obj[["percent.rb"]] <<- PercentageFeatureSet(the.obj, features = ribo.percent)
    #
    the.obj <<- NormalizeData(the.obj)
    the.obj <<- FindVariableFeatures(the.obj, nfeatures = varGenes)
    if(use_all_genes){
      all.genes <- rownames(the.obj)
      the.obj <<- ScaleData(the.obj, features = all.genes, vars.to.regress = vars2regress)
    }else{
      the.obj <<- ScaleData(the.obj, vars.to.regress = vars2regress)
    }
    the.obj <<- RunPCA(the.obj, npcs = npcas)
    #
    the.obj <<- IntegrateLayers(object = the.obj, method = CCAIntegration,  orig.reduction = "pca", new.reduction = "integrated.cca",verbose = FALSE)
    the.obj <<- FindNeighbors(the.obj, reduction = "integrated.cca", dims = 1:ndims)
    #
    for(resL in resLs) the.obj <<- FindClusters(the.obj, resolution = resL, cluster.name = paste0("integrated_snn_res.",resL))
    #
    the.obj <<- RunUMAP(the.obj, reduction = "integrated.cca", dims = 1:ndims, reduction.name = "umap.cca", verbose = F)
    the.obj <<- JoinLayers(the.obj)
    Idents(the.obj) = the.obj$integrated_snn_res.1
    the.obj$cca_clusters = the.obj$integrated_snn_res.1
    return(the.obj)
}
getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}
readMatrix = function(thePath, trt){
    fls = dir(thePath)
    #
    if("matrix.mtx" %in% fls){
      matrix_file = "matrix.mtx"
    }else if("matrix.mtx.gz" %in% fls){
      matrix_file = "matrix.mtx.gz"
    }else{
      matrix_file = NULL
    }
    #
    if("barcodes.tsv" %in% fls){
      barcodes_file = "barcodes.tsv"
    }else if("barcodes.tsv.gz" %in% fls){
      barcodes_file = "barcodes.tsv.gz"
    }else{
      barcodes_file = NULL
    }
    #
    if("features.tsv" %in% fls){
      features_file = "features.tsv"
    }else if("features.tsv.gz" %in% fls){
      features_file = "features.tsv.gz"
    }else{
      barcodes_file = NULL
    }
    #
    if(sum(is.null(c(matrix_file, barcodes_file, features_file))) != 0){
      print("one of the following files is missing: ")
      print(c("Matrix File:   ", matrix_file))
      print(c("Barcodes File: ", barcodes_file))
      print(c("Features File: ", features_file))
    }else{    
    print(fls)
    mtx = readMM(paste0(thePath,"/",matrix_file))
    brc = read.delim(paste0(thePath,"/",barcodes_file), header = F)
    brc = brc$V1
    brc = gsub("-1$","",brc)
    brc = gsub("_1$","",brc)
    brc = paste(trt, brc, sep = "_")
    print(head(brc))
    ftr = read.delim(paste0(thePath,"/",features_file), header = F)
    rownames(mtx) = getUniqs(ftr$V2) #rwn
    colnames(mtx) = brc
    return(mtx)
  }
}
# use scDblFinder
find_doublets_by_scDblFinder = function(the.seu.obj, the.clusters = "seurat_clusters", the.dbr = 0.1, nreps = 3, min.reps = 0, n.iter = 3, n.dims = NULL){
  # the.clusters; seurat_clusters (default), or use another type of cluster (e.g., main.cells)
  if(is.null(n.dims)) n.dims = the.seu.obj@commands$RunPCA.RNA$npcs
  the.counts = as(the.seu.obj[["RNA"]]$counts, "sparseMatrix")
  sce = SingleCellExperiment(list(counts = the.counts), colData = the.seu.obj@meta.data) #, rowData=DataFrame() ' one can add ensebembl and gene names'
  the.dfs = vector("list", nreps)
  names(the.dfs) = paste("rep", seq(1, nreps))
  for(rplcs in names(the.dfs)){
    the.doublets = scDblFinder(sce, dbr = the.dbr, clusters = the.clusters, iter = n.iter, dims = n.dims)
    the.doublets = the.doublets$scDblFinder.class
    the.dfs[[rplcs]] = the.doublets
  }
  the.counts = NULL
  sce = NULL
  gc()
  the.dfs = as.data.frame(the.dfs)
  the.dfs = rowSums(the.dfs == "singlet")
  the.dfs = the.dfs >= nreps-min.reps
  the.dfs[the.dfs == TRUE] = "Singlet"
  the.dfs[the.dfs != "Singlet"] = "Doublet"
  return(the.dfs)
}
# use doubletFinder
find_doublets_by_DoubletFinder = function(the.seu.obj, npcas = NULL, ann2use = "seurat_clusters", the.perc = 0.075){
  if(is.null(npcas)) npcas = the.seu.obj@commands$RunPCA.RNA$npcs
  sweep.res.list <- paramSweep(the.seu.obj, PCs = 1:npcas, sct = FALSE)
  sweep.stats_list <- summarizeSweep(sweep.res.list, GT = FALSE)
  pk_values <- find.pK(sweep.stats_list) # TO DO: how to use pk_values and use as input for doubletFinder
  homotypic.prop <- modelHomotypic(the.seu.obj@meta.data[, ann2use])
  nExp_poi <- round(the.perc*nrow(the.seu.obj@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  the.seu.obj <- doubletFinder(the.seu.obj, PCs = 1:npcas, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(the.seu.obj@meta.data)[ncol(the.seu.obj@meta.data)] = "DF"
  return(the.seu.obj$DF)
}


####### VENN DIAGRAM
ggvenn <- function(data, columns = NULL,
                   show_elements = FALSE,
                   show_percentage = TRUE,
                   digits = 1,
                   fill_color = c("blue", "yellow", "green", "red"),
                   fill_alpha = .5,
                   stroke_color = "black",
                   stroke_alpha = 1,
                   stroke_size = 1,
                   stroke_linetype = "solid",
                   set_name_color = "black",
                   set_name_size = 6,
                   text_color = "black",
                   text_size = 4,
                   label_sep = ",",
                   count_column = NULL,
                   show_outside = c("auto", "none", "always"),
                   auto_scale = FALSE) {
  show_outside <- match.arg(show_outside)
  venn <- prepare_venn_data(data, columns, show_elements, show_percentage, digits,
                            label_sep, count_column, show_outside, auto_scale)
  g <- venn$shapes %>%
    mutate(group = LETTERS[group]) %>%
    ggplot() +
    geom_polygon(aes(x = x, y = y, group = group, fill = group),
                 alpha = fill_alpha) +
    geom_polygon(aes(x = x, y = y, group = group),
                 fill = NA,
                 color = stroke_color,
                 size = stroke_size,
                 alpha = stroke_alpha,
                 linetype = stroke_linetype)
  if (nrow(venn$labels) > 0) {
    g <- g +
      geom_text(data = venn$labels,
                aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
                color = set_name_color,
                size = set_name_size)
  }
  if (nrow(venn$texts) > 0) {
    g <- g +
      geom_text(data = venn$texts,
                aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
                color = text_color,
                size = text_size)
  }
  if (nrow(venn$segs) > 0) {
    g <- g +
      geom_segment(data = venn$segs,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = text_color,
                   size = 0.5)
  }
  g <- g +
    scale_fill_manual(values = fill_color) +
    guides(fill = "none") +
    coord_fixed() +
    theme_void()
  return(g)
}
gen_element_df_2 <- function() {
  df <- tribble(~name, ~A,    ~B,
                "A",   TRUE,  FALSE,
                "B",   FALSE, TRUE,
                "AB",  TRUE,  TRUE,
                "-",   FALSE, FALSE)
  stopifnot(all((df %>% dplyr::count(A, B) %>% with(n)) == 1))
  return(df %>% mutate(n = 0, text = ""))
}
gen_element_df_3 <- function() {
  df <- tribble(~name, ~A,    ~B,    ~C,
                "A",   TRUE,  FALSE, FALSE,
                "B",   FALSE, TRUE,  FALSE,
                "C",   FALSE, FALSE, TRUE,
                "AB",  TRUE,  TRUE,  FALSE,
                "AC",  TRUE,  FALSE, TRUE,
                "BC",  FALSE, TRUE,  TRUE,
                "ABC", TRUE,  TRUE,  TRUE,
                "-",   FALSE, FALSE, FALSE)
  stopifnot(all((df %>% dplyr::count(A, B, C) %>% with(n)) == 1))
  return(df %>% mutate(n = 0, text = ""))
}
gen_element_df_4 <- function() {
  df <- tribble(~name, ~A,    ~B,    ~C,    ~D,
                "A",   TRUE,  FALSE, FALSE, FALSE,
                "B",   FALSE, TRUE,  FALSE, FALSE,
                "C",   FALSE, FALSE, TRUE,  FALSE,
                "D",   FALSE, FALSE, FALSE, TRUE,
                "AB",  TRUE,  TRUE,  FALSE, FALSE,
                "BC",  FALSE, TRUE,  TRUE,  FALSE,
                "CD",  FALSE, FALSE, TRUE,  TRUE,
                "AC",  TRUE,  FALSE, TRUE,  FALSE,
                "BD",  FALSE, TRUE,  FALSE, TRUE,
                "AD",  TRUE,  FALSE, FALSE, TRUE,
                "ABC", TRUE,  TRUE,  TRUE,  FALSE,
                "BCD", FALSE, TRUE,  TRUE,  TRUE,
                "ACD", TRUE,  FALSE, TRUE,  TRUE,
                "ABD", TRUE,  TRUE,  FALSE, TRUE,
                "ABCD",TRUE,  TRUE,  TRUE,  TRUE,
                "-",   FALSE, FALSE, FALSE, FALSE)
  stopifnot(all((df %>% dplyr::count(A, B, C, D) %>% with(n)) == 1))
  return(df %>% mutate(n = 0, text = ""))
}
gen_circle <- function(group, x_offset = 0, y_offset = 0, radius = 1,
                       radius_b = radius, theta_offset = 0, length.out = 100) {
  tibble(group = group,
         theta = seq(0, 2 * pi, length.out = length.out)) %>%
    mutate(x_raw = radius * cos(theta),
           y_raw = radius_b * sin(theta),
           x = x_offset + x_raw * cos(theta_offset) - y_raw * sin(theta_offset),
           y = y_offset + x_raw * sin(theta_offset) + y_raw * cos(theta_offset))
}
calc_scale_info_2 <- function(auto_scale, n_sets, max_scale_diff = 5) {
  if (auto_scale) {
    stopifnot(length(n_sets) == 4)
    if (n_sets[[1]] == 0 && n_sets[[2]] == 0 && n_sets[[3]] == 0) { # both sets are empty
      a_radius <- 1
      b_radius <- 1
      overlap_size <- -0.2
    } else if (n_sets[[1]] + n_sets[[3]] == 0) { # set A is empty
      a_radius <- 1 / max_scale_diff
      b_radius <- 1
      overlap_size <- -0.2
    } else if (n_sets[[2]] + n_sets[[3]] == 0) { # set B is empty
      a_radius <- 1
      b_radius <- 1 / max_scale_diff
      overlap_size <- -0.2
    } else if (n_sets[[1]] >= n_sets[[2]]) { # set A is larger than or equal to set B
      a_radius <- 1
      b_radius <- (n_sets[[2]] + n_sets[[3]]) / (n_sets[[1]] + n_sets[[3]])
      overlap_size <- ifelse(n_sets[[3]] == 0, -0.2, n_sets[[3]] / (n_sets[[1]] + n_sets[[3]]))
      if (b_radius < 1 / max_scale_diff) {
        b_radius <- 1 / max_scale_diff
        if (overlap_size > 0) {
          overlap_size <- b_radius * (n_sets[[3]] / (n_sets[[2]] + n_sets[[3]]))
        }
      }
    } else { # set A is smaller than set B
      a_radius <- (n_sets[[1]] + n_sets[[3]]) / (n_sets[[2]] + n_sets[[3]])
      b_radius <- 1
      overlap_size <- ifelse(n_sets[[3]] == 0, -0.2, n_sets[[3]] / (n_sets[[2]] + n_sets[[3]]))
      if (a_radius < 1 / max_scale_diff) {
        a_radius <- 1 / max_scale_diff
        if (overlap_size > 0) {
          overlap_size <- a_radius * (n_sets[[3]] / (n_sets[[1]] + n_sets[[3]]))
        }
      }
    }
  } else {
    a_radius = 1
    b_radius = 1
    overlap_size = 1/3
  }
  return(c(auto_scale = auto_scale,
           a_radius = a_radius,
           b_radius = b_radius,
           overlap_size = overlap_size))
}
calc_scale_info_3 <- function(auto_scale, n_sets, max_scale_diff = 5) {
  if (auto_scale) {
    stop("Error: 'auto_scale' parameter is supported for only two set venn so far.")
  }
  return(NULL)
}
calc_scale_info_4 <- function(auto_scale, n_sets, max_scale_diff = 5) {
  if (auto_scale) {
    stop("Error: 'auto_scale' parameter is supported for only two set venn so far.")
  }
  return(NULL)
}

min_overlap_for_text <- 0.2

gen_circle_2 <- function(scale_info) {
  x_dist <- (scale_info['a_radius'] + scale_info['b_radius'] - scale_info['overlap_size'] * 2) / 2
  rbind(gen_circle(1L, -x_dist, 0, scale_info['a_radius']),
        gen_circle(2L, x_dist, 0, scale_info['b_radius']))
}
gen_text_pos_2 <- function(scale_info) {
  df <- tribble(~name, ~x,    ~y,  ~hjust, ~vjust,
                "A",   -0.8,  0,   0.5,    0.5,
                "B",    0.8,  0,   0.5,    0.5,
                "AB",   0,    0,   0.5,    0.5,
                "-",    0,   -1.2, 0.5,    0.5)
  if (scale_info['auto_scale']) {
    x_dist <- (scale_info['a_radius'] + scale_info['b_radius'] - scale_info['overlap_size'] * 2) / 2
    if (scale_info['overlap_size'] <= 0) {
      df$x[[1]] <- -x_dist
      df$x[[2]] <- x_dist
      df <- df %>% filter(name != "AB")
    } else {
      if (scale_info['overlap_size'] < min_overlap_for_text) {
        df$x[[1]] <- -x_dist - scale_info['overlap_size']
        df$x[[2]] <- x_dist + scale_info['overlap_size']
        if (scale_info['a_radius'] < min_overlap_for_text) {
          df$x[[3]] <- -x_dist + (scale_info['a_radius'] - scale_info['overlap_size']) / 2
          df$y[[3]] <- -1.5 * scale_info['a_radius']
        } else if (scale_info['b_radius'] < min_overlap_for_text) {
          df$x[[3]] <- x_dist - (scale_info['a_radius'] - scale_info['overlap_size']) / 2
          df$y[[3]] <- -1.5 * scale_info['b_radius']
        } else {
          df$x[[3]] <- -x_dist + scale_info['a_radius'] - scale_info['overlap_size']
          df$y[[3]] <- -1.2
        }
        df$x[[4]] <- -x_dist - scale_info['a_radius']
        df$y[[4]] <- -1.6
        df$hjust[[4]] <- 0
      } else {
        df$x[[1]] <- -x_dist - scale_info['overlap_size']
        df$x[[2]] <- x_dist + scale_info['overlap_size']
        df$x[[3]] <- -x_dist + scale_info['a_radius'] - scale_info['overlap_size']
      }
      if (scale_info['a_radius'] <= scale_info['overlap_size']) {
        df <- df %>% filter(name != "A")
      } else if (scale_info['b_radius'] <= scale_info['overlap_size']) {
        df <- df %>% filter(name != "B")
      }
    }
  }
  return(df)
}
gen_seg_pos_2 <- function(scale_info) {
  df <- tibble(x = 0, y = 0, xend = 0, yend = 0)[-1,]
  if (scale_info['overlap_size'] > 0 && scale_info['auto_scale']) {
    x_dist <- (scale_info['a_radius'] + scale_info['b_radius'] - scale_info['overlap_size'] * 2) / 2
    if (scale_info['overlap_size'] < min_overlap_for_text) {
      x_pos <- -x_dist + scale_info['a_radius'] - scale_info['overlap_size']
      if (scale_info['a_radius'] < min_overlap_for_text) {
        x2_pos <- -x_dist + 1.2 * (scale_info['a_radius'] - scale_info['overlap_size']) / 2
        df <- tibble(x = x_pos, y = 0, xend = x2_pos, yend = -1.2 * scale_info['a_radius'])
      } else if (scale_info['b_radius'] < min_overlap_for_text) {
        x2_pos <- x_dist - 1.2 * (scale_info['a_radius'] - scale_info['overlap_size']) / 2
        df <- tibble(x = x_pos, y = 0, xend = x2_pos, yend = -1.2 * scale_info['a_radius'])
      } else {
        df <- tibble(x = x_pos, y = 0, xend = x_pos, yend = -1)
      }
    }
  }
  return(df)
}
gen_label_pos_2 <- function(scale_info) {
  df <- tribble(~name, ~x,   ~y,  ~hjust, ~vjust,
                "A",   -0.8, 1.2, 0.5,    0,
                "B",    0.8, 1.2, 0.5,    0)
  if (scale_info['auto_scale']) {
  }
  return(df)
}
gen_circle_3 <- function() {
  rbind(gen_circle(1L, -2/3, (sqrt(3) + 2) / 6, 1),
        gen_circle(2L, 2/3,(sqrt(3) + 2) / 6, 1),
        gen_circle(3L, 0, -(sqrt(3) + 2) / 6, 1))
}
gen_text_pos_3 <- function() {
  tribble(~name, ~x,    ~y,   ~hjust, ~vjust,
          "A",   -0.8,  0.62, 0.5,    0.5,
          "B",    0.8,  0.62, 0.5,    0.5,
          "C",    0,   -0.62, 0.5,    0.5,
          "AB",   0,    0.8,  0.5,    0.5,
          "AC",  -0.5,  0,    0.5,    0.5,
          "BC",   0.5,  0,    0.5,    0.5,
          "ABC",  0,    0.2,  0.5,    0.5,
          "-",    1.2, -0.8,  0,      0.5)
}
gen_seg_pos_3 <- function(scale_info) {
  df <- tibble(x = 0, y = 0, xend = 0, yend = 0)[-1,]
  return(df)
}
gen_label_pos_3 <- function() {
  tribble(~name, ~x,    ~y,  ~hjust, ~vjust,
          "A",   -0.8,  1.8, 0.5,    0,
          "B",    0.8,  1.8, 0.5,    0,
          "C",    0,   -1.8, 0.5,    1)
}
gen_circle_4 <- function() {
  rbind(gen_circle(1L, -.7, -1/2, .75, 1.5, pi/4),
        gen_circle(2L, -.72+2/3, -1/6, .75, 1.5, pi/4),
        gen_circle(3L, .72-2/3, -1/6, .75, 1.5, -pi/4),
        gen_circle(4L, .7, -1/2, .75, 1.5, -pi/4))
}
gen_text_pos_4 <- function() {
  tribble(~name, ~x,    ~y,  ~hjust, ~vjust,
          "A",   -1.5,  0,   0.5,    0.5,
          "B",   -0.6,  0.7, 0.5,    0.5,
          "C",    0.6,  0.7, 0.5,    0.5,
          "D",    1.5,  0,   0.5,    0.5,
          "AB",  -0.9,  0.3, 0.5,    0.5,
          "BC",   0,    0.4, 0.5,    0.5,
          "CD",   0.9,  0.3, 0.5,    0.5,
          "AC",  -0.8, -0.9, 0.5,    0.5,
          "BD",   0.8, -0.9, 0.5,    0.5,
          "AD",   0,   -1.4, 0.5,    0.5,
          "ABC", -0.5, -0.2, 0.5,    0.5,
          "BCD",  0.5, -0.2, 0.5,    0.5,
          "ACD", -0.3, -1.1, 0.5,    0.5,
          "ABD",  0.3, -1.1, 0.5,    0.5,
          "ABCD", 0,   -0.7, 0.5,    0.5,
          "-",    0,   -1.9, 0.5,    0.5)
}
gen_seg_pos_4 <- function(scale_info) {
  df <- tibble(x = 0, y = 0, xend = 0, yend = 0)[-1,]
  return(df)
}
gen_label_pos_4 <- function() {
  tribble(~name, ~x,   ~y,   ~hjust, ~vjust,
          "A",   -1.5, -1.3, 1,      1,
          "B",   -0.8,  1.2, 0.5,    0,
          "C",    0.8,  1.2, 0.5,    0,
          "D",    1.5, -1.3, 0,      1)
}
prepare_venn_data <- function(data, columns = NULL,
                              show_elements = FALSE, show_percentage = TRUE, digits = 1,
                              label_sep = ",", count_column = NULL,
                              show_outside = c("auto", "none", "always"),
                              auto_scale = FALSE) {
  show_outside <- match.arg(show_outside)
  if (is.data.frame(data)) {
    if (is.null(columns)) {
      columns = data %>% select_if(is.logical) %>% names
    }
    if (!identical(show_elements, FALSE)) {
      if (!{
        if (is.character(show_elements)) {
          show_elements <- show_elements[[1]]
          show_elements %in% names(data)
          } else { FALSE }}) {
        stop("Value ", deparse(show_elements), 
             " in `show_elements` does not correspond to any column name of the data frame.",
             call. = FALSE)
      }
    }
    if (length(columns) == 2) {
      stopifnot(is.logical(as_tibble(data)[,columns[[1]], drop = TRUE]))
      stopifnot(is.logical(as_tibble(data)[,columns[[2]], drop = TRUE]))
      df_element <- gen_element_df_2()
      for (i in 1:nrow(df_element)) {
        idx <- ((!xor(df_element$A[[i]], as_tibble(data)[,columns[[1]]])) &
                  (!xor(df_element$B[[i]], as_tibble(data)[,columns[[2]]])))
        if (is.null(count_column)) {
          df_element$n[[i]] <- sum(idx)
        } else {
          df_element$n[[i]] <- sum(as_tibble(data)[,count_column][idx,])
        }
        if (!identical(show_elements, FALSE)) {
          df_element$text[[i]] <- paste(unlist(as_tibble(data)[idx,show_elements]), collapse = label_sep)
        }
      }
      scale_info <- calc_scale_info_2(auto_scale, df_element$n)
      df_shape <- gen_circle_2(scale_info)
      df_text <- gen_text_pos_2(scale_info) %>% inner_join(df_element, by = "name")
      df_label <- gen_label_pos_2(scale_info)
      df_seg <- gen_seg_pos_2(scale_info)
    } else if (length(columns) == 3) {
      stopifnot(is.logical(as_tibble(data)[,columns[[1]], drop = TRUE]))
      stopifnot(is.logical(as_tibble(data)[,columns[[2]], drop = TRUE]))
      stopifnot(is.logical(as_tibble(data)[,columns[[3]], drop = TRUE]))
      df_element <- gen_element_df_3()
      for (i in 1:nrow(df_element)) {
        idx <- ((!xor(df_element$A[[i]], as_tibble(data)[,columns[[1]]])) &
                  (!xor(df_element$B[[i]], as_tibble(data)[,columns[[2]]])) &
                  (!xor(df_element$C[[i]], as_tibble(data)[,columns[[3]]])))
        if (is.null(count_column)) {
          df_element$n[[i]] <- sum(idx)
        } else {
          df_element$n[[i]] <- sum(as_tibble(data)[,count_column][idx,])
        }
        if (!identical(show_elements, FALSE)) {
          df_element$text[[i]] <- paste(unlist(as_tibble(data)[idx,show_elements]), collapse = label_sep)
        }
      }
      scale_info <- calc_scale_info_3(auto_scale, df_element$n)
      df_shape <- gen_circle_3()
      df_text <- gen_text_pos_3() %>% inner_join(df_element, by = "name")
      df_label <- gen_label_pos_3()
      df_seg <- gen_seg_pos_3(scale_info)
    } else if (length(columns) == 4) {
      stopifnot(is.logical(as_tibble(data)[,columns[[1]], drop = TRUE]))
      stopifnot(is.logical(as_tibble(data)[,columns[[2]], drop = TRUE]))
      stopifnot(is.logical(as_tibble(data)[,columns[[3]], drop = TRUE]))
      stopifnot(is.logical(as_tibble(data)[,columns[[4]], drop = TRUE]))
      df_element <- gen_element_df_4()
      for (i in 1:nrow(df_element)) {
        idx <- ((df_element$A[[i]] == as_tibble(data)[,columns[[1]], drop = TRUE]) &
                  (df_element$B[[i]] == as_tibble(data)[,columns[[2]], drop = TRUE]) &
                  (df_element$C[[i]] == as_tibble(data)[,columns[[3]], drop = TRUE]) &
                  (df_element$D[[i]] == as_tibble(data)[,columns[[4]], drop = TRUE]))
        if (is.null(count_column)) {
          df_element$n[[i]] <- sum(idx)
        } else {
          df_element$n[[i]] <- sum(as_tibble(data)[,count_column][idx,])
        }
        if (!identical(show_elements, FALSE)) {
          df_element$text[[i]] <- paste(unlist(as_tibble(data)[idx,show_elements]), collapse = label_sep)
        }
      }
      scale_info <- calc_scale_info_4(auto_scale, df_element$n)
      df_shape <- gen_circle_4()
      df_text <- gen_text_pos_4() %>% inner_join(df_element, by = "name")
      df_label <- gen_label_pos_4()
      df_seg <- gen_seg_pos_4(scale_info)
    } else {
      stop("logical columns in data.frame `data` or vector `columns` should be length between 2 and 4")
    }
    df_label <- df_label %>% mutate(text = columns)
    show_elements <- !identical(show_elements, FALSE)
  } else if (is.list(data)) {
    if (is.null(columns)) {
      columns <- names(data) %>% head(4)
    }
    a2 <- na.omit(unique(unlist(data[columns])))
    if (length(columns) == 2) {
      df_element <- gen_element_df_2()
      for (i in 1:nrow(df_element)) {
        idx <- ((!xor(df_element$A[[i]], a2 %in% data[[columns[[1]]]])) &
                  (!xor(df_element$B[[i]], a2 %in% data[[columns[[2]]]])))
        df_element$n[[i]] <- sum(idx)
        df_element$text[[i]] <- paste(a2[idx], collapse = label_sep)
      }
      scale_info <- calc_scale_info_2(auto_scale, df_element$n)
      df_shape <- gen_circle_2(scale_info)
      df_text <- gen_text_pos_2(scale_info) %>% inner_join(df_element, by = "name")
      df_label <- gen_label_pos_2(scale_info)
      df_seg <- gen_seg_pos_2(scale_info)
    } else if (length(columns) == 3) {
      df_element <- gen_element_df_3()
      for (i in 1:nrow(df_element)) {
        idx <- ((!xor(df_element$A[[i]], a2 %in% data[[columns[[1]]]])) &
                  (!xor(df_element$B[[i]], a2 %in% data[[columns[[2]]]])) &
                  (!xor(df_element$C[[i]], a2 %in% data[[columns[[3]]]])))
        df_element$n[[i]] <- sum(idx)
        df_element$text[[i]] <- paste(a2[idx], collapse = label_sep)
      }
      scale_info <- calc_scale_info_3(auto_scale, df_element$n)
      df_shape <- gen_circle_3()
      df_text <- gen_text_pos_3() %>% inner_join(df_element, by = "name")
      df_label <- gen_label_pos_3()
      df_seg <- gen_seg_pos_3(scale_info)
    } else if (length(columns) == 4) {
      df_element <- gen_element_df_4()
      for (i in 1:nrow(df_element)) {
        idx <- ((!xor(df_element$A[[i]], a2 %in% data[[columns[[1]]]])) &
                  (!xor(df_element$B[[i]], a2 %in% data[[columns[[2]]]])) &
                  (!xor(df_element$C[[i]], a2 %in% data[[columns[[3]]]])) &
                  (!xor(df_element$D[[i]], a2 %in% data[[columns[[4]]]])))
        df_element$n[[i]] <- sum(idx)
        df_element$text[[i]] <- paste(a2[idx], collapse = label_sep)
      }
      scale_info <- calc_scale_info_4(auto_scale, df_element$n)
      df_shape <- gen_circle_4()
      df_text <- gen_text_pos_4() %>% inner_join(df_element, by = "name")
      df_label <- gen_label_pos_4()
      df_seg <- gen_seg_pos_4(scale_info)
    } else {
      stop("list `data` or vector `column` should be length between 2 and 4")
    }
    df_label <- df_label %>% mutate(text = columns)
  } else {
    stop("`data` should be either a list or a data.frame")
  }
  if ((show_outside == "none") || (show_outside == "auto" && df_text$n[[nrow(df_text)]] == 0)) {
    if (df_text$n[[nrow(df_text)]] > 0)
      warning("Although not display in plot, outside elements are still count in percentages.")
    df_text <- df_text[-nrow(df_text), ]
  }
  if (!show_elements) {
    fmt <- sprintf("%%d\n(%%.%df%%%%)", digits)
    if (show_percentage) {
      df_text <- df_text %>% mutate(text = sprintf(fmt, n, 100 * n / sum(n)))
    } else {
      df_text <- df_text %>% mutate(text = sprintf("%d", n))
    }
  }
  list(shapes = df_shape, texts = df_text, labels = df_label, segs = df_seg)
}
my.ggvenn <- function(data, data.plus, columns = NULL,
                   show_elements = FALSE,
                   show_percentage = TRUE,
                   digits = 1,
                   fill_color = c("blue", "yellow", "green", "red"),
                   fill_alpha = .5,
                   stroke_color = "black",
                   stroke_alpha = 1,
                   stroke_size = 1,
                   stroke_linetype = "solid",
                   set_name_color = "black",
                   set_name_size = 6,
                   text_color = "black",
###This is the color of the Y numbers underneath
                   text_color.plus = "blue",
                   text_size = 4,
###This is the vertical justification
                   vjust.plus = 1.2,
                   label_sep = ",",
                   count_column = NULL,
                   show_outside = c("auto", "none", "always"),
                   auto_scale = FALSE)
{
  show_outside <- match.arg(show_outside)

### we run this function 2 times, one for X (data) and one for Y (data.plus).  
### The objective is to perform the identical diagrams but adjust the Y one lower with vjust.

  venn <- prepare_venn_data(data, columns, show_elements, show_percentage, digits,
                            label_sep, count_column, show_outside, auto_scale)

### this was added by me; you can check out what it does but it creates the formatting

  venn.plus <- prepare_venn_data(data.plus, columns, show_elements, show_percentage, digits,
                            label_sep, count_column, show_outside, auto_scale)
  
  g <- venn$shapes %>%
    mutate(group = LETTERS[group]) %>%
    ggplot() +
    geom_polygon(aes(x = x, y = y, group = group, fill = group),
                 alpha = fill_alpha) +
    geom_polygon(aes(x = x, y = y, group = group),
                 fill = NA,
                 color = stroke_color,
                 size = stroke_size,
                 alpha = stroke_alpha,
                 linetype = stroke_linetype)
  if (nrow(venn$labels) > 0) {
    g <- g +
      geom_text(data = venn$labels,
                aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
                color = set_name_color,
                size = set_name_size)
  }
  if (nrow(venn$texts) > 0) {
    g <- g +
      geom_text(data = venn$texts,
                aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
                color = text_color,
                size = text_size)
  }
### this adjusts all of the Y values down by vjust + vjust.plus and colors them.
  if (nrow(venn.plus$texts) > 0) {
    g <- g +
      geom_text(data = venn.plus$texts,
                aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust+vjust.plus),
                color = text_color.plus,
                size = text_size)
  }
  if (nrow(venn$segs) > 0) {
    g <- g +
      geom_segment(data = venn$segs,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = text_color,
                   size = 0.5)
  }
  g <- g +
    scale_fill_manual(values = fill_color) +
    guides(fill = "none") +
    coord_fixed() +
    theme_void()
  return(g)
}

#my.ggvenn(X,Y, show_percentage = F)
#my.ggvenn(X,Y, show_percentage = F, set_name_color=c("red","blue","green","violet"), fill_color = c("red","blue","green","violet"), text_color = c("red"), text_color.plus = "blue")
