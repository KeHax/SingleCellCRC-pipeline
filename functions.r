Correctnormdata <- function(seur,
                            samples="orig.ident",
                            assay="RNA",
                            normalized.samples){
  cat("correcting counts matrix...\n")
  #subset samples whose data have been normalized
  corss_sampels <- intersect(normalized.samples, unique(seur@meta.data[[samples]]))
  if(length(corss_sampels)==0){
    stop("no sample was normalized, skipping.\n")
  }
  
  
  eval(parse(text = paste("datanrom <- subset(seur, ", samples," %in% corss_sampels)", sep = "")))
  
  #subset samples whose data have not been normalized
  counts.samples <- setdiff(unique(seur@meta.data[[samples]]), corss_sampels)
  
  eval(parse(text = paste("dataraw <- subset(seur, ", samples," %in% counts.samples)", sep = "")))
  
  if(dim(datanrom)[2] + dim(dataraw)[2] !=dim(seur)[2]){
    stop("Splitting samples falied.\n")
  }
  cat(dim(datanrom)[2], "cells will be corrected...\n")
  
  #transformation by log2p1 was altered to be logep1
  
  a=datanrom@assays[[assay]]@counts@x
  b=log(2^a, base = exp(1))
  datanrom@assays[[assay]]@data@x <- b
  
  new_seur <- merge(x=datanrom,y=dataraw)
  
  cat("correcting hvg and dr...\n")
  #add hvg and DR
  new_seur@assays[[assay]]@var.features <- seur@assays[[assay]]@var.features
  new_seur@assays[[assay]]@meta.features <- seur@assays[[assay]]@meta.features 
  
  rlist <- seur@reductions
  for (i in 1:length(rlist)) {
    rlist[[i]]@cell.embeddings <- rlist[[i]]@cell.embeddings[colnames(new_seur),]
  }
  
  new_seur@reductions <- rlist
  
  return(new_seur)
}

Filtervariablefeatures <- function(object,
                                   pattern=NULL){
  #filter non-essential variable features of lineages.
  if(!pattern %in% c("T_cell", "B_cell", "NBT")){
    cat("pattern must be one of T_cell, B_cell and NBT.\n")
  }
  
  tcr.gene <- read.csv("/data/khaoxian/project/metagene/TCR_gene.csv")
  tcr.gene <- tcr.gene$Approved.symbol
  
  if(pattern=="T_cell"){
    tcr.gene <- tcr.gene[-grep("^TRA$|^TRAC$|^TRBC1$|^TRBC2$|^TRD$|^TRDC$|^TRG$|^TRGC1$|^TRGC2$", tcr.gene)]
    #retain TRA, TRAC, TRBC1, TRBC2, TRD, TRDC, TRG, TRGC1, TRGC2 
    grepl <- "^IG[HJKL]|^RNA|^MT-|^RPS|^RPL"
  }
  if(pattern=="B_cell"){
    grepl <- "^IG[JKL]|^IGH[VDJ]|^RNA|^MT-|^RPS|^RPL"
  }
  if(pattern=="NBT"){
    grepl <- "^IG[HJKL]|^RNA|^MT-|^RPS|^RPL"
  }
  
  hvg <- VariableFeatures(object)
  hvg <- hvg[-grep(grepl,hvg)]
  hvg <- hvg[!hvg %in% tcr.gene]
  VariableFeatures(object) <- hvg
  return(object)
}

FormsurviResult <- function(DEList,
                            cohort="all",  
                            group.by=cluster[1],
                            survival.message=c("os", "dfs"),
                            group.by.quantile=0.5,
                            auto=T,
                            autofit=NULL,
                            return_ggplot=F,
                            line.size=0.4,
                            censor.size=0.5,
                            axis.title.size=5,
                            axis.line.size=0.2,
                            axis.ticks.size=0.2,
                            axis.text.size=4,
                            legend.position = "right", 
                            legend.line.size=0.4,
                            legend.title.size=3,
                            legend.text.size=3,
                            legend.key.size=2,
                            plot.title.size=5,
                            p.val.text.size=2,
                            ncol=5){
  #survival analysis
  #1. select cohort and remove normal sample
  #sorting DEList according to cohort parameter
  if(cohort!="all"){
    if(cohort %in% names(DEList)){
      DEList=DEList[cohort] 
      
      
    }else{
      stop("cohort did not indicated specific names of inputed list.\n")
    }
  }
  
  sum_gene <- purrr::map(DEList, function(x){
    x <- rownames(x[["counts"]])
    return(x)
  })
  sum_meta <- purrr::map(DEList, function(x){
    x <- colnames(x[["clinic"]])
    return(x)
  })
  
  gene_detect <- purrr::map(sum_gene, function(x){
    x <- c(group.by %in% x)
    return(x)
  }) %>% unlist()
  
  meta_detect <- purrr::map(sum_meta, function(x){
    x <- c(group.by %in% x)
    return(x)
  }) %>% unlist()
  
  if(any(gene_detect)){
    cat("detect",group.by,"in gene names.\n")
    
    #sorting DEList according to gene_detect
    DEList <- DEList[gene_detect]
    module <- "gene"
    
  }else{
    if(any(meta_detect)){
      cat("detect",group.by,"in clinic message.\n")
      
      #sorting DEList according to meta_detect
      DEList <- DEList[meta_detect]
      module <- "meta"
      
    }else{
      stop("group.by can not be found in gene names of clinic message")
    }
  }
  

  #remove normal sample in RNA-seq datatsets
  DEList <- purrr::map(DEList, function(x){
    clin <- x[["clinic"]]
    x[["clinic"]] <- dplyr::filter(x[["clinic"]], tissue == "adenocarcinoma")
    x[["counts"]] <- x[["counts"]][,x[["clinic"]]$barcode]
    return(x)
  })
  
  #2. select survival cohort
  sum_meta <- purrr::map(DEList, function(x){
    x <- colnames(x[["clinic"]])
    return(x)
  })
  survival.message=survival.message[1]
  sur.names.event <- paste(survival.message,"event",sep = "_")
  sur.names.time <- paste(survival.message,"time",sep = "_")
  survname <- c(sur.names.time, sur.names.event)
  
  meta_detect <- purrr::map(sum_meta, function(x){
    x <- c(survname %in% x)
    x <- any(x)
    return(x)
  }) %>% unlist()
  
  #sorting DEList according to survival message
  if(any(meta_detect)){
    DEList <- DEList[meta_detect]
  }else{
    stop("please check survival message.\n")
  }
  
  #3. group by 
  #firstly filter 0 and na
  #filter time==0 and event==NA
  DEList.group <- purrr::map(DEList,function(x){
    
    a <- !is.na(x[["clinic"]][[survname[2]]])
    x[["counts"]] <- x[["counts"]][,a]
    x[["clinic"]] <- x[["clinic"]][a,]
    
    eval(parse(text = paste("x[['clinic']] <- filter(x[['clinic']], ",survname[1], " > 0)", sep = "")))
    x[["counts"]] <- x[["counts"]][,x[["clinic"]]$barcode]
    
    return(x)
  })
  
  if(module == "gene"){
    if(!auto){
      group_key <- paste(group.by,"_",sep = "")
      DEList.group <- purrr::map(DEList.group,function(x){
        a=x[["counts"]][group.by,]  %>% as.numeric()
        thres <- quantile(a,group.by.quantile)
        groupx <- ifelse(a > thres, "high", "low")
        x[["clinic"]][[group_key]] <- groupx
        return(x)
      })
    }else{
      cat("auto-detecting threshold.\n")
      if(is.null(autofit)){
        check_list <- list()
        for (i in 20:80) {
          check_list[[i-19]] <- c(i/100, i/100) 
        }
        check_list[[62]] <- c(0.25,0.75)
        autofit <- check_list
      }
      

      poolss <- c(paste(letters,letters, sep = ""),
                  paste(letters,letters,letters, sep = ""),
                  paste(letters,letters,letters,letters, sep = ""),
                  paste(letters,letters,letters,letters,letters, sep = ""))
      autofit_name <- poolss[1:length(autofit)]
      names(autofit) <- autofit_name
      DEList.group_pval <- DEList.group
      for (i in 1:length(autofit)) {
        cutof_low <- autofit[[i]][1]
        cutof_hight <- autofit[[i]][2]
        
        DEList.group_pval <- purrr::map(DEList.group_pval,function(x){
          a=x[["counts"]][group.by,]  %>% as.numeric()
          thres_low <- quantile(a,cutof_low)
          thres_high <- quantile(a,cutof_hight)
          groupx <- ifelse(a < thres_low, "low", 
                           ifelse(a >= thres_high, "high", "removed"))
          #filter only one group because threshold could not group them
          if(length(setdiff(unique(groupx), "removed")) != 1){
            x[["clinic"]][[autofit_name[i]]] <- groupx
          }
          return(x)
        })
      }
      
      #globally filter DEList.group_pval, if it did not have autofit_name
      check_colnames <- purrr::map(DEList.group_pval, list("clinic"))
      check_colnames <- purrr::map(check_colnames, function(x){
        x <- colnames(x)
        if(any(x %in% autofit_name)){
          x <- TRUE
        }else{
          x <- FALSE
        }
        return(x)
      })
      
      check_colnames <- unlist(check_colnames)
      
      #cat("threshold1:",check_colnames)
      if(all(check_colnames==F)){
        stop("threshold1: no threshold can make groups.\n")
      }
      
      DEList.group_pval <- DEList.group_pval[check_colnames]
      DEList.group <- DEList.group[check_colnames]
      
      #check p_val
      p_val_list <- list()
      pass_fitname.list <- list()
      for (i in 1:length(DEList.group_pval)) {
        metadata=DEList.group_pval[[i]][["clinic"]]
        pass_fitname <- colnames(metadata)[colnames(metadata) %in% autofit_name]
        pass_fitname.list[[i]] <- pass_fitname
        p_val <- list()
        for (id in 1:length(pass_fitname)) {
          eval(parse(text = paste("metadata_filtered <- filter(metadata, ",
                                  pass_fitname[id], " != 'removed')", sep = "")))
          
          eval(parse(text = paste("fit <- survfit(Surv(", survname[1], ", ", survname[2], ") ~ ", pass_fitname[id], ", data = metadata_filtered)", sep = "")))
          
          p_val[[id]] <- surv_pvalue(fit, data = metadata_filtered)$pval
          
        }
        p_val <- unlist(p_val)
        p_val_list[[i]] <- p_val
        
      }
      names(p_val_list) <- names(DEList.group_pval)
      names(pass_fitname.list) <- names(DEList.group_pval)
      
      #select minimal p value 
      p_val_list_select <- list()
      for (i in 1:length(p_val_list)) {
        if(length(which.min(p_val_list[[i]])) != 0){
          p_val_list_select[[i]] <- pass_fitname.list[[i]][which.min(p_val_list[[i]])]
        }else{
          p_val_list_select[[i]] < -NA
        }
      }
      p_val_list_select <- as.character(unlist(p_val_list_select))
      
      #filter cohort if all p_val is NA
      if(all(is.na(p_val_list_select))){
        stop("threshold2: no threshold can return p values.\n")
      }
      
      DEList.group <- DEList.group[!is.na(p_val_list_select)]
      DEList.group_pval <- DEList.group_pval[!is.na(p_val_list_select)]
      p_val_list_select <- p_val_list_select[!is.na(p_val_list_select)]
      
      group_key <- paste(group.by,"_",sep = "")
      for (i in 1:length(p_val_list_select)) {
        DEList.group[[i]]$clinic[[group_key]] <- DEList.group_pval[[i]]$clinic[[p_val_list_select[i]]]
        a <- DEList.group[[i]]$clinic
        eval(parse(text = paste("a <- dplyr::filter(a, ", group_key,  "!= 'removed')", sep = "")))
        DEList.group[[i]]$clinic <- a
      }
    }
    
  }else{
    group_key <- group.by
    
  }
  
  #4. plot
  plot_theme <- theme_classic()+
    theme(aspect.ratio = 0.68,
          axis.title = element_text(size = axis.title.size),
          axis.line = element_line(size = axis.line.size),
          axis.ticks = element_line(size=axis.ticks.size),
          axis.text = element_text(size=axis.text.size),
          legend.text = element_text(colour = "black",size = legend.text.size),
          legend.position = legend.position,
          legend.key.size = unit(legend.line.size, "lines"),
          legend.title = element_text(size = legend.title.size),
          plot.title = element_text(size=plot.title.size, vjust = 0.5))
  
  plot.list <- purrr::map(names(DEList.group),function(x){
    metadata=DEList.group[[x]][["clinic"]]
    n_heigh <- sum(metadata[[group_key]] == "high")
    n_low <- sum(metadata[[group_key]] == "low")
    eval(parse(text = paste("fit <- survfit(Surv(", survname[1], ", ", survname[2], ") ~ ", group_key, ", data = metadata)", sep = "")))
    ggsurvplot(fit, data = metadata, 
               pval = TRUE, conf.int = F, pval.method = F,#palette = c("#E7B800", "#2E9FDF"),
               size=line.size,
               censor.shape="|", 
               censor.size = censor.size,
               legend.labs = c(paste(group.by," high (n = ", n_heigh, ")", sep = ""), 
                               paste(group.by," low (n = ", n_low, ")", sep = "")),
               ggtheme = plot_theme) + 
      ggtitle(x)
    
  })
  names(plot.list) <- names(DEList.group)
  plot.list2 <- purrr::map(plot.list, list("plot"))
  plot.list2 <- purrr::map(plot.list2,function(x){
    x[["layers"]][[4]][["aes_params"]][["size"]] <- p.val.text.size
    return(x)
  })
  
  if(return_ggplot){
    return(plot.list)
  }else{
    p <- cowplot::plot_grid(plotlist = plot.list2, ncol = ncol)
    return(p)
  }
}

Addsignature2DEList <- function(DEList,
                                signature.list){
  DEList_added <- DEList
  for (i in 1:length(signature.list)) {
    genes <- signature.list[[i]]
    genes_name <- names(signature.list)[i]
    DEList_added <- purrr::map(DEList_added, function(x){
      genesframe <- x[["counts"]][intersect(rownames(x[["counts"]]), genes),,drop=F]
      if(nrow(genesframe)==0){
        x <- NULL
        return(x)
      }else{
        genesmean <- colSums(genesframe)/nrow(genesframe)
        x[["counts"]][genes_name, ] <- genesmean
        return(x)
      }
    })
  }
  return(DEList_added)
}

CytoTRACEz <- function(mat, 
                       enableFast = TRUE, 
                       ncores = 1, 
                       subsamplesize = 1000){
  cat("input nromalized matrix which has been batch-corrected.\n")  #like mnn output
  cat("this function was formated from CytoTRACE.\n")
  
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  a1 <- mat
  a2 <- NULL #batch=NULL
  batch=NULL
  if (ncol(mat) < 3000) {
    enableFast = FALSE
    message("The number of cells in your dataset is less than 3,000. Fast mode has been disabled.")
  }
  else {
    message("The number of cells in your dataset exceeds 3,000. CytoTRACE will now be run in fast mode (see documentation). You can multi-thread this run using the 'ncores' flag. To disable fast mode, please indicate 'enableFast = FALSE'.")
  }
  
  pqgenes <- is.na(rowSums(mat > 0)) | apply(mat, 1, var) == 0
  num_pqgenes <- length(which(pqgenes == TRUE))
  mat <- mat[!pqgenes, ]
  if (num_pqgenes > 0) {
    warning(paste(num_pqgenes, "genes have zero expression in the matrix and were filtered"))
  }
  
  
  if (enableFast == FALSE) {
    size <- ncol(mat)
  }
  else if (enableFast == TRUE & subsamplesize < ncol(mat)) {
    size <- subsamplesize
  }
  else if (enableFast == TRUE & subsamplesize >= ncol(mat)) {
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }
  
  #sampling
  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat)%%chunk)))
  
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately", 
                round(mean(unlist(lapply(subsamples, length)))), "cells each using", 
                min(chunk, ncores), "/", ncores, "core(s)"))
  message(paste("Pre-processing data and generating similarity matrix..."))
  
  batches <- parallel::mclapply(subsamples,
                                mc.cores = min(chunk,ncores), 
                                function(subsample) {
                                  mat <- mat[, subsample]
                                  batch <- batch[subsample]
                                  
                                  #skpping cell QC: na or nFeature_RNA < 10
                                  #do not normalize data again, because input was normalized data
                                  counts <- apply(mat > 0, 2, sum)
                                  mat2 <- mat
                                  
                                  #select 1000 hvg
                                  mvg <- function(matn) {
                                    A <- matn
                                    n_expr <- rowSums(A > 0) 
                                    
                                    #filter gene express in less than 5% cells
                                    A_filt <- A[n_expr >= 0.05 * ncol(A), ]
                                    
                                    #calculated disp delta/miu, retain only 1000 genes
                                    vars <- apply(A_filt, 1, var)
                                    means <- apply(A_filt, 1, mean)
                                    disp <- vars/means
                                    last_disp <- tail(sort(disp), 1000)[1]
                                    A_filt <- A_filt[disp >= last_disp, ]
                                    
                                    return(A_filt)
                                  }
                                  
                                  mat2.mvg <- mvg(mat2)
                                  rm1 <- colSums(mat2.mvg) == 0
                                  mat2 <- mat2[, !rm1]
                                  counts <- counts[!rm1]
                                  
                                  similarity_matrix_cleaned <- function(similarity_matrix) {
                                    D <- similarity_matrix
                                    cutoff <- mean(as.vector(D))
                                    diag(D) <- 0
                                    D[which(D < 0)] <- 0
                                    D[which(D <= cutoff)] <- 0
                                    Ds <- D
                                    D <- D/rowSums(D)
                                    D[which(rowSums(Ds) == 0), ] <- 0
                                    return(D)
                                  }
                                  
                                  D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))
                                  return(list(mat2 = mat2, counts = counts, D = D))
                                })
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  
  filter <- colnames(a1)[-which(colnames(a1) %in% colnames(mat2))]
  
  if (length(filter) > 0) {
    warning(paste(length(filter), "poor quality cells were filtered based on low or no expression. See 'filteredCells' in returned object for names of filtered cells."))
  }
  
  message("Calculating gene counts signature...")
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x, 
  ], counts))
  names(ds2) <- rownames(mat2)
  
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])), 
  ], 2, mean)
  
  samplesize <- unlist(lapply(lapply(batches, function(x) x$counts), 
                              length))
  
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))
  
  D2 <- lapply(batches, function(x) x$D)
  
  regressed <- function(similarity_matrix_cleaned, score) {
    out <- nnls::nnls(similarity_matrix_cleaned, score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }
  
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9) {
    vals <- score
    v_prev <- rep(vals)
    v_curr <- rep(vals)
    for (i in 1:10000) {
      v_prev <- rep(v_curr)
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% 
                           v_curr) + (1 - ALPHA) * vals
      diff <- mean(abs(v_curr - v_prev))
      if (diff <= 1e-06) {
        break
      }
    }
    return(v_curr)
  }
  
  message("Smoothing values with NNLS regression and diffusion...")
  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, 
                                  function(i) {
                                    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
                                    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
                                    cytotrace <- rank(gcs_diffused)
                                  })
  
  #calculate cytotrace score
  cytotrace <- cytotrace_ranked <- unlist(cytotrace)
  cytotrace <- range01(cytotrace)
  
  cytogenes <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x, ], cytotrace))
  names(cytogenes) <- rownames(mat2)
  
  message("Calculating genes associated with CytoTRACE...")
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2)
  
  cytotrace <- cytotrace[colnames(a1)]
  cytotrace_ranked <- cytotrace_ranked[colnames(a1)]
  gcs <- gcs[colnames(a1)]
  counts <- counts[colnames(a1)]
  mat2 <- t(data.frame(t(mat2))[colnames(a1), ])
  
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2) <- colnames(a1)
  
  message("Done")
  return(list(CytoTRACE = cytotrace, CytoTRACErank = cytotrace_ranked, 
              cytoGenes = sort(cytogenes, decreasing = T), GCS = gcs, 
              gcsGenes = sort(ds2, decreasing = T), Counts = counts, 
              filteredCells = filter, exprMatrix = mat2))
}

seurat2monocle3 <- function(seur, 
                            CID=NULL,
                            assay="RNA",
                            slot="data"){
  if(is.null(CID)){CID=colnames(seur)}
  seur <- seur[,match(CID, colnames(seur))]
  data <- GetAssayData(seur, assay = assay, slot = slot)
  cell_metadata <- seur@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  return(cds)
}

compute_gene_trends <- function(seur.palantir,
                                gene_exprs=NULL,
                                branch,
                                markers,
                                branch.probability=0.7){
  if(is.null(gene_exprs)){
    gene_exprs <- seur.palantir$X
  }
  #raw pseudo-time: x
  pseudotimes <- seur.palantir$obs$pseudotime
  branch.prob <- list()
  
  #raw weights: weights
  for (i in 1:length(branch)) {
    branch.prob[[branch[i]]] <- seur.palantir$obs[[branch[i]]]
  }
  
  #raw gene_exprs: y
  gene_exprs <- gene_exprs[,markers,drop=F]
  gene_exprs <- as.matrix(gene_exprs)
  
  #bins: 500
  bins.list <- list()
  for (i in 1:length(branch)) {
    bins.list[[branch[i]]] <- seq(0,max(pseudotimes[branch.prob[[i]] > branch.probability]),length.out=500)
  }
  
  #making models for each gene for each branch
  r_data <- as.data.frame(gene_exprs)
  r_data$x=pseudotimes
  
  branch_bodel <- list()
  for (id in 1:length(branch)) {
    branch_bodel[[id]] <- purrr::map(markers, function(xmarker){
      r_model_data <- r_data[,c("x",xmarker)]
      colnames(r_model_data) <- c("x", "y")
      models <- gam::gam(formula = y~s(x),
                         data = r_model_data[branch.prob[[id]] > branch.probability,],
                         weights = branch.prob[[id]][branch.prob[[id]] > branch.probability])
      xmarker <- models
      return(xmarker)
    })
    names(branch_bodel[[id]]) <- markers
    
  }
  names(branch_bodel) <- branch
  
  #predict
  #bins.list x predict.y + std
  predict.y.list <- list()
  for (i in 1:length(branch)) {
    predict.y.gene.list <- list()
    for (id in 1:length(markers)) {
      predict.y.gene.list[[id]] <- predict(branch_bodel[[branch[i]]][[markers[id]]], 
                                           data.frame(x=bins.list[[branch[i]]]))
      
    }
    names(predict.y.gene.list) <- markers
    predict.y.list[[i]] <- predict.y.gene.list
    
  }
  names(predict.y.list) <- branch
  
  #predict all y by all x
  p.list <- list()
  for (i in 1:length(branch)) {
    p.gene.list <- list()
    for (id in 1:length(markers)) {
      p.gene.list[[id]] <- predict(branch_bodel[[branch[i]]][[markers[id]]], 
                                   data.frame(x=r_data$x[branch.prob[[branch[i]]] > branch.probability]))
    }
    names(p.gene.list) <- markers
    p.list[[i]] <- p.gene.list
  }
  names(p.list) <- branch
  
  #std
  std.list <- list()
  for (i in 1:length(branch)) {
    std.gene.list <- list()
    for (id in 1:length(markers)) {
      
      use_inds=c(1:length(branch.prob[[branch[i]]]))[branch.prob[[branch[i]]] > branch.probability]
      
      n = length(use_inds)
      
      sigma = sqrt(sum((gene_exprs[,markers[id]][use_inds] - p.list[[branch[i]]][[markers[id]]]) ** 2) / (n - 2))
      
      stds = (
        sqrt(1 + 1 / n + (bins.list[[branch[i]]] - mean(r_data$x)) ** 2 /sum((r_data$x - mean(r_data$x)) ** 2) )* sigma / 2
      )
      
      std.gene.list[[id]] <- stds
      
    }
    names(std.gene.list) <- markers
    std.list[[branch[i]]] <- std.gene.list
  }
  names(std.list) <- branch
  
  newdata.list <- list()
  for (i in 1:length(branch)) {
    newdata.gene.list <- list()
    for (id in 1:length(markers)) {
      newdatax <- data.frame(
        x=bins.list[[branch[i]]],
        y=as.numeric(predict.y.list[[branch[i]]][[markers[id]]]),
        std=std.list[[branch[i]]][[markers[id]]],
        branch=rep(branch[i],length(bins.list[[branch[i]]]))
      )
      newdata.gene.list[[id]] <- newdatax
    }
    names(newdata.gene.list) <- markers
    newdata.list[[i]] <- newdata.gene.list
  }
  names(newdata.list) <- branch
  
  markers.data.list <- list()
  for (id in 1:length(markers)) {
    markers.data.list[[id]] <- purrr::map(newdata.list, function(x){
      x <- x[[markers[id]]]
    })
    markers.data.list[[id]] <- bind_rows(markers.data.list[[id]])
    
    
  }
  names(markers.data.list) <- markers 
  return(markers.data.list)
  
}

plot_gene_trends <- function(gene_trends,
                             xlab="pseudotime",
                             axis.title.size=6,
                             axis.line.size=0.25,
                             axis.ticks.size=0.25,
                             axis.text.size=5.2,
                             legend.text.size=7,
                             legend.key.size=3,
                             plot.title.size=7,
                             legend.line.size=0.9,
                             legend.position = "top",
                             aspect.ratio = 0.65,
                             ncol=5){
  
  plot_theme <- theme(aspect.ratio = aspect.ratio,
                      axis.title = element_text(size = axis.title.size),
                      axis.line = element_line(size = axis.line.size),
                      axis.ticks = element_line(size=axis.ticks.size),
                      axis.text = element_text(size=axis.text.size),
                      legend.text = element_text(colour = "black",size = legend.text.size),
                      legend.position = legend.position,
                      legend.key.size = unit(legend.line.size, "lines"),
                      legend.title = element_blank(),
                      plot.title = element_text(size=plot.title.size))
  
  p <- purrr::map(gene_trends, function(xplot){
    if(length(unique(xplot$branch)) != 1){
      ggplot(xplot, aes(x = x, y = y)) +
        geom_ribbon(aes(ymin = y-std, ymax = y+std, fill=branch), alpha=0.1) +
        geom_line(aes(color=branch))+
        scale_color_manual(values = RColorBrewer::brewer.pal(8,"Set2"))+
        scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set2"))+
        theme_classic()+
        plot_theme
    }else{
      ggplot(xplot, aes(x = x, y = y)) +
        geom_ribbon(aes(ymin = y-std, ymax = y+std, fill=branch), alpha=0.1) +
        geom_line(aes(color=branch))+
        scale_color_manual(values = RColorBrewer::brewer.pal(8,"Set2"))+
        scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set2"))+
        theme_classic()+
        plot_theme+
        theme(legend.position = "none")
    }
  })
  
  p <- purrr::map(names(p), function(xplot){
    p[[xplot]] <- p[[xplot]]+
      labs(x=xlab, y=xplot)
    xplot <- p[[xplot]] 
    xplot <- xplot+ theme(plot.margin = unit(c(0,0,0,0), "cm"))
    return(xplot)
  })
  p <- cowplot::plot_grid(plotlist = p,
                          ncol = ncol)
  return(p)
}

plot_gene_trend_heatmaps <- function(gene_trends, 
                                     genes=NULL,
                                     lineages=NULL,
                                     Annotategenecluster=NULL,
                                     return_data=F,
                                     show_row_names=F,
                                     cluster_rows=T,
                                     right_annotation=NULL,
                                     clust_by_Annotategenecluster=F,
                                     clust_by_Annotategenecluster_rankded_clust=F,
                                     clust_by_Annotategenecluster_rankded_clust.algorithmn=c("weight", "hclust"), 
                                     clust_by_Annotategenecluster_rankded_clust.algorithmn.scale=c("exp","ind2"),
                                     filter_NA=F,
                                     fontsize=5,
                                     link_width.size=5,
                                     link_height.size=5,
                                     padding.size=1){
  if(is.null(lineages)){
    lineages=as.character(unique(gene_trends[[1]]$branch))
  }
  cat("The number of lineages is", length(lineages),".\n")

  if(!is.null(genes)){
    gene_trends <- gene_trends[genes]
  }
  
  lineages.list <- list()
  for (i in 1:length(lineages)) {
    lineages.list[[i]] <- purrr::map(gene_trends,function(xtrend){
      xtrend <- filter(xtrend, branch == lineages[i])
      xtrend <- select(xtrend,x,y)
      rownames(xtrend) <- xtrend$x
      xtrend <- xtrend[,"y", drop=F]
      return(xtrend)
    })
  }
  names(lineages.list) <- lineages
  markers=names(gene_trends)
  for (i in 1:length(lineages)) {
    for (id in 1:length(markers)) {
      colnames(lineages.list[[lineages[i]]][[markers[id]]]) <- markers[id]
    }
    
  }
  lineages.list <- purrr::map(lineages.list,function(x){
    x <- bind_cols(x)

    return(x)
  })
  lineages.list.scaled <- purrr::map(lineages.list,function(x){
    x_train <- scale(x)
    x <-  scale(x, center=attr(x_train, "scaled:center"), 
                scale=attr(x_train, "scaled:scale"))
    return(x)
  })
  #return_data
  if(return_data){
    for (i in 1:length(lineages.list.scaled)) {
      lineages.list.scaled[[i]] <- t(as.data.frame(lineages.list.scaled[[i]]))
    }
    return(lineages.list.scaled)
  }
  
  myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  colorRampin=circlize::colorRamp2(seq(-3,3,length.out=100),
                                   myPalette(100))
  heatmap.list <- list()
  
  for (i in 1:length(lineages)) {
    
    #Heatmap
    #set data
    mat <- t(as.matrix(lineages.list.scaled[[i]]))
    
    #filter NA
    if(filter_NA==T){
      mat <- mat[!is.nan(apply(mat,1,sum)),]
    }
    
    #set Annotategenecluster data
    if(!is.null(Annotategenecluster)){
      
      #check Annotategenecluster
      if(is.list(Annotategenecluster) &!is.data.frame(Annotategenecluster)){
        if(length(Annotategenecluster)!=length(heatmap.list)){
          stop("the length of Annotategenecluster is not equal to branch.\n")
        }
        annotaterowdata <- Annotategenecluster[i]
      }else{
        annotaterowdata <- Annotategenecluster
      }
      
      if(!is.null(genes)){
        annotaterowdata <- annotaterowdata[annotaterowdata$gene %in% genes,]
      }
      
      if(nrow(annotaterowdata) != nrow(mat)){
        cat("the row of Annotategenecluster is not equal to gene_trends.\n")
        cat("forcing the row of Annotategenecluster equal to gene_trends.\n")
        annotaterowdata <- annotaterowdata[annotaterowdata$gene %in% rownames(mat), ]
      } 
    }else{
      rowAnno=NULL
    }
    
    if(!is.null(Annotategenecluster)){
      if(!clust_by_Annotategenecluster){
        #match annotaterowdata
        annotaterowdata <- annotaterowdata[match(rownames(mat), annotaterowdata$gene),]
        
        name_color <- scCustomize::scCustomize_Palette(num_groups = 36, 
                                                       ggplot_default_colors = F, 
                                                       color_seed = 123)[1:length(unique(annotaterowdata$cluster))]
        
        if(!is.null(levels(annotaterowdata$cluster))){
          names(name_color) <- levels(annotaterowdata$cluster)
        }else{
          annotaterowdata$cluster <- factor(annotaterowdata$cluster)
          names(name_color) <- levels(annotaterowdata$cluster)
        }
        
        rowAnno <- ComplexHeatmap::rowAnnotation(GeneCluster=annotaterowdata$cluster,
                                                 col = list(GeneCluster=name_color))
        rowAnno@anno_list[["GeneCluster"]]@label <- ""
      }else{
        #when clust_by_Annotategenecluster, cluster_rows was set to F
        cluster_rows=F
        
        #ranking annotaterowdata by factor levels firstly.
        if(!is.null(levels(annotaterowdata$cluster))){
          annotaterowdata.levels <- levels(annotaterowdata$cluster)
        }else{
          annotaterowdata$cluster <- factor(annotaterowdata$cluster)
          annotaterowdata.levels <- levels(annotaterowdata$cluster)
          #cat("the length of annotaterowdata.levels is", length(annotaterowdata.levels), "\n")
        }
        
        annotaterowdata.list <- split(annotaterowdata,f=annotaterowdata$cluster)
        annotaterowdata.list <- annotaterowdata.list[annotaterowdata.levels]
        
        if(clust_by_Annotategenecluster_rankded_clust){
          clust_by_Annotategenecluster_rankded_clust.algorithmn=clust_by_Annotategenecluster_rankded_clust.algorithmn[1]
          
          #hclust
          if(clust_by_Annotategenecluster_rankded_clust.algorithmn=="hclust"){
            #in each gene cluster, ranked by hclust again secondly.
            heatmap_ranked <- ComplexHeatmap::Heatmap(mat)
            heatmap_ranked <- rownames(mat)[ComplexHeatmap::row_order(heatmap_ranked)]
            annotaterowdata.list <- purrr::map(annotaterowdata.list, function(x){
              index <- match(heatmap_ranked,x$gene)
              index <- index[!is.na(index)]
              x <- x[index,]
              return(x)
            })
          }
          
          #weight
          if(clust_by_Annotategenecluster_rankded_clust.algorithmn=="weight"){
            #return(annotaterowdata.list)
            annotaterowdata.list <- purrr::map(annotaterowdata.list, function(x){
              matx=as.matrix(mat)
              matx <- mat[x$gene, ,drop=F]
              matx=t(apply(matx,1,function(x){as.numeric(x)*(as.numeric(colnames(matx)))^2}))  #^2 enhance larger pseudotime
              matx <- apply(matx,1,sum)
              matx<- match(sort(matx),matx)
              x <- x[matx,]
              return(x)
            })
            
            
          }
        }
        
        #end ranking
        annotaterowdata <- bind_rows(annotaterowdata.list)
        
        mat <- mat[match(annotaterowdata$gene, rownames(mat)),]
        #
        if(nrow(annotaterowdata) != nrow(mat)){
          stop("the row of Annotategenecluster is not equal to that of gene_trends after clust_by_Annotategenecluster_rankded_clust.\n")
        }
        
        #add Rowannotation
        name_color <- scCustomize::scCustomize_Palette(num_groups = 36, 
                                                       ggplot_default_colors = F, 
                                                       color_seed = 123)[1:length(unique(annotaterowdata$cluster))]
        
        names(name_color) <- annotaterowdata.levels
        rowAnno <- ComplexHeatmap::rowAnnotation(GeneCluster=annotaterowdata$cluster,
                                                 col = list(GeneCluster=name_color))
        rowAnno@anno_list[["GeneCluster"]]@label <- ""
        
      }
    }
    
    if(clust_by_Annotategenecluster==F & clust_by_Annotategenecluster_rankded_clust){
      clust_by_Annotategenecluster_rankded_clust.algorithmn=clust_by_Annotategenecluster_rankded_clust.algorithmn[1]
      if(clust_by_Annotategenecluster_rankded_clust.algorithmn=="hclust"){
        cluster_rows=T
      }
      
      if(clust_by_Annotategenecluster_rankded_clust.algorithmn=="weight"){
        matx=as.matrix(mat)
        if(clust_by_Annotategenecluster_rankded_clust.algorithmn.scale=="exp"){
          matx=t(apply(matx,1,function(x){as.numeric(x)*exp((as.numeric(colnames(matx))))}))
        }else{
          indes <- gsub("ind","",clust_by_Annotategenecluster_rankded_clust.algorithmn.scale) %>% as.numeric()
          matx=t(apply(matx,1,function(x){as.numeric(x)*(as.numeric(colnames(matx)))^indes}))  #^2 enhance larger pseudotime
        }
        matx <- apply(matx,1,sum)
        matx<- match(sort(matx),matx)
        mat <- mat[matx,]
        cluster_rows=F
      }
      
    }
    
    #set right_annotation
    if(!is.null(right_annotation)){
      right_annotation_select <- intersect(right_annotation, rownames(mat))
      loci <- match(right_annotation_select,rownames(mat))
      
      geneMark = ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = loci, 
                                                                                labels = right_annotation_select, 
                                                                                labels_gp = grid::gpar(fontface = "italic", fontsize = fontsize),
                                                                                link_width = unit(link_width.size, "mm"),  #left_right
                                                                                padding = unit(padding.size, "mm"), #space between texts
                                                                                link_height  = unit(link_height.size, "mm"))) #up_down
      show_row_names=F
    }else{
      geneMark=NULL
    }
    
    #plot
    heatmap.list[[i]] <- ComplexHeatmap::Heatmap(mat, 
                                                 name = "Z-score",
                                                 left_annotation =rowAnno,
                                                 col = colorRampin, 
                                                 cluster_columns = F, 
                                                 cluster_rows = cluster_rows,
                                                 show_column_names = F,
                                                 show_row_names = show_row_names,
                                                 right_annotation = geneMark)
  }
  names(heatmap.list) <- lineages
  return(heatmap.list)
}

Samplingseurat <- function(seur, 
                           clustering_to_use="seurat_clusters", 
                           sample_number=500,
                           IdentedCell=T){
  metas <- seur@meta.data 
  metas <- split(metas, f=metas[[clustering_to_use]])
  metas <- purrr::map(metas, function(x){
    x <- rownames(x)[sample(c(1:nrow(x)), ifelse(nrow(x) > (sample_number-1), sample_number, nrow(x)))]
    return(x)
  })
  metas <- unlist(metas)
  seur <- seur[,metas]
  if(IdentedCell){
    Idents(seur) <- clustering_to_use 
  }
  
  return(seur)
}



