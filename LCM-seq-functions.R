# load libraries 
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(plyr)
  library(ggpubr)
  library(limma)
  library(edgeR)
  library(WGCNA)
  library(variancePartition)
  library(DESeq2)
  library(parallel)
  library(BiocParallel)
  library(readxl)
  library(genefilter)
  library(FactoMineR)
  library(umap)
  library(Rtsne)
  library(gridExtra)
  library(cowplot)
  library(RColorBrewer)
  library(dendextend)
  library(biomaRt)
})

# for parallel processing
ncores = detectCores()
param = SnowParam(ncores, "SOCK", progressbar = TRUE)

#' Gene Filtering
#'
#' Source: https://seqqc.wordpress.com/2020/02/17/removing-low-count-genes-for-rna-seq-downstream-analysis/
#' Keep genes that have at least min.count reads in a minimum proportion of samples given by N,
#' evaluates to ‘TRUE’ if >=25% of the samples have count-per-million (CPM) above k, 
#' where k is determined by the default value of min.count=10 and by the sample library
#' 
#' @param counts count matrix 
#' @param min.count minimum number of counts 
#' @param N minimum proportion of samples 
#'
#' @return Genes to keep 
selectGenes = function(counts, min.count = 10, N = 0.25){
  
  lib.size = colSums(counts)
  MedianLibSize = median(lib.size)
  CPM.Cutoff = min.count / MedianLibSize*1e6
  CPM = edgeR::cpm(counts, lib.size=lib.size)
  
  min.samples = round(N * ncol(counts))
  
  keep = apply(CPM, 1, function(x, n = min.samples){
    t = sum(x >= CPM.Cutoff) >= n
    t
  })
  print(paste("Equivalent CPM cutoff:", CPM.Cutoff))
  print(paste("Minimum samples:", min.samples))
  return(keep)
}

#' Plot mean vs. standard-deviation relationship
#'
#' @param data data matrix where columns represent samples to be clustered
#' 
#' @return ggplot 
plotMeanSD = function(data, stabilized){ 

  p_mean_sd =
    p_df = data %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
    ungroup()
    if(stabilized == TRUE) {
      p = ggplot(p_df, aes(x = gene_average, y = gene_stdev)) 
    } else {
      p =  ggplot(p_df, aes(x = log10(gene_average), y = log10(gene_stdev))) 
    }
    p = p + geom_point(alpha = 0.5, fill = "grey", colour = "black") +
    labs(x = "Gene count average", 
         y = "Gene count standard deviation") +
    ggtitle("Mean - Standard deviation relationship")
  
  return(p)
}

#' Hierarchical clustering for removing outlier samples 
#'
#' @param data data matrix where columns represent samples to be clustered
#' @param metadata metadata where the row names match the columns of `data`
#' @param cluster.number the expected number of clusters, used for `cutree`
#' @param method the agglomeration method to be used, one of "ward.D", "ward.D2", "single", "complete", 
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#'
#' @return dendrogram with leafs colored  by DX, and dataframe of cluster membership 
plotClusters = function(data, metadata, cluster.number, method){  
  
  if(ncol(data) < nrow(data)) data = t(data)
  groups = metadata$DX
  color_codes = c(CTRL="red", MDD="green", BD="blue", SCZ="orange")
  
  dist.df = data %>% na.omit %>% dist(method = "euclidean")
  
  hc = hclust(dist.df, method = method)
  dend = as.dendrogram(hc)
  labels_colors(dend) = color_codes[groups][order.dendrogram(dend)]
  labels_cex(dend) = 0.3
  
  par(cex = 1)
  plot(dend, horiz = T)
  
  split.tree = cutree(dend, k = cluster.number)
  unname(split.tree)

  return(split.tree)
}

#' Cluster samples using method of choice
#'
#' @param data data matrix where columns represent samples to be clustered  
#' @param metadata metadata where the row names match the columns of `data`
#' @param variable name of the metadata variable to identify samples by
#' @param method one of `UMAP`, `PCA`, or `TSNE`
#' @param ellipses should ellipses be drawn around samples of the same identity `TRUE` or `FALSE`
#' @param legend should the legend be displayed `TRUE` or `FALSE`
#'
#' @return ggplot of clusters 
plotDims = function(data, metadata, variable, method = "UMAP", ellipses = FALSE, legend = TRUE){
  
  if(ncol(data) < nrow(data)) data = t(data)
  sample_ids = rownames(data)
  
  if(method == "UMAP"){
    res = umap(data)
    res = data.frame(ID = sample_ids, X = res$layout[,1], Y = res$layout[,2])
  } 
  if(method == "PCA"){
    res = PCA(data, scale.unit = T, graph = F)
    res = data.frame(ID = sample_ids, X = res$ind$coord[,1], Y = res$ind$coord[,2])
  }
  if(method == "TSNE"){
    res = Rtsne(data)
    res = data.frame(ID = sample_ids, X = res$Y[,1], Y = res$Y[,2])
  } 
  plotdata = left_join(metadata %>% dplyr::select(c(variable,"ID")), res, 'ID')
  names(plotdata)[1] = "var"
  
  p = ggplot(plotdata, aes(x = X, y = Y)) +
      geom_point(aes(color = var), shape = 16) + 
      theme_minimal() 
  if(ellipses == TRUE) p + stat_ellipse(aes(color = var), level = 0.3)
  if(legend == FALSE) p + theme(legend.position = "none")
  
  return(p)
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density = function(x, y, ...) {
  dens = MASS::kde2d(x, y, ...)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

#' DESeq2 within each cell type 
#'
#' @param counts matrix of counts data 
#' @param metadata metadata where the row names match the columns of `data`
#' @param ind.filt whether to use independent filtering as described by the `results()` function in DESeq2, `TRUE` or `FALSE`
#' @param p.adjust.method method used to adjust the p-values for multiple testing. Options, in increasing conservatism, include "none", "BH", "BY" and "holm"
#' 
#' @return list of differential expression testing results 
DESeq2 = function(counts, metadata, ind.filt = FALSE, alpha = 0.05, p.adjust.method = "BH"){
  
  res_lst = dlply(metadata, .(CT), .fun = function(x){
    
    dds = DESeqDataSetFromMatrix(countData = counts[,x$ID], colData = x, design = ~ Age_scaled + Sex + DX)
    dds = DESeq(dds, parallel = TRUE, BPPARAM = param, quiet = TRUE)
    
    contrasts = resultsNames(dds)[-c(1:3)]
    
    tmp_lst = lapply(contrasts, function(z){
      
      res = results(dds, name = z, independentFiltering = ind.filt, alpha = alpha, pAdjustMethod = p.adjust.method)
      res = as.data.frame(res) %>% rownames_to_column(var = "gene_ensembl")
      res$DX = z
      
      return(res)
    })
    res = do.call(rbind, tmp_lst)
    return(res)
  })
  names(res_lst) = unique(metadata$CT)
  return(res_lst)
}

#' Limma-voom within each cell type 
#'
#' @param counts matrix of counts data 
#' @param metadata metadata where the row names match the columns of `data`
#' @param p.adjust.method method used to adjust the p-values for multiple testing. Options, in increasing conservatism, include "none", "BH", "BY" and "holm"
#' 
#' @return list of differential expression testing results 
LimmaVoom = function(counts, metadata, p.adjust.method = "BH", robust.ebayes = TRUE, quantile.norm = TRUE){
  
  res_lst = dlply(metadata, .(CT), .fun = function(x){
    
    print(as.character(unique(x$CT)))
    
    dge = DGEList(counts[x$ID], group = x$DX)
    if(quantile.norm == TRUE) { 
      dge = calcNormFactors(dge, method = "none")
    } else {
      dge = calcNormFactors(dge, method = "TMM")
    }
    
    design = model.matrix(~ 0 + DX + Age_scaled + Sex, x)
    contrs = makeContrasts(DXMDD = DXMDD-DXCTRL,
                           DXBD = DXBD-DXCTRL,
                           DXSCZ = DXSCZ-DXCTRL,
                           levels = colnames(design))
        par(mfrow=c(1,2))
    
    if(quantile.norm == TRUE){ 
      vm = voom(dge, design, plot = TRUE, normalize.method = "quantile")
    } else {
      vm = voom(dge, design, plot = TRUE)
    }
        
    fit = lmFit(vm, design)
    fit = contrasts.fit(fit, contrs)
    fit = eBayes(fit, robust = robust.ebayes)
    
        plotSA(fit, main = "Final model: Mean-variance trend")
      
        lcpm = edgeR::cpm(dge, log = TRUE)
        boxplot(lcpm, las = 2, main = "")
        title(main = paste0("calNormFactors\n",as.character(unique(x$CT))), ylab = "Log-cpm")
        
        boxplot(vm$E, las = 2, main = "")
        title(main = paste0("voom\n",as.character(unique(x$CT))), ylab = "Log-cpm")
        
        print(summary(decideTests(fit, p.value = 0.1, adjust.method = p.adjust.method)))
    
    contrasts = colnames(contrs)
    tmp_lst = lapply(contrasts, function(x){
      res = topTable(fit, coef = x, sort = "none", n = Inf, adjust.method = p.adjust.method) %>% rownames_to_column(var = "gene_ensembl")
      res$DX = x
      return(res)
    })
    res = do.call(rbind, tmp_lst)
    return(res)
  })
  names(res_lst) = unique(metadata$CT)
  return(res_lst)
}

#' edgeR-LRT within each cell type 
#'
#' @param counts matrix of counts data 
#' @param metadata metadata where the row names match the columns of `data`
#' @param p.adjust.method method used to adjust the p-values for multiple testing. Options, in increasing conservatism, include "none", "BH", "BY" and "holm"
#' 
#' @return list of differential expression testing results 
EdgeRLRT = function(counts, metadata, p.adjust.method = "BH"){
  
  i = 0
  res_lst = dlply(metadata, .(CT), .fun = function(x){
    
    print(paste0(i+1, ": ", as.character(unique(x$CT))))
    
    dge = DGEList(counts[x$ID], group = x$DX)
    dge = calcNormFactors(dge, method = "TMM")
    
    design = model.matrix(~ 0 + DX + Age_scaled + Sex, x)
    contrs = makeContrasts(DXMDD = DXMDD-DXCTRL,
                           DXBD = DXBD-DXCTRL,
                           DXSCZ = DXSCZ-DXCTRL,
                           levels = colnames(design))
    
    dge = estimateDisp(dge, design = design)
    print(dge$common.dispersion)
    plotBCV(dge)
    
    fit = glmFit(dge, design = design)
    
    contrasts = colnames(contrs)
    tmp_lst = lapply(contrasts, function(z){
      lrt = glmLRT(fit, contrast = contrs[,z])
      res = topTags(lrt, n = Inf, adjust.method = p.adjust.method, sort.by = "none") %>% as.data.frame() %>% rownames_to_column(var = "gene_ensembl")
      res$DX = z
      return(res)
    })
    res = do.call(rbind, tmp_lst)
    return(res)
  })
  names(res_lst) = unique(metadata$CT)
  return(res_lst)
}




























































