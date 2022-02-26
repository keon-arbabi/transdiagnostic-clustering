suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(plyr)
  library(ggpubr)
  library(limma)
  library(edgeR)
  library(WGCNA)
  library(variancePartition)
})

# for parallel processing
param = SnowParam(4, "SOCK", progressbar=TRUE)
# for plotting voom trend 
setMethod("plot", signature(x="EList", y="missing"),
          function(x,  y, ...){
            if( is.null(x$voom.xy) || is.null(x$voom.line)){
              stop("x does not contain the trend information.\nvoom() must be run with save.plot=TRUE")
            }
            # points
            df = data.frame(x = x$voom.xy$x,
                            y = x$voom.xy$y)
            # trend line
            df.line = data.frame(x = x$voom.line$x,
                                 y = x$voom.line$y)
            ggplot(df, aes(x,y)) + geom_point(size=.1) + theme_bw(15) + 
              theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + 
              geom_line(data=df.line, aes(x,y), color="red") + xlab(bquote(log[2](count + 0.5))) + 
              ylab(expression( sqrt("standard deviation"))) + ggtitle("Voom: Mean-variance trend") + 
              ylim(0, max(df$y))
          }
)

load("C:/Users/arbabik/Downloads/PITT Tetrad (Control, MDD, BPD, SCZ) SCT-RNAseq/seFullData_5p-3bp_3p-15bp.rData")

metadata = seFullData@colData@listData %>% as.data.frame() %>% column_to_rownames(var = "X")
# make proper names 
metadata$Cell.Type = gsub("-", "", metadata$Cell.Type)
names(metadata)[2:3] = c("HU","CT")

countMatrix = seFullData@assays$data@listData$counts %>% as.data.frame()
rownames(countMatrix) = seFullData@rowRanges@partitioning@NAMES

# fix sample naming issue 
names(countMatrix)[names(countMatrix) %in% c('PV-Hu1195.bam','PV-Hu1222.bam')] = c("PVALB-Hu1195.bam", "PVALB-Hu1222.bam")
countMatrix = countMatrix[, match(rownames(metadata), colnames(countMatrix))]
# check
all(rownames(metadata) == colnames(countMatrix)) 

# filter genes by expression 
isexpr = rowSums(edgeR::cpm(countMatrix)>0.1) >= 130
# normalize 
geneExpr = DGEList(countMatrix[isexpr,] )
geneExpr = calcNormFactors(geneExpr, "none") # TMM

# convert factors 
contVars = c("Age","PMI","pH","RNA.Ratio","RIN")
catVars = c("HU","Tetrad","CT","Subject.Group","MOD","Sex","Race")
metadata[,catVars] = lapply(metadata[,catVars], factor)

form = ~ CT + Subject.Group + CT*Subject.Group + scale(Age) + Sex + (1|HU)

vobj = voomWithDreamWeights(geneExpr, form, metadata, save.plot = T, BPPARAM = param)
plot(vobj)
fitExtractVarPartModel(vobj_dream, form, METADATA)

# compare each cell-type to all others 
L = makeContrastsDream(form, metadata, 
                       contrasts = c(PVALB = "CTPVALB - (CTPyrL2n3 + CTPyrL5n6 + CTSST + CTVIP)/4",
                                     PyrL2n3 = "CTPyrL2n3 - (CTPVALB + CTPyrL5n6 + CTSST + CTVIP)/4",
                                     PyrL5n6 = "CTPyrL5n6 - (CTPVALB + CTPyrL2n3 + CTSST + CTVIP)/4",
                                     SST = "CTSST - (CTPVALB + CTPyrL2n3 + CTPyrL5n6 + CTVIP)/4",
                                     VIP = "CTVIP - (CTPVALB + CTPyrL2n3 + CTPyrL5n6 + CTSST)/4"))
plotContrasts(L)
# fit model for each gene 
fit = dream(vobj, form, metadata, L)
fit = eBayes(fit)

# Get background genes 
backgroundGenes = data.frame(ensembl_gene_id = rownames(countMatrix))
# Define biomart object
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        host = "uswest.ensembl.org",
                        dataset = "hsapiens_gene_ensembl")
# Query biomart
Ensemble2HGNC = biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype", "chromosome_name"),
                               filters = "ensembl_gene_id", 
                               values = backgroundGenes$ensembl_gene_id,
                               mart = mart)
# get DE genes 
resLst = list()
cts = c("PVALB", "PyrL2n3", "PyrL5n6", "SST", "VIP")

for(i in 1:length(cts)){
  resLst[[i]] = topTable(fit, coef = cts[i], number = Inf) %>%
    rownames_to_column(var = 'ensembl_gene_id') %>% 
    left_join(backgroundGenes) %>% 
    left_join(Ensemble2HGNC) 
}
names(resLst) = cts

saveRDS(resLst, file = "C:/Users/arbabik/Downloads/resLst.rds")


















































# form = ~ 0 + CT + Subject.Group + Sex + scale(Age) + scale(PMI) + scale(RIN) 
# # ideally we would used a mixed-effects model...
# #form = ~ 0 + CT + Subject.Group + Sex + scale(Age) + scale(PMI) + scale(RIN) + (1|HU)
# 
# vobj = voomWithDreamWeights(geneExpr, form, metadata, save.plot = T, BPPARAM = param)
# plot(vobj)
# # compare each cell-type to all others 
# L = makeContrastsDream(form, metadata, 
#                        contrasts = c(PVALB = "CTPVALB - (CTPyrL2n3 + CTPyrL5n6 + CTSST + CTVIP)/4",
#                                      PyrL2n3 = "CTPyrL2n3 - (CTPVALB + CTPyrL5n6 + CTSST + CTVIP)/4",
#                                      PyrL5n6 = "CTPyrL5n6 - (CTPVALB + CTPyrL2n3 + CTSST + CTVIP)/4",
#                                      SST = "CTSST - (CTPVALB + CTPyrL2n3 + CTPyrL5n6 + CTVIP)/4",
#                                      VIP = "CTVIP - (CTPVALB + CTPyrL2n3 + CTPyrL5n6 + CTSST)/4"))
# plotContrasts(L)
# # fit model for each gene 
# fit = dream(vobj, form, metadata, L)
# fit = eBayes(fit)
# 
# # Get background genes 
# backgroundGenes = data.frame(ensembl_gene_id = rownames(countMatrix))
# # Define biomart object
# mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         host = "uswest.ensembl.org",
#                         dataset = "hsapiens_gene_ensembl")
# # Query biomart
# Ensemble2HGNC = biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype", "chromosome_name"),
#                                filters = "ensembl_gene_id", 
#                                values = backgroundGenes$ensembl_gene_id,
#                                mart = mart)
# # get DE genes 
# resLst = list()
# cts = c("PVALB", "PyrL2n3", "PyrL5n6", "SST", "VIP")
# 
# for(i in 1:length(cts)){
#   resLst[[i]] = topTable(fit, coef = cts[i], number = Inf) %>%
#     rownames_to_column(var = 'ensembl_gene_id') %>% 
#     left_join(backgroundGenes) %>% 
#     left_join(Ensemble2HGNC) 
# }
# names(resLst) = cts
# 
# saveRDS(resLst, file = "C:/Users/arbabik/Downloads/resLst.rds")

