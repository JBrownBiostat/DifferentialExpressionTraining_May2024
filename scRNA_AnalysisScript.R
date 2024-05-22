##############################################################################
## Tutorial script for single-cell RNA-seq differential expression analysis ##
##############################################################################


###############
## Libraries ##
###############
library(DESeq2)
library(ggplot2)
library(TENxPBMCData)
library(Matrix)
library(irlba)
library(matrixStats)
library(scran)
library(BiocParallel)
library(scater)


##########
## Data ##
##########
set.seed(93051421)
pbmc68k <- TENxPBMCData(dataset = "pbmc68k")
pbmc68k <- pbmc68k[, sample.int(ncol(pbmc68k), 5e3)]
counts(pbmc68k) <- as(counts(pbmc68k), "sparseMatrix")
pbmc68k <- pbmc68k[rowSums(counts(pbmc68k) > 0) >= 5e-3 * ncol(pbmc68k), ]
print(pbmc68k)

prllWorkers <- 4 # number of cores (strictly threads) to use for parallel computation


################################
## Normalization size factors ##
################################
summary(rowMeans(counts(pbmc68k) > 0))
pbmc68k <- computeSumFactors(
    pbmc68k, BPPARAM = SnowParam(workers = prllWorkers)
)
pbmc68k <- logNormCounts(pbmc68k)


###################
## Quick cluster ##
###################
# Top PCs #
top_pbmc68k <- getTopHVGs(pbmc68k, n = 2500)
pbmc68k <- fixedPCA(pbmc68k, subset.row = top_pbmc68k, rank = 20) 
plotReducedDim(pbmc68k, dimred="PCA")

# umap dimension reduction #
pbmc68k <- runUMAP(pbmc68k, dimred="PCA")
plotReducedDim(pbmc68k, dimred="UMAP")

# clustering #
graphClusters <- clusterCells(pbmc68k, use.dimred="PCA")
graphClusters <- LETTERS[graphClusters]
table(graphClusters)
colLabels(pbmc68k) <- graphClusters
plotReducedDim(pbmc68k, "UMAP", colour_by="label")


#############################
## Differential expression ##
#############################
metadat <- data.frame(
    cluster = factor(graphClusters, levels = sort(unique(graphClusters)))
)
dds <- DESeqDataSetFromMatrix(
    countData = counts(pbmc68k),
    colData = metadat,
    design = ~ cluster - 1
)
sizeFactors(dds) <- pbmc68k$sizeFactor
dds <- DESeq(
    dds,
    minReplicatesForReplace = Inf,
    minmu = 1e-6,
    parallel = TRUE,
    BPPARAM = SnowParam(workers = prllWorkers)
)
plotDispEsts(dds)
resultsNames(dds)
shrunk_AvsB <- lfcShrink(
    dds,
    contrast = c("cluster", "A", "B"),
    type = "ashr",
    svalue = TRUE,
    parallel = TRUE,
    BPPARAM = SnowParam(workers = prllWorkers)
)
head(shrunk_AvsB)
table(shrunk_AvsB$svalue <= 1e-3)
plot(
    x = log2(shrunk_AvsB$baseMean),
    y = shrunk_AvsB$log2FoldChange,
    col = ifelse(shrunk_AvsB$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25
)



