#######################################################################
## Tutorial script for bulk RNA-seq differential expression analysis ##
#######################################################################


###############
## Libraries ##
###############
library(pasilla)
library(DESeq2)


####################################################
## Data: scripts adapted from the DESeq2 vignette ##
####################################################
# Import Data #
cntFilePath <- system.file(
    "extdata",
    "pasilla_gene_counts.tsv",
    package="pasilla", 
    mustWork=TRUE
)
metaFilePath <- system.file(
    "extdata",
    "pasilla_sample_annotation.csv",
    package="pasilla", 
    mustWork=TRUE
)
cntMat <- as.matrix(read.csv(cntFilePath, sep = "\t", row.names = "gene_id"))
metaDat <- read.csv(metaFilePath, row.names = 1)
metaDat$condition <- factor(metaDat$condition)
metaDat$type <- factor(metaDat$type)

# Check sample structure #
print(metaDat)


#################################
## Quick-start simple analysis ##
#################################
# Re-format condition as a factor type #
metaDat$type <- factor(
    x = as.character(metaDat$type),
    levels = c("untreated", "treated")
)

# Match rows of meta data against columns of cntMat #
rownames(metaDat) <- gsub("fb", "", rownames(metaDat))
all(rownames(metaDat) %in% colnames(cntMat))
all(colnames(cntMat) %in% rownames(metaDat))
cntMat <- cntMat[, rownames(metaDat)]
all(colnames(cntMat) == rownames(metaDat))

# Setup DESeqDataSet #
dds_quickStart <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)

# Run analysis and view results #
dds_quickStart <- DESeq(dds_quickStart)
head(results(dds_quickStart))
plotMA(dds_quickStart)






