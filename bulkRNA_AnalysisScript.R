#######################################################################
## Tutorial script for bulk RNA-seq differential expression analysis ##
#######################################################################


###############
## Libraries ##
###############
library(pasilla)
library(DESeq2)
library(ggplot2)


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


##################################
## Custom normalization factors ##
##################################
# Factors by library size: Reads per million #
rpmSF <- colSums(cntMat) / 1e6

# Factors by Median-of-Ratios #
filterCnts <- cntMat[rowSums(cntMat <= 5) == 0, ]
dim(filterCnts)
filterGeoMeans <- exp(rowMeans(log(filterCnts)))
mrSF <- as.numeric(apply(filterCnts / filterGeoMeans, 2, median))

# Factors calculated by default #
defaultSF <- sizeFactors(dds_quickStart)

# Compare custom to default size factors #
plotOrd <- order(defaultSF)
plot(
    x = log2(defaultSF[plotOrd]),
    y = log2(rpmSF[plotOrd]),
    type = "b", col = "red", pch = 19, 
    ylim = range(log2(c(rpmSF, mrSF))) + c(-0.1, 0.1),
    xlim = range(log2(defaultSF)) + c(-0.1, 0.1),
    xlab = "Default Factors (log2)", ylab = "Custom Factors (log2)"
)
lines(
    x = log2(defaultSF[plotOrd]),
    y = log2(mrSF[plotOrd]),
    type = "b", col = "blue", pch = 19, lty = 2
)
abline(0, 1, col = "black", lty = 2)
legend("topleft", legend = c("RPM", "MR"), col = c("red", "blue"), lty = 1:2)

# Visualize normalization quality #
normDat <- data.frame(
    normExpr = c(
        as.numeric(t(t(cntMat) / rpmSF)),
        as.numeric(t(t(cntMat) / mrSF))
    ),
    normMethod = rep(
        c("RPM", "MR"),
        rep(prod(dim(cntMat)), 2)
    ),
    sample = rep(
        rep(
            colnames(cntMat),
            rep(nrow(cntMat), ncol(cntMat))
        ),
        2
    )
)
normQualityPlot <- ggplot(
    normDat, 
    aes(x = sample, y = normExpr, color = sample, fill = sample)
) + 
    theme_classic() + 
    facet_grid(. ~ normMethod) + 
    geom_boxplot(alpha = 0.5) + 
    scale_color_brewer(palette = "Dark2") + 
    scale_fill_brewer(palette = "Dark2") + 
    scale_y_continuous(
        trans = "log1p",
        breaks = 10 ^ (0:6)
    ) + 
    labs(
        x = "Sample", y = "Normalized expression", color = "", fill = ""
    ) + 
    theme(
        axis.text.x = element_text(angle = -25, hjust = 0)
    )
normQualityPlot

# DESeq Analysis: RPM #
dds_RPM <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)
sizeFactors(dds_RPM) <- rpmSF
dds_RPM <- DESeq(dds_RPM)
head(results(dds_RPM))
plotMA(dds_RPM)

# DESeq Analysis: MR #
dds_MR <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)
sizeFactors(dds_MR) <- mrSF
dds_MR <- DESeq(dds_MR)
head(results(dds_MR))
plotMA(dds_MR)


###############################################
## Exploring estimated dispersion parameters ##
###############################################
# Dispersion plot: MR #
plotDispEsts(dds_MR)

# Dispersion plot: RPM #
plotDispEsts(dds_RPM)



