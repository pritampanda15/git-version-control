# HGEN 473 - Genomics
# Spring 2017
# Tuesday, May 9 & Thursday, May 11
#
# RNA-seq analysis with R/Bioconductor
#
# John Blischak
#
# Last updated: 2020-04-08

# Introduction -------------------------------------------------------

# The goal of this tutorial is to introduce you to the analysis of
# RNA-seq data using some of the powerful, open source software
# packages provides by R, and specifically the Bioconductor project.
#
# Before the lecture, please make sure you have R and RStudio
# installed, and also run the code in the sections "Installation" and
# "Download data" below.
#
# The code below is adapted from the paper "RNA-seq analysis is easy
# as 1-2-3 with limma, Glimma and edgeR" by Charity et al., 2017. The
# original paper and this derivative work are freely available under
# the CC-BY license.
#
# Publication: https://f1000research.com/articles/5-1408/v2
#
# Source code: https://bioconductor.org/packages/release/workflows/html/RNAseq123.html
#
# CC-BY: https://creativecommons.org/licenses/by/4.0/
#
# Full citation:
#
# Law CW, Alhamdoosh M, Su S et al. RNA-seq analysis is easy as 1-2-3
# with limma, Glimma and edgeR [version 2; referees: 3 approved].
# F1000Research 2016, 5:1408 (doi: 10.12688/f1000research.9005.2)

# Installation -------------------------------------------------------

# To install the Bioconductor and CRAN packages used in this tutorial, run the
# lines below to install the RNAseq123 workflow package. If it asks if you would
# like to install the packages in a personal directory, confirm yes (Y).

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RNAseq123")

# Download data ------------------------------------------------------

# The example data used in this tutorial is RNA-seq of 3 cell
# populations in the mouse mammary gland collected by Sheridan et al.,
# 2015. The code below downloads the count data from the Gene
# Expression Omnibus (GEO), an NCBI repository of genomics data.
#
# If you receive an error message from the function `download.file`,
# try setting the argument `method`. For example try `method =
# "wget"`` or `method = "curl"`.

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile = "GSE63310_RAW.tar", mode = "wb")
utils::untar("GSE63310_RAW.tar", exdir = ".")
files_gz <- Sys.glob("GSM*txt.gz")
for(f in files_gz)
  R.utils::gunzip(f, overwrite = TRUE)

# Citation for data set:
#
# Sheridan JM, Ritchie ME, Best SA, et al.: A pooled shRNA screen for
# regulators of primary mammary stem and progenitor cells identifies
# roles for Asap1 and Prox1. BMC Cancer. 2015; 15(1): 221.

# Setup --------------------------------------------------------------

# Load the packages we will use.

library("limma")   # Linear models for differential expression
library("Glimma")  # Interactive plots for exploration
library("edgeR")   # Process count data from NGS experiments
library("Mus.musculus") # Gene annotations for the Mus musculus genome

# Import the data.

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",
           "GSM1545545_JMS9-P8c.txt")

read.delim(files[1], nrow = 5)

x <- readDGE(files, columns = c(1, 3))
class(x)
dim(x)
names(x)
str(x)

# Annotate the samples.

x$samples
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP",
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004", "L006", "L008"), c(3, 4, 2)))
x$samples$lane <- lane
x$samples

# Annotate the genes.

head(x$counts)
dim(x$counts)
geneid <- rownames(x)
genes <- select(Mus.musculus, keys = geneid,
                columns = c("SYMBOL", "TXCHROM"),
                keytype = "ENTREZID")
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]

x$genes <- genes
x

# Filter genes -------------------------------------------------------

# Calculate (log) counts-per-million (cpm)

cpm <- cpm(x)
lcpm <- cpm(x, log = TRUE)

# Detect genes with zero counts across all 9 samples

table(rowSums(x$counts == 0) == 9)

# Visualize distribution of gene expression levels

plotDensities(lcpm, legend = FALSE, main = "Before filtering")
abline(v = 0, lty = 3)

# Only keep genes which have cpm greater than 1 in at least 3 samples.

keep.exprs <- rowSums(cpm > 1) >= 3
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
dim(x)
lcpm <- cpm(x, log=TRUE)

# Visualize distribution of gene expression levels after filtering

plotDensities(lcpm, legend = FALSE, main = "After filtering")
abline(v = 0, lty = 3)

# Discuss: What information are we ignoring when we use
# counts-per-million? Is it a valid measurement for performing a
# differential expression analysis?

# Normalization ------------------------------------------------------

# The default normalization provided by edgeR is TMM (trimmed mean of
# M-values), which prevents differences in highly expressed genes from
# biasing the entire distribution of gene expression. It often has a
# modest effect, as observed here.

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# But here is a extreme toy example that demonstrates it will work if
# necessary.

x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[, 1] * 0.05)
x2$counts[,2] <- x2$counts[, 2] * 5

lcpm2 <- cpm(x2, log = TRUE)
boxplot(lcpm2, las = 2, main = "Before normalization")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors
lcpm2 <- cpm(x2, log = TRUE)
boxplot(lcpm2, las=2, main = "After normalization")

# Exploration --------------------------------------------------------

# Visualize sample relationships with multidimensional scaling (MDS).

library("RColorBrewer")
group
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
lane
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels = group, col = col.group,
        main = "group")
plotMDS(lcpm, labels = lane, col = col.lane, dim = c(3, 4),
        main = "lane")

# Interactive MDS plot

glMDSPlot(lcpm, labels = paste(group, lane, sep = "_"),
          groups = x$samples[, c(2, 5)], launch = TRUE)

# Construct linear model ---------------------------------------------

# The linear model we construct will not have an intercept. This is
# referred to as the group-means parametrization in the limma User's
# Guide. The advantage of having each coefficient (beta) model the
# mean expression level for that group is that it makes it more
# straightforward to test specific hypotheses.

design <- model.matrix(~0 + group + lane)
colnames(design) <- gsub("group", "", colnames(design))
design

# Here is the model specified in LaTeX, where Y is the expression
# level of a particular gene in a particular sample:
#
# Y=\beta_{basal}+\beta_{LP}+\beta_{ML}+\beta_{L06}+\beta_{L08}+\epsilon

# We create the following 3 contrasts to test the following 3
# hypotheses:
#
# BasalvsLP: Which genes are DE between Basal and LP cells?
#
# BasalvsML: Which genes are DE between Basal and ML cells?
#
# LPvsML: Which genes are DE between LP and ML cells?

contr.matrix <- makeContrasts(BasalvsLP = Basal - LP,
                              BasalvsML = Basal - ML,
                              LPvsML = LP - ML,
                              levels = colnames(design))
contr.matrix

# Convert counts to be used in linear model --------------------------

# The RNA-seq counts cannot be used directly in the linear model
# because they violoate its assumptions, specifically that the
# variance should not depend on the mean. One option is to perform a
# test that directly models the counts (e.g. such models are provided
# by edgeR and DESeq2). However, the linear modelling framework is
# generally more flexible, and limma has many nice downstream
# functions for further testing the data (e.g. testing for enrichment
# of functional categories).
#
# Even after converting to log-cpm, the RNA-seq data still has the
# mean-variance relationship. Thus the function `voom` calculates
# weights to offset this relationship.

v <- voom(x, design, plot = TRUE)
v

# Note that for convenience `voom` also calculates the log-cpm and
# optionally applies additional normalization. These side-effects
# should not be mistaken as the primary purpose of `voom` since
# standardization to log-cpm and normalization can easily be done by
# other functions (which is what `voom` does under the hood). The
# primary purpose of `voom` is to calculate the weights to be used in
# the linear model.

# Test for differential expression (DE) ------------------------------

# Fit a linear model per gene.

vfit <- lmFit(v, design)

# Calculate the statistics for our specific contrasts of interest.

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

# Calculate levels of significance. Uses an empirical Bayes algorithm
# to shrink the gene-specific variances towards the average variance
# across all genes.

efit <- eBayes(vfit)

# As seen in this diagnostic plot of the residual variation versus the
# mean gene expression level, `voom` successfully removed the
# relationship between the mean and the variance.

plotSA(efit, main = "Final model: Mean−variance trend")

# Discuss: What is the value of sharing information across genes? In
# other words, why not just use the original variance calculated per
# gene?

# Explore the results ------------------------------------------------

# Tabulate the results

summary(decideTests(efit))

# If the magnitude of the effect size is important for your downstream
# analysis (e.g. perhaps you are trying to prioritize genes for futher
# experiments), you can specify a minimal log-fold-change with the
# function `treat` isntead of using `eBayes`.

tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
summary(dt)

# Create a venn diagram of the results.

head(dt)
de.common <- which(dt[, 1] != 0 & dt[, 2] != 0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n = 20)
vennDiagram(dt[, 1:2], circle.col = c("turquoise", "salmon"))
write.fit(tfit, dt, file = "results.txt")

# Identify top DE genes.

basal.vs.lp <- topTreat(tfit, coef = 1, n = Inf)
basal.vs.ml <- topTreat(tfit, coef = 2, n = Inf)
head(basal.vs.lp)
head(basal.vs.ml)

# Visualize DE genes.

plotMD(tfit, column = 1, status = dt[, 1], main = colnames(tfit)[1],
       xlim = c(-8, 13))

glMDPlot(tfit, coef = 1, status = dt, main = colnames(tfit)[1],
         side.main = "ENTREZID", counts = x$counts, groups = group,
         launch = TRUE)

# View heatmap of top 100 DE genes between Basal and LP cells.

library("gplots")
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v$E[i, ], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group,
          col = mycol, trace = "none", density.info = "none")

# Discuss: Why did we scale by row (per gene) in the heatmap?

# Test for enrichment of gene sets -----------------------------------

# A common first-step post-DE testing is to test for enrichment of
# specific gene sets, e.g. Gene Ontology (GO) categories or KEGG
# pathways. Here we will use the gene set collection MSigDB c2 from
# the Broad Institute.

load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
class(Mm.c2)
Mm.c2$KEGG_GLYCOLYSIS_GLUCONEOGENESIS
idx <- ids2indices(Mm.c2, id = rownames(v))

# The limma function `camera` tests for enrichment for a specific
# contrast of interest by comparing all the gene sets against one
# another.

cam.BasalvsLP <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.BasalvsLP, 5)
cam.BasalvsML <- camera(v, idx, design, contrast = contr.matrix[, 2])
head(cam.BasalvsML, 5)
cam.LPvsML <- camera(v, idx, design, contrast = contr.matrix[, 3])
head(cam.LPvsML, 5)

# Visualize gene set enrichment with a barcode plot.

barcodeplot(efit$t[, 3], index = idx$LIM_MAMMARY_LUMINAL_MATURE_UP,
            index2 = idx$LIM_MAMMARY_LUMINAL_MATURE_DN,
            main = "LPvsML")

# Note that the directions of the effect are reversed from the results
# in the original study because we tested the contrast LP - ML,
# whereas the original analysis tested ML - LP.

# Report session information -----------------------------------------

# It's always a good idea to report the versions of the software you
# used to generate results.

sessionInfo()

# The code in this tutorial was tested with the following setup:

# > sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
#
# Matrix products: default
#
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252
#
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] gplots_3.0.3                              RColorBrewer_1.1-2
# [3] Mus.musculus_1.3.1                        TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
# [5] org.Mm.eg.db_3.10.0                       GO.db_3.10.0
# [7] OrganismDbi_1.28.0                        GenomicFeatures_1.38.2
# [9] GenomicRanges_1.38.0                      GenomeInfoDb_1.22.1
# [11] AnnotationDbi_1.48.0                      IRanges_2.20.2
# [13] S4Vectors_0.24.3                          Biobase_2.46.0
# [15] BiocGenerics_0.32.0                       edgeR_3.28.0
# [17] Glimma_1.14.0                             limma_3.42.0
#
# loaded via a namespace (and not attached):
#   [1] httr_1.4.1                  bit64_0.9-7                 jsonlite_1.6.1              R.utils_2.9.2
# [5] gtools_3.8.2                assertthat_0.2.1            askpass_1.1                 BiocManager_1.30.10
# [9] BiocFileCache_1.10.2        RBGL_1.62.1                 blob_1.2.1                  GenomeInfoDbData_1.2.2
# [13] Rsamtools_2.2.3             progress_1.2.2              pillar_1.4.3                RSQLite_2.2.0
# [17] lattice_0.20-40             glue_1.3.1                  digest_0.6.25               XVector_0.26.0
# [21] Matrix_1.2-18               R.oo_1.23.0                 XML_3.99-0.3                pkgconfig_2.0.3
# [25] biomaRt_2.42.1              zlibbioc_1.32.0             purrr_0.3.3                 gdata_2.18.0
# [29] BiocParallel_1.20.1         tibble_2.1.3                openssl_1.4.1               SummarizedExperiment_1.16.1
# [33] magrittr_1.5                crayon_1.3.4                memoise_1.1.0               R.methodsS3_1.8.0
# [37] graph_1.64.0                tools_3.6.3                 prettyunits_1.1.1           hms_0.5.3
# [41] matrixStats_0.56.0          stringr_1.4.0               locfit_1.5-9.1              DelayedArray_0.12.2
# [45] Biostrings_2.54.0           compiler_3.6.3              caTools_1.18.0              rlang_0.4.5
# [49] grid_3.6.3                  RCurl_1.98-1.1              rstudioapi_0.11             rappdirs_0.3.1
# [53] bitops_1.0-6                DBI_1.1.0                   curl_4.3                    R6_2.4.1
# [57] GenomicAlignments_1.22.1    dplyr_0.8.5                 rtracklayer_1.46.0          bit_1.1-15.2
# [61] KernSmooth_2.23-16          stringi_1.4.6               Rcpp_1.0.3                  vctrs_0.2.4
# [65] dbplyr_1.4.2                tidyselect_1.0.0

# Further reading ----------------------------------------------------

# The limma User's Guide (run `limmaUsersGuide()` in the R console) is
# very useful. Especially see Section 2.1 for citations to the primary
# publications that describe the methods and Chapter 9 for how to
# contruct the model for different types of study designs.

# The Bioconductor support site (https://support.bioconductor.org/)
# has years worth of questions and answers about RNA-seq analysis and
# other topics in bioinformatics.

# Overview of RNA-seq analysis
#
# A. Oshlack, M. D. Robinson and M. D. Young. “From RNA-seq reads to
# differential expression results”. In: _Genome Biology_ 11.12 (2010),
# p. 220. DOI: 10.1186/gb-2010-11-12-220.

# R/Bioconductor tutorial starting from fastq files
#
# Chen Y, Lun ATL and Smyth GK. From reads to genes to pathways:
# differential expression analysis of RNA-Seq experiments using
# Rsubread and the edgeR quasi-likelihood pipeline [version 2;
# referees: 5 approved]. F1000Research 2016, 5:1438 (doi:
# 10.12688/f1000research.8987.2)

# Comparisons of RNA-seq methods for differential expression testing
#
# F. Rapaport, R. Khanin, Y. Liang, et al. “Comprehensive evaluation
# of differential gene expression analysis methods for RNA-seq data”.
# In: _Genome Biology_ 14.9 (2013), p. R95. DOI:
# 10.1186/gb-2013-14-9-r95.
#
# C. Soneson and M. Delorenzi. “A comparison of methods for
# differential expression analysis of RNA-seq data”. In: _BMC
# Bioinformatics_ 14.1 (2013), p. 91. DOI: 10.1186/1471-2105-14-91.

# limma
#
# Ritchie ME, Phipson B, Wu D, et al.: limma powers differential
# expression analyses for RNA-sequencing and microarray studies.
# Nucleic Acids Res. 2015; 43(7): e47.

# Glimma
#
# Su S, Ritchie ME: Glimma: Interactive HTML graphics for RNA-seq
# data. 2016; R package version 1.1.1.

# edgeR
#
# Robinson MD, McCarthy DJ, Smyth GK: edgeR: a Bioconductor package
# for differential expression analysis of digital gene expression
# data. Bioinformatics. 2010; 26(1): 139–140.

# Bioconductor project
#
# Huber W, Carey VJ, Gentleman R, et al.: Orchestrating
# high-throughput genomic analysis with Bioconductor. Nat Methods.
# 2015; 12(2): 115–121.

# TMM normalization
#
# Robinson MD, Oshlack A: A scaling normalization method for
# differential expression analysis of RNA-seq data. Genome Biol. 2010;
# 11(3): R25.

# limma+voom
#
# Law CW, Chen Y, Shi W, et al.: voom: Precision weights unlock linear
# model analysis tools for RNA-seq read counts. Genome Biol. 2014;
# 15(2): R29

# Empirical Bayes to estimate gene expression variance
#
# Smyth GK: Linear models and empirical bayes methods for assessing
# differential expression in microarray experiments. Stat Appl Genet
# Mol Biol. 2004; 3(1): Article3.
