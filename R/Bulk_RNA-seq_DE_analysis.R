# ==============================================================================
# Bulk RNA-seq DE analysis
# statistical inference, creating figures/plots of the data, interpreting the results
# Dataset (from dbGa): phs000007.v23.p8, phs000280.v2.p1, and phs000209.v10.p2.
# Author: Giuseppe Castiglione
# Creation date: 10/06/2026
# ==============================================================================

# 1. LIBRARIES # ==================================================================

#install.packages("tidyverse")
library(tidyverse)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(BiocManager)

#BiocManager::install("limma") #package for analyzing bulk RNA-seq data
#BiocManager::install("edgeR") #package for analyzing bulk RNA-seq data
#BiocManager::install("ggrepel")
library(limma)
library(edgeR)
library(ggrepel)

library(here)
library(writexl)

# 2. PATH # =======================================================================
path <- here("Data", "salmon.merged.gene_counts.tsv")
path
# dr_here()


# 3. PROCESSING DATA # =========================================================

# LOAD GENE-LEVEL ABUNDANCE DATA
# salmon.merged.gene_counts.tsv <- "gene_id", "gene_name" (gene symbol)
# Remove id and name from the data table in order to manipulate the data into a matrix format
expression_data <-read.table(path,header=TRUE,row.names="gene_name", quote="")
expression_data$gene_id <- NULL
head(expression_data)


# LOAD SAMPLE INFORMATION
# pop_tempo <- concatenation of the two categorical variables to filter out lowly expressed genes for sufficient counts
sample_info <- read.table(here("Data", "dme_elev_samples.tsv"),header = TRUE, stringsAsFactors=FALSE) #separate columns "population" and "temp"
sample_info$pop_temp <- paste(sample_info$population, sample_info$temp,sep="_")

# Verify the new object:
#head(sample_info)
#View(sample_info)
# Correct output: in "pop_temp" column, result as "maine_low" or "maine_high"


# CREATE DGE OBJECT
# Verify that rows order in sample_info = cols order in expression_data
#colnames(expression_data)
#sample_info$sample
# Output: the order of the two outputs should be exactly the same

DGE <- DGEList(expression_data,samples=sample_info$sample,group=sample_info$pop_temp)
DGE$samples # To verify group and lib.size for each sample are ok
sort(DGE$samples$lib.size) # Output: lib.size (seq depth) should be > 1x10^7,  


# FILTERING LOWLY EXPRESSED GENES
# filterByExpr() in `edgeR` to identify genes have a minimum count for a minimum number of samples within each group defined DGE.  
keep.exprs <- filterByExpr(DGE)
DGE <- DGE[keep.exprs,, keep.lib.sizes=FALSE]
# Output: boolean output of that function (TRUE=minimum count criterion is satisfied) to filter the DGE.


# TMM NORMALIZATION
# Important to normalize highly activated genes among genes with same expressions between conditions
DGE <- calcNormFactors(DGE,method =c("TMM"))


# QUICK DATA VIEW FOR OUTLIERS
tempvals <- sample_info$temp
popvals <- sample_info$population
mds <- plotMDS(DGE,top=500,plot=FALSE,gene.selection="common")
mds_dataframe <- cbind(sample_info$sample,sample_info$pop_temp,mds$x,mds$y)
colnames(mds_dataframe) <- c("sample","poptemp","pc1","pc2")
mds_dataframe <- as_tibble(mds_dataframe) %>% mutate(across(c(pc1, pc2), as.numeric))

mds_plot <- mds_dataframe %>% ggplot(aes(x=pc1,y=pc2,color=poptemp)) +
  geom_point(size=3) +
  xlab("Principal coordinate 1") +
  ylab("Principal coordinate 2") +
  geom_text_repel(aes(label = sample), size=3)
print(mds_plot)
#ggsave(mds_plot,
#       path = "Outputs",
#       filename = "mds_plot_outliers.pdf",
#       device = "pdf",
#       height = 6, width = 8, units = "in")
# Output: clear separation between temperature regimes and geographic locations, but no outliers indicating bad samples or potential label swaps.
# SRR1576458 only seems a bit separated, hard to know exactly which is the cause

# DESIGN MATRIX
# For DE analysis, limma needs design matrix in order to fit the linear model with boolean variables
# The model looks at the effects of low and high temperature treatments, ignoring geography
design_temp <- model.matrix(~temp, data=sample_info)
design_temp
# Output: row number=the row in the sample table (associated sample id), templow=samples have received the low temperature treatment (=1)


# 4. RUNNING LIMMA # ===========================================================


# RUNNING LIMMA
# Goal: fitting linear models, and involves transforming count data to log-transformed 
# counts per million (CPM), estimating the gene-wise relationship between mean and variance 
# of expression, and calculating "observation-level precision weights"

pdf(here::here("Outputs", "voom_plot.pdf"))
v <- voom(DGE, design=design_temp, plot=TRUE)
dev.off()

# NOTE: in RNA-seq, there is variation in sample quality --> instead of discard outliers,
# apply weights to samples such that outlier samples are down-weighted during DE calculations


# RUNNING LIMMA WITH SAMPLE QUALITY WEIGHTS
pdf(here::here("Outputs", "voom_qualityweights_plot.pdf"), width = 20, height = 10)
vwts <- voomWithQualityWeights(DGE, design = design_temp, normalize.method = "none", plot = TRUE)
dev.off()

# normalize.method = "none" because already TMM normalization


# RUN LINEAR MODEL
# fit the linear model specified by the design matrix (design_temp) + estimate fold changes and standard errors
fit <- lmFit(vwts,design_temp)


# EMPIRICAL BAYES PROCEDURES
# Problem: noise because small number of biological replicates in a bulk RNA-experiment and estimation of gene dispersion
# Solution: standard errors estimates are adjusted (assumption: genes with similar expression levels should have similar variances)
fit <- eBayes(fit,robust=TRUE)

# robust=TRUE because we are using voomWithQualityWeights


# SUMMARY TABLE AND VOLCANO PLOT
# Threshold=fdr, BH method=Benjamini Hochberg adjusted p-values
summary(decideTests(fit,adjust.method="fdr",p.value = 0.05))
# Output: tot genes=13307, DE genes= 3317 (~ 24.9%) for the temperature treatment
# Down = the low temperature condition is down-regulated relative to the high-temperature condition (i.e. LFC <0)
# Up = the low temperature condition is up-regulated relative to the high-temperature condition (i.e. LFC >0)
# NotSig = the null hypothesis that expression does not differ between the low and high temperature conditions was not rejected

pdf(here::here("Outputs", "dirty_volcano_plot.pdf"))
volcanoplot(fit,coeff=1)
dev.off()
# Output: there are far more up-regulated genes in the low-temp condition than up-regulated genes in the high-temp


# EXPLORING DATA - TOP 10 DE GENES (sorted by p-value)
topTable(fit, adjust="BH", coef="templow", resort.by="P")
# Output: gene name (row name), log-FC, Average expression, t-stat for t-test, raw P.value, adjusted p-value, B-stat (log-odds that the gene is differentially expressed)
help(topTable)

# RESULTS TABLE FOR ALL GENES
# Useful for querying specific hypothesized candidate genes. topTable() with p.value=1
all_genes <- topTable(fit, adjust="BH", coef="templow", p.value=1, number=Inf ,resort.by="P")
# Output: coeff = the coefficient or contrast you want to extract; number = the max number of genes to list;
# adjust = the P value adjustment method; resort.by determines what criteria with which to sort the table
head(all_genes)

# saving the all_genes table
# CSV
write.csv(all_genes, here::here("Outputs", "topTable_all_genes_templow.csv"), row.names = TRUE)

# Excel
all_genes_export <- tibble::rownames_to_column(all_genes, var = "gene_id")
write_xlsx(all_genes_export, here::here("Outputs", "topTable_all_genes_templow.xlsx"))

save.image(file = here("Outputs", "limma_processing_24062026.Rdata")
#load("limma_processing_24062026.Rdata")

# 5. DE ANALYSIS: 2-FACTOR DESIGN # ===============================================
# 2 factors: temperature and population

# Build design matrix for 2-factor model
population <- factor(sample_info$population,levels=c("maine","panama"))
temperature <- factor(sample_info$temp, levels=c("high","low"))
design_2factor <- model.matrix(~population+temperature)
head(design_2factor)


# Run limma with quality weights with 2-factor design matrix
pdf(here::here("Outputs", "voom_qualityweights__2factor_plot.pdf"), width = 20, height = 10)
vwts_2factor <- voomWithQualityWeights(DGE, design=design_2factor,normalize.method="none", plot=TRUE)
dev.off()
fit_2factor <- lmFit(vwts_2factor,design_2factor)
fit_2factor <- eBayes(fit_2factor,robust=TRUE)
summary(decideTests(fit_2factor,adjust.method="fdr",p.value = 0.05))
# Thanks to the covariate (population) we have a more accurate estimation of the DE:
## 1-factor (low_temp): tot genes=13307, DE genes= 3317 (~ 24.9%) == more noise
## 2-factor (low_temp): DE genes = 4736 (~ 36.4%) == higher precision


# RESULTS TABLE FOR ALL GENES
# Fitting models with limma with a multi-factor design, two possible coefficients to select.
# 1. Interest in "effects of temperature": specify *temperaturelow*
# 2. Interest in "effects of population": specify *populationpanama*
all_genes2 <- topTable(fit_2factor, adjust="BH",coef="temperaturelow", p.value=1, number=Inf ,resort.by="P")
all_genes2$geneid <- row.names(all_genes2)
head(all_genes2)

# saving the all_genes2 table
# CSV
write.csv(all_genes2, here::here("Outputs", "topTable_all_genes2_templow.csv"), row.names = TRUE)

# Excel
all_genes_export2 <- tibble::rownames_to_column(all_genes2, var = "gene_id")
write_xlsx(all_genes_export2, here::here("Outputs", "topTable_all_genes2_templow.xlsx"))

save.image(file = here("Outputs", "limma_processing_2factor_24062026.Rdata")
#load("limma_processing__2factor_24062026.Rdata")


# 6. DATA SUBSETTING #=============================================================

# If we are interested in a specific comparison with a particular condition between two groups (e.g. high, med, low temp)

# As said in 5., two possible coefficients to select:
# 1. Interest in "effects of temperature": specify *temperaturelow*
# 2. Interest in "effects of population": specify *populationpanama*

# How to subset:
# 1. Subset the DGE to only include samples from the conditions of interest
# 2. Create a design matrix from this subsetted data
# 3. Run the DE testing with limma


# high vs low temp only within the Panama population
panama_samples <- sample_info$sample[sample_info$population=="panama"]
panama_DGE <- DGE[,panama_samples]


# Subset the sample info table
# Selecting rows with "panama" for the value of the population factor + include all columns
# limma and edgeR work with data frames (not tibbles) --> use base R code for doing this
panama_sample_info <-sample_info[sample_info$population=="panama",]
panama_sample_info$temp <- factor(panama_sample_info$temp, levels=c("high","low"))
head(panama_sample_info)


# Create design matrix for subsetted data
panama_design_temp <- model.matrix(~temp, data=panama_sample_info)
panama_design_temp


# DE testing on Panama samples
pdf(here::here("Outputs", "voom_qualityweights__subset-panama_plot.pdf"), width = 20, height = 10)
vwts_panama <- voomWithQualityWeights(panama_DGE, design=panama_design_temp,normalize.method="none", plot=TRUE)
dev.off()

panama_fit <- lmFit(vwts_panama,panama_design_temp)
panama_fit <- eBayes(panama_fit,robust=TRUE)
summary(decideTests(panama_fit,adjust.method="fdr",p.value = 0.05))
topTable(panama_fit, adjust="BH",resort.by="P")

save.image(file = here("Outputs", "limma_processing_subset_panama_24062026.Rdata")
#load("limma_processing_subset_panama_24062026.Rdata")


topTable(panama_fit, adjust="BH",resort.by="P")


all_genes3 <- topTable(fit_2factor, adjust="BH",coef="temperaturelow", p.value=1, number=Inf ,resort.by="P")
all_genes3$geneid <- row.names(all_genes3)

# saving the all_genes3 table
# CSV
write.csv(all_genes3, here::here("Outputs", "topTable_all_genes3_templow.csv"), row.names = TRUE)
# Excel
all_genes_export3 <- tibble::rownames_to_column(all_genes3, var = "gene_id")
write_xlsx(all_genes_export3, here::here("Outputs", "topTable_all_genes3_templow.xlsx"))


# Building a contrast design matrix
# Variable concatenation: instead of modeling population and temperature as separate factors, 
# we combine them into a single factor to gain full flexibility in specifying any biological comparison we want
poptemp2 <- factor(sample_info$pop_temp,levels=c("maine_low","maine_high","panama_low","panama_high"))

contrast_design <- model.matrix(~0+poptemp2)
colnames(contrast_design) <- levels(poptemp2)
contrast_design
# Output: each value for poptemp gets a column, with ones indicating if the sample is assigned to that value


# Naming contrasts in the contrast matrix
# Mlowhi = contrast (and model coefficient) for the difference in expression between low and high temp expression in Maine
# Plowhi = contrast for the difference in expression between low and high temperature in Panama, and
# Diff = interaction term for evaluating whether the nature of the difference between low and high temp depends upon geographic location
cont.matrix <- makeContrasts(Mlowhi=maine_low-maine_high,Plowhi=panama_low-panama_high,Diff=(maine_low-maine_high)-(panama_low-panama_high),levels=contrast_design)
cont.matrix
# Output: if Diff is significant for a gene, means temp has a different effect in Maine compared to Panama (difference vary by geography)


# Build a voom object for the contrast design
pdf(here::here("Outputs", "voom_qualityweights__contrast_plot.pdf"), width = 20, height = 10)
vwts_contrast <- voomWithQualityWeights(DGE, design=contrast_design,normalize.method="none", plot=TRUE)
dev.off()


# DE TESTING #

# Run linear model fitting
fit_2factor_contrast <- lmFit(vwts_contrast,contrast_design)

# fit contrast
fit_2factor_contrast <- contrasts.fit(fit_2factor_contrast, cont.matrix)

# Run empirical Bayes procedure
fit_2factor_contrast <- eBayes(fit_2factor_contrast,robust=TRUE)

# Summary table
contrast_summary <- decideTests(fit_2factor_contrast,adjust.method="fdr",p.value = 0.05)
summary(decideTests(fit_2factor_contrast,adjust.method="fdr",p.value = 0.05))
pdf(here::here("Outputs", "DEtest__contrast_summary_venn.pdf"))
vennDiagram(contrast_summary)
dev.off()
# Output:
## Diff (9851): many genes exhibit a temperature response that differs quantitatively 
## between Maine and Panama (statistical power for the individual contrasts is lower than that of the interaction test)
## Mlowhi-Plowhi interaction (1382): gene temp response independent from population


# Create merged results table
Mlowhi_results <- topTable(fit_2factor_contrast, adjust="BH",coef=1, p.value=1, number=Inf ,resort.by="P")
colnames(Mlowhi_results) <- paste("Mlowhi",colnames(Mlowhi_results),sep=".")
Mlowhi_results$geneid <- row.names(Mlowhi_results)
Plowhi_results <- topTable(fit_2factor_contrast, adjust="BH",coef=2, p.value=1, number=Inf ,resort.by="P")
colnames(Plowhi_results) <- paste("Plowhi",colnames(Plowhi_results),sep=".")
Plowhi_results$geneid <- row.names(Plowhi_results)
Diff_results <- topTable(fit_2factor_contrast, adjust="BH",coef=3, p.value=1, number=Inf ,resort.by="P")
colnames(Diff_results) <- paste("Diff",colnames(Diff_results),sep=".")
Diff_results$geneid <- row.names(Diff_results)
results_merge <- merge(Mlowhi_results,Plowhi_results,by=c("geneid"))
results_merge <- as_tibble(merge(results_merge,Diff_results,by=c("geneid")))
head(results_merge)

# saving the results_merge table
# CSV
write.csv(results_merge, here::here("Outputs", "results_merge_contrast.csv"), row.names = TRUE)
# Excel
write_xlsx(results_merge, here::here("Outputs", "results_merge_contrast.xlsx"))


# 7. OTHER QUERY # =============================================================

# For example, which genes are up-reg in high temp conditions in a parallel fashion between Panama and Maine?
# identify genes where there is a non-significant interaction term, log-FC is negative in both regions and significant for both regions, assuming an FDR level of 0.05

parallel_high_upreg_genes <- results_merge %>% filter(Diff.adj.P.Val > 0.05) %>%
  filter(Plowhi.adj.P.Val <=0.05) %>%
  filter(Mlowhi.adj.P.Val <=0.05) %>%
  filter(Plowhi.logFC < 0) %>%
  filter(Mlowhi.logFC < 0)
head(parallel_high_upreg_genes)

# Output: a total of 378 genes fit the pattern of parallel, significant upreg with high temp