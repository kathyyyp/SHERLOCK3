# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"
# .libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")

library("sva")
library("readxl")
library("stringr")
library("limma")
library("rstatix")
library("tibble")
library("ggvenn")
library("ggplot2")
library("ggrepel")
library("ggfortify")
library("EnsDb.Hsapiens.v79")
library("ggpubr")
library("edgeR")
library("DESeq2")
library("tidyverse")
library("PCAtools")

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory #my SHERLOCK3 directory

#Data directory
data.dir <- file.path(main.dir,"data")
if(!exists(data.dir)) dir.create(data.dir, recursive = TRUE)

processed.data.dir <- file.path(data.dir,"processed")
if(!exists(processed.data.dir)) dir.create(processed.data.dir)

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
if(!exists(postQC.data.dir)) dir.create(postQC.data.dir)

#Output directory
output.dir <- file.path(main.dir,"output")
if(!exists(output.dir)) dir.create(output.dir, recursive = TRUE)

#SHERLOCK2 dir
sherlock2.dir <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK2"

setwd(file.path(main.dir))

# ================================================================================== #
# 1. COMBAT-SEQ TO COMBINE SHERLOCK1, SHERLOCK2 AND SHERLOCK3 =================================
# ================================================================================== #
cat("Starting 1. COMBAT-SEQ TO COMBINE SHERLOCK2 AND SHERLOCK3", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#Load in SHERLOCK3 data #66,138 genes
# counts_sk3 <- readRDS(file.path(processed.data.dir, "datawrangling_qc_sk3_only", "counts_sk3_postQC.rds")) #225 samples, 66138 genes
counts_sk3_brush <- readRDS(file.path(data.dir, "processed", "datawrangling_qc_sk3_only","counts_sk3_brush.rds")) #112 brushes

#Load in SHERLOCK2 data #69,972 genes
# counts_sk2 <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts.rds")) #327 samples, 69972 genes
counts_sk2_brush <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts_brush.rds")) #162 brushes

#Load in SHERLOCK1 data #60,683 genes #all are brushes
counts_sk1 <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds")) #169 brush samples, 60683 genes,

#FilterByExpr() R code
#https://rdrr.io/bioc/edgeR/src/R/filterByExpr.R


# Code copied from source code - annotated w comments to explain
# In summary
# FilterByExpr()
# keeps genes that have count-per-million (CPM) above k in n samples,
#k is determined by min.count and by the sample library sizes
# n is determined by the design matrix.
# n is essentially the smallest group sample size or, more generally, the minimum inverse leverage of any fitted value. 
# If all the group sizes are larger than large.n, then this is relaxed slightly, 
# but with n always greater than min.prop of the smallest group size (70% by default).

# In addition, each kept gene is required to have at least min.total.count reads across all the samples.


filterByExpr.default <- function(y, 
                                 design=NULL, 
                                 group=NULL, 
                                 lib.size=NULL, #if NULL it uses colSums(y) (total counts per sample).
                                 min.count=10, #minimum count required in at least some sample/s to consider “expressed” (default = 10).
                                 min.total.count=15, #minimum total count across all samples required (default = 15).
                                 large.n=10, #Number of samples per group that is considered to be “large”.
                                 min.prop=0.7, ...) #Minimum proportion of samples in the smallest group that express the gene.
  
  #	Filter low expressed genes given count matrix
  #	Computes TRUE/FALSE index vector indicating which rows to keep
  #	Gordon Smyth
  #	Created 13 Nov 2017. Last revised 26 Jan 2020.

  #Start of function
  {
  y <- as.matrix(y)
  if(mode(y) != "numeric") stop("y is not a numeric matrix")
  if(is.null(lib.size)) lib.size <- colSums(y) #use colSums(y) as libsize. ie. sum down the columns (total counts per sample)
  
  # 1: Determine minimum effect sample size
  ## a) if no group or design specified (will keep genes that have sufficient counts in at least a minimum number of samples overall)
  if(is.null(group)) {
    if(is.null(design)) {
      message("No group or design set. Assuming all samples belong to one group.")
      MinSampleSize <- ncol(y)   #this is number of samples (112 for sherlock3 brush)
      
  ## b) if design is specified
    } else {
      h <- hat(design)
      MinSampleSize <- 1/max(h)
    }
  } 
  ## c) if group specified (will keep genes that are expressed enough, in enough samples per group)
  else {
    group <- as.factor(group)
    n <- tabulate(group)
    MinSampleSize <- min(n[n > 0L])
  }
  
  
  # this step is to avoid requiring overly large numbers of samples if you have many replicates: i.e., after you pass large.n replicates, you only require a proportion of them.
  if(MinSampleSize > large.n) MinSampleSize <- large.n + (MinSampleSize-large.n)*min.prop #if there are more samples than your specified large.n
  
  #	CPM cutoff
  MedianLibSize <- median(lib.size) #median libsize (ie. total counts per sample)
  CPM.Cutoff <- min.count/MedianLibSize*1e6 #min.count default is 10, CPM.cutoff is the CPM equivalent of 10 counts in a library of median size
  CPM <- cpm(y,lib.size=lib.size) # convert all count sto cpm. eg y[,1] = 301, CPM[,y] = 304/sum(y[,1]) X 1,000,000 = 134.60
  tol <- 1e-14 #tolerance for floating point. (rounding R might cause errors, value '3' might be stored in R as 2.99999999)
  keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol) #TRUE/FALSE for cpm > cpm cutoff. 
  #add up how many samples per gene have TRUE. keep only genes where cpm > cpm cutoff in enough sample. in this case enough = 112 (total no. of samples for sherlock3 brush). Keep only genes that are expressed above CPM cutoff in every sample. 
  # Note: above is very strict, genes that are lowly expressed in even one sample will be dropped
  # If we specified group, then MinSampleSize is the smallest group's size, and it keeps genes that are expressed in enough samples within at least one group.
  
  #	Total count cutoff
  keep.TotalCount <- (rowSums(y) >= min.total.count - tol) #min.total.count default is 15
  
  keep.CPM & keep.TotalCount #keep = TRUE when counts exceed CPM cutoff and total count cutoff
}

## when running the code abve manually on sherlock3 counts
# original y = 66,138 genes for 112 samples
# manual_filt <- y[keep.CPM & keep.TotalCount, ] 
# dim(manual_filt)
#> 19,734 remaining for 112 samples

#run function
DGEL<- DGEList(counts=y) 
keep <- filterByExpr(DGEL) 
DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] 
#note this is the same as counts_sk3_brush[filterByExpr(counts_sk3_brush),]
# dim(DGEL)
# > 19,734
# all(manual_filt == DGEL$counts)
# > TRUE





# --- Dummy RNA-seq count matrix --- to explain filterByExpr()
set.seed(123)

# Realistic gene names (ENSEMBL-like or common gene symbols)
genes <- c("ACTB", "GAPDH", "TP53", "MYC", "CD19", "IL2", "FOXP3", "NANOG", "HBB", "INS")

# Suppose we have 6 samples (3 controls, 3 treated)
samples <- paste0("Sample", 1:6)

# Simulate raw counts (some genes low, some high)
y <- matrix(
  c( # manually tuned so you see variation
    45000, 52000, 47000, 41000, 38000, 42000,   # ACTB: high, housekeeping
    38000, 40000, 42000, 35000, 37000, 36000,   # GAPDH: high, housekeeping
    50, 80, 60, 200, 150, 220,                  # TP53: moderate, DE gene
    300, 250, 400, 100, 120, 80,                # MYC: moderate, DE gene
    10, 8, 12, 0, 0, 5,                         # CD19: expressed in controls only
    0, 0, 5, 90, 100, 110,                      # IL2: expressed in treated only
    0, 0, 0, 0, 0, 0,                           # FOXP3: never expressed
    200, 250, 220, 240, 180, 190,               # NANOG: always low–medium
    5000, 6000, 5500, 5200, 4800, 5000,         # HBB: medium
    0, 20, 0, 0, 10, 0                          # INS: rare
  ),
  nrow = length(genes), byrow = TRUE,
  dimnames = list(genes, samples)
)

# Check
y

DGEL <- DGEList(y, group = c("Control", "Control", "Mild", "Mild", "Severe", "Severe"))
keep <- filterByExpr(DGEL)
y_dgelgroup <- DGEL[keep, ,keep.lib.sizes = FALSE]
y_dgelgroup$counts


# Parameters
design=NULL
group=NULL
lib.size=NULL #if NULL it uses colSums(y) (total counts per sample).
min.count=10 #minimum count required in at least some sample/s to consider “expressed” (default = 10).
min.total.count=15 #minimum total count across all samples required (default = 15).
large.n=10 #Number of samples per group that is considered to be “large”.
min.prop=0.7 #Minimum proportion of samples in the smallest group that express the gene.
  
if(is.null(lib.size)) lib.size <- colSums(y)
MinSampleSize <- ncol(y) #no design or group specified

if(MinSampleSize > large.n) MinSampleSize <- large.n + (MinSampleSize-large.n)*min.prop


MedianLibSize <- median(lib.size) #libsize is colsums (ie. total counts per sample)
CPM.Cutoff <- min.count/MedianLibSize*1e6 #min.count default is 10, CPM.cutoff is the CPM equivalent of 10 counts in a library of median size
CPM <- cpm(y,lib.size=lib.size) # convert all count sto cpm. eg y[,1] = 301, CPM[,y] = 304/sum(y[,1]) X 1,000,000 = 134.60
tol <- 1e-14 #tolerance for floating point. (rounding R might cause errors, value '3' might be stored in R as 2.99999999)
keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol) #TRUE/FALSE for cpm > cpm cutoff. 
#add up how many samples per gene have TRUE. keep only genes where cpm > cpm cutoff in enough sample. in this case enough = 112 (total no. of samples for sherlock3 brush). Keep only genes that are expressed above CPM cutoff in every sample. 
# Note: above is very strict, genes that are lowly expressed in even one sample will be dropped
# If we specified group, then MinSampleSize is the smallest group's size, and it keeps genes that are expressed in enough samples within at least one group.

#	Total count cutoff
keep.TotalCount <- (rowSums(y) >= min.total.count - tol) #min.total.count default is 15

keep.CPM & keep.TotalCount #keep = TRUE when counts exceed CPM cutoff and total count cutoff

y[keep.CPM & keep.TotalCount,]



# Method 1)Get common genes, filter, get common genes again, integrate --------------------------------------------------------------------------------------------------------------------------
#SHERLOCK1, 2 and 3
common_genes <- intersect(
  intersect(row.names(counts_sk1),
            row.names(counts_sk2_brush)),
  row.names(counts_sk3_brush)) #48,039 common genes

# SHERLOCK1 = 169 samples
counts_sk1_common <- counts_sk1[common_genes,]
counts_sk1_keep <- filterByExpr(counts_sk1_common) 
counts_sk1_filtered <- counts_sk1_common[counts_sk1_keep, ] 
#60,683 genes --> #18,144 genes


# SHERLOCK2 Brush = 162 samples
counts_sk2_brush_common <- counts_sk2_brush[common_genes,] 
counts_sk2_brush_keep <- filterByExpr(counts_sk2_brush_common)
counts_sk2_brush_filtered <- counts_sk2_brush_common[counts_sk2_brush_keep, ] 
#69,972 genes -->  #15,649 genes


# SHERLOCK3 Brush = 112 samples
counts_sk3_brush_common <- counts_sk3_brush[common_genes,] 
counts_sk3_brush_keep <- filterByExpr(counts_sk3_brush_common)
counts_sk3_brush_filtered <- counts_sk3_brush_common[counts_sk3_brush_keep, ] 
#66,138 genes  --> #18,116 genes 

### COMBATSEQ ----------------------------------------------------------------------------------------------------------------------
counts_sk1_sk2_all <- merge(counts_sk1_filtered, counts_sk2_brush_filtered, by = "row.names", all = FALSE) # all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
row.names(counts_sk1_sk2_all) <- counts_sk1_sk2_all[,1]
counts_sk1_sk2_all <- as.matrix(counts_sk1_sk2_all[,-1])

# Merge the sk1_sk2 with sk3 (keep only the matched Ensembl IDS)
counts_merged <- merge(counts_sk1_sk2_all,  counts_sk3_brush_filtered,  by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1]) #15,286 genes total 

medina.data.dir <- file.path(processed.data.dir, "datawrangling_qc_allbrush_medina")
if(!exists(medina.data.dir))dir.create(medina.data.dir, recursive = TRUE)

combat.processed.data.dir <- file.path(medina.data.dir, "combat_results_filterbyexprbeforecombat")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir, recursive = TRUE)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.csv")) #48,039 genes and 721 samples
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.rds"))

#combat
batch <- c(rep(1, ncol(counts_sk1_filtered)), rep(2, ncol(counts_sk2_brush_filtered)), rep(3, ncol(counts_sk3_brush_filtered)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.rds"))



# plot settings
fontsize <- theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14)
)

counts_merged <- readRDS(file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.rds"))
counts_combat <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.rds"))
#PCA plot - PC1 AND 2 -------------------------------------------------------------------------------------------
clinical_brush <- readRDS(file.path(medina.data.dir, "clinical_sherlock_brush_simple.rds"))
clinical_subset <- clinical_brush[colnames(counts_merged),]
clinical_subset$batch <- as.factor(clinical_subset$batch )


qc.dir <- file.path(output.dir, "qc", "medina")
filterbyexpr.dir <- file.path(qc.dir, "filterbyExpr")
if(!exists(filterbyexpr.dir)) dir.create(filterbyexpr.dir)

#PRE COMBAT  -------------------------------------------------------------------------------------------
#PCA plot - individual
variable = "batch"
counts_norm <- voom(counts_merged) #voom is log2cpm values
counts_norm <- counts_norm$E
pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
ggsave(autoplot(pca_res,data=clinical_subset, colour = variable)+ fontsize + guides(color=guide_legend(variable)),
       file = file.path(filterbyexpr.dir,paste0("pca_batch_precombat_filtered.png")),
       width = 2300,
       height = 1500,
       units = "px")

#PCA plot - pairs
png(filename = file.path(filterbyexpr.dir,paste0("pca_pairs_batch_precombat_filtered.png")), width = 26, height = 23, units = "cm", res = 1200)
clinical_pca <- cbind(pca_res$x, clinical_subset)

variable <- as.factor(clinical_pca[,"batch"])

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,3,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "batch")
dev.off()

#POST COMBAT -------------------------------------------------------------------------------------------
#PCA plot - individual
variable = "batch"
counts_norm <- voom(counts_combat) #voom is log2cpm values
counts_norm <- counts_norm$E
pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
ggsave(autoplot(pca_res,data=clinical_subset, colour = variable)+ fontsize + guides(color=guide_legend(variable)),
       file = file.path(filterbyexpr.dir,paste0("pca_batch_postcombat_filtered.png")),
       width = 2300,
       height = 1500,
       units = "px")

#PCA plot - pairs
png(filename = file.path(filterbyexpr.dir,paste0("pca_pairs_batch_postcombat_filtered.png")), width = 26, height = 23, units = "cm", res = 1200)
clinical_pca <- cbind(pca_res$x, clinical_subset)

variable <- as.factor(clinical_pca[,"batch"])

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,3,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "batch")
dev.off()


# Method 2) Filter, get common genes, integrate --------------------------------------------------------------------------------------------------------------------------

# SHERLOCK1 = 169 samples
counts_sk1_keep <- filterByExpr(counts_sk1) 
counts_sk1_filtered <- counts_sk1[counts_sk1_keep, ] 
#60,683 genes --> #18,732 genes


# SHERLOCK2 Brush = 162 samples
counts_sk2_brush_keep <- filterByExpr(counts_sk2_brush)
counts_sk2_brush_filtered <- counts_sk2_brush[counts_sk2_brush_keep, ] 
#69,972 genes -->  #16,612 genes


# SHERLOCK3 Brush = 112 samples
counts_sk3_brush_keep <- filterByExpr(counts_sk3_brush)
counts_sk3_brush_filtered <- counts_sk3_brush[counts_sk3_brush_keep, ] 
#66,138 genes  --> #19,734 genes 

### COMBATSEQ ----------------------------------------------------------------------------------------------------------------------
counts_sk1_sk2_all <- merge(counts_sk1_filtered, counts_sk2_brush_filtered, by = "row.names", all = FALSE) # all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
row.names(counts_sk1_sk2_all) <- counts_sk1_sk2_all[,1]
counts_sk1_sk2_all <- as.matrix(counts_sk1_sk2_all[,-1])

# Merge the sk1_sk2 with sk3 (keep only the matched Ensembl IDS)
counts_merged <- merge(counts_sk1_sk2_all,  counts_sk3_brush_filtered,  by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1]) #15,292 genes total 

medina.data.dir <- file.path(processed.data.dir, "datawrangling_qc_allbrush_medina")
if(!exists(medina.data.dir))dir.create(medina.data.dir, recursive = TRUE)

combat.processed.data.dir <- file.path(medina.data.dir, "combat_results_filterbyexprbeforecombat")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir, recursive = TRUE)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr_method2.csv")) #48,039 genes and 721 samples
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr_method2.rds"))

#combat
batch <- c(rep(1, ncol(counts_sk1_filtered)), rep(2, ncol(counts_sk2_brush_filtered)), rep(3, ncol(counts_sk3_brush_filtered)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr_method2.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr_method2.rds"))


#PCA plot - PC1 AND 2 -------------------------------------------------------------------------------------------
clinical_brush <- readRDS(file.path(medina.data.dir, "clinical_sherlock_brush_simple.rds"))
clinical_subset <- clinical_brush[colnames(counts_merged),]
clinical_subset$batch <- as.factor(clinical_subset$batch )

variable = "batch"

qc.dir <- file.path(output.dir, "qc", "medina")
filterbyexpr.dir <- file.path(qc.dir, "filterbyExpr")
if(!exists(filterbyexpr.dir)) dir.create(filterbyexpr.dir)

#PRE COMBAT  -------------------------------------------------------------------------------------------
#PCA plot - individual
counts_norm <- voom(counts_merged) #voom is log2cpm values
counts_norm <- counts_norm$E
pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
ggsave(autoplot(pca_res,data=clinical_subset, colour = variable)+ fontsize + guides(color=guide_legend(variable)),
       file = file.path(filterbyexpr.dir,paste0("pca_batch_precombat_filtered_method2.png")),
       width = 2300,
       height = 1500,
       units = "px")

#PCA plot - pairs
png(filename = file.path(filterbyexpr.dir,paste0("pca_pairs_batch_precombat_filtered_method2.png")), width = 26, height = 23, units = "cm", res = 1200)
clinical_pca <- cbind(pca_res$x, clinical_subset)

variable <- as.factor(clinical_pca[,"batch"])

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,3,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "batch")
dev.off()

#POST COMBAT -------------------------------------------------------------------------------------------
#PCA plot - individual
counts_norm <- voom(counts_combat) #voom is log2cpm values
counts_norm <- counts_norm$E
pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
ggsave(autoplot(pca_res,data=clinical_subset, colour = variable)+ fontsize + guides(color=guide_legend(variable)),
       file = file.path(filterbyexpr.dir,paste0("pca_batch_postcombat_filtered_method2.png")),
       width = 2300,
       height = 1500,
       units = "px")

#PCA plot - pairs
png(filename = file.path(filterbyexpr.dir,paste0("pca_pairs_batch_postcombat_filtered_method2.png")), width = 26, height = 23, units = "cm", res = 1200)
clinical_pca <- cbind(pca_res$x, clinical_subset)

variable <- as.factor(clinical_pca[,"batch"])

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,3,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "batch")
dev.off()






dim(counts_sk1_filtered)
dim(counts_sk2_brush_filtered)
dim(counts_sk3_brush_filtered)

# Method 3) Filter with grouping, get common genes, integrate --------------------------------------------------------------------------------------------------------------------------
clinical_brush <- readRDS(file.path(medina.data.dir, "clinical_sherlock_brush_simple.rds"))
clinical_brush[colnames(counts_sk1), "classification"]

# SHERLOCK1 = 169 samples
counts_sk1_keep <- filterByExpr(counts_sk1, group = clinical_brush[colnames(counts_sk1), "classification"]) 
counts_sk1_filtered <- counts_sk1[counts_sk1_keep, ] 
#60,683 genes --> #23,7134 genes


# SHERLOCK2 Brush = 162 samples
counts_sk2_brush_keep <- filterByExpr(counts_sk2_brush, group = clinical_brush[colnames(counts_sk2_brush), "classification"] )
counts_sk2_brush_filtered <- counts_sk2_brush[counts_sk2_brush_keep, ] 
#69,972 genes -->  #20,675 genes


# SHERLOCK3 Brush = 112 samples
counts_sk3_brush_keep <- filterByExpr(counts_sk3_brush, group = clinical_brush[colnames(counts_sk3_brush), "classification"])
counts_sk3_brush_filtered <- counts_sk3_brush[counts_sk3_brush_keep, ] 
#66,138 genes  --> #24,038 genes 

### COMBATSEQ ----------------------------------------------------------------------------------------------------------------------
counts_sk1_sk2_all <- merge(counts_sk1_filtered, counts_sk2_brush_filtered, by = "row.names", all = FALSE) # all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
row.names(counts_sk1_sk2_all) <- counts_sk1_sk2_all[,1]
counts_sk1_sk2_all <- as.matrix(counts_sk1_sk2_all[,-1])

# Merge the sk1_sk2 with sk3 (keep only the matched Ensembl IDS)
counts_merged <- merge(counts_sk1_sk2_all,  counts_sk3_brush_filtered,  by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1]) #18,491 genes total 

medina.data.dir <- file.path(processed.data.dir, "datawrangling_qc_allbrush_medina")
if(!exists(medina.data.dir))dir.create(medina.data.dir, recursive = TRUE)

combat.processed.data.dir <- file.path(medina.data.dir, "combat_results_filterbyexprbeforecombat")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir, recursive = TRUE)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr_method3.csv")) #48,039 genes and 721 samples
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr_method3.rds"))

#combat
batch <- c(rep(1, ncol(counts_sk1_filtered)), rep(2, ncol(counts_sk2_brush_filtered)), rep(3, ncol(counts_sk3_brush_filtered)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr_method3.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr_method3.rds"))


#PCA plot - PC1 AND 2 -------------------------------------------------------------------------------------------
clinical_brush <- readRDS(file.path(medina.data.dir, "clinical_sherlock_brush_simple.rds"))
clinical_subset <- clinical_brush[colnames(counts_merged),]
clinical_subset$batch <- as.factor(clinical_subset$batch )

variable = "batch"

qc.dir <- file.path(output.dir, "qc", "medina")
filterbyexpr.dir <- file.path(qc.dir, "filterbyExpr")
if(!exists(filterbyexpr.dir)) dir.create(filterbyexpr.dir)

#PRE COMBAT  -------------------------------------------------------------------------------------------
#PCA plot - individual
counts_norm <- voom(counts_merged) #voom is log2cpm values
counts_norm <- counts_norm$E
pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
ggsave(autoplot(pca_res,data=clinical_subset, colour = variable)+ fontsize + guides(color=guide_legend(variable)),
       file = file.path(filterbyexpr.dir,paste0("pca_batch_precombat_filtered_method3.png")),
       width = 2300,
       height = 1500,
       units = "px")

#PCA plot - pairs
png(filename = file.path(filterbyexpr.dir,paste0("pca_pairs_batch_precombat_filtered_method3.png")), width = 26, height = 23, units = "cm", res = 1200)
clinical_pca <- cbind(pca_res$x, clinical_subset)

variable <- as.factor(clinical_pca[,"batch"])

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,3,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "batch")
dev.off()

#POST COMBAT -------------------------------------------------------------------------------------------
#PCA plot - individual
counts_norm <- voom(counts_combat) #voom is log2cpm values
counts_norm <- counts_norm$E
pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
ggsave(autoplot(pca_res,data=clinical_subset, colour = variable)+ fontsize + guides(color=guide_legend(variable)),
       file = file.path(filterbyexpr.dir,paste0("pca_batch_postcombat_filtered_method3.png")),
       width = 2300,
       height = 1500,
       units = "px")

#PCA plot - pairs
png(filename = file.path(filterbyexpr.dir,paste0("pca_pairs_batch_postcombat_filtered_method3.png")), width = 26, height = 23, units = "cm", res = 1200)
clinical_pca <- cbind(pca_res$x, clinical_subset)

variable <- as.factor(clinical_pca[,"batch"])

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,3,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "batch")
dev.off()

