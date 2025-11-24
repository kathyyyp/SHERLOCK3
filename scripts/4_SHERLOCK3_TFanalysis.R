#SHERLOCK3 - TF Analysis

options(error = function() { traceback(2); quit(status = 1) })
# Infer TF activity using bulk RNASeq data

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

# Load extra functions
source(file.path(my_directory, "functions","process-utils.R"))
source(file.path(my_directory, "functions","cluster-functions.R"))
source(file.path(my_directory, "functions","viper-utils.R"))

# Load packages
library(Matrix)
library(dplyr)
# library(DCGL) #?
library(scales)
library(viper) #viper 1.34
library(PISCES)
library(gplots)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(limma)
library(ggvenn)
library(edgeR)

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir, "data")
processed.data.dir <- file.path(data.dir,"processed")
postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir, "output")

#TF output directory
tf.dir <- file.path(output.dir, "tf_analysis")
if(!exists(tf.dir)) dir.create(tf.dir)

brush.tf.dir <- file.path(tf.dir, "brush")
if(!exists(brush.tf.dir)) dir.create(brush.tf.dir)

biopt.tf.dir <- file.path(tf.dir, "biopt")
if(!exists(biopt.tf.dir)) dir.create(biopt.tf.dir)


# ================================================================================== #
# C. SET UP FUNCTIONS ===============================================
# ================================================================================== #
# Created Ranking and CPM transformation functions
RankTransform <- function(dat.mat) {
  rank.mat <- apply(dat.mat, 2, rank)
  median <- apply(rank.mat, 1, median)
  mad <- apply(rank.mat, 1, mad)
  rank.mat <- (rank.mat - median) / mad
  return(rank.mat)
}

CPM_Transform <- function(dat.mat, l2 = FALSE) {
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  if (l2) {
    cpm.mat <- log2(cpm.mat + 1)
  }
  return(cpm.mat)
}


# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
# ##-- Post batch correction
all_counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))

all_clinical<- readRDS(file.path(postQC.data.dir, "clinical_brushbiopt_simple.rds")) #552 samples
clinical_brush <- readRDS(file.path(postQC.data.dir, "clinical_brush_simple.rds"))
clinical_biopt <- readRDS(file.path(postQC.data.dir, "clinical_biopt_simple.rds"))

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


tf_analysis_func <- function(sampletype, this.tf.dir){
  
#TF output directory BRUSH
  if(sampletype == "brush"){
    
    clinical <- all_clinical[which(all_clinical$sampletype == "Brush"),]
    counts<- all_counts[, row.names(clinical)]
  }
  
  #TF output directory BIOPT
  if(sampletype == "biopt"){
    
    clinical<- all_clinical[which(all_clinical$sampletype == "Biopt"),]
    counts <- all_counts[, row.names(clinical)]
  }
    
    

  tf.results.dir <- file.path(this.tf.dir, "results")
  if(!exists(tf.results.dir)) dir.create(tf.results.dir)
    
  tf.figures.dir <- file.path(this.tf.dir, "figures")
  if(!exists(tf.figures.dir)) dir.create(tf.figures.dir)
    

my_ensembl_gene_ids <- rownames(counts) #

library(EnsDb.Hsapiens.v79)
##### #with library(EnsDb.Hsapiens.v79) ##### 
hgnc_symbols_db <- ensembldb::select(EnsDb.Hsapiens.v79, keys= my_ensembl_gene_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

row.names(hgnc_symbols_db) <- hgnc_symbols_db$GENEID

n_occur <- data.frame(table(hgnc_symbols_db$SYMBOL)) #Get frequencies of each gene/hgnc_symbol
n_occur[n_occur$Freq > 1,]

counts_mapped <- merge(hgnc_symbols_db, data.frame(counts), by.x = "GENEID", by.y = "row.names")
counts_mapped[,-c(1,2)] <- sapply(counts_mapped[, -c(1,2)], as.numeric)
rownames(counts_mapped) <- make.unique(counts_mapped$SYMBOL)
counts_mapped <- counts_mapped[,-c(1,2)]
colnames(counts_mapped) <- colnames(counts)


# ================================================================================== #
# 2. NORMALISATION AND RANKING =====================================================
# ================================================================================== #
cat("Starting 2. NORMALISATION AND RANKING", Sys.time(), "\n")

raw.mat = as.matrix(counts_mapped)
cpm.mat <- CPM_Transform(raw.mat);nrow(cpm.mat) #Apply edgeR's Counts Per Million transformation to raw counts


sample.names <- colnames(cpm.mat)
gene.ids <- rownames(cpm.mat)
m <- cpm.mat
mm <- rbind(c("gene", sample.names), cbind(gene.ids, m))
write.table(x = mm, file = file.path(this.tf.dir,"counts_cpm.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = F)


# }# stop function here when saving cpm counts files to make the network file

rank.mat <- RankTransform(cpm.mat);nrow(rank.mat)

#run script to make network file

cat("Load TF list and network", Sys.time(), "\n")
# Load TF list, could be human or mouse [should have only one column  without header]
TF_list = read.table(file.path(tf.dir, "humanTF.tsv"), sep = "\t")

network.dir <- file.path(this.tf.dir, "network_res")

network<- read.table(file.path(network.dir,"network.txt"))
network <- network[-1,]

cat("Save network as .tsv", Sys.time(), "\n")
write.table(network[,c(1:3)], file.path(this.tf.dir,"network.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names  = FALSE)

cat("Reload network.tsv", Sys.time(), "\n")
network <-  read.table(file.path(this.tf.dir,"network.tsv"))
colnames(network) <- c("Regulator", "Target", "MI")
head(network)
str(network)


# ================================================================================== #
# 3. NETWORK GENERATION ============================================================
# ================================================================================== #
cat("Starting 3. NETWORK GENERATION", Sys.time(), "\n")

network.file.path <- file.path(this.tf.dir,"network.tsv")

# #Processes ARACNe results into a regulon object compatible with VIPER.
# RegProcess(network.file.path, cpm.mat, out.dir = tf.data.dir, out.name = 'data_r1-net-')
# 
# #Loads the final processed regulatory network 
# r1.net <- readRDS(file.path(tf.data.dir,'data_r1-net-pruned.rds'))

# Processes the ARACNe-generated network (network.tsv) and applies it to cpm.mat
r1.net <- aracne2regulon(network.file.path, cpm.mat, format = "3col") #3col means the input has 3 columns
saveRDS(r1.net, file = file.path(this.tf.dir, "data_r1-net.rds"))

# ================================================================================== #
# 4. INFER TF ACTIVITY =============================================================
# ================================================================================== #
cat("Starting 4. INFER TF ACTIVITY", Sys.time(), "\n")

#This function performs Virtual Inference of Protein-activity by Enriched Regulon analysis
#VIPER transforms gene expression data into TF activity score
# output of viper is a DF where rows = TFs, columns = samples and values are estimated TF activity scores
#This means that for each TF, VIPER computes an activity score across all samples, based on how its target genes behave in rank.mat
r1.pAct <- viper(rank.mat, r1.net, method = 'none') 

# Computes similarity between TF activities using VIPER similarity metric, converts it into a distance mtric
r1.viperDist <- as.dist(viperSimilarity(r1.pAct))

#Clusters TF activity profiles using partitioning around medoids (PAM).
#not sure why this part is even needed - r1.cluts is not used in subsequent steps??
# r1.clusts <- PamKRange(r1.viperDist, kmin = 2, kmax = 10)
library(fpc)
r1.clusts <- pamk(r1.viperDist, krange = 2:10) #tests different numbers of clusters (k) to test the best

# ================================================================================== #
# 5. IDENTIFY MOST REPRESENTATIVE PROTEINS =========================================
# ================================================================================== #
cat("Starting 5. IDENTIFY MOST REPRESENTATIVE PROTEINS ", Sys.time(), "\n")

r2.cbcMRs <- CBCMRs(r1.pAct) #a vector of proteins


# ================================================================================== #
# 6.  "FILTER THE PROTEIN ACTIVITY MATRIX" =========================================
# ================================================================================== #
cat("Starting 6.  FILTER THE PROTEIN ACTIVITY MATRIX", Sys.time(), "\n")

r2.pAct.cbc <- r1.pAct[ r2.cbcMRs ,] #subset to only include the most representative proteins
head(r2.pAct.cbc)
nrow(r2.pAct.cbc)
TFs_data <- as.data.frame(r2.pAct.cbc)
TFs_names <- unique(rownames(TFs_data))

#Save output: Activity scores for each TF
#If a TF activates a gene, and that gene is highly expressed, the TF is likely active.
#################################################
write.csv(TFs_data, file.path(tf.results.dir, "TF_analysis_results.csv"))
write.csv(TFs_names, file.path(tf.results.dir, "TF_analysis_TFnames.csv"))
saveRDS(TFs_data, file.path(tf.results.dir,"TF_analysis_results.rds"))

# Differential expression +----------------------------------------------------------------------------------------------------------
# TF_analysis_results <- readRDS(file.path("TF_analysis", "TF_analysis_results.rds"))
# TF_analysis_results <- as.matrix(TF_analysis_results)

# TF_analysis_results <-readRDS(file.path(output.dir,"tf_analysis","sherlock1","results","TF_analysis_results.rds"))
TF_analysis_results <- as.matrix(TFs_data)

## Normality Tests ----------------------------------------------------------------------------------------------------------

if(!exists(file.path(tf.figures.dir, "normality"))) dir.create(file.path(tf.figures.dir, "normality"))
pdf(file= file.path(tf.figures.dir, "normality", "normality_tests.pdf"), width = 8, height = 8)

par(mar = c(7, 4, 4, 2))  
hist(TF_analysis_results) 
# Add label of ks.test result
ktest <- ks.test(TF_analysis_results, "pnorm", mean(TF_analysis_results), sd(TF_analysis_results))
mtext(paste(ktest$method, ktest$p.value), side = 1, line = 4, cex = 1)

# Add label of ad.test result
library("nortest")
adtest <- ad.test(TF_analysis_results)
mtext(paste(adtest$method, adtest$p.value), side = 1, line = 5.5, cex = 1)

qqnorm(TF_analysis_results)
# qqline(TF_analysis_results, col="red")

plot(density(TF_analysis_results))

dev.off()


#### Limma  -----------------------------------------------------------------------------------------------------------------

# dir.create(result.dir, "tT")
# dir.create(result.dir, "tT2")

design <- model.matrix(
  ~0 + classification + age + sex + smoking_status,
  data = clinical)  # clinical$classificationation is your grouping factor

colnames(design)[1:3] = c("Control", "Mild.moderate.COPD", "Severe.COPD")

# Fit the linear model to TF activity scores
fit <- lmFit(TF_analysis_results, 
             design)


cont.matrix <- makeContrasts(
  test1 = Mild.moderate.COPD - Control,
  test2 = Severe.COPD - Control,
  test3 = Severe.COPD - Mild.moderate.COPD,
  levels = design) 

# Apply empirical Bayes moderation
fit2_1 <- contrasts.fit(fit, contrast=cont.matrix[,'test1'])
fit2_2 <- contrasts.fit(fit, contrast=cont.matrix[,'test2'])
fit2_3 <- contrasts.fit(fit, contrast=cont.matrix[,'test3'])

fit2_1 <- eBayes(fit2_1)
fit2_2 <- eBayes(fit2_2)
fit2_3 <- eBayes(fit2_3)


#topTableF ranks cpg sites on the basis of moderated F-statistics for one or more coefficients. FDR calculated using Benjamin-Hochberg method
tT_1 <- topTable(fit2_1, adjust="BH", sort.by="P", number=nrow(fit2_1))
tT_2 <- topTable(fit2_2, adjust="BH", sort.by="P", number=nrow(fit2_2))
tT_3 <- topTable(fit2_3, adjust="BH", sort.by="P", number=nrow(fit2_3))


listoftT <- list("Mild.moderate.COPDvsControl" = tT_1, 
                 "Severe.COPDvsControl" = tT_2, 
                 "Severe.COPDvsMild.moderate.COPD" = tT_3)

listoftT2 <- list()
listofvolcano <- list()

for (i in 1:length(listoftT)){
tT <- listoftT[[i]]
tT$Legend <- ifelse(
  tT$adj.P.Val < 0.05 & tT$logFC > 0, "Upregulated", 
  ifelse(
    tT$adj.P.Val < 0.05 & tT$logFC < -0, "Downregulated",
    "Not Significant"))

tT$genename <- row.names(tT)
selection <- which(tT$adj.P.Val<0.05)

tT2 <- tT[selection,]

listoftT2[[names(listoftT[i])]] <- tT2

write.csv(tT2, file.path(tf.results.dir, paste0("tT2_", names(listoftT[i]), ".csv")))


volcano <- ggplot(tT, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Legend)) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"))+
  geom_hline(yintercept =-log10(max(tT2$P.Value)),colour="black", linetype="dashed")+
  geom_text_repel(data = subset(tT2[1:20,]),
                  aes(label= genename),size = 4, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines") ) +
  theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                   legend.text = element_text(size = 14),
                                   legend.title = element_text(size = 16)) +
  labs(title = names(listoftT)[i])

listofvolcano[[i]] <- volcano

ggsave(volcano, filename = file.path(tf.figures.dir, paste0("volcano_", names(listoftT[i]), ".png")),
       width = 22,
       height = 22,
       units = "cm")

}

saveRDS(listoftT, file.path(tf.results.dir, "listoftT.rds"))
saveRDS(listoftT2, file.path(tf.results.dir, "listoftT2.rds"))

} #close function

tf_analysis_func(sampletype = "brush", this.tf.dir = brush.tf.dir)
tf_analysis_func(sampletype = "biopt", this.tf.dir = biopt.tf.dir)

# 
# 
# kruskal.test(TF_analysis_results ~ group_vector)
# pairwise.wilcox.test(TF_activity_vector, group_vector, p.adjust.method="BH")
# 
# 
# 
# wilcox_test(TF_analysis_results,
#             p.adjust.method = BH)


# ================================================================================== #
# 7.  VALIDATION WITH SHERLOCK1 ====================================================
# ================================================================================== #
cat("Starting 7.  VALIDATION WITH SHERLOCK1", Sys.time(), "\n")

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
  
  validation.tf.dir <- file.path(output.dir, "tf_analysis", "sherlock1")
  if(!exists(validation.tf.dir)) dir.create(validation.tf.dir, recursive = TRUE)


counts <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
clinical <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_simple.rds"))

tf.results.dir <- file.path(validation.tf.dir, "results")
if(!exists(tf.results.dir)) dir.create(tf.results.dir)

tf.figures.dir <- file.path(validation.tf.dir, "figures")
if(!exists(tf.figures.dir)) dir.create(tf.figures.dir)


my_ensembl_gene_ids <- rownames(counts) #

library(EnsDb.Hsapiens.v79)
##### #with library(EnsDb.Hsapiens.v79) ##### 
hgnc_symbols_db <- ensembldb::select(EnsDb.Hsapiens.v79, keys= my_ensembl_gene_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

row.names(hgnc_symbols_db) <- hgnc_symbols_db$GENEID

n_occur <- data.frame(table(hgnc_symbols_db$SYMBOL)) #Get frequencies of each gene/hgnc_symbol
n_occur[n_occur$Freq > 1,]

counts_mapped <- merge(hgnc_symbols_db, data.frame(counts), by.x = "GENEID", by.y = "row.names")
counts_mapped[,-c(1,2)] <- sapply(counts_mapped[, -c(1,2)], as.numeric)
rownames(counts_mapped) <- make.unique(counts_mapped$SYMBOL)
counts_mapped <- counts_mapped[,-c(1,2)]
colnames(counts_mapped) <- colnames(counts)


# ================================================================================== #
# 2. NORMALISATION AND RANKING =====================================================
# ================================================================================== #
cat("Starting 2. NORMALISATION AND RANKING", Sys.time(), "\n")


raw.mat = as.matrix(counts_mapped)
cpm.mat <- CPM_Transform(raw.mat);nrow(cpm.mat) #Apply edgeR's Counts Per Million transformation to raw counts


sample.names <- colnames(cpm.mat)
gene.ids <- rownames(cpm.mat)
m <- cpm.mat
mm <- rbind(c("gene", sample.names), cbind(gene.ids, m))

write.table(x = mm, file = file.path(validation.tf.dir,"counts_cpm.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = F)

#run script to make network file

rank.mat <- RankTransform(cpm.mat);nrow(rank.mat)

cat("Load TF list and network", Sys.time(), "\n")
# Load TF list, could be human or mouse [should have only one column  without header]
TF_list = read.table(file.path(tf.dir, "humanTF.tsv"), sep = "\t")

network.dir <- file.path(validation.tf.dir, "network_res")

network<- read.table(file.path(network.dir,"network.txt"))
network <- network[-1,]

# network$Regulator <- toupper(network$Regulator) ###NOTE some count row names also aren't capitalised
# network$Target <- toupper(network$Target)

cat("Save network as .tsv", Sys.time(), "\n")
write.table(network[,c(1:3)], file.path(validation.tf.dir,"network.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names  = FALSE)

cat("Reload network.tsv", Sys.time(), "\n")
network <-  read.table(file.path(validation.tf.dir,"network.tsv"))
head(network)
str(network)


# ================================================================================== #
# 3. NETWORK GENERATION ============================================================
# ================================================================================== #
cat("Starting 3. NETWORK GENERATION", Sys.time(), "\n")

network.file.path <- file.path(validation.tf.dir,"network.tsv")

# #Processes ARACNe results into a regulon object compatible with VIPER.
# RegProcess(network.file.path, cpm.mat, out.dir = tf.data.dir, out.name = 'data_r1-net-')
#
# #Loads the final processed regulatory network
# r1.net <- readRDS(file.path(tf.data.dir,'data_r1-net-pruned.rds'))

# Processes the ARACNe-generated network (network.tsv) and applies it to cpm.mat
r1.net <- aracne2regulon(network.file.path, cpm.mat, format = "3col") #3col means the input has 3 columns
saveRDS(r1.net, file = file.path(validation.tf.dir, "data_r1-net.rds"))

# ================================================================================== #
# 4. INFER TF ACTIVITY =============================================================
# ================================================================================== #
cat("Starting 4. INFER TF ACTIVITY", Sys.time(), "\n")

#This function performs Virtual Inference of Protein-activity by Enriched Regulon analysis
#VIPER transforms gene expression data into TF activity score
# output of viper is a DF where rows = TFs, columns = samples and values are estimated TF activity scores
#This means that for each TF, VIPER computes an activity score across all samples, based on how its target genes behave in rank.mat
r1.pAct <- viper(rank.mat, r1.net, method = 'none')

# Computes similarity between TF activities using VIPER similarity metric, converts it into a distance mtric
r1.viperDist <- as.dist(viperSimilarity(r1.pAct))

#Clusters TF activity profiles using partitioning around medoids (PAM).
#not sure why this part is even needed - r1.cluts is not used in subsequent steps??
# r1.clusts <- PamKRange(r1.viperDist, kmin = 2, kmax = 10)
library(fpc)
r1.clusts <- pamk(r1.viperDist, krange = 2:10) #tests different numbers of clusters (k) to test the best

# ================================================================================== #
# 5. IDENTIFY MOST REPRESENTATIVE PROTEINS =========================================
# ================================================================================== #
cat("Starting 5. IDENTIFY MOST REPRESENTATIVE PROTEINS ", Sys.time(), "\n")

r2.cbcMRs <- CBCMRs(r1.pAct) #a vector of proteins


# ================================================================================== #
# 6.  "FILTER THE PROTEIN ACTIVITY MATRIX" =========================================
# ================================================================================== #
cat("Starting 6.  FILTER THE PROTEIN ACTIVITY MATRIX", Sys.time(), "\n")

r2.pAct.cbc <- r1.pAct[ r2.cbcMRs ,] #subset to only include the most representative proteins
head(r2.pAct.cbc)
nrow(r2.pAct.cbc)
TFs_data <- as.data.frame(r2.pAct.cbc)
TFs_names <- unique(rownames(TFs_data))

#Save output: Activity scores for each TF in each sample (normalised enrichment-like scores)
#If a TF activates a gene, and that gene is highly expressed, the TF is likely active.
# positive = activated
# negative = repressed 
# 0 = no strong signal
#################################################
write.csv(TFs_data, file.path(tf.results.dir, "TF_analysis_results.csv"))
write.csv(TFs_names, file.path(tf.results.dir, "TF_analysis_TFnames.csv"))
saveRDS(TFs_data, file.path(tf.results.dir,"TF_analysis_results.rds"))

# Differential expression +----------------------------------------------------------------------------------------------------------
# TF_analysis_results <- readRDS(file.path("TF_analysis", "TF_analysis_results.rds"))
# TF_analysis_results <- as.matrix(TF_analysis_results)

TF_analysis_results <- as.matrix(TFs_data)

## Normality Tests ----------------------------------------------------------------------------------------------------------

if(!exists(file.path(tf.figures.dir, "normality"))) dir.create(file.path(tf.figures.dir, "normality"))
pdf(file= file.path(tf.figures.dir, "normality", "normality_tests.pdf"), width = 8, height = 8)

par(mar = c(7, 4, 4, 2))
hist(TF_analysis_results)
# Add label of ks.test result
ktest <- ks.test(TF_analysis_results, "pnorm", mean(TF_analysis_results), sd(TF_analysis_results))
mtext(paste(ktest$method, ktest$p.value), side = 1, line = 4, cex = 1)

# Add label of ad.test result
library("nortest")
adtest <- ad.test(TF_analysis_results)
mtext(paste(adtest$method, adtest$p.value), side = 1, line = 5.5, cex = 1)

qqnorm(TF_analysis_results)
# qqline(TF_analysis_results, col="red")

plot(density(TF_analysis_results))

dev.off()


#### Limma  -----------------------------------------------------------------------------------------------------------------

# dir.create(result.dir, "tT")
# dir.create(result.dir, "tT2")


design <- model.matrix(
  ~0 + classification + age + sex + packyears, ##all sherlock1 are ex smokers
  data = clinical)  # clinical$classificationation is your grouping factor

colnames(design)[1:3] = c("Control", "Mild.moderate.COPD", "Severe.COPD")

# Fit the linear model to TF activity scores
fit <- lmFit(TF_analysis_results,
             design)


cont.matrix <- makeContrasts(
  test1 = Mild.moderate.COPD - Control,
  test2 = Severe.COPD - Control,
  test3 = Severe.COPD - Mild.moderate.COPD,
  levels = design)

# Apply empirical Bayes moderation
fit2_1 <- contrasts.fit(fit, contrast=cont.matrix[,'test1'])
fit2_2 <- contrasts.fit(fit, contrast=cont.matrix[,'test2'])
fit2_3 <- contrasts.fit(fit, contrast=cont.matrix[,'test3'])

fit2_1 <- eBayes(fit2_1)
fit2_2 <- eBayes(fit2_2)
fit2_3 <- eBayes(fit2_3)


#topTableF ranks cpg sites on the basis of moderated F-statistics for one or more coefficients. FDR calculated using Benjamin-Hochberg method
tT_1 <- topTable(fit2_1, adjust="BH", sort.by="P", number=nrow(fit2_1))
tT_2 <- topTable(fit2_2, adjust="BH", sort.by="P", number=nrow(fit2_2))
tT_3 <- topTable(fit2_3, adjust="BH", sort.by="P", number=nrow(fit2_3))



listoftT <- list("Mild.moderate.COPDvsControl" = tT_1, 
                 "Severe.COPDvsControl" = tT_2, 
                 "Severe.COPDvsMild.moderate.COPD" = tT_3)

listoftT2 <- list()
listofvolcano <- list()


for (i in 1:length(listoftT)){
  tT <- listoftT[[i]]
  tT$Legend <- ifelse(
    tT$adj.P.Val < 0.05 & tT$logFC > 0, "Upregulated", 
    ifelse(
      tT$adj.P.Val < 0.05 & tT$logFC < -0, "Downregulated",
      "Not Significant"))
  
  
  tT$genename <- row.names(tT)
  selection <- which(tT$adj.P.Val<0.05)

  tT2 <- tT[selection,]

  listoftT2[[names(listoftT[i])]] <- tT2

  write.csv(tT2, file.path(tf.results.dir, paste0("tT2_", names(listoftT[i]), ".csv")))
  
  #Volcano Plot
  volcano <- ggplot(tT, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = Legend)) +
    # scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey"))+
    scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"))+
    geom_hline(yintercept =-log10(max(tT2$P.Value)),colour="black", linetype="dashed")+
    geom_text_repel(data = subset(tT2[1:20,]),
                    aes(label= genename),size = 3, box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines") ) +
    theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                     legend.text = element_text(size = 14),
                                     legend.title = element_text(size = 16)) +
    labs(title = names(listoftT)[i])

  listofvolcano[[i]] <- volcano

  ggsave(volcano, filename = file.path(tf.figures.dir, paste0("volcano_", names(listoftT[i]), ".png")),
         width = 22,
         height = 22,
         units = "cm")

}

saveRDS(listoftT, file.path(tf.results.dir, "listoftT.rds"))
saveRDS(listoftT2, file.path(tf.results.dir, "listoftT2.rds"))

cat("END OF THIS JOB", Sys.time(), "\n")


# VIPER output is typically a normalized enrichment score (NES) for each TF per sample (or group), based on how the expression of its target genes (its regulon) behaves.
# A positive NES indicates the TF’s targets are activated as a group, suggesting TF activation.
# A negative NES means the TF’s targets are repressed as a group, suggesting TF inhibition.


# ================================================================================== #
# 8. SHERLOCK1 VS SHERLOCK3 TF RESULTS =============================================
# ================================================================================== #
## Load SHERLOCK1 TF and SHERLOCK2 TF results
Sherlock1_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "listoftT2.rds"))


Sherlock3_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "brush","results", "listoftT2.rds"))
# Sherlock3_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "brush_exsmoker","results", "listoftT2.rds"))



# ================================================================================== #
# 8.1. VENN DIAGRAMS ===============================================================
# ================================================================================== #

tf.venn.dir <- file.path(output.dir, "tf_analysis", "sherlock1", "figures", "venn")
if(!exists(tf.venn.dir)) dir.create(tf.venn.dir)
# 
# tf.venn.dir <- file.path(output.dir, "tf_analysis", "sherlock1", "figures", "venn_exsmoker")
# if(!exists(tf.venn.dir)) dir.create(tf.venn.dir)


for (direction in c("Upregulated", "Downregulated")){
  
pdf(file= file.path(tf.venn.dir, paste0("tf_validation_venn_",direction,".pdf")), width = 8, height = 8)


# Venn diagram
x <- list(SHERLOCK1 = row.names(Sherlock1_tT2[["mCOPDvsNon_COPD"]][which(Sherlock1_tT2[["mCOPDvsNon_COPD"]]$Legend == direction),]),
          SHERLOCK2 = row.names(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]][which(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]]$Legend == direction),]))

if(all(sapply(x,function(i) length(i))) != 0){
ggvenn(x, stroke_size = 1.5) + ggtitle(paste("Mild.moderate.COPD - Control:", direction))}

x <- list(SHERLOCK1 = row.names(Sherlock1_tT2[["sCOPDvsNon_COPD"]][which(Sherlock1_tT2[["sCOPDvsNon_COPD"]]$Legend == direction),]),
          SHERLOCK2 = row.names(Sherlock3_tT2[["Severe.COPDvsControl"]][which(Sherlock3_tT2[["Severe.COPDvsControl"]]$Legend == direction),]))

if(all(sapply(x,function(i) length(i))) != 0){
print(ggvenn(x, stroke_size = 1.5) + ggtitle(paste("Severe.COPD - Control:", direction)))}


x <- list(SHERLOCK1 = row.names(Sherlock1_tT2[["sCOPDvsmCOPD"]][which(Sherlock1_tT2[["sCOPDvsmCOPD"]]$Legend == direction),]),
          SHERLOCK2 = row.names(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]][which(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]]$Legend == direction),]))

if(all(sapply(x,function(i) length(i))) != 0){
print(ggvenn(x, stroke_size = 1.5) + ggtitle(paste("Severe.COPD - Mild.moderate.COPD:", direction)))}

dev.off()


# ================================================================================== #
# 8.2. COMMON LIST =================================================================
# ================================================================================== #
write.csv(intersect(row.names(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]][which(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]]$Legend == direction),]),
                    row.names(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]][which(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]]$Legend == direction),])),
          file.path(tf.venn.dir,paste0("Common_Mild.moderate.COPDvsControl_",direction,".csv")))

write.csv(intersect(row.names(Sherlock1_tT2[["Severe.COPDvsControl"]][which(Sherlock1_tT2[["Severe.COPDvsControl"]]$Legend == direction),]),
                    row.names(Sherlock3_tT2[["Severe.COPDvsControl"]][which(Sherlock3_tT2[["Severe.COPDvsControl"]]$Legend == direction),])),
          file.path(tf.venn.dir,paste0("Common_Severe.COPDvsControl_",direction,".csv")))

write.csv(intersect(row.names(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]][which(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]]$Legend == direction),]),
                    row.names(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]][which(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]]$Legend == direction),])),
          file.path(tf.venn.dir,paste0("Common_Severe.COPDvsMild.moderate.COPD_",direction,".csv")))
}

# Common - not seperated by upregulated and downregulated
common1 <- intersect(row.names(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]]), row.names(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]]))
if(length(common1) >0){
write.csv(cbind(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]][common1,],Sherlock3_tT2[["Mild.moderate.COPDvsControl"]][common1,]),
                file.path(tf.venn.dir,"Common_Mild.moderate.COPDvsControl.csv"))
}

common2 <- intersect(row.names(Sherlock1_tT2[["Severe.COPDvsControl"]]), row.names(Sherlock3_tT2[["Severe.COPDvsControl"]]))
if(length(common2) >0){
write.csv(cbind(Sherlock1_tT2[["Severe.COPDvsControl"]][common2,],Sherlock3_tT2[["Severe.COPDvsControl"]][common2,]),
          file.path(tf.venn.dir,"Common_Severe.COPDvsControl.csv"))
}

common3 <- intersect(row.names(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]]), row.names(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]]))
if(length(common3) >0){
write.csv(cbind(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]][common3,],Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]][common3,]),
          file.path(tf.venn.dir,"Common_Severe.COPDvsMild.moderate.COPD.csv"))
}


# ================================================================================== #
# 9. BOXPLOTS of TF Activity Scores ================================================
# ================================================================================== #

#LOAD Clinical files
sherlock1_clinical <- readRDS(file.path(data.dir, "processed", "SHERLOCK1", "clinical_sk1_simple.rds"))
all_clinical<- readRDS(file.path(postQC.data.dir, "clinical_brushbiopt_simple.rds"))
sherlock3_brush_clinical <- all_clinical[which(all_clinical$sampletype == "Brush"),]

#LOAD TF activity scores for all genes
sherlock1_TF_results <- readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "TF_analysis_results.rds"))
sherlock3_TF_results <- readRDS(file.path(output.dir, "tf_analysis", "brush", "results", "TF_analysis_results.rds"))
# sherlock3_TF_results <- readRDS(file.path(output.dir, "tf_analysis", "brush_exsmoker", "results", "TF_analysis_results.rds"))

#LOAD tT2 genes
Sherlock1_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "listoftT2.rds"))
Sherlock3_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "brush","results", "listoftT2.rds"))
# Sherlock3_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "brush_exsmoker","results", "listoftT2.rds"))

#LOAD Intersecting tT2 genes
common1 <- intersect(row.names(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]]), row.names(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]]))
common2 <- intersect(row.names(Sherlock1_tT2[["Severe.COPDvsControl"]]), row.names(Sherlock3_tT2[["Severe.COPDvsControl"]]))
common3 <- intersect(row.names(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]]), row.names(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]]))

listofboxplots <- list()
boxplot_func <- function(study, comparison){
  

  listofboxplots_study <- list()
  
  
  #Set up for different options for study and comparison arguments in the function
  if(study == "SHERLOCK1"){
    TF_results <- sherlock1_TF_results
    clinical <- sherlock1_clinical
    x_order <- c("Non_COPD","mCOPD","sCOPD")
    
  }
  
  if(study == "SHERLOCK2"){
    TF_results <- sherlock3_TF_results
    clinical <- sherlock2_brush_clinical
    # clinical <- clinical[which(clinical$Smoking.status == "Ex.smoker"),]
    x_order <- c("Control","Mild.moderate.COPD","Severe.COPD")
    
  }
  
  if(comparison == "MildModvsControl"){
    common_tfs <- common1}
  
  if(comparison == "SeverevsControl"){
    common_tfs <- common2}
  
  if(comparison == "SeverevsMildMod"){
    common_tfs <- common3}
  

  
  #Check that sample IDs are in same order
  colnames(TF_results) == row.names(clinical)
  
  if(all(colnames(TF_results) == row.names(clinical)) == FALSE){
    stop("all(colnames(TF_results) == row.names(clinical)) == FALSE") }
  
for (tf in common_tfs){
#bind together TF score and sample disease classification
boxplot_data <- as.data.frame(cbind(tf_activity = t(TF_results[tf,]),
                                    group = clinical$classification))

colnames(boxplot_data)[1] <- "tf_activity"

boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None")

my_comparisons <- combn(unique(clinical$classification), 2, simplify = FALSE)

boxplot_data$tf_activity <- as.numeric(boxplot_data$tf_activity)
boxplot_data$group <- factor(boxplot_data$group, levels = x_order)

stat.table <- boxplot_data  %>%
  t_test(tf_activity ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group")

stat.table <- stat.table[which(stat.table$p < 0.05),]
lowest_bracket <- max(boxplot_data$tf_activity) + 0.05*(max(boxplot_data$tf_activity))
stat.table$y.position <- seq(lowest_bracket, by= 0.1*max(boxplot_data$tf_activity), length.out = nrow(stat.table))

boxplotfinal <- ggplot(boxplot_data, aes(
  x = factor(group, level = x_order),
  y = tf_activity,
  group = group)) +
  
  theme_bw()+
  
  boxplot_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.2) +
  
  {if(nrow(stat.table) >0 )
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       size = 5)
  } +
  stat_summary(fun = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste(study, "TF Activity:", tf)
       # ,
  #      caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1B"
       # caption = "Signature:S100A8"
       # caption = paste0("Signature: ", paste0(gene_set_list[[1]], collapse = " "))
       
       
  ) +
  ylab (label = "TF Activity Score") +
  xlab (label = "Disease")




listofboxplots_study[[tf]] <- boxplotfinal


} #close loop
  return(listofboxplots_study)
  
} #close function


#Run the function for each comparison and study - this returns a list of boxplots
for (j in c("SeverevsControl", "MildModvsControl", "SeverevsMildMod")){
  
listofboxplots[[j]][["SHERLOCK1"]] <- boxplot_func(study = "SHERLOCK1", comparison = j)
listofboxplots[[j]][["SHERLOCK3"]] <- boxplot_func(study = "SHERLOCK3", comparison = j)
# }

#Save all boxplots
this.boxplot.dir <- file.path(output.dir, "tf_analysis", "sherlock1", "figures", "boxplot")
if(!exists(this.boxplot.dir)) dir.create(this.boxplot.dir)


# for (comparison in c("SeverevsControl", "MildModvsControl", "SeverevsMildMod")){
  
# pdf(file= file.path(this.boxplot.dir, paste0(j,"_tf_boxplots.pdf")), width = 9, height = 12)
pdf(file= file.path(this.boxplot.dir, paste0(j,"_tf_boxplots.pdf")), width = 9, height = 12)

for (i in names(listofboxplots[[j]][["SHERLOCK1"]])){
  sherlock1_plot <- listofboxplots[[j]][["SHERLOCK1"]][[i]]
  sherlock3_plot <- listofboxplots[[j]][["SHERLOCK3"]][[i]]

  print(ggarrange(plotlist=list(sherlock1_plot, sherlock3_plot),
            nrow = 2, ncol =1))
}
  dev.off()
  
}










# ================================================================================== #
# 10. BOXPLOTS of Gene Expression  ================================================
# ================================================================================== #

#LOAD Clinical files
sherlock1_clinical <- readRDS(file.path(data.dir, "processed", "SHERLOCK1", "clinical.rds"))
all_clinical <- readRDS(file.path(data.dir, "processed", "datawrangling_qc", "clinical.rds"))
sherlock2_brush_clinical <- all_clinical[which(all_clinical$sampletype == "Brush"),]
hgnc_symbols_db <- readRDS(file.path(data.dir,"processed", "datawrangling_qc", "hgnc_symbols_db.rds"))


#LOAD Counts for all genes
sherlock1_counts <- readRDS(file.path(data.dir, "processed", "SHERLOCK1", "counts.rds"))
all_counts <- readRDS(file.path(data.dir, "processed", "datawrangling_qc", "counts.rds"))
sherlock2_brush_counts<- all_counts[, row.names(sherlock2_brush_clinical)]


listofboxplots <- list()
boxplot_func <- function(study){
  
  
  listofboxplots_study <- list()
  
  
  #Set up for different options for study and comparison arguments in the function
  if(study == "SHERLOCK1"){
    counts <- sherlock1_counts
    counts_voom<- voom(counts)
    clinical <- sherlock1_clinical
    x_order <- c("Non_COPD","mCOPD","sCOPD")
    
  }
  
  if(study == "SHERLOCK2"){
    counts <- sherlock2_brush_counts
    counts_voom<- voom(counts)
    clinical <- sherlock2_brush_clinical
    # clinical <- clinical[which(clinical$Smoking.status == "Ex.smoker"),]
    x_order <- c("Control","Mild.moderate.COPD","Severe.COPD")
    
  }
  
  #Check that sample IDs are in same order
  colnames(counts_voom) == row.names(clinical)
  
  if(all(colnames(counts_voom) == row.names(clinical)) == FALSE){
    stop("all(colnames(counts_voom) == row.names(clinical)) == FALSE") }
  
  for (gene in c("NR3C1", "STAT3", "NFE2L1", "FKBP5", "KEL")){
    
    geneofinterestid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == gene), "GENEID"]
    
    #bind together TF score and sample disease classification
    boxplot_data <- as.data.frame(cbind(gene = counts_voom$E[geneofinterestid,],
                                        group = clinical$classification))
    
    colnames(boxplot_data)[1] <- "gene"
    
    boxplot_theme <- theme(axis.title = element_text(size = 24),
                           axis.text = element_text(size = 24),
                           title = element_text(size = 20),
                           legend.position = "None")
    
    # my_comparisons <- combn(unique(clinical$classificationation), 2, simplify = FALSE)
    
    boxplot_data$gene <- as.numeric(boxplot_data$gene)
    boxplot_data$group <- factor(boxplot_data$group, levels = x_order)
    
    stat.table <- boxplot_data  %>%
      t_test(gene ~ group,
             paired = FALSE) %>%
      add_xy_position(x = "group")
    
    stat.table <- stat.table[which(stat.table$p < 0.05),]
    lowest_bracket <- max(boxplot_data$gene) + 0.01*(max(boxplot_data$gene))
    stat.table$y.position <- seq(lowest_bracket, by= 0.01*max(boxplot_data$gene), length.out = nrow(stat.table))
    
    boxplotfinal <- ggplot(boxplot_data, aes(
      x = factor(group, level = x_order),
      y = gene,
      group = group)) +
      
      theme_bw()+
      
      boxplot_theme +
      
      geom_boxplot(position = position_dodge(1)) +
      
      geom_jitter(aes(color = group),
                  alpha = 0.5,
                  size = 2.5, 
                  width = 0.2) +
      
      {if(nrow(stat.table) >0 )
        stat_pvalue_manual(stat.table,
                           label = "p",
                           tip.length = 0.01,
                           size = 5)
      } +
      stat_summary(fun = mean, fill = "red",
                   geom = "point", shape = 21, size =4,
                   show.legend = TRUE) +
      
      theme(axis.text.x = element_text(size = 18))+
      labs(title = paste(study, "Gene Expression:", gene)
           
      ) +
      ylab (label = "Gene Expression") +
      xlab (label = "Disease")
    
    
    
    
    listofboxplots_study[[gene]] <- boxplotfinal
    
    
  } #close loop
  return(listofboxplots_study)
  
} #close function


#Run the function for each comparison and study - this returns a list of boxplots

  
  listofboxplots[["SHERLOCK1"]] <- boxplot_func(study = "SHERLOCK1")
  listofboxplots[["SHERLOCK2"]] <- boxplot_func(study = "SHERLOCK2")

  
  #Save all boxplots
  this.boxplot.dir <- file.path(output.dir, "tf_analysis", "sherlock1", "figures", "boxplot")
  if(!exists(this.boxplot.dir)) dir.create(this.boxplot.dir)
  
  pdf(file= file.path(this.boxplot.dir, paste0("gene_expression_boxplots.pdf")), width = 9, height = 12)
  
  for (i in names(listofboxplots[["SHERLOCK1"]])){
    sherlock1_plot <- listofboxplots[["SHERLOCK1"]][[i]]
    sherlock2_plot <- listofboxplots[["SHERLOCK2"]][[i]]
    
    print(ggarrange(plotlist=list(sherlock1_plot, sherlock2_plot),
                    nrow = 2, ncol =1))
  }
  dev.off()
  


# ================================================================================== #
# 9. HEATMAPS ======================================================================
# ================================================================================== #

