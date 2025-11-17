#Microplastics differential expression
#use robust= TRUE 
#DESeq2 better for reducing false positives (more conservative)
#can i even correct for sex if theres only one female and two males in the highpolyure group 


# Project with Daan & Barbro Melgert measuring the levels of specific microplastic particles in bronchial wash samples of the Sherlock study.
# Aim: see if their results are correlated with specific gene expression signals (differential gene expression analyses with the RNA-seq data from the bronchial brushes)

#In the attached file you will find the SEO numbers of the subjects in which microplastics were measured and the results of the microplastics measurements.
#Here you will find 3 tables with subjects listed in green and subjects listed in orange, can you perform a differential gene expression analysis between the green and orange subjects (3 different analyses for the 3 different tables).
#I know the power will probably be too low to generate genome-wide significant hits, but we would still like to see the table with the results. Could you send the excel file with results for all genes?
  


options(error = function() { traceback(2); quit(status = 1) })

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("ggrepel")
library("readxl")
library(ggfortify)
library(ggplot2)
library(heatmap3)
library(gplots)
library(biomaRt)
library(edgeR)
library(limma)
library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(DESeq2)
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

#MP output directory
mp.dir <- file.path(output.dir, "microplastics")
if(!exists(mp.dir)) dir.create(mp.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(data.dir, "raw", "microplastics"))

total_mp <- as.data.frame(read_xlsx("microplastics_patients.xlsx", sheet = "total_microplastics")) %>%  column_to_rownames(var = "Number")
nylon <- as.data.frame(read_xlsx("microplastics_patients.xlsx", sheet = "nylon"))  %>%  column_to_rownames(var = "Number")
polyu <- as.data.frame(read_xlsx("microplastics_patients.xlsx", sheet = "polyurethane"))  %>%  column_to_rownames(var = "Number")
  
setwd(file.path(main.dir))

clinical_brush <- readRDS(file.path(postQC.data.dir, "clinical_brush_simple.rds"))
counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


#Match up patients
#SEO250, SEO131 and SEO304 only have biopsy sample, not brush

matching_patients <- intersect(row.names(total_mp),clinical_brush$Study.ID)
clinical_total_mp <- clinical_brush[match(matching_patients,clinical_brush$Study.ID),] %>% 
  mutate(quantity = total_mp[matching_patients,'Total nr. of microplastics']) %>% 
  mutate(exposure = ifelse(quantity < 8.45, "low", "high"))

matching_patients <- intersect(row.names(nylon),clinical_brush$Study.ID)
clinical_nylon <- clinical_brush[match(matching_patients,clinical_brush$Study.ID),] %>% 
mutate(quantity = nylon[matching_patients,'Nylon']) %>% 
  mutate(exposure = ifelse(quantity < 22.88, "low", "high"))

matching_patients <- intersect(row.names(polyu),clinical_brush$Study.ID)
clinical_polyu <- clinical_brush[match(matching_patients,clinical_brush$Study.ID),] %>% 
  mutate(quantity = polyu[matching_patients,'Polyurethane']) %>% 
  mutate(exposure = ifelse(quantity < 67.44, "low", "high"))



diffexp_edgeR <- function(microplastic) {

  
  #Make output directory for total_mp or nylon or polyu
  this.mp.dir <- file.path(mp.dir, microplastic)
  if(!exists(this.mp.dir)) dir.create(this.mp.dir, recursive = TRUE)
  
  diffexp.dir <- file.path(this.mp.dir, "diffexp")
  if(!exists(diffexp.dir)) dir.create(diffexp.dir, recursive = TRUE)
  
  
  #Set results paths
  diffexp.results.dir <- file.path(diffexp.dir, "results")
  if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)
  
  
  diffexp.figures.dir <- file.path(diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir)) dir.create(diffexp.figures.dir, recursive = TRUE)
  
  

  if(microplastic == "total_mp"){
  clinical <- clinical_total_mp}
  
  if(microplastic == "nylon"){
  clinical <- clinical_nylon}
  
  
  if(microplastic == "polyu"){
  clinical <- clinical_polyu}
  
  
  counts <- counts_brush[,row.names(clinical)]
  


  design <- model.matrix(~ 0 + exposure + age + sex + smoking_status, 
                          data = clinical) #Design matrix
  
  colnames(design)[1:2] <- c(levels(as.factor(clinical$exposure)))
  
  DGEL<- DGEList(counts=counts) 
  
  # FILTER
  keep <- filterByExpr(DGEL, group = clinical$classification) 
  # The filterBy Expr function keeps rows that have worthwhile counts in a minimum number of samples.
  # keep <- which(rowMedians(as.matrix(DGEL))>10) 
  # keep <- rowSums(cpm(expression)>100) >= 2
  
  
  DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
  # expression <- DGEL$counts
  
  # NORMALISE
  DGEL<- calcNormFactors(DGEL,method = "TMM")
  
  # ESTIMATE DISPERSON
  DGEL <- estimateDisp(DGEL, design)
  
  # FIT MODEL
  fit <- glmQLFit(DGEL, design, robust = TRUE) #fit the GLM (design) to the DGEL(the DGEL object, which contains the counts data that has been filtered,normalised and dispersons estimtated)
  
  my.contrasts <- makeContrasts(contrast1 = high - low,
                                levels = design)
  
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,1], robust)
  
  tT <- topTags(qlf,n=nrow(DGEL))$table #topTags gets top genes, here we want all of the genes, edgeR's default p.adjust method is BH
  tT$Legend <- ifelse(
    tT$FDR < 0.05 & tT$logFC > 0, "Upregulated",
    ifelse(
      tT$FDR < 0.05 & tT$logFC < 0, "Downregulated",
      "Not Significant"))
  
  tT$Legend[is.na(tT$Legend)]="Not Significant"
  
  tT$Legend <- factor(tT$Legend, levels = c("Downregulated", "Upregulated", "Not Significant"))
  
  tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols
  
  # if(showEnsemblID == TRUE){
    #for those with no hgnc symbol, label with ensembl id
    tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl
  # }
  # 
  # else{
  #   #for those with no hgnc symbol, remove
  #   tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
  # }
  
  selection <-which(tT$FDR<0.05)
  
  tT2 <- tT[selection,]
  
  write.csv(tT, file = paste0(diffexp.results.dir, "tT_", microplastic, ".csv"))
  write.csv(head(tT, n = 100) , file = paste0(diffexp.results.dir, "tT_top100_", microplastic, ".csv"))
  
  
  # ================================================================================== #
  # 2.1. VOLCANO PLOT ================================================================
  # ================================================================================== #
  cat(paste("Starting 2.1. VOLCANO PLOT", microplastic), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
    ggtitle(paste(microplastic)) +
    geom_point(aes(color = Legend)) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
    geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
    geom_text_repel(data = subset(tT2[1:30,]),
                    aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines") ) +
    theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                     legend.text = element_text(size = 14),
                                     legend.title = element_text(size = 16)) 
  
  ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(microplastic,"_volcano_plot.png")),
         width = 25, height = 25,
         units = "cm")
  
  
  listofresults <- list(tT = tT, tT2 = tT2, volcano = volcano)
  #Save all results
  saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))
  

} #close function


diffexp_edgeR(microplastic = "total_mp")
diffexp_edgeR(microplastic = "nylon")
diffexp_edgeR(microplastic = "polyu")


#### DESEQ2####


diffexp_deseq <- function(microplastic) {
  
  
  #Make output directory for total_mp or nylon or polyu
  this.mp.dir <- file.path(mp.dir, microplastic)
  if(!exists(this.mp.dir)) dir.create(this.mp.dir, recursive = TRUE)
  
  diffexp.dir <- file.path(this.mp.dir, "diffexp_deseq")
  if(!exists(diffexp.dir)) dir.create(diffexp.dir, recursive = TRUE)
  
  
  #Set results paths
  diffexp.results.dir <- file.path(diffexp.dir, "results")
  if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)
  
  
  diffexp.figures.dir <- file.path(diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir)) dir.create(diffexp.figures.dir, recursive = TRUE)
  
  
  
  if(microplastic == "total_mp"){
    clinical <- clinical_total_mp}
  
  if(microplastic == "nylon"){
    clinical <- clinical_nylon}
  
  
  if(microplastic == "polyu"){
    clinical <- clinical_polyu}
  
  
  counts <- counts_brush[,row.names(clinical)]
  
##### Q1a DGE #####
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = clinical,
                              design = ~ 0 + exposure + age + sex + smoking_status)

# Filter the genes that are lowly expressed and normalize
# Low exp genes affects statistics - can make p value very significant even if only a few samples are lowly expressed amongst other samples with no expression.  also affects multiple testing.
# One method = keep row medians that are greater than 10 (ie. half of the samples for a gene must have a minimum number of 10 counts)
keep <- rowMedians(counts(dds)) >= 10

dds <- dds[keep,]
dds <- DESeq(dds)

# results extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values;
resultsNames(dds)

results <- results(dds, contrast = c("exposure", "high", "low")) #This means we have set high as control and we define fold change based on high as baseline (log fold change = high - low). switch the order of low and high to make low the baseline

tT <- as.data.frame(results)

# is this done auto
# tT$padj=p.adjust(tT$pvalue,method="BH")

tT <- as.data.frame(results)
tT <-tT[order(tT$pvalue),]
# baseMean is the average of the normalized count values, dividing by size factors, taken over all samples.


tT$Legend <- ifelse(
  tT$padj < 0.05 & tT$log2FoldChange > 0, "Upregulated",
  ifelse(
    tT$padj < 0.05 & tT$log2FoldChange < 0, "Downregulated",
    "Not Significant"))

tT$Legend[is.na(tT$Legend)]="Not Significant"

tT$Legend <- factor(tT$Legend, levels = c("Downregulated", "Upregulated", "Not Significant"))

tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols

# if(showEnsemblID == TRUE){
#for those with no hgnc symbol, label with ensembl id
tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl
# }
# 
# else{
#   #for those with no hgnc symbol, remove
#   tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
# }

selection <-which(tT$padj<0.05)

tT2 <- tT[selection,]

write.csv(tT, file = paste0(diffexp.results.dir, "tT_", microplastic, ".csv"))
write.csv(head(tT, n = 100) , file = paste0(diffexp.results.dir, "tT_top100_", microplastic, ".csv"))


# ================================================================================== #
# 2.1. VOLCANO PLOT ================================================================
# ================================================================================== #
cat(paste("Starting 2.1. VOLCANO PLOT", microplastic), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

volcano <- ggplot(tT, aes(x = log2FoldChange, y = -log10(pvalue))) +
  ggtitle(paste(microplastic)) +
  geom_point(aes(color = Legend)) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
  geom_hline(yintercept =-log10(max(tT2$pvalue)),colour="black", linetype="dashed")+
  geom_text_repel(data = subset(tT2[1:30,]),
                  aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines") ) +
  theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                   legend.text = element_text(size = 14),
                                   legend.title = element_text(size = 16)) 

ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(microplastic,"_volcano_plot.png")),
       width = 25, height = 25,
       units = "cm")


listofresults <- list(tT = tT, tT2 = tT2, volcano = volcano)
#Save all results
saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))


} #close function


diffexp_deseq(microplastic = "total_mp")
diffexp_deseq(microplastic = "nylon")
diffexp_deseq(microplastic = "polyu")

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
