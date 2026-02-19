#Microplastics differential expression
#use robust= TRUE 
#DESeq2 better for reducing false positives (more conservative)
#can i even correct for sex if theres only one female and two males in the highpolyure group ?


# Project with Daan & Barbro Melgert measuring the levels of specific microplastic particles in bronchial wash samples of the Sherlock study.
# Aim: see if their results are correlated with specific gene expression signals (differential gene expression analyses with the RNA-seq data from the bronchial brushes)

# November 2025 (Batch 1)
# "In the attached file you will find the SEO numbers of the subjects in which microplastics were measured and the results of the microplastics measurements.
# Here you will find 3 tables with subjects listed in green and subjects listed in orange, can you perform a differential gene expression analysis between the green and orange subjects (3 different analyses for the 3 different tables).
# I know the power will probably be too low to generate genome-wide significant hits, but we would still like to see the table with the results.
# Could you send the excel file with results for all genes?"


# 03/02/2026 (Batch 2)
# Daan sent 260124 Summary results and plan correlation with RNAseq signatures PyGCMS
# "we just received new microplastic data. These were generated using a different method to quantify microplastics.
# Kathy, would you be able to also run the same analyses with these microplastic data (See attached).
# It will be interesting to see if we get similar or very different results."



options(error = function() { traceback(2); quit(status = 1) })

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library(ggrepel)
library(readxl)
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
library(rstatix)
library(ggpubr)

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
clinical123_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
clinical_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush" & clinical123_master$batch != "1"),] 
clinical_biopt <- clinical123_master[which(clinical123_master$sampletype == "Biopt" & clinical123_master$batch != "1"),] 

counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds")) #only includes sherlock2and3
counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds")) #only includes sherlock2and3

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


setwd(file.path(data.dir, "raw", "microplastics"))

# Batch 1 data (November 2025)
total_mp1 <- as.data.frame(read_xlsx("microplastics_patients.xlsx", sheet = "total_microplastics")) %>%  column_to_rownames(var = "Number")
nylon1 <- as.data.frame(read_xlsx("microplastics_patients.xlsx", sheet = "nylon"))  %>%  column_to_rownames(var = "Number")
poly1 <- as.data.frame(read_xlsx("microplastics_patients.xlsx", sheet = "polyurethane"))  %>%  column_to_rownames(var = "Number")
  
# Batch 2 data (03/02/36)
total_mp2 <- as.data.frame(read_xlsx("260124 Summary results and plan correlation with RNAseq signatures PyGCMS.xlsx", sheet = "Total")) %>%  column_to_rownames(var = "Number")
nylon2 <- as.data.frame(read_xlsx("260124 Summary results and plan correlation with RNAseq signatures PyGCMS.xlsx", sheet = "Nylon"))  %>%  column_to_rownames(var = "Number")
poly2 <- as.data.frame(read_xlsx("260124 Summary results and plan correlation with RNAseq signatures PyGCMS.xlsx", sheet = "Polypropylene"))  %>%  column_to_rownames(var = "Number")
  
setwd(file.path(main.dir))


# ================================================================================== #
# 2. Diff Exp(define function) ==================================================
# ================================================================================== #

#### Define DESEQ2 function ####
diffexp_deseq_func <- function(microplastic) {
  
  #batch 1 has polyurethane data, batch 2 has polypropylene
  microplastic_fullname <- microplastic
  if (microplastic_fullname == "poly"){
    
    if(batch == 1){
      microplastic_fullname <- "polyurethane"}
    
    else if(batch == 2){
      microplastic_fullname <- "polypropylene"}
  }
  
  cat(paste("Starting 2. Differential expression", microplastic_fullname, sample_type), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  
  #Make output directory for total_mp or nylon or poly
  this.mp.dir <- file.path(mp.dir, microplastic)
  if(!exists(this.mp.dir)) dir.create(this.mp.dir, recursive = TRUE)
  
  diffexp.dir <- file.path(this.mp.dir, "diffexp_deseq")
  if(!exists(diffexp.dir)) dir.create(diffexp.dir, recursive = TRUE)
  
  
  #Set results paths
  diffexp.results.dir <- file.path(diffexp.dir, "results")
  if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)
  
  
  diffexp.figures.dir <- file.path(diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir)) dir.create(diffexp.figures.dir, recursive = TRUE)
  
  if(batch == 1){ #batch is the index in loop after deseq and boxplot function
    total_mp <- total_mp1
    nylon <- nylon1
    poly <- poly1
    
    #Match up patients
    #Batch 1 : SEO250, SEO131 and SEO304 only have biopsy sample, not brush.  print sample names that donti ntersect
    print(paste0("Batch ",batch," samples do not have ", sample_type, " data", paste0(row.names(total_mp)[!row.names(total_mp) %in% clinical$Study.ID], collapse = ",")))
    
    if(microplastic == "total_mp"){
    matching_patients <- intersect(row.names(total_mp),clinical$Study.ID)
    clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
      mutate(quantity = total_mp[matching_patients,'Total nr. of microplastics']) %>% 
      mutate(exposure = ifelse(quantity < 8.45, "low", "high"))}
    
    if(microplastic == "nylon"){
    matching_patients <- intersect(row.names(nylon),clinical$Study.ID)
    clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
      mutate(quantity = nylon[matching_patients,'Nylon']) %>% 
      mutate(exposure = ifelse(quantity < 22.88, "low", "high"))}
    
    if(microplastic == "poly"){
    matching_patients <- intersect(row.names(poly),clinical$Study.ID)
    clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
      mutate(quantity = poly[matching_patients,'Polyurethane']) %>% 
      mutate(exposure = ifelse(quantity < 67.44, "low", "high"))}
    
  }
  
  
  if(batch == 2){
    total_mp <- total_mp2
    nylon <- nylon2
    poly <- poly2
    
    
    #Match up patients
    if(microplastic == "total_mp"){
      
      #SEO532 - no expression data , print sample names that donti ntersect
      print(paste0("Batch ",batch," samples not in clinical file: ", row.names(total_mp)[!row.names(total_mp) %in% clinical$Study.ID]))
      matching_patients <- intersect(row.names(total_mp),clinical$Study.ID)
      clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
        mutate(quantity = total_mp[matching_patients,'Total ng/ml of microplastics']) %>% 
        mutate(exposure = total_mp[matching_patients,'Category'])}
    
    if(microplastic == "nylon"){
      matching_patients <- intersect(row.names(nylon),clinical$Study.ID)
      clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
        mutate(quantity = nylon[matching_patients,'Nylon ng/ml']) %>% 
        mutate(exposure = nylon[matching_patients,'Category'])}
    
    
    if(microplastic == "poly"){
      matching_patients <- intersect(row.names(poly),clinical$Study.ID)
      clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
        mutate(quantity = poly[matching_patients,'Polypropylene ng/ml']) %>% 
        mutate(exposure = poly[matching_patients,'Category'])}
    
    clinical$exposure <- tolower(clinical$exposure) #make lowercase to match batch 1
  }
  
  counts <- counts[,row.names(clinical)]
  
##### Q1a DGE #####
  
if(batch == 2 && microplastic == nylon){
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = clinical,
                                design = ~ 0 + exposure + age + smoking_status) # all are high sampples are males, no sex correction
  } else{

    dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = clinical,
                              design = ~ 0 + exposure + age + sex + smoking_status) #note the order of the variables here doesn't matter. exposure at the end would be the same
}
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

write.csv(tT, file = file.path(diffexp.results.dir, paste0("batch", batch,"_tT_", microplastic, ".csv")))
write.csv(head(tT, n = 100) , file.path(diffexp.results.dir, paste0("batch", batch,"_tT_top100_", microplastic, ".csv")))




# ================================================================================== #
# 2.1. VOLCANO PLOT ================================================================
# ================================================================================== #
cat(paste("Starting 2.1. VOLCANO PLOT", microplastic_fullname, sample_type), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

volcano <- ggplot(tT, aes(x = log2FoldChange, y = -log10(pvalue))) +
  ggtitle(paste(microplastic_fullname)) +
  geom_point(aes(color = Legend)) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
  geom_hline(yintercept =-log10(max(tT2$pvalue)),colour="black", linetype="dashed")+
  geom_text_repel(data = subset(tT2[1:30,]),
                  aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines") ) +
  theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                   legend.text = element_text(size = 14),
                                   legend.title = element_text(size = 16)) 

ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0("batch", batch,"_",microplastic_fullname,"_volcano_plot.png")),
       width = 25, height = 25,
       units = "cm")


listofresults <- list(tT = tT, tT2 = tT2, volcano = volcano)
#Save all results
saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))


} #close diffexp function




# ================================================================================== #
# 2.2. BOXPLOTS(definme function) ==================================================
# ================================================================================== #

boxplot_func <- function(microplastic){
  
  #batch 1 has polyurethane data, batch 2 has polypropylene
  microplastic_fullname <- microplastic
  if (microplastic_fullname == "poly"){
    
    if(batch == 1){
      microplastic_fullname <- "polyurethane"}
    
    else if(batch == 2){
      microplastic_fullname <- "polypropylene"}
  }
cat("Starting plots", microplastic_fullname, sample_type, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  this.mp.dir <- file.path(mp.dir, microplastic)
diffexp.dir <- file.path(this.mp.dir, "diffexp_deseq")
diffexp.results.dir <- file.path(diffexp.dir, "results")
diffexp.figures.dir <- file.path(diffexp.dir, "figures")

listofresults <- readRDS(file.path(diffexp.results.dir, "listofresults.RDS"))

if(batch == 1){ #batch is the index in loop after deseq and boxplot function
  total_mp <- total_mp1
  nylon <- nylon1
  poly <- poly1
  
  #Match up patients
  #SEO250, SEO131 and SEO304 only have biopsy sample, not brush
  
  if(microplastic == "total_mp"){
  matching_patients <- intersect(row.names(total_mp),clinical$Study.ID)
  clinical<- clinical[match(matching_patients,clinical$Study.ID),] %>% 
    mutate(quantity = total_mp[matching_patients,'Total nr. of microplastics']) %>% 
    mutate(exposure = ifelse(quantity < 8.45, "low", "high"))}
  
  if(microplastic == "nylon"){
  matching_patients <- intersect(row.names(nylon),clinical$Study.ID)
  clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
    mutate(quantity = nylon[matching_patients,'Nylon']) %>% 
    mutate(exposure = ifelse(quantity < 22.88, "low", "high"))}
  
  if(microplastic == "poly"){
  matching_patients <- intersect(row.names(poly),clinical$Study.ID)
  clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
    mutate(quantity = poly[matching_patients,'Polyurethane']) %>% 
    mutate(exposure = ifelse(quantity < 67.44, "low", "high"))}
  
}


if(batch == 2){
  total_mp <- total_mp2
  nylon <- nylon2
  poly <- poly2
  
  
  #Match up patients
  if(microplastic == "total_mp"){
    
    #SEO532 - no expression data
    matching_patients <- intersect(row.names(total_mp),clinical$Study.ID)
    clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
      mutate(quantity = total_mp[matching_patients,'Total ng/ml of microplastics']) %>% 
      mutate(exposure = total_mp[matching_patients,'Category'])}
  
  if(microplastic == "nylon"){
    
    matching_patients <- intersect(row.names(nylon),clinical$Study.ID)
    clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
      mutate(quantity = nylon[matching_patients,'Nylon ng/ml']) %>% 
      mutate(exposure = nylon[matching_patients,'Category'])}
  
  
  if(microplastic == "poly"){
    matching_patients <- intersect(row.names(poly),clinical$Study.ID)
    clinical <- clinical[match(matching_patients,clinical$Study.ID),] %>% 
      mutate(quantity = poly[matching_patients,'Polypropylene ng/ml']) %>% 
      mutate(exposure = poly[matching_patients,'Category'])}
  
  clinical$exposure <- tolower(clinical$exposure) #make lowercase to match batch 1
}

boxplotcounts <- counts[,row.names(clinical)]
counts_voom <- voom(boxplotcounts)

boxplotdata <- as.data.frame(t(counts_voom$E))

boxplot <- cbind(boxplotdata,
                 quantity = clinical$quantity,
                 exposure = clinical$exposure)

tT2_genes<- row.names(listofresults[["tT2"]])

top10 <- if (length(tT2_genes) >= 10) {
  tT2_genes[1:10]
} else {
  tT2_genes
}


microplastic_fullname <- microplastic
if (microplastic_fullname == "poly"){
  
  if(batch == 1){
    microplastic_fullname <- "polyurethane"}
  
  else if(batch == 2){
    microplastic_fullname <- "polypropylene"}
}



pdf(file = file.path(diffexp.figures.dir, paste0("batch", batch,"_",microplastic_fullname,"_boxplot",".pdf")),
    height = 8,
    width= 6)
cat("Starting plots", microplastic_fullname, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


for (gene in c(top10)){
  
  gene_hgnc <- ifelse(is.na(hgnc_symbols_db[gene, "SYMBOL"]),
                      gene,
                hgnc_symbols_db[gene, "SYMBOL"])


plot <- boxplot[,c(gene,
                   "quantity",
                   "exposure")]

colnames(plot)[1] <- "gene"

plot <- as.data.frame(plot)


## Get P-Values --------------------------------------------------------------------------------------

# THE TABLES IN STEP 1 AND 2 ARE TO GET X POSITIONS, THE T-TEST VALUES WILL NOT BE USED #
#### STEP 1) THIS IS  A TABLE OF THE COMPARISONS I ACTUALLY WANT --------------------------------------------------------------------------------------
stat.table <- plot %>%
  t_test(gene ~ exposure)

stat.table<- stat.table %>%
  add_xy_position(x = "exposure", dodge = 0.8)

# stat.table$contrast <- paste0(stat.table$group1, "-" ,stat.table$group2)



#### STEP 2) MODIFY STAT TABLE TO INCLUDE COMPARISONS OF INTEREST --------------------------------------------------------------------------------------
stat.table3 <- stat.table


stat.table3[,"p"] <- listofresults[["tT"]][gene, "pvalue"]
stat.table3$p <- signif(as.numeric(stat.table3$p), digits = 4)
# stat.table3$y.position <- max(plot[,"gene"]) + 0.025*(max(plot[,"gene"]))

logFC_res <- listofresults[["tT"]][gene, "log2FoldChange"]
pval_res <- listofresults[["tT"]][gene, "pvalue"]

cat("   Making boxplot", microplastic_fullname, gene_hgnc, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


boxplotimage <- ggplot(plot, aes(
  x = as.factor(exposure),
  y = gene
  # fill = smokingstatus
)) +
  theme_bw()+
  
  geom_boxplot(aes(color = exposure),position = position_dodge(1)) +
  
  geom_jitter(aes(color = exposure),
              alpha = 0.5,
              size = 2.5,
              width = 0.3) +
  
  labs(title = paste(gene_hgnc, "expression:", microplastic_fullname, "exposure" )) +
  ylab (label =  paste(gene_hgnc, "expression")) +
  xlab (label = paste(microplastic_fullname, "exposure")) +
  
  
  
  stat_pvalue_manual(stat.table3,
                     label = "p",
                     tip.length = 0.01,
                     size = 3) 



scatterplotimage <- ggplot(plot, aes(
  x = quantity,
  y = gene,
  color = exposure
)) +
  theme_bw()+
  
  
  geom_point(aes(color = exposure)) +
  
  
  labs(title = paste(gene_hgnc, "expression vs", microplastic_fullname, "quantity" )) +
  ylab (label = paste(gene_hgnc, "expression")) +
  xlab (label = paste("Quantity of", microplastic_fullname, "exposure")) +
  
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  
  annotate(
    "text",
    x = Inf,# adjust horizontally
    y = -Inf,  # adjust vertically
    hjust = 1.1,
    vjust = -0.1,
    label = paste0(
      "High - Low exposure", "\n",
      "logFC = ", signif(logFC_res, 3), "\n", 
      "p = ", signif(pval_res, 3)),
    size = 3
  )

image <- ggarrange(plotlist = list(boxplotimage,
                            scatterplotimage),
                   nrow = 2,
                   ncol = 1)

print(image)


} #close loop of top10genes
dev.off()

} # close the boxplot function 






# ================================================================================== #
# 3. RUN THE FUNCTIONS ================================================================
# ================================================================================== #

#set sampletype of interest !! brush or biopsy
for(sample_type in c("Brush", "Biopt")){
  if(sample_type == "Brush"){
    clinical <- clinical_brush
    counts <- counts_brush
  }
  
  if(sample_type == "Biopt"){
    clinical <- clinical_biopt
    counts <- counts_biopt
  }
  
  print(sample_type)
  
  
  
  for (batch in 1:2){
  
  mp.dir <- file.path(output.dir, "microplastics", paste0("data_batch", batch), sample_type)
  if(!exists(mp.dir)) dir.create(mp.dir, recursive = TRUE)

  diffexp_deseq_func(microplastic = "total_mp")
  diffexp_deseq_func(microplastic = "nylon")
  diffexp_deseq_func(microplastic = "poly")

  boxplot_func(microplastic = "total_mp")
  boxplot_func(microplastic = "nylon")
  boxplot_func(microplastic = "poly")
  } # close batch loop
} #close sample_type  loop

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
stop()


# ================================================================================== #
# 4. CELLULAR DECONVLUTION =========================================================
# ================================================================================== #

# Mixture/Count Matrix - CPM normalise
counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
clinical<- readRDS(file.path(postQC.data.dir, "clinical_brushbiopt_simple.rds")) #552 samples
hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))

expression_cpm <- cpm(counts)
expression_cpm <- expression_cpm[hgnc_symbols_db$GENEID,] #552 samples

hidesample <- as.data.frame(cbind(colnames(expression_cpm), seq(1:ncol(expression_cpm))))
colnames(expression_cpm) <- hidesample$V2


## RUN CIBERSORT with ref1_sigmatrix from fia (from Nov 2024, "ref1_1_1411") - check if there is a new one - this is for bronchial
## 1. Create signature matrix was done already - created Ref1_sigmatrix
## Go to 2. Impute cell fractions
## Batch correction, S-mode, 100 permutations
## refmatrix = Ref_1_1411 (single cell reference matrix)
## sigmatrix = Ref1_sigmatrix (signature matrix)
## mixture = expression_cpm


#Upload results from cibersort (results.txt) to processed data folder and run
cibersort_results <- read.delim(file.path(data.dir, "processed", "celldecon", "results_20251114.txt"))
row.names(cibersort_results) <- hidesample$V1
cibersort_results <- cibersort_results[,-1]

row.names(cibersort_results) == rownames(clinical)
cibersort_results <- cibersort_results[row.names(clinical),]

clinical_celldecon_all <- as.data.frame(cbind(clinical,cibersort_results))


# #split into biopt and brush
## Cell Decon vs Exposure Plot --------------------------------------------------------------------------------------

  #### Define DESEQ2 function ####
  celldecon_func <- function(microplastic) {
    
    clinical_celldecon <- clinical_celldecon_all[which(clinical_celldecon_all$sampletype == sample_type),]
    
    #batch 1 has polyurethane data, batch 2 has polypropylene
    microplastic_fullname <- microplastic
    if (microplastic_fullname == "poly"){
      
      if(batch == 1){
        microplastic_fullname <- "polyurethane"}
      
      else if(batch == 2){
        microplastic_fullname <- "polypropylene"}
    }
    
    cat(paste("Starting 4. Cell Decon", microplastic_fullname, sample_type), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    
    #Make output directory for total_mp or nylon or poly
    this.mp.dir <- file.path(mp.dir, microplastic)
    if(!exists(this.mp.dir)) dir.create(this.mp.dir, recursive = TRUE)
    
    celldecon.dir <- file.path(this.mp.dir, "cell_decon")
    if(!exists(celldecon.dir)) dir.create(celldecon.dir, recursive = TRUE)
    
    celldecon.results.dir <- file.path(celldecon.dir, "results")
    if(!exists(celldecon.results.dir)) dir.create(celldecon.results.dir, recursive = TRUE)
    
    celldecon.figures.dir <- file.path(celldecon.dir, "figures")
    if(!exists(celldecon.figures.dir)) dir.create(celldecon.figures.dir, recursive = TRUE)
    
    
  if(batch == 1){ #batch is the index in loop after deseq and boxplot function
    total_mp <- total_mp1
    nylon <- nylon1
    poly <- poly1
    
    #Match up patients
    #Batch 1 : SEO250, SEO131 and SEO304 only have biopsy sample, not brush.  print sample names that donti ntersect
    print(paste0("Batch ",batch," samples do not have ", sample_type, " data", paste0(row.names(total_mp)[!row.names(total_mp) %in% clinical_celldecon$Study.ID], collapse = ",")))
    
    if(microplastic == "total_mp"){
      matching_patients <- intersect(row.names(total_mp),clinical_celldecon$Study.ID)
      clinical_celldecon <- clinical_celldecon[match(matching_patients,clinical_celldecon$Study.ID),] %>% 
        mutate(quantity = total_mp[matching_patients,'Total nr. of microplastics']) %>% 
        mutate(exposure = ifelse(quantity < 8.45, "low", "high"))}
    
    if(microplastic == "nylon"){
      matching_patients <- intersect(row.names(nylon),clinical_celldecon$Study.ID)
      clinical_celldecon <- clinical_celldecon[match(matching_patients,clinical_celldecon$Study.ID),] %>% 
        mutate(quantity = nylon[matching_patients,'Nylon']) %>% 
        mutate(exposure = ifelse(quantity < 22.88, "low", "high"))}
    
    if(microplastic == "poly"){
      matching_patients <- intersect(row.names(poly),clinical_celldecon$Study.ID)
      clinical_celldecon <- clinical_celldecon[match(matching_patients,clinical_celldecon$Study.ID),] %>% 
        mutate(quantity = poly[matching_patients,'Polyurethane']) %>% 
        mutate(exposure = ifelse(quantity < 67.44, "low", "high"))}
    
  }
  
  
  if(batch == 2){
    total_mp <- total_mp2
    nylon <- nylon2
    poly <- poly2
    
    
    #Match up patients
    if(microplastic == "total_mp"){
      
      #SEO532 - no expression data , print sample names that donti ntersect
      print(paste0("Batch ",batch," samples do not have ", sample_type, " data", paste0(row.names(total_mp)[!row.names(total_mp) %in% clinical_celldecon$Study.ID], collapse = ",")))
            matching_patients <- intersect(row.names(total_mp),clinical_celldecon$Study.ID)
      clinical_celldecon <- clinical_celldecon[match(matching_patients,clinical_celldecon$Study.ID),] %>% 
        mutate(quantity = total_mp[matching_patients,'Total ng/ml of microplastics']) %>% 
        mutate(exposure = total_mp[matching_patients,'Category'])}
    
    if(microplastic == "nylon"){
      matching_patients <- intersect(row.names(nylon),clinical_celldecon$Study.ID)
      clinical_celldecon <- clinical_celldecon[match(matching_patients,clinical_celldecon$Study.ID),] %>% 
        mutate(quantity = nylon[matching_patients,'Nylon ng/ml']) %>% 
        mutate(exposure = nylon[matching_patients,'Category'])}
    
    
    if(microplastic == "poly"){
      matching_patients <- intersect(row.names(poly),clinical_celldecon$Study.ID)
      clinical_celldecon <- clinical_celldecon[match(matching_patients,clinical_celldecon$Study.ID),] %>% 
        mutate(quantity = poly[matching_patients,'Polypropylene ng/ml']) %>% 
        mutate(exposure = poly[matching_patients,'Category'])}
    
    clinical_celldecon$exposure <- tolower(clinical_celldecon$exposure) #make lowercase to match batch 1
  }
  
  row.names(clinical_celldecon) <- clinical_celldecon$Study.ID
boxplot_celldecon= cbind(id = row.names(clinical_celldecon),
                         sampletype = clinical_celldecon$sampletype,
                         classification = clinical_celldecon$classification,
                         clinical_celldecon[,which(colnames(clinical_celldecon) == "Alveolar.macrophages"):ncol(clinical_celldecon)])



boxplot_celldecon <- as.data.frame(boxplot_celldecon)

boxplot_celldecon_long <- pivot_longer(boxplot_celldecon, cols = 4:21)

ciber <-boxplot_celldecon_long
colnames(ciber)[9] <- "celltype"
colnames(ciber)[10] <- "proportion"


write.csv(ciber,file.path(celldecon.results.dir,paste0(microplastic,"_batch",batch,"_",sample_type, "celldecon_ciber.csv")))

#Normality test
# ggqqplot(residuals(lm(proportion ~ classification, data = ciber_biopt)))

### Get P Values --------------------------------------------------------------------------------
#Step 1) Get list of comparisons
celldecon_stats <- ciber %>%
  group_by(celltype) %>%
  wilcox_test(proportion ~ exposure,
              paired = FALSE
              # p.adjust.method = "BH" #dont need if we're just testing high vs low (one comparison per cell type group)
              )%>%
  add_significance("p") %>% 
  add_xy_position(x = "celltype")


#Manually change y positions because right now they are unecessarily staggered
#Don't want all the p-value brackets to be adjacent in the plot because the boxplots that are far down are then very far from the brackets up top
#Get max proportion for each cell type and stagger the brackets for each cell type
maxRows <- by(data = ciber, INDICES = ciber$celltype, FUN = function(X) {X[which.max(X$proportion),]})
#by() applies a funcion to subsets of a dataframe, grouped by a factor or variable (in this case cell type)
#celltype is the grouping factor, for each unique celltype, the function will be executed
#function(X) defines a function that is applied to every susbet of data (here, X is all rows corresponding to a specific celltype), those rows become X and the function looks for the max proportion in X
max_per_celltype <- as.data.frame(do.call("rbind", maxRows))

celldecon_stats$y.position <- as.numeric(max_per_celltype[match(celldecon_stats$celltype, max_per_celltype$celltype), "proportion"])
celldecon_stats$y.position <- celldecon_stats$y.position + (celldecon_stats$y.position)*0.05

celldecon_stats$contrast <- paste0(celldecon_stats$group1, "_" ,celldecon_stats$group2)

celldecon_stats <- as.data.frame(celldecon_stats)

#Only label the significant comparisons
celldecon_stats2 <- celldecon_stats[which(celldecon_stats$p.adj.signif != "ns"),]

png(file.path(celldecon.figures.dir,paste0(microplastic,"_batch",batch,"_",sample_type, "_celldecon_ciber.png")), width = 60, height = 30, units = "cm",res = 1200)

ggplot(ciber,
       aes(x=celltype,
           y=proportion))+
  geom_boxplot(aes(fill=exposure))+ #have to put fill here, because if we put it in the ggplot() aes, it becomes a global aes and needs 'exposure' in each subsequent step. theres no 'exposure' column in celldecon_stats so it doesnt work if its global
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5,
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour =  "white"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "grey"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.margin = margin(0,0,1,2, "cm"))+
  
  stat_pvalue_manual(celldecon_stats2,
                     label = "p.signif",
                     tip.length = 0.01,
                     # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                     #  bracket.shorten = 0.1,
                     size = 3)+
  scale_fill_manual(values=c("low" = "#619CFF" ,
                             "high" = "#F8766D")) +
  
  xlab("Cell Type") +
  ylab("Proportion")+
  labs(fill = paste(sample_type, ":",microplastic,"Exposure"))+
  
  # stat_compare_means(aes(group = classification),
  #                    comparisons = my_comparisons,
  #                    method = "wilcox.test")+
  
  ggtitle(paste0("Cellular Proportions (",sample_type,")"))

dev.off()
  } #close cell decon function




# ================================================================================== #
# 4.1 RUN THE FUNCTIONS ================================================================
# ================================================================================== #


for(sample_type in c("Brush", "Biopt")){
  print(sample_type)
  
for (batch in 1:2){
  
  mp.dir <- file.path(output.dir, "microplastics", paste0("data_batch", batch), sample_type)
  if(!exists(mp.dir)) dir.create(mp.dir, recursive = TRUE)
  
  celldecon_func(microplastic = "total_mp")
  celldecon_func(microplastic = "nylon")
  celldecon_func(microplastic = "poly")
  
} # close batch loop
  } #close sampletype  loop



cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


# #### edgeR ####
# 
# diffexp_edgeR <- function(microplastic) {
# 
#   
#   #Make output directory for total_mp or nylon or poly
#   this.mp.dir <- file.path(mp.dir, microplastic)
#   if(!exists(this.mp.dir)) dir.create(this.mp.dir, recursive = TRUE)
#   
#   diffexp.dir <- file.path(this.mp.dir, "diffexp")
#   if(!exists(diffexp.dir)) dir.create(diffexp.dir, recursive = TRUE)
#   
#   
#   #Set results paths
#   diffexp.results.dir <- file.path(diffexp.dir, "results")
#   if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)
#   
#   
#   diffexp.figures.dir <- file.path(diffexp.dir, "figures")
#   if(!exists(diffexp.figures.dir)) dir.create(diffexp.figures.dir, recursive = TRUE)
#   
#   
# 
#   if(microplastic == "total_mp"){
#   clinical <- clinical_total_mp}
#   
#   if(microplastic == "nylon"){
#   clinical <- clinical_nylon}
#   
#   
#   if(microplastic == "poly"){
#   clinical <- clinical_poly}
#   
#   
#   counts <- counts[,row.names(clinical)]
#   
# 
# 
#   design <- model.matrix(~ 0 + exposure + age + sex + smoking_status, 
#                           data = clinical) #Design matrix
#   
#   colnames(design)[1:2] <- c(levels(as.factor(clinical$exposure)))
#   
#   DGEL<- DGEList(counts=counts) 
#   
#   # FILTER
#   keep <- filterByExpr(DGEL, group = clinical$classification) 
#   # The filterBy Expr function keeps rows that have worthwhile counts in a minimum number of samples.
#   # keep <- which(rowMedians(as.matrix(DGEL))>10) 
#   # keep <- rowSums(cpm(expression)>100) >= 2
#   
#   
#   DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
#   # expression <- DGEL$counts
#   
#   # NORMALISE
#   DGEL<- calcNormFactors(DGEL,method = "TMM")
#   
#   # ESTIMATE DISPERSON
#   DGEL <- estimateDisp(DGEL, design)
#   
#   # FIT MODEL
#   fit <- glmQLFit(DGEL, design, robust = TRUE) #fit the GLM (design) to the DGEL(the DGEL object, which contains the counts data that has been filtered,normalised and dispersons estimtated)
#   
#   my.contrasts <- makeContrasts(contrast1 = high - low,
#                                 levels = design)
#   
#   qlf <- glmQLFTest(fit, contrast=my.contrasts[,1], robust)
#   
#   tT <- topTags(qlf,n=nrow(DGEL))$table #topTags gets top genes, here we want all of the genes, edgeR's default p.adjust method is BH
#   tT$Legend <- ifelse(
#     tT$FDR < 0.05 & tT$logFC > 0, "Upregulated",
#     ifelse(
#       tT$FDR < 0.05 & tT$logFC < 0, "Downregulated",
#       "Not Significant"))
#   
#   tT$Legend[is.na(tT$Legend)]="Not Significant"
#   
#   tT$Legend <- factor(tT$Legend, levels = c("Downregulated", "Upregulated", "Not Significant"))
#   
#   tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols
#   
#   # if(showEnsemblID == TRUE){
#     #for those with no hgnc symbol, label with ensembl id
#     tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl
#   # }
#   # 
#   # else{
#   #   #for those with no hgnc symbol, remove
#   #   tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
#   # }
#   
#   selection <-which(tT$FDR<0.05)
#   
#   tT2 <- tT[selection,]
#   
#   write.csv(tT, file = paste0(diffexp.results.dir, "tT_", microplastic, ".csv"))
#   write.csv(head(tT, n = 100) , file = paste0(diffexp.results.dir, "tT_top100_", microplastic, ".csv"))
#   
#   
#   # ================================================================================== #
#   # 2.1. VOLCANO PLOT ================================================================
#   # ================================================================================== #
#   cat(paste("Starting 2.1. VOLCANO PLOT", microplastic), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
#   
#   volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
#     ggtitle(paste(microplastic)) +
#     geom_point(aes(color = Legend)) +
#     scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
#     geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
#     geom_text_repel(data = subset(tT2[1:30,]),
#                     aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
#                     point.padding = unit(0.3, "lines") ) +
#     theme_bw(base_size = 18) + theme(legend.position = "bottom",
#                                      legend.text = element_text(size = 14),
#                                      legend.title = element_text(size = 16)) 
#   
#   ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(microplastic,"_volcano_plot.png")),
#          width = 25, height = 25,
#          units = "cm")
#   
#   
#   listofresults <- list(tT = tT, tT2 = tT2, volcano = volcano)
#   #Save all results
#   saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))
#   
# 
# } #close function
# 
# 
# diffexp_edgeR(microplastic = "total_mp")
# diffexp_edgeR(microplastic = "nylon")
# diffexp_edgeR(microplastic = "poly")
