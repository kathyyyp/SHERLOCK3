# SHERLOCK3 - QC on merged SHERLOCK2&3 data and batch corrected with Combat_seq 
options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"


library("readxl")
library("limma")
library("rstatix")
library("tibble")
library("ggvenn")
library("ggplot2")
library("ggrepel")
library("ggfortify")
library("stringr")
library("EnsDb.Hsapiens.v79")
library("ggpubr")
library("edgeR")
library("DESeq2")
library("tidyverse")
library("PCAtools")
library("sva")


# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir,"data")
if(!exists(data.dir)) dir.create(data.dir, recursive = TRUE)

processed.data.dir <- file.path(data.dir,"processed")
if(!exists(processed.data.dir)) dir.create(processed.data.dir)

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")

#Output directory
output.dir <- file.path(main.dir,"output")
if(!exists(output.dir)) dir.create(output.dir, recursive = TRUE)

#SHERLOCK2 dir
sherlock2.dir <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK2"

#QC directory
qc.dir <- file.path(output.dir, "qc")
if(!exists(qc.dir)) dir.create(qc.dir, recursive = TRUE)

setwd(file.path(main.dir))

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

# ##-- Pre batch correction
counts_brush_merged <- readRDS(file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat.rds"))
counts_biopt_merged <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_merged_pre_combat.rds"))

# ##-- Post batch correction
counts_brush_combat <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat_raw.rds"))
counts_biopt_combat <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat_raw.rds"))

clinical_brushbiopt <- readRDS(file.path(postQC.data.dir, "clinical_brushbiopt.rds")) #552 samples
clinical_brush <- readRDS(file.path(postQC.data.dir, "clinical_brush.rds"))
clinical_biopt <- readRDS(file.path(postQC.data.dir, "clinical_biopt.rds"))

#Combine brush and biopt for first PCA
counts_brushbiopt_merged <- as.matrix(cbind(counts_brush_merged, counts_biopt_merged))
counts_brushbiopt_combat <- as.matrix(cbind(counts_brush_combat, counts_biopt_combat))
saveRDS(counts_brushbiopt_combat, file.path(combat.processed.data.dir, "counts_brushbiopt_combat.rds"))
write.csv(counts_brushbiopt_combat, file.path(combat.processed.data.dir, "counts_brushbiopt_combat.csv"))



clinical_brushbiopt$batch <- as.factor(clinical_brushbiopt$batch)
clinical_brush$batch <- as.factor(clinical_brush$batch)
clinical_biopt$batch <- as.factor(clinical_biopt$batch)


# ================================================================================== #
# 2. PRINCIPLE COMPONENT ANALYSIS =======================================
# ================================================================================== #
### Create FUNCTIONS

cat("Starting 2. PRINCIPLE COMPONENT ANALYSIS",format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

pca_pipeline <- function(clinical_file, counts_file, directory_name){
  cat("Starting 2. PRINCIPLE COMPONENT ANALYSIS - ", directory_name, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  
  this.qc.dir <- file.path(qc.dir, directory_name)
  if(!exists(this.qc.dir)) dir.create(this.qc.dir, recursive = TRUE)
  
  pca.dir <- file.path(this.qc.dir, "pca")
  if(!exists(pca.dir)) dir.create(pca.dir, recursive = TRUE)
  if(!exists(file.path(pca.dir,"pca_5"))) dir.create(file.path(pca.dir,"pca_5"), recursive = TRUE)
  
  pca.top500.dir <- file.path(this.qc.dir, "pca_top500")
  if(!exists(pca.top500.dir)) dir.create(pca.top500.dir)
  if(!exists(file.path(pca.top500.dir,"pca_5"))) dir.create(file.path(pca.top500.dir,"pca_5"), recursive = TRUE)
  
  clinical <- clinical_file
  counts <- counts_file
  
  if(all(row.names(clinical) == colnames(counts)) == FALSE){
    stop("all(row.names(clinical) == colnames(counts)) == FALSE") }
  
  ## Voom normalise counts ---------------------------------------------------------
  counts_norm <- voom(counts) #voom is log2cpm values
  counts_norm <- counts_norm$E
  
  ## Top 500 highly variable genes -----------------------------------------------------
  gene_variance <- apply(counts_norm, 1, var) #1 means apply over rows, 2 means columns #variance calculated by taking the differences between each number in the data set and the mean, squaring the differences to make them positive, and then dividing the sum of the squares by the number of values in the data set.
  top_genes <- order(gene_variance, decreasing = TRUE)[1:500]
  top_counts_norm <- as.matrix(counts_norm[top_genes, ])
  
  
  
  ## Plot 1: EIGENCORPLOT -----------------------------------------------
  cat("Starting Plot 1: Eigencorplot", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # main variables
  clinical_numeric <- sapply(clinical[,c("batch",
                                         "sampletype",
                                         "sex",
                                         "smoking_status",
                                         "classification",
                                         "packyears",
                                         "ics_use")],
                             function(col){
                               as.numeric(as.factor(col))
                             }
  )
  
  clinical_numeric <- as.matrix(cbind(clinical_numeric, age = clinical$age))
  #batch
 
  if(length(unique(clinical_numeric[,"sampletype"])) == 1){
    clinical_numeric <- clinical_numeric[,-which(colnames(clinical_numeric) == "sampletype")]
    }
  
  row.names(clinical_numeric) <- row.names(clinical)
  
  
  ## All genes
  eigen.pca <- pca(counts_norm,
                   metadata = clinical_numeric,
                   center = TRUE,
                   scale = TRUE)
  
  png(filename = file.path(pca.dir, paste0("eigencorplot_all.png")),
      width = 3000, height = 1500, units = "px", res = 300)
  
  eigencorplot(eigen.pca, metavars = c(colnames(clinical_numeric)))
  
  dev.off()
  
  
  # Top 500 variable genes
  eigen.pca <- pca(top_counts_norm,
                   metadata = clinical_numeric,
                   center = TRUE,
                   scale = TRUE)
  
  png(filename = file.path(pca.top500.dir, paste0("eigencorplot_top500.png")),
      width = 3000, height = 1500, units = "px", res = 300)
  
  eigencorplot(eigen.pca, metavars = c(colnames(clinical_numeric)))
  
  dev.off()
  
  
  
  
  ### PCA plots ---------------------------------------------------------------------------------------
  listofcounts <- list(all = counts_norm, 
                       top500 = top_counts_norm)
  listofdir <- list(all = pca.dir, top500 = pca.top500.dir)
  
  ## LOOP PCA for i) all genes and ii)top 500 genes ===============================================
  for (i in names(listofcounts)){
    
    counts_norm <- listofcounts[[i]]
    this.dir <- listofdir[[i]]
    
    
    ## PCA Analysis  ===================================================================================================================
    pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
    
    
    # plot settings
    fontsize <- theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
    
    ## Plot 2: Percentage variance explained ============================================================================================
    cat("Starting Plot 2: Percentage variance explained", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    eigenvalues <- pca_res$sdev^2
    var_explained = pca_res$sdev^2 / sum(pca_res$sdev^2)
    
    png(filename = paste0(this.dir,"/var_explained.png"),
        width = 2400, height = 1500, units = "px", res = 300)
    
    plot(var_explained[1:20], type = "b",
         xlab = "Principal Component",
         ylab = "Percentage of Variance Explained",
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
    
    dev.off()
    
    
    ## Plot 3: Individual PCA plots (1 PC) coloured by clinical variable ============================================================================================
    cat("Starting Plot 3:  Individual PCA plot (1 PC)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    library("viridis")
    
    # main variables
    for (variable in c("smoking_status", "sex", "classification", "sampletype", "batch", "ics_use")){
      ggsave(autoplot(pca_res,data=clinical, colour = variable)+ fontsize  + guides(color=guide_legend(variable)),
             file = file.path(this.dir,paste0("pca_",variable,".png")),
             width = 2300,
             height = 1500,
             units = "px")
    }
    
    #age and packyears (colour gradient)
    for (variable in c("age", "packyears")){

    ggsave(
      autoplot(pca_res, data = clinical, colour = variable) +  fontsize + guides(color=guide_legend(variable))+
        scale_color_viridis_c(option = "plasma") ,
      
      file = file.path(this.dir,paste0("pca_",variable,".png")),
      width = 2200,
      height = 1500,
      units = "px")
    }
    
    ## labelled sampletype
    ggsave(autoplot(pca_res,data=clinical, colour = "sampletype", label = TRUE)+ fontsize + guides(color=guide_legend("Sample Type")),
           file =paste0(this.dir,"/pca_sampletype_labelled.png"),
           width = 2200,
           height = 1500,
           units = "px")

    

    
    
    ## Plot 4: Pairs plot (5 PCs) coloured by clinical variable ============================================================================================
    cat("Starting Plot4: Pairs plot (5 PCs) ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    clinical_pca <- cbind(pca_res$x, clinical)
    # clinical_pca <- cbind(eigen.pca$rotated, clinical) #pca() does the same thing as prcomp()
    
    #main variables
    for (i in c("smoking_status", "sex", "classification", "sampletype", "batch", "packyears", "ics_use")){
      png(filename = paste0(this.dir,"/pca_5/", i, ".png"), width = 26, height = 23, units = "cm", res = 1200)
      
      variable <- as.factor(clinical_pca[,i])
      
      pcaplot <- pairs(clinical_pca[,c(1:5)],
                       pch = 19,
                       cex = 0.75,
                       oma=c(3,3,3,20),
                       col = variable)
      
      par(xpd=TRUE)
      legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = i)
      dev.off()
    }
    
    # Age - need to break into categories
    png(filename = paste0(this.dir,"/pca_5/", "age", ".png"), width = 26, height = 23, units = "cm", res = 1200 )
    
    variable <- clinical_pca[,"age"]
    
    colors <- colorRampPalette(c("#F0F921", "#FC8C3F", "#F8766D", "#C51B7D","#5E4FA2"))(100)
    
    pcaplot <- pairs(clinical_pca[,c(1:5)],
                     pch = 19,
                     cex = 0.75,
                     oma=c(3,3,3,20),
                     col =  colors[as.numeric(cut(clinical_pca[,"age"], breaks = 100))])
      
    breaks <- seq(min(clinical_pca[,"age"]), max(clinical_pca[,"age"]), length.out = 9)
    
    # Add custom color legend
    par(xpd=TRUE)
    
    legend("bottomright",
           legend = round(breaks, 1),
           fill = colorRampPalette(c("#F0F921", "#FC8C3F", "#F8766D", "#C51B7D","#5E4FA2"))(9),
           title = "age",
           border = "white")
    
    dev.off()

    # packyears - need to break into categories
    png(filename = paste0(this.dir,"/pca_5/", "packyears", ".png"), width = 26, height = 23, units = "cm", res = 1200 )
    
    variable <- clinical_pca[,"packyears"]
    
    
    # Custom bins
    breaks <- c(0, 10, 40, 60, 80, Inf)       
    labels <- c("<10", "10-39", "40-59", "60-79", "80+")
    
    bins <- cut(variable, breaks = breaks, labels = labels, right = FALSE)
    
    colors <- c("#F0F921", "#FC8C3F", "#F8766D", "#C51B7D", "#5E4FA2")
    cols <- colors[as.numeric(bins)]
    
    pcaplot <- pairs(clinical_pca[,c(1:5)],
                     pch = 19,
                     cex = 0.75,
                     oma=c(3,3,3,20),
                     col =  cols )
    
   
    # Add custom color legend
    par(xpd=TRUE)
    
    legend("bottomright",
           legend = labels,
           fill = colors,
           title = "packyears",
           border = "white")
    
    dev.off()
  } #close loop
  
} #close function


# PCA on the pre-batch corrected data brush + biopt to show the sampletype split
## clinical = clinical_brushbiopt
## counts = counts_brushbiopt_merged

# PCA on merged data - brush
## clinical_brush
## counts_brush_merged

# PCA on merged data - biopt
## clinical_biopt
## counts_biopt_merged

# PCA on combat data - brush
## clinical_brush
## counts_brush_combat

# PCA on combat data - biopt
## clinical_biopt
## counts_biopt_combat

pca_pipeline(clinical_file = clinical_brushbiopt, counts_file = counts_brushbiopt_merged, directory_name = "precombat_brushbiopt")
pca_pipeline(clinical_file = clinical_brushbiopt, counts_file = counts_brushbiopt_combat, directory_name = "postcombat_brushbiopt")
pca_pipeline(clinical_file = clinical_brush, counts_file = counts_brush_merged, directory_name = "brush_merged")
pca_pipeline(clinical_file = clinical_biopt, counts_file = counts_biopt_merged, directory_name = "biopt_merged")
pca_pipeline(clinical_file = clinical_brush, counts_file = counts_brush_combat, directory_name = "brush_combat")
pca_pipeline(clinical_file = clinical_biopt, counts_file = counts_biopt_combat, directory_name = "biopt_combat")




# ================================================================================== #
# 3. CONVERT ENSEMBL ID TO HGNC SYMBOL =====================
# ================================================================================== #
cat("Starting 3. CONVERT ENSEMBL ID TO HGNC SYMBOL", Sys.time(), "\n")

# ##### with gProfiler online ##### 36,239
# hgnc_symbols_gprofiler <- read.csv("processed/gprofiler_hgnc.csv")
# length(which(hgnc_symbols_gprofiler$converted_alias != "None"))
#
# ##### with syngo ##### 38,371
# hgnc_symbols_syngo <- read.table("processed/target_symbols.txt")

my_ensembl_gene_ids <- row.names(counts_brushbiopt_merged) #35,111 biomart
#
##### #with library(EnsDb.Hsapiens.v79) ##### 46,666
hgnc_symbols_db <- ensembldb::select(EnsDb.Hsapiens.v79, keys= my_ensembl_gene_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
row.names(hgnc_symbols_db) <- hgnc_symbols_db$GENEID

n_occur <- data.frame(table(hgnc_symbols_db$SYMBOL)) #Get frequencies of each gene/hgnc_symbol
n_occur[n_occur$Freq > 1,]

saveRDS(hgnc_symbols_db, file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
write.csv(hgnc_symbols_db, file.path(postQC.data.dir,"hgnc_symbols_db.csv"))



# ================================================================================== #
# 4. PATIENT DEMOGRAPHICS (SHERLOCK 2 and 3 combined) ==============================
# ================================================================================== #
cat("Starting 4. PATIENT DEMOGRAPHICS (SHERLOCK 2 and 3 combined)", Sys.time(), "\n")

# Patient demographics
unique_patients_all <- clinical_brushbiopt[match(unique(clinical_brushbiopt$Study.ID), clinical_brushbiopt$Study.ID),]

patient_demographics <- t(unique_patients_all %>% 
                            mutate(packyears = as.numeric(as.character(packyears))) %>% 
                            
                            group_by(classification) %>% 
                            
                            summarise(
                              
                              #Total patients
                              total_patients = n(),
                              
                              #sex
                              male_patients = sum(sex == "Male"), 
                              male_percentage = (male_patients / total_patients) * 100,
                              sex = paste0(male_patients, "(", round(male_percentage, digits=3), ")"),
                              
                              #age
                              mean_age = median(age),
                              range_age =paste0(min(age, na.rm = T), "-", max(age, na.rm = T)),
                              age = paste0(mean_age,"(",range_age, ")"),
                              
                              #Smoking status
                              currentsmoker = sum(smoking_status == "Current.smoker"),
                              smokerpercentage = currentsmoker/total_patients *100,
                              smoke = paste0(currentsmoker, "(", round(smokerpercentage, digits = 3), ")"),
                              
                              #packyears
                              median_packyears = median(packyears, na.rm = TRUE),
                              range_packyears = paste0(min(packyears, na.rm = T), "-", max(packyears, na.rm = T)),
                              packyears = paste0(median_packyears, "(", range_packyears, ")"),
                              
                              
                              #FEV1
                              median_FEV1_percent_pred =  median(FEV1_percent_pred, na.rm = T),
                              range_postfev1percpred = paste0(round(min(FEV1_percent_pred, na.rm = T),digits=3), "-", round(max(FEV1_percent_pred, na.rm = T), digits = 3)),
                              fev1 = paste0(round(median_FEV1_percent_pred, digits = 3), "(", range_postfev1percpred, ")"),
                              
                              #FEV/FVC
                              median_postfev1fvcpercpred =  median(FEV1_FVC_post, na.rm = TRUE),
                              range_postfev1fvcpercpred = paste0(round(min(FEV1_FVC_post, na.rm = T)), "-", round(max(FEV1_FVC_post, na.rm = T), digits = 3 )),
                              fev_fvc = paste0(round(median_postfev1fvcpercpred, digits = 3), "(", range_postfev1fvcpercpred, ")"),
                            )
                          %>%
                            dplyr::select(
                              classification,
                              total_patients,
                              sex,
                              age,
                              smoke,
                              packyears,
                              fev1,
                              fev_fvc
                            )
)

colnames(patient_demographics) <- patient_demographics[1,]
patient_demographics <- patient_demographics[-1,]
patient_demographics <- rbind(as.numeric(table(clinical_brushbiopt$classification)),patient_demographics)

row.names(patient_demographics) <- c(
  "Samples, n",
  "Patients, n",
  "Sex Male, n (%)",
  "Age, median (Range)",
  "Current smoker, n (%)",
  "Packyears, median (Range)",
  "FEV1 % pred. (post-bronchodilater), median (Range)",
  "FEV1/FVC % (post-bronchodilater), median (Range)"
)

write.csv(patient_demographics, file = file.path(qc.dir, "patient_demographics_postcombat.csv"))
saveRDS(patient_demographics, file = file.path(qc.dir, "patient_demographics_postcombat.rds"))


# ================================================================================== #
# 5. QC - Library size, gene counts and sex check ==================================
# ================================================================================== #
cat("Starting 5.  QC - Library size, gene counts and sex check", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

## Library size: colsums/total reads per sample  ------------------------------------------------------------------------------
colsums_persample=as.matrix(colSums(counts_brushbiopt_combat))
write.csv(colsums_persample, file = file.path(qc.dir,"sk2sk3_colsums.csv"))
# (lowest is 954,000)

## Number of detected genes ------------------------------------------------------------------------------
genes_per_sample = colSums(counts_brushbiopt_combat >0) # same as sapply(counts_brushbiopt_combat, function(x) sum(x >0)) 
write.csv(genes_per_sample, file = file.path(qc.dir,"sk2sk3_colsums_genes_detec_per_sample.csv"))

## Sex checks (this step is probably not useful?? sex checks were done for each batch seperately, as mixups would have happened per batch ------------------------------------------------------------------------------
# cat("Starting sex check", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# 
# sex.dir <- file.path(qc.dir, "sex", "s2k_sk3")
# if(!exists(sex.dir)) dir.create(sex.dir, recursive = TRUE)
# 
# # male - DDX3Y (ENSG00000067048), XIST (ENSG00000229807)
# # female - XIST
# counts_cpm <- cpm(counts_brushbiopt_combat)
# xist_counts <- counts_cpm["ENSG00000229807",] #XIST
# ddx3y_counts <- counts_cpm["ENSG00000067048",] #DDX3Y
# 
# ddx3y_samples <- colnames(counts_cpm)[counts_cpm["ENSG00000067048",] != 0]
# clinical_brushbiopt[ddx3y_samples, "sex"]
# 
# 
# #Raw counts
# gendercheck_counts <- as.data.frame(cbind(t(counts_brushbiopt_combat[c("ENSG00000229807","ENSG00000067048"),]), clinical_brushbiopt$sex))
# colnames(gendercheck_counts) <- c("XIST", "DDX3Y", "Sex")
# gendercheck_counts$Sample <- row.names(gendercheck_counts)
# 
# png(filename = file.path(sex.dir,"XIST_vs_DDX3Y_rawcounts.png"), width = 15, height = 12, unit = "cm", res = 1200)
# ggplot(data = gendercheck_counts,
#        aes(x = as.numeric(DDX3Y),
#            y = as.numeric(XIST)
#        )) + 
#   theme_bw() + 
#   geom_point(aes(color = Sex)) +
#   # geom_text(aes(label = Sample), angle = 90, vjust = 0.5, hjust = -0.2, size = 0.8) +   # add labels
#   ylab(lab = "XIST Counts") +
#   xlab(lab = "DDX3Y Counts")+
#   ggtitle("Sex Check")
# 
# dev.off()
# 
# 
# #CPM
# gendercheck_cpm <- as.data.frame(cbind(t(counts_cpm[c("ENSG00000229807","ENSG00000067048"),]), clinical_brushbiopt$sex))
# colnames(gendercheck_cpm) <- c("XIST", "DDX3Y", "Sex")
# gendercheck_cpm$Sample <- row.names(gendercheck_cpm)
# 
# png(filename = file.path(sex.dir,"XIST_vs_DDX3Y_cpm.png"), width = 15, height = 12, unit = "cm", res = 1200)
# ggplot(data = gendercheck_cpm,
#        aes(x = as.numeric(DDX3Y),
#            y = as.numeric(XIST))) + 
#   theme_bw() + 
#   geom_point(aes(color = Sex)) +
#   # geom_text(aes(label = Sample), angle = 90, vjust = 0.5, hjust = -0.2, size = 0.8) +   # add labels
#   ylab(lab = "XIST CPM") +
#   xlab(lab = "DDX3Y CPM")+
#   ggtitle("Sex Check")
# 
# dev.off()
# 
# 
# gendercheck <- merge(gendercheck_counts, gendercheck_cpm, by = "Sample")
# gendercheck <- gendercheck[,-which(colnames(gendercheck) == "Sex.x")]
# colnames(gendercheck) <- c("Sample", "XIST_count", "DDX3Y_count", "XIST_cpm", "DDX3Y_cpm", "Sex")
# gendercheck <- gendercheck[order(gendercheck$DDX3Y_cpm, decreasing = TRUE),]
# write.csv(gendercheck, file.path(sex.dir, "gendercheck_table.csv"))
# 
# # 107165-001-146 cpm = 3.13 rawcount = 58 DDX37
# # 107165-001-184 cpm = 2.82 rawcount = 5 DDX3Y ? will not remove for now as raw count is low, cpm only high because total reads are low for the sample?
# 
# 
# #XIST
# png(filename = file.path(sex.dir,"XIST_cpm.png"), width = 15, height = 12, unit = "cm", res = 1200)
# ggplot(data = gendercheck,
#        aes(x = as.factor(Sex),
#            y = as.numeric(XIST_cpm),
#            fill = Sex)) + 
#   theme_bw() + 
#   geom_boxplot() +
#   ylab(lab = "XIST CPM") +
#   xlab(lab = "Sex")+
#   ggtitle("XIST")
# 
# dev.off()
# 
# #DDX3Y
# png(filename = file.path(sex.dir,"DDX3Y_cpm.png"), width = 15, height = 12, unit = "cm", res = 1200)
# ggplot(data = gendercheck,
#        aes(x = as.factor(Sex),
#            y = as.numeric(DDX3Y_cpm),
#            fill = Sex)) + 
#   theme_bw() + 
#   geom_boxplot()+
#   # geom_dotplot(binaxis = "y",
#   #              stackdir = "center",
#   #              stackratio = 0.3,
#   #              dotsize = 0.55) +
#   ylab(lab = "DDX3Y CPM") +
#   xlab(lab = "Sex") +
#   ggtitle("DDX3Y")
# dev.off()



cat("END OF THIS JOB", Sys.time(), "\n")

