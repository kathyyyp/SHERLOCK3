# SHERLOCK3 - integration of new data with SHERLOCK2
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

#Output directory
output.dir <- file.path(main.dir,"output")
if(!exists(output.dir)) dir.create(output.dir, recursive = TRUE)

qc.dir <- file.path(output.dir, "qc")
if(!exists(qc.dir)) dir.create(qc.dir, recursive = TRUE)

#SHERLOCK2 dir
sherlock2.dir <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK2"

setwd(file.path(main.dir))

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(data.dir))

#counts
raw_expression <- read.table(file.path("raw","20250902_Sherlock3_readcount_UMI_dedup.txt"), header = TRUE, check.names = FALSE, row.names = 1)
raw_clinical_sk3 <- read_xlsx(file.path("raw","25_08_20_SHERLOCk_1_2 mvdb.xlsx"), .name_repair = "universal") 
clinical_sk3_ids <- read_xlsx(file.path("raw", "Copy of Sample Submission Form _107165 (1).xlsx"), sheet = "sample_submission_info", col_names = TRUE, .name_repair = "universal") # unique SEO/patient IDs
# 'sample_submission_info' is copied from 'Isolation samples' sheet, last 58 samples had wrong prefix, Daan sent 'Kathy_GenomeScan_Codes.xlsx' IDs with correct prefixes, which I replaced in 'sample_submission_info'

setwd(file.path(main.dir))



# ================================================================================== #
# 2. CLEAN UP DATA ======================================================
# ================================================================================== #
cat("Starting 2. CLEAN UP DATA", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# COPD: FEV1/FVC < 70
# GOLD 1 (Mild): FEV1percpred >= 80%
# GOLD 2 (Moderate): 50% <= FEV1percpred < 80%
# GOLD 3 (Severe): 30% <= FEV1percpred < 50%
# GOLD 4 (Very severe): FEV1percpred < 30%

# For our classifications
# Mild/Moderate COPD =  FEV1percpred > 50%
# Severe COPD = FEV1percpred < 50%


# Remove uneccessary rows -------------------------
clinical_sk3_ids <- as.data.frame(clinical_sk3_ids[,-c(1,2)])
raw_clinical_sk3 <- as.data.frame(raw_clinical_sk3)

# Check data ----------------------------------------
colnames(raw_expression)
head(as.data.frame(raw_clinical_sk3))
tail(as.data.frame(raw_clinical_sk3))
head(raw_expression)
tail(raw_expression)

sapply(raw_clinical_sk3, function(x) sum(is.na(x)))
# There are 5 samples with no gold classifcation

sapply(raw_expression, function(x) sum(is.na(x)))
any(raw_expression == "NA")
any(raw_expression == " ")

sapply(clinical_sk3_ids, function(x) sum(is.na(x)))


# Add row for sample type - biopsy or brushing
clinical_sk3_ids$sampletype <- word(clinical_sk3_ids$Remarks, 2)
clinical_sk3_ids$sampletype <- sub("BRONCHOSCOPIE-","",clinical_sk3_ids$sampletype)
clinical_sk3_ids[which(clinical_sk3_ids$sampletype == "Broch.brush"), "sampletype"] <- "Brush"

# Remve unnecessary columns
clinical_sk3_ids <- clinical_sk3_ids[,c("GS_ID", "Customer.ID", "sampletype")]

dim(raw_expression) #66,138 genes #226 samples
dim(raw_clinical_sk3) #626 samples
dim(clinical_sk3_ids) #226 samples

counts <- raw_expression 
colnames(raw_clinical_sk3)[colnames(raw_clinical_sk3) == "class_incl_study_id"] <- "Study.ID"

# Sanity checks
all(colnames(counts) %in% clinical_sk3_ids$GS_ID) #False
intersect(colnames(counts),clinical_sk3_ids$GS_ID) #331 samples total

# Check any mismatch between GenomeScan IDs and sample ID names in counts file
setdiff(colnames(counts),clinical_sk3_ids$GS_ID)
setdiff(clinical_sk3_ids$GS_ID, colnames(counts))

# Check any IDs in clinical_sk3 file that are missing in clinical_sk3_ids
setdiff(clinical_sk3_ids$GS_ID, row.names(raw_clinical_sk3))
setdiff(row.names(raw_clinical_sk3), clinical_sk3_ids$GS_ID)

all(clinical_sk3_ids$Customer.ID %in% raw_clinical_sk3$Study.ID)

# ================================================================================== #
# 2.1 Resolve duplicated patients in raw_clinical_sk3 =====================================
# ================================================================================== #
cat("Starting 2.1. Resolve duplicated patients in raw_clinical_sk3", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# Resolve duplicate patients in SHERLOCK3 clinical_sk3 fle
dup_ids <- raw_clinical_sk3[duplicated(raw_clinical_sk3$Study.ID), "Study.ID"]

# Contains all duplicate rows (all instances)
raw_clinical_sk3_dups <- raw_clinical_sk3[raw_clinical_sk3$Study.ID %in% dup_ids,]
dim(raw_clinical_sk3_dups) #50 rows = 25 duplicates
# Merge duplicate rows
merged_clinical_sk3_dups <- raw_clinical_sk3_dups %>%
  group_by(Study.ID) %>% #create a group for each patient
  # In each group,
  summarise(#across applies same function (that we define) to all columns in the dataset
    across(everything(), #summarise/collapse each group to one row (without this, output will still be 50 rows, with 2 idenical rows for each patient)
                   function(column_values) { #create a function. eg. column_values would be c(21,21) for the age column for patient SEO23, think abt this like sapply(df, function(x) x+1), where x is every value in the column of df)
                     vals <- unique(na.omit(column_values)) #remove NA values from the column values for that patient and remove duplicates among remaining values (take only unique values) eg. 21 (for previous example)
                     if (length(vals) == 0) NA #if all values (vals) were empty, return NA
                     else paste(vals, collapse = "; ") #Otherwise, concatenate into a string seperated by ; (eg if column_values was c(21,22), output would be 21;22)
                     } #close function
                   ), 
            .groups = "drop") #remove group_by


#to visualise
# write.csv(merged_clinical_sk3_dups, file.path(data.dir, "processed","merged_clinical_sk3_dups.csv"))

# Remove all duplicate rows
raw_clinical_sk3_nondups <- raw_clinical_sk3[-which(raw_clinical_sk3$Study.ID %in% dup_ids),]
dim(raw_clinical_sk3_nondups) #576 616 #50 dups removed

# Add the merged duplicate row back in
raw_clinical_sherlock_full <- rbind(raw_clinical_sk3_nondups,  merged_clinical_sk3_dups)
raw_clinical_sherlock_full <- raw_clinical_sherlock_full[order(raw_clinical_sherlock_full$Study.ID),]
dim(raw_clinical_sherlock_full) #601 616 (25 added back in)

# SAVE MASTER SHERLOCK CLINICAL FILE ===============================================================
write.csv(raw_clinical_sherlock_full, file.path(data.dir, "processed","clinical_sherlock_full_master_preQC.csv"))

# ================================================================================== #
# 2.2 Subset for SHERLOCK3 patients only ===========================================
# ================================================================================== #
cat("Starting 2.2. Subsetting for SHERLOCK3 patients only", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# Include only SHERLOCK3 patients
clinical_sk3_ids <- clinical_sk3_ids[order(clinical_sk3_ids$Customer.ID),]

# > table(clinical_sk3_ids$sampletype)
# 
# Biopt Brush
# 114   112

row.names(raw_clinical_sherlock_full) <- raw_clinical_sherlock_full$Study.ID
all(colnames(counts) %in% clinical_sk3_ids$GS_ID)

#Reorder counts to match clinical_sk3_ids GS_ID
counts <- counts[clinical_sk3_ids$GS_ID]
counts_sk3 <- counts

# Pivot the clinical_sk3 table long to include all samples (clinical_sk3_ids Customer.ID has all samples with their matching patient IDs)
clinical_sk3_master <- raw_clinical_sherlock_full[clinical_sk3_ids$Customer.ID,]
row.names(clinical_sk3_master) <- clinical_sk3_ids$GS_ID
clinical_sk3_master[,1:3] #226 616
clinical_sk3_master <- cbind(clinical_sk3_master, sampletype = clinical_sk3_ids$sampletype)

row.names(clinical_sk3_master) == colnames(counts_sk3)

clinical_sk3 <- clinical_sk3_master #incase I want to subset columns (currently have 600+ clinical_sk3 variables)


# ================================================================================== #
# 2.3. Save SHERLOCK3 counts and clinical_sk3 ===========================================
# ================================================================================== #
cat("2.3. Saving SHERLOCK3 counts and clinical_sk3 info", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

write.csv(counts_sk3, file.path(processed.data.dir, "counts_sk3_preQC.csv"))
saveRDS(counts_sk3, file.path(processed.data.dir, "counts_sk3_preQC.rds"))


write.csv(clinical_sk3_master, file.path(processed.data.dir, "clinical_sk3_master_preQC.csv"))
saveRDS(clinical_sk3_master, file.path(processed.data.dir, "clinical_sk3_master_preQC.rds"))


# ================================================================================== #
# 3.1 Subset main clinical_sk3 file ================================================
# ================================================================================== #
cat("Starting 3. Subset main clinical_sk3 file", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
length(unique(clinical_sk3_master$Study.ID))

clinical_sk3 <- clinical_sk3_master %>% 
  dplyr::select(
  Study.ID,
  age,
  crf_gender,
  crf_smoking,
  crf_packyears,
  crf_corticosteroid,
  ics_use,
  ics_name_factor,
  lama_use,
  laba_use,
  sex_numeric,
  FEV1, #note this is the same as postbodybox_fev1_post
  FEV1_pred,
  FEV1_percent_pred,
  postbodybox_fvc_post,
  postbodybox_fev1_fvc_post, #this is postbodybox_fev1_post/postbodybox_fvc_post
  gold_classification,
  classification,
  sampletype

  ) %>% 
  dplyr::rename(
  sex = crf_gender,
  smoking_status = crf_smoking,
  packyears = crf_packyears,
  corticosteroid = crf_corticosteroid,
  FVC_post= postbodybox_fvc_post ,
  FEV1_FVC_post = postbodybox_fev1_fvc_post
  )
  

clinical_sk3$classification <- make.names(clinical_sk3$classification)
clinical_sk3$smoking_status <- make.names(clinical_sk3$smoking_status)
clinical_sk3[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(clinical_sk3[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], function(x) as.numeric(x))


## Save pre_QC clinical_sk3 (main patient info, smoking, ics and lung function data)
write.csv(clinical_sk3, file.path(processed.data.dir, "clinical_sk3_preQC.csv"))
saveRDS(clinical_sk3, file.path(processed.data.dir, "clinical_sk3_preQC.rds"))


# ================================================================================== #
# 3. PATIENT DEMOGRAPHICS (SK3 Pre-QC) =============================================
# ================================================================================== #

cat("Starting 3. PATIENT DEMOGRAPHICS (SK3 Pre-QC) ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

unique_patients_all <- clinical_sk3[match(unique(clinical_sk3$Study.ID), clinical_sk3$Study.ID),]

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
patient_demographics <- rbind(as.numeric(table(clinical_sk3$classification)),patient_demographics)

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

write.csv(patient_demographics, file = file.path(qc.dir, "patient_demographics_sk3_preQC.csv"))
saveRDS(patient_demographics, file = file.path(qc.dir, "patient_demographics_sk3_preQC.rds"))


# ================================================================================== #
# 4. QUALITY CONTROL ===============================================================
# ================================================================================== #

#Load in data#
clinical_sk3 <- readRDS(file.path(processed.data.dir, "clinical_sk3_preQC.rds"))
counts_sk3 <- readRDS(file.path(processed.data.dir, "counts_sk3_preQC.rds"))

# ================================================================================== #
# 4.1 PRINCIPLE COMPONENT ANALYSIS =================================================
# ================================================================================== #
cat("Starting 4.1 PRINCIPLE COMPONENT ANALYSIS",format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

pca.dir <- file.path(qc.dir, "pca")
if(!exists(pca.dir)) dir.create(pca.dir, recursive = TRUE)
if(!exists(file.path(pca.dir,"pca_5"))) dir.create(file.path(pca.dir,"pca_5"), recursive = TRUE)

pca.top500.dir <- file.path(qc.dir, "pca_top500")
if(!exists(pca.top500.dir)) dir.create(pca.top500.dir)
if(!exists(file.path(pca.top500.dir,"pca_5"))) dir.create(file.path(pca.top500.dir,"pca_5"), recursive = TRUE)

row.names(clinical_sk3) == colnames(counts_sk3)

## Voom normalise counts ---------------------------------------------------------
counts_sk3_norm <- voom(counts_sk3) #voom is log2cpm values
# hist(counts_norm$E, xlab = "Counts (log2CPM)")


## Top 500 highly variable genes -----------------------------------------------------
gene_variance <- apply(counts_sk3_norm$E, 1, var) #1 means apply over rows, 2 means columns #variance calculated by taking the differences between each number in the data set and the mean, squaring the differences to make them positive, and then dividing the sum of the squares by the number of values in the data set.
top_genes <- order(gene_variance, decreasing = TRUE)[1:500]
top_counts_sk3_norm <- counts_sk3_norm$E[top_genes, ]
pca_res_top <- prcomp(t(top_counts_sk3_norm), scale. = TRUE, center = TRUE)





## EIGENCORPLOT -----------------------------------------------
cat("Starting eigencorplots", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
clinical_sk3_numeric <- sapply(clinical_sk3[,c("sampletype",
                                               "sex",
                                               "smoking_status",
                                               "classification")],
                               function(col){
                                 as.numeric(as.factor(col))
                               }
)

clinical_sk3_numeric <- cbind(clinical_sk3_numeric, age = clinical_sk3$age)
row.names(clinical_sk3_numeric) <- row.names(clinical_sk3)

## All genes
eigen.pca <- pca(counts_sk3_norm$E,
                 metadata = clinical_sk3_numeric,
                 center = TRUE,
                 scale = TRUE)

png(filename = file.path(pca.dir, paste0("eigencorplot_",i,".png")),
    width = 3000, height = 1500, units = "px", res = 300)

eigencorplot(eigen.pca, metavars = c(colnames(clinical_sk3_numeric)))

dev.off()


# Top 500 variable genes
eigen.pca <- pca(top_counts_sk3_norm,
                 metadata = clinical_sk3_numeric,
                 center = TRUE,
                 scale = TRUE)

png(filename = file.path(pca.top500.dir, paste0("eigencorplot_",i,".png")),
    width = 3000, height = 1500, units = "px", res = 300)

eigencorplot(eigen.pca, metavars = c(colnames(clinical_sk3_numeric)))

dev.off()

### PCA plots ---------------------------------------------------------------------------------------
listofcounts <- list(all = counts_sk3_norm$E, top500 = top_counts_sk3_norm)
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
  
  ## Plot 1: Percentage variance explained ============================================================================================
  cat("Starting Percentage variance explained plots", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  eigenvalues <- pca_res$sdev^2
  var_explained = pca_res$sdev^2 / sum(pca_res$sdev^2)
  
  png(filename = paste0(this.dir,"/var_explained.png"),
      width = 2400, height = 1500, units = "px", res = 300)
  
  plot(var_explained[1:20], type = "b",
       xlab = "Principal Component",
       ylab = "Percentage of Variance Explained",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  
  dev.off()
  
  
  ## Plot 2: Individual PCA plots (1 PC) coloured by clinical_sk3 cariable ============================================================================================
  cat("Starting individual pca plots", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  library("viridis")
  
  for (variable in c("smoking_status", "sex", "classification", "age", "sampletype")){
    ggsave(autoplot(pca_res,data=clinical_sk3, colour = variable)+ fontsize  + guides(color=guide_legend(variable)),  
           file = file.path(this.dir,paste0("pca_",variable,".png")),
           width = 2300,
           height = 1500,
           units = "px")
  }
  

  ## Plot 3: Pairs plot (5 PCs) coloured by clinical_sk3 cariable ============================================================================================
  cat("Starting pairs plot (5 PCs)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  clinical_sk3_pca <- cbind(pca_res$x, clinical_sk3)
  # clinical_sk3_pca <- cbind(eigen.pca$rotated, clinical_sk3) #pca() does the same thing as prcomp()
  
  for (i in c("smoking_status", "sex", "classification", "sampletype")){
    png(filename = paste0(this.dir,"/pca_5/", i, ".png"), width = 26, height = 23, units = "cm", res = 1200)
    
    variable <- as.factor(clinical_sk3_pca[,i])
    
    pcaplot <- pairs(clinical_sk3_pca[,c(1:5)],
                     pch = 19, 
                     cex = 0.75,
                     oma=c(3,3,3,20),
                     col = variable)
    
    par(xpd=TRUE)
    legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = i)
    dev.off()
  }
  
  # Age - need to break age into categories
  png(filename = paste0(this.dir,"/pca_5/Age.png"), width = 26, height = 23, units = "cm", res = 1200 )
  
  variable <- clinical_sk3_pca[,"age"]
  
  colors <- colorRampPalette(c("yellow2", "purple"))(100)
  
  pcaplot <- pairs(clinical_sk3_pca[,c(1:5)],
                   pch = 19, 
                   cex = 0.75,
                   oma=c(3,3,3,20),
                   col =  colors[as.numeric(cut(clinical_sk3_pca$age, breaks = 100))])
  
  breaks <- seq(min(clinical_sk3_pca$age), max(clinical_sk3_pca$age), length.out = 9)
  
  # Add custom color legend
  par(xpd=TRUE)
  
  legend("bottomright", 
         legend = round(breaks, 1), 
         fill = colorRampPalette(c("yellow2", "purple"))(9), 
         title = "Age", 
         border = "white")
  
  dev.off()
  
  
  ggsave(autoplot(pca_res,data=clinical_sk3, colour = "sampletype", label = TRUE)+ fontsize + guides(color=guide_legend("Sample Type")),  
         file =paste0(this.dir,"/pca_sampletype_labelled.png"),
         width = 2200,
         height = 1500,
         units = "px")
  
}

# ================================================================================== #
# 5. QC - Library size, gene counts and sex check ==================================
# ================================================================================== #
cat("Starting 5.  QC - Library size, gene counts and sex check", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

## Library size: colsums/total reads per sample  ------------------------------------------------------------------------------
colsums_persample=as.matrix(colSums(counts_sk3))
write.csv(colsums_persample, file = file.path(qc.dir,"colsums.csv"))
# (lowest is 954,000)

## Number of detected genes ------------------------------------------------------------------------------
genes_per_sample = colSums(counts_sk3 >0) # same as sapply(counts_sk3, function(x) sum(x >0)) 
write.csv(genes_per_sample, file = file.path(qc.dir,"colsums_genes_detec_per_sample.csv"))

## Sex checks ------------------------------------------------------------------------------
cat("Starting sex check", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

sex.dir <- file.path(qc.dir, "sex")
if(!exists(sex.dir)) dir.create(sex.dir, recursive = TRUE)

# male - DDX3Y (ENSG00000067048), XIST (ENSG00000229807)
# female - XIST
counts_cpm <- cpm(counts_sk3)
xist_counts <- counts_cpm["ENSG00000229807",] #XIST
ddx3y_counts <- counts_cpm["ENSG00000067048",] #DDX3Y

ddx3y_samples <- colnames(counts_cpm)[counts_cpm["ENSG00000067048",] != 0]
clinical_sk3[ddx3y_samples, "sex"]


#Raw counts
gendercheck_counts <- as.data.frame(cbind(t(counts_sk3[c("ENSG00000229807","ENSG00000067048"),]), clinical_sk3$sex))
colnames(gendercheck_counts) <- c("XIST", "DDX3Y", "Sex")
gendercheck_counts$Sample <- row.names(gendercheck_counts)

png(filename = file.path(sex.dir,"XIST_vs_DDX3Y_rawcounts.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggplot(data = gendercheck_counts,
       aes(x = as.numeric(DDX3Y),
           y = as.numeric(XIST)
           )) + 
  theme_bw() + 
  geom_point(aes(color = Sex)) +
  # geom_text(aes(label = Sample), angle = 90, vjust = 0.5, hjust = -0.2, size = 0.8) +   # add labels
  ylab(lab = "XIST Counts") +
  xlab(lab = "DDX3Y Counts")+
  ggtitle("Sex Check")

dev.off()


#CPM
gendercheck_cpm <- as.data.frame(cbind(t(counts_cpm[c("ENSG00000229807","ENSG00000067048"),]), clinical_sk3$sex))
colnames(gendercheck_cpm) <- c("XIST", "DDX3Y", "Sex")
gendercheck_cpm$Sample <- row.names(gendercheck_cpm)

png(filename = file.path(sex.dir,"XIST_vs_DDX3Y_cpm.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggplot(data = gendercheck_cpm,
       aes(x = as.numeric(DDX3Y),
           y = as.numeric(XIST))) + 
  theme_bw() + 
  geom_point(aes(color = Sex)) +
  # geom_text(aes(label = Sample), angle = 90, vjust = 0.5, hjust = -0.2, size = 0.8) +   # add labels
  ylab(lab = "XIST CPM") +
  xlab(lab = "DDX3Y CPM")+
  ggtitle("Sex Check")

dev.off()


gendercheck <- merge(gendercheck_counts, gendercheck_cpm, by = "Sample")
gendercheck <- gendercheck[,-which(colnames(gendercheck) == "Sex.x")]
colnames(gendercheck) <- c("Sample", "XIST_count", "DDX3Y_count", "XIST_cpm", "DDX3Y_cpm", "Sex")
gendercheck <- gendercheck[order(gendercheck$DDX3Y_cpm, decreasing = TRUE),]
write.csv(gendercheck, file.path(sex.dir, "gendercheck_table.csv"))

# 107165-001-146 cpm = 3.13 rawcount = 58 DDX37
# 107165-001-184 cpm = 2.82 rawcount = 5 DDX3Y ? will not remove for now as raw count is low, cpm only high because total reads are low for the sample?


#XIST
png(filename = file.path(sex.dir,"XIST_cpm.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggplot(data = gendercheck,
       aes(x = as.factor(Sex),
           y = as.numeric(XIST_cpm),
           fill = Sex)) + 
  theme_bw() + 
  geom_boxplot() +
  ylab(lab = "XIST CPM") +
  xlab(lab = "Sex")+
  ggtitle("XIST")

dev.off()

#DDX3Y
png(filename = file.path(sex.dir,"DDX3Y_cpm.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggplot(data = gendercheck,
       aes(x = as.factor(Sex),
           y = as.numeric(DDX3Y_cpm),
           fill = Sex)) + 
  theme_bw() + 
  geom_boxplot()+
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              stackratio = 0.3,
  #              dotsize = 0.55) +
  ylab(lab = "DDX3Y CPM") +
  xlab(lab = "Sex") +
  ggtitle("DDX3Y")
dev.off()

# ggplot(data = gendercheck, aes(x = as.numeric(XIST), y= as.numeric(DDX3Y), fill = as.factor(Sex))) + geom_point()



# ================================================================================== #
# 5.1 Remove samples/outliers ======================================================
# ================================================================================== #
cat("Starting 5.  QC - Library size, gene counts and sex check", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

clinical_sk3 <- clinical_sk3[-which(row.names(clinical_sk3) == "107165-001-146"),]
counts_sk3 <- counts_sk3[,-which(colnames(counts_sk3) == "107165-001-146"),]


this.processed.data.dir <- file.path(data.dir, "processed", "datawrangling_qc_sk3_only")
if(!exists(this.processed.data.dir))dir.create(this.processed.data.dir)

saveRDS(clinical_sk3, file.path(this.processed.data.dir, "clinical_sk3_postQC.rds")) #225 samples 66138 genes
saveRDS(counts_sk3, file.path(this.processed.data.dir, "counts_sk3_postQC.rds"))

# ================================================================================== #
# 6. Split clinical_sk3 and counts into brush and biopsy ===========================
# ================================================================================== #
cat("Starting 6. Split clinical_sk3 and counts into brush and biopsy", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

clinical_sk3_biopt <- clinical_sk3[which(clinical_sk3$sampletype == "Biopt"),] #116 samples
counts_sk3_biopt <- counts_sk3[, row.names(clinical_sk3_biopt)]

clinical_sk3_brush <- clinical_sk3[which(clinical_sk3$sampletype == "Brush"),] #114 samples
counts_sk3_brush <- counts_sk3[, row.names(clinical_sk3_brush)]


saveRDS(clinical_sk3_biopt, file.path(this.processed.data.dir, "clinical_sk3_biopt.rds"))
saveRDS(counts_sk3_biopt, file.path(this.processed.data.dir,"counts_sk3_biopt.rds"))

saveRDS(clinical_sk3_brush, file.path(this.processed.data.dir,"clinical_sk3_brush.rds"))
saveRDS(counts_sk3_brush, file.path(this.processed.data.dir,"counts_sk3_brush.rds"))



# > setdiff(clinical_sk3_brush$Study.ID, clinical_sk3_biopt$Study.ID) #has brush but not biopt
# [1] "SEO077" "SEO095" "SEO185" "SEO267" "SEO310" "SEO411" "SEO413" "SEO422"
# [9] "SEO423" "SEO424" "SEO537" "SEO553" "SEO556" "SEO563" "SEO574" "SEO585"
# [17] "SEO587" "SEO590"
# > setdiff(clinical_sk3_biopt$Study.ID, clinical_sk3_brush$Study.ID) #has biopt but not brush
# [1] "SEO069" "SEO070" "SEO071" "SEO084" "SEO085" "SEO096" "SEO131" "SEO142"
# [9] "SEO215" "SEO252" "SEO265" "SEO290" "SEO304" "SEO308" "SEO317" "SEO472"
# [17] "SEO580" "SEO592" "SEO603" "SEO610"

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
