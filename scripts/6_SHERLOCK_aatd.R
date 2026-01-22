# SHERLOCK analysis of genotyping data (sepcifically SERPINA1 mutations) in relation to gene expression
# Aim: comparing gene expression for MZ and ZZ patients (in COPD)

# SETTING FOR HPC BATCH JOB (run script from 2.DIFFERENTIAL EXPRESSION as a batch job, everything before that was run in interactive)
hpc_batch_job = TRUE


# ### BASH ####
# # note the .vcf.gz file needs the accompanying .tbi file
# # View header
# module load BCFtools
# bcftools view -h data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz
# 
# #see one variant
# bcftools view -H data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz | head -n 1
# 
# # Subset chr1
# bcftools view -r chr1:1-5e6 data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz -Oz -o data/processed/snp_array/chr1_subset.vcf.gz
# 
# # -Oz means Output format is bgZipped vcf
# 
# Subset chr14
# bcftools view -r chr14 data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz -Oz -o data/processed/snp_array/chr14_array.vcf.gz
# bcftools view -r chr14:94376747-94388602 data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz -Oz -o data/processed/snp_array/chr14_serpina1_array.vcf.gz


#Check size
# bcftools stats data/processed/snp_array/chr14_serpina1_array.vcf.gz | grep "^SN"

# should be hg38 ! according to alen
# Subset samples of interest?
# bcftools view -S samples.txt data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gzz -Oz -o data/processed/snp_array/samples_subset.vcf.gz

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT 
# chr1    
# 10390   chr1:10390_CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA_C       CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA    
# C       
# .       
# LOWIMPQUAL   
# DR2=0.16;AF=0.0009;INFO=0.163865;MAF=0.0009;AF_EUR=0.00085644;AF_ADMIXED=0.001;AF_EAS=0;AF_AFR=0;IMPUTED;RSID=rs1557426845      
# GT:DS:GP    
# 0|0:0:1,0,0 

# Z rs28929474		
# S rs17580		
# I rs28931570		
# M(Herleen) rs199422209	
# Q0amersfoort rs199422210		

if(hpc_batch_job == FALSE){ #don't need this part for hpc


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"


library("readxl")
# library("limma")
# library("rstatix")
# library("tibble")
# library("ggplot2")
# library("ggrepel")
# library("ggfortify")
library("stringr")
# library("ggpubr")
# library("edgeR")
# library("DESeq2")
# library("tidyverse")
# library("PCAtools")
library("vcfR")


# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir,"data")

processed.data.dir <- file.path(data.dir,"processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir,"output")
output.dir <- file.path(output.dir, "aatd")
if(!exists(output.dir))dir.create(output.dir)

# ================================================================= #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))

#master clinical file - output from 1a_SHERLOCK3_qc_combat.R, after removing outliers
# this file still has unecessary columns we won't use (like survey questions etc. also emphysema values are weird)
# these issues were fixed in Sherlock_database_07_25_Final.xlsx so match to that
clinical_master <- read.csv("/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/clinical_master.csv", row.names = 1)


raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #598 SEO/patient IDs
# Note that this file has less columns than previous version, survey columns have been removed
# Also, this file has one row per patient, clinical_sk_all has one row per sample (same patients can have more than one sample bc brush and biopt)
# IMPORTANT - radiomics data in this file is updated too. the clinical_brushbiopt_master file has weird emphysema % values 


clinical_sk_all <- as.data.frame(raw_clinical[match(clinical_master$Study.ID, raw_clinical$class_incl_study_id),])
# clinical_sk_all <- clinical_sk_all[match(clinical_master$Study.ID, clinical_sk_all$Study.ID),]

clinical_sk_all <- cbind(clinical_sk_all, sampletype = clinical_master$sampletype, 
                             batch = clinical_master$batch,
                             classification = clinical_master$classification)
row.names(clinical_sk_all) <- row.names(clinical_master)

clinical_sk_all <- clinical_sk_all %>% 
  dplyr::rename(
    Study.ID = class_incl_study_id,
    age = crf_age,
    sex = crf_gender,
    smoking_status = crf_smoking,
    packyears = crf_packyears,
    corticosteroid = crf_corticosteroid,
    FVC_post= postbodybox_fvc_post ,
    FEV1 = postbodybox_fev1_post,
    FEV1_FVC_post = postbodybox_fev1_fvc_post
  ) 


clinical_sk_all[,c("age", "packyears", "FEV1", 
                       "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(
                         clinical_sk_all[,c("age", "packyears", "FEV1", 
                                                "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], 
                         function(x) as.numeric(x))


# LIB5426634_SAM24375606 
#this patient had no data in raw_clinical but is in clinical_sk1_master.rds(generated from mastertable_sherlock3.csv ie. the SHERLOCK1 data alen gave me in 2024 
#(mastertable_sherlock3.csv does not mean sherlock3 study, just the name of alen's sherlock1 file version 3))
clinical_sk1_master <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_master.rds"))
clinical_sk1_master["LIB5426634_SAM24375606",]
clinical_sk1_master <- clinical_sk1_master %>% 
  dplyr::rename(
    sex = crf_gender,
    smoking_status = crf_smoking,
    packyears = crf_packyears,
    corticosteroid = crf_corticosteroid,
    FVC_post= postbodybox_fvc_post ,
    FEV1_FVC_post = postbodybox_fev1_fvc_post
  ) 
clinical_sk1_master$age <- as.numeric(clinical_sk1_master$age)
common_cols <- intersect(colnames(clinical_sk1_master), colnames(clinical_sk_all))
clinical_sk_all["LIB5426634_SAM24375606", common_cols] <- clinical_sk1_master["LIB5426634_SAM24375606", common_cols]


## Get SERPINA1 Z mutation genotypes ------------------------------------------------------------------------------------------------------------
# Pull out chromosome 14
vcf_serpina1 <- read.vcfR(
  "data/processed/snp_array/chr14_serpina1_array.vcf.gz",
  verbose = FALSE
)

head(vcf_serpina1, n= 10)

fix <- as.data.frame(getFIX(vcf_serpina1))
serpina1_snps_df <- cbind(fix[, c("CHROM","POS", "ID","REF","ALT", "QUAL", "FILTER")], 
                          extract.gt(vcf_serpina1, element = "DS", as.numeric = TRUE)
                          # extract.gt(vcf_serpina1, element = "GT", as.numeric = FALSE)
                          ) #626 samples
 
# Current naming "RES0225_SEO568_GSAv3+". Convert to SEO and A number alone ("SEO568")
colnames(serpina1_snps_df)[8:ncol(serpina1_snps_df)] <- str_extract(colnames(serpina1_snps_df)[8:ncol(serpina1_snps_df)], "(?<=_)[^_]+(?=_)")

# Pull out Z mutation (ch14 pos94378610 for hg38)
# This should be rs28929474 (SERPINA1 Pi*Z mutation) but RS IDs aren't labelled
serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610),] #Z mutation
# chr14 94378610   ref=C   alt=T 

#Remove the first 7 column (genotyping meta) so we are left with just dosage (ie. genotype as decimal)
serpina1_z_snp <- as.data.frame(cbind(Study.ID = colnames(serpina1_snps_df)[-c(1:7)], 
                        z_mutation = as.numeric(serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610), -c(1:7)])))
# write.csv(serpina1_z_snp, file.path(processed.data.dir, "snp_array","serpina1_z_snp_dosage.csv"), row.names = FALSE)

# standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
# MM: DS < 0.2
# MZ: 0.8 ≤ DS ≤ 1.2
# ZZ: DS ≥ 1.8

#Convert names in clinical file to new names provided by Daan (the A numbers in genotyping file are different to clinical file)
A_number_conversion_file <- as.data.frame(read_excel(file.path(data.dir, "raw", "A_number_conversion_file.xlsx")))
A_number_conversion_file$genotyping_ID <- paste0("A", A_number_conversion_file$'A-numbers Genotyping')
A_number_conversion_file$clinical_ID <- paste0("A", A_number_conversion_file$'A-numbers Sherlock database')



#Add mutation to clinical file

A_IDs <- serpina1_z_snp[grep("^A", serpina1_z_snp$Study.ID),]
not_A_IDs <- serpina1_z_snp[!serpina1_z_snp$Study.ID %in% A_IDs$Study.ID,]


A_IDs$Study.ID %in% A_number_conversion_file$genotyping_ID
A_IDs$Study.ID.corrected <- A_number_conversion_file[match(A_IDs$Study.ID, A_number_conversion_file$genotyping_ID), "clinical_ID"]


A_IDs$Study.ID.corrected <- sub("^A", "A_", A_IDs$Study.ID.corrected )


serpina1_z_snp <- rbind(cbind(Study.ID = A_IDs$Study.ID.corrected, 
                                   z_mutation = A_IDs$z_mutation),
                             not_A_IDs)

row.names(serpina1_z_snp) <- serpina1_z_snp$Study.ID #There are 626 patients with genotyping data
length(unique(serpina1_z_snp$Study.ID))

matching_ids <- intersect(clinical_sk_all$Study.ID, row.names(serpina1_z_snp)) #432 patients that have genotyping data AND gene expresion data. 194 have genotyping but not gene expression

clinical2 <- clinical_sk_all[clinical_sk_all$Study.ID %in% matching_ids,] #(674 samples from 432 patients)
clinical2 <- cbind(clinical2, serpina1_z_snp_GRCh38_ch14_pos94378610 = serpina1_z_snp[match(clinical2$Study.ID, serpina1_z_snp$Study.ID),
                                                          "z_mutation"])



# standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
# MM: DS < 0.2
# MZ: 0.8 ≤ DS ≤ 1.2
# ZZ: DS ≥ 1.8

#For all unique PATIENTS -------------------------------------------------------------------------------------------------------#
table(clinical2[!duplicated(clinical2$Study.ID),"classification"], clinical2[!duplicated(clinical2$Study.ID),"serpina1_z_snp_GRCh38_ch14_pos94378610"])


#                      0   1   2
# Control             94   4   0
# Mild-moderate COPD 124   9   0
# Severe COPD        173  19   9



# For samples ------------------------------------------------------------------#
clinical_brush <-  clinical2[which(clinical2$sampletype == "Brush"),] #143 (sherlock1, 2 and 3)
table(clinical_brush$classification, clinical_brush$serpina1_z_snp_GRCh38_ch14_pos94378610)


#                      0   1   2
# Control             88   3   0
# Mild-moderate COPD 116   9   0
# Severe COPD        171  17   9



clinical_biopt <-  clinical2[which(clinical2$sampletype == "Biopt"),] #185 (sherlock 2 and 3)
table(clinical_biopt$classification, clinical_biopt$serpina1_z_snp_GRCh38_ch14_pos94378610)


#                     0  1  2
# Control            69  2  0
# Mild-moderate COPD 98  6  0
# Severe COPD        75  6  5



# Make master table!
clinical123_master <- clinical_sk_all
clinical123_master$serpina1_z_snp_GRCh38_ch14_pos94378610 <- serpina1_z_snp[match(clinical_sk_all$Study.ID, serpina1_z_snp$Study.ID), "z_mutation"]



### START EXTRA WRANGLING ------------------------------------------------------------------------------------------------------
# 15/01/26: Came back to add this section after running 5_SHERLOCK_radiomics and SHERLOCK_sysmex_diffexp.R scripts,  before starting SHERLOCK_sysmex_multivariate.R script
# The below wrangling is unrelated to this script. editing radiomics and sysmex variables columns so that clinical_sherlock123_master.rds IS THE MOST UPDATED FILE THAT CAN BE USED FOR EVERYTHING ! ###


## FOR RADIOMICS COLUMNS
#Make the names valid for R
radiomics_index_start <-which(colnames(clinical123_master) =="RL_insp_vol_ml")
radiomics_index_end <- which(colnames(clinical123_master) =="LLL_airtrapping_emphysema_%")

colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] <- gsub("%", "perc", colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] )
colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] <- gsub(">", "over", colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] )
colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] <- gsub("-", ".", colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] )
clinical123_master[, radiomics_index_start:radiomics_index_end] <- lapply(
  clinical123_master[, radiomics_index_start:radiomics_index_end],
  function(x) {
    if (is.list(x)) {
      as.numeric(unlist(x))
    } else {
      as.numeric(x)
    }
  }
)
sapply(clinical123_master[radiomics_index_start:radiomics_index_end],class)


# FOR SYSMEX COLUMNS
# Some rows for sysmex variable have "----" which turn it into a character. make these NA
# Sysmex first column = "Mono", last column = "EO-Z", 53 total columns
sysmex_index_start <- which(colnames(clinical123_master) == "Mono")
sysmex_index_end <- which(colnames(clinical123_master) == "EO-Z")

colnames(clinical123_master)[sysmex_index_start:sysmex_index_end] <- make.names(colnames(clinical123_master)[sysmex_index_start:sysmex_index_end])

sapply(clinical123_master[sysmex_index_start:sysmex_index_end],class)

clinical123_master[which(clinical123_master$RELYMP.103uL == "----"),"RELYMP.103uL"] <- NA
clinical123_master$RELYMP.103uL <- as.numeric(clinical123_master$RELYMP.103uL)

clinical123_master[which(clinical123_master$RELYMP == "----"),"RELYMP"] <- NA
clinical123_master$RELYMP <- as.numeric(clinical123_master$RELYMP)

sapply(clinical123_master[sysmex_index_start:sysmex_index_end],class)

### END EXTRA WRANGLING ---------------------------------------------------------------------------

saveRDS(clinical123_master, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
write.csv(clinical123_master, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.csv"))

### SUMMARISING NUMBERS ----------------------------------------------------------------------------------


sherlock1_ids <- unique(clinical123_master[which(clinical123_master$batch == 1), "Study.ID"])
sherlock2_ids <- unique(clinical123_master[which(clinical123_master$batch == 2), "Study.ID"])
sherlock3_ids <- unique(clinical123_master[which(clinical123_master$batch == 3), "Study.ID"])


### 626 Patients have genotyping data
### 598 patients have clinical data in master file
### 541 patients have expression data
### 24 PATIENTS THAT HAVE EXPRESSION DATA AND CLINICAL DATA BUT NO GENOTYPING data  ------------------------------------------------- #
setdiff(unique(clinical123_master$Study.ID), row.names(serpina1_z_snp))
cat(setdiff(unique(clinical123_master$Study.ID), row.names(serpina1_z_snp)), sep = "\n")

# [1] "A_1494" "A_922"  "A_978"  "A_1655" "A_1472" "A_1542" "A_1680" "A_1541"
# [9] "A_2700" "A_2642" "A_3063" "A_2890" "A_960"  "SEO066" "SEO069" "SEO070"
# [17] "SEO075" "SEO077" "SEO078" "SEO087" "SEO185" "SEO196" "SEO507" "SEO418"

# Breakdown------------------------------------------#
# SHERLOCK1
sherlock1_ids[!(sherlock1_ids %in% row.names(serpina1_z_snp))]
# [1] "A_1494" "A_922"  "A_978"  "A_1655" "A_1472" "A_1542" "A_1680" "A_1541"
# [9] "A_2700" "A_2642" "A_3063" "A_2890" "A_960"


# SHERLOCK2 
sherlock2_ids[!(sherlock2_ids %in% row.names(serpina1_z_snp))]
# [1] "SEO066" "SEO069" "SEO070" "SEO075" "SEO077" "SEO078" "SEO087" "SEO185"
# [9] "SEO196" "SEO507"


# SHERLOCK3 
sherlock3_ids[!(sherlock3_ids %in% row.names(serpina1_z_snp))]
# [1] "SEO069" "SEO070" "SEO077" "SEO185" "SEO418"

### 432 PATIENTS THAT HAVE GENOTYPING DATA BUT NO EXPRESSION data ### -------------------------------------------------#
cat(intersect(clinical_sk_all$Study.ID, row.names(serpina1_z_snp)), sep = "\n")

### 194 PATIENTS THAT HAVE GENOTYPING DATA BUT NO EXPRESSION data ### -------------------------------------------------#
setdiff(row.names(serpina1_z_snp), clinical123_master$Study.ID)
cat(setdiff(row.names(serpina1_z_snp), clinical123_master$Study.ID), sep = "\n")


### 179 PATIENTS THAT HAVE GENOTYPING DATA BUT NOT CLINICAL DATA  ### -------------------------------------------------#
setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id))
cat(setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id)), sep = "\n")


### 74 PATIENTS THAT HAVE GENOTYPING DATA AND CLINICAL DATA BUT NO EXPRESSION DATA ### ---------------------------------#
setdiff(setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id)), matching_ids)
cat(setdiff(setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id)), matching_ids), sep = "\n")


# ##-- Post batch correction
counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
sherlock1_counts<- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))

unique(clinical123_master[intersect(colnames(counts_brush), row.names(clinical123_master)), "Study.ID"]) #270 patients have brushes (sherlock 2,3)
unique(clinical123_master[intersect(colnames(counts_biopt), row.names(clinical123_master)), "Study.ID"]) #271 patients have biopsies (sherlock 2,3)
unique(clinical123_master[intersect(colnames(sherlock1_counts), row.names(clinical123_master)), "Study.ID"]) #167 patients have brushes (sherlock 1)


intersect(row.names(serpina1_z_snp), raw_clinical$class_incl_study_id)
### --------------------------------------------------------------------------------------------------------------------------------------#

} #close the hpc_batch_job. script from here on will be run 

# ================================================================================== #
# 2. DIFFERENTIAL EXPRESSION =======================================================
# ================================================================================== #

options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED

# ================================================================================== #
## 2A. SCRIPT SET UP ===============================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("readxl")
library("ggplot2")
library("DESeq2")
library("ggrepel")
library("ggfortify")
library("stringr")
library("tidyverse")

# ================================================================================== #
## 2B. SET UP DIRECTORY & OUTPUT PATHS =============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir,"data")

processed.data.dir <- file.path(data.dir,"processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir, "output","aatd")

setwd(file.path(main.dir))

# ================================================================================== #
## 2.1. LOAD IN DATA ===============================================================
# ================================================================================== #

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
clinical123_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
colnames(clinical123_master)[which(colnames(clinical123_master) == "serpina1_z_snp_GRCh38_ch14_pos94378610")] <- "serpina1_snp"

# Brush SHERLOCK1, 2 and 3 counts
counts123_brush <- read.csv(file.path(combat.processed.data.dir, "sherlock1_2_3_counts_brush_integrated.csv"), check.names = FALSE, row.names = 1) #from SHERLOCK_SOP_integration.R script (copied from  "/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/sherlock1_2_3_combat)
clinical_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush"),] 

# Biopsy SHERLOCK 2 and 3 counts
counts23_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
clinical_biopt <- clinical123_master[which(clinical123_master$sampletype == "Biopt"),]
counts23_biopt <- counts23_biopt[,row.names(clinical_biopt)]

#General COPD vs Severe COPD
# MM = 0
# MZ = 1
# ZZ = 2

table(clinical_brush$classification, clinical_brush$serpina1_snp)
#                      0   1   2
# Control             88   3   0
# Mild-moderate COPD 116   9   0
# Severe COPD        171  17   9

table(clinical_biopt$classification, clinical_biopt$serpina1_snp)
#                     0  1  2
# Control            69  2  0
# Mild-moderate COPD 98  6  0
# Severe COPD        75  6  5

#make valid names

# ================================================================================== #
# 3. DIFFEERENTIAL EXPRESSION (DESeq2) =============================================
# General COPD (Subset only for Mild/Mod and Severe COPD. Compare the genotypes)
# Severe COPD (Subset for only Severe COPD. Compare the genotypes)
# ================================================================================== #
#DESEQ DIRECTORY
diffexp.dir <- file.path(output.dir, "diffexp_serpina1_deseq")
if(!exists(diffexp.dir))dir.create(diffexp.dir)


diffexp_deseq_func <- function(sampletype, disease_group){
  cat(paste("Starting 3. DIFFERENTIAL EXPRESSION (DESeq)", sampletype), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  if(disease_group == "general_copd"){
  this.diffexp.dir <- file.path(diffexp.dir, "general_copd")}
  
  if(disease_group == "severe_copd"){
  this.diffexp.dir <- file.path(diffexp.dir, "severe_copd")}
  
  if(!exists(this.diffexp.dir))dir.create(this.diffexp.dir)
  
  #Results and figures directory
  diffexp.results.dir <- file.path(this.diffexp.dir, "results")
  if(!exists(diffexp.results.dir))dir.create(diffexp.results.dir)
  
  diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir))dir.create(diffexp.figures.dir)
  
  
  if(sampletype == "brush"){
      clinical <- clinical_brush
      counts <- counts123_brush
  }
  
  if(sampletype == "biopt"){
    clinical <- clinical_biopt
    counts <- counts23_biopt
  }

  
  if(all(colnames(counts) == row.names(clinical)) == FALSE){ 
    stop("all(colnames(counts) == row.names(clinical) = FALSE)") }
  
  #Remove samples with NA serpina1 snp data
  clinical <- clinical[-which(is.na(clinical$serpina1_snp)),]
  counts <- counts[,row.names(clinical)]
  
  #make names valid
  clinical$smoking_status<- make.names(clinical$smoking_status)
  
if(disease_group == "general_copd"){
#Only include COPD samples
clinical <- clinical[-which(clinical$classification == "Control"),]
counts <- counts[,row.names(clinical)]

# DGE
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = clinical,
                              design = ~ 0 + serpina1_snp + age + sex + smoking_status)
}



if(disease_group == "severe_copd"){
  #Only include COPD samples
  clinical <- clinical[which(clinical$classification == "Severe COPD"),]
  counts <- counts[,row.names(clinical)]
  
  # DGE
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = clinical,
                                design = ~ 0 + serpina1_snp + age + sex + packyears) #no smokers in severe group
  
}



# Filter the genes that are lowly expressed and normalize
# Low exp genes affects statistics - can make p value very significant even if only a few samples are lowly expressed amongst other samples with no expression.  also affects multiple testing.
# One method = keep row medians that are greater than 10 (ie. half of the samples for a gene must have a minimum number of 10 counts)
keep <- rowMedians(counts(dds)) >= 10

dds <- dds[keep,]
dds <- DESeq(dds)

# results extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values;
resultsNames(dds)

listofcontrasts <- list(
  MZ_MM = c("serpina1_snp", "1", "0"),#This means we have set snp1 (MZ) as control and we define fold change based on snp1 as baseline (log fold change = MZ - MMM)
  ZZ_MM = c("serpina1_snp", "2", "0"), 
  ZZ_MZ = c("serpina1_snp", "2", "1")

)  

listoftT <- list()
listoftT2 <- list()
listofvolcano <- list()

for (contrast in names(listofcontrasts)){
  cat(paste(sampletype, contrast), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
results <- results(dds, contrast = c(listofcontrasts[[contrast]][1], #variable name
                                     listofcontrasts[[contrast]][2], #reference 
                                     listofcontrasts[[contrast]][3])) #other comparison

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

listoftT[[contrast]] <- tT 
listoftT2[[contrast]] <- tT2

write.csv(tT2, file = file.path(diffexp.results.dir, paste0(sampletype, "_", contrast, "_tT2.csv")))


# ================================================================================== #
# 3.1. VOLCANO PLOT ================================================================
# ================================================================================== #
cat(paste("Starting 3.1. VOLCANO PLOT", contrast), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

volcano <- ggplot(tT, aes(x = log2FoldChange, y = -log10(pvalue))) +
  ggtitle(paste0(sampletype,": ",contrast)) +
  geom_point(aes(color = Legend)) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
  geom_hline(yintercept =-log10(max(tT2$pvalue)),colour="black", linetype="dashed")+
  geom_text_repel(data = subset(tT2[1:30,]),
                  aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines") ) +
  theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                   legend.text = element_text(size = 14),
                                   legend.title = element_text(size = 16)) 

listofvolcano[[contrast]] <- volcano

ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(sampletype, "_", contrast, "_111volcano.png")),
       width = 25, height = 25,
       units = "cm")

} #close listofcontrasts loop

listofresults[[disease_group]][[sampletype]] <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)

} #end diffexp deseq function



listofresults <- list(
  general_copd = list(),
  severe_copd  = list()
)

diffexp_deseq_func(disease_group = "general_copd", sampletype = "brush")
diffexp_deseq_func(disease_group = "general_copd", sampletype = "biopt")

diffexp_deseq_func(disease_group = "severe_copd", sampletype = "brush")
diffexp_deseq_func(disease_group = "severe_copd", sampletype = "biopt")



#Save all results
saveRDS(listofresults[["general_copd"]], file = file.path(diffexp.dir, "general_copd", "results", "listofresults.rds"))
saveRDS(listofresults[["severe_copd"]], file = file.path(diffexp.dir, "severe_copd", "results","listofresults.rds"))





# ================================================================================== #
# 4. DIFFEERENTIAL EXPRESSION (edgeR) ================================================
# General COPD (Subset only for Mild/Mod and Severe COPD. Compare the genotypes)
# Severe COPD (Subset for only Severe COPD. Compare the genotypes)
# ================================================================================== #
library(edgeR)
diffexp.dir <- file.path(output.dir, "diffexp_serpina1_edgeR")
if(!exists(diffexp.dir)) dir.create(diffexp.dir)

diffexp_edgeR_func <- function(sampletype, disease_group){
  cat(paste("Starting 4. DIFFERENTIAL EXPRESSION (edgeR)", sampletype), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  if(disease_group == "general_copd"){
    this.diffexp.dir <- file.path(diffexp.dir, "general_copd")}
  
  if(disease_group == "severe_copd"){
    this.diffexp.dir <- file.path(diffexp.dir, "severe_copd")}
  
  if(!exists(this.diffexp.dir))dir.create(this.diffexp.dir)
  
  #Results and figures directory
  diffexp.results.dir <- file.path(this.diffexp.dir, "results")
  if(!exists(diffexp.results.dir))dir.create(diffexp.results.dir)
  
  diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir))dir.create(diffexp.figures.dir)
  
  
  if(sampletype == "brush"){
    clinical <- clinical_brush
    counts <- counts123_brush
  }
  
  if(sampletype == "biopt"){
    clinical <- clinical_biopt
    counts <- counts23_biopt
  }
  
  
  if(all(colnames(counts) == row.names(clinical)) == FALSE){ 
    stop("all(colnames(counts) == row.names(clinical) = FALSE)") }
  
  #Remove samples with NA serpina1 snp data
  clinical <- clinical[-which(is.na(clinical$serpina1_snp)),]
  counts <- counts[,row.names(clinical)]
  
  
  #make names valid
  clinical$smoking_status<- make.names(clinical$smoking_status)
  
  if(disease_group == "general_copd"){
    #Only include COPD samples
    clinical <- clinical[-which(clinical$classification == "Control"),]
    counts <- counts[,row.names(clinical)]
    
    design <- model.matrix(~ 0 + serpina1_snp + age + sex + smoking_status,
                           data = clinical) 
  }
  
  
  
  if(disease_group == "severe_copd"){
    #Only include COPD samples
    clinical <- clinical[which(clinical$classification == "Severe COPD"),]
    counts <- counts[,row.names(clinical)]
    
    design <- model.matrix(~ 0 + serpina1_snp + age + sex + packyears,
                           data = clinical) 
    
    
  }
  
  # colnames(design)[1:3] <- c(levels(as.factor(clinical$serpina1_snp))) #cant do this as names will be invalid for making contrasts ("0","1","2")
  
  DGEL<- DGEList(counts=counts, group = clinical$classification) 
  
  #FILTER
  keep <- filterByExpr(DGEL) 
  # keep <- which(rowMedians(as.matrix(DGEL))>10) 
  # keep <- rowSums(cpm(expression)>100) >= 2
  
  DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
  
  # NORMALISE
  DGEL<- calcNormFactors(DGEL,method = "TMM")
  
  # ESTIMATE DISPERSON
  DGEL <- estimateDisp(DGEL, design)
  
  # FIT MODEL
  fit <- glmQLFit(DGEL, design)  #fit the GLM (design) to the DGEL(the DGEL object, which contains the counts data that has been filtered,normalised and dispersons estimtated)
  
  
  my.contrasts <- makeContrasts(
    MZ_MM = serpina1_snp1 - serpina1_snp0,
    ZZ_MM = serpina1_snp2 - serpina1_snp0,
    ZZ_MZ = serpina1_snp2 - serpina1_snp1,
    levels = design) 
  
  
  listoftT <- list()
  listoftT2 <- list()
  listofvolcano <- list()  
  
  for (contrast in colnames(my.contrasts)){
    cat(paste(sampletype, contrast), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    qlf <- glmQLFTest(fit, contrast=my.contrasts[,contrast]) #after fitting GLM to counts data (glmQLFit), glmQLFTest perform quasi likelihood f-tests to test for differential expression. ie. hypothesis testing for differential expression. (make inferences or draw conclusions about the data)
    tT <- topTags(qlf,n=nrow(DGEL))$table #topTags gets top genes, here we want all of the genes, edgeR's default p.adjust method is BH
    tT$Legend <- ifelse(
      tT$FDR < 0.05 & tT$logFC > 0, "Upregulated", # try >1
      ifelse(
        tT$FDR < 0.05 & tT$logFC < 0, "Downregulated",# try <-1
        "Not Significant"))
    
    tT$Legend[is.na(tT$Legend)]="Not significant"
    
    tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols
    
    #for those with no hgnc symbol, label with ensembl id
    tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl

    
    selection <-which(tT$FDR<0.05)
    # selection <-which((tT$logFC>1|tT$logFC< -1)&tT$FDR<0.05)
    
    tT2 <- tT[selection,]
    
    listoftT[[contrast]] <- tT 
    listoftT2[[contrast]] <- tT2
    
    write.csv(tT2, file = file.path(diffexp.results.dir, paste0(sampletype, "_", contrast, "_tT2.csv")))
    
    
# ================================================================================== #
# 4.1. VOLCANO PLOT ================================================================
# ================================================================================== #
    cat(paste("Starting 4.1. VOLCANO PLOT", contrast), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
      ggtitle(paste0(sampletype,": ",contrast)) +
      geom_point(aes(color = Legend)) +
      scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
      geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
      geom_text_repel(data = subset(tT2[1:30,]),
                      aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines") ) +
      theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                       legend.text = element_text(size = 14),
                                       legend.title = element_text(size = 16)) 
    
    listofvolcano[[contrast]] <- volcano
    
    ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(sampletype, "_", contrast, "_volcano.png")),
           width = 25, height = 25,
           units = "cm")
    
  } #close listofcontrasts loop
  
  listofresults[[disease_group]][[sampletype]] <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)
  
} #end diffexp deseq function



listofresults <- list(
  general_copd = list(),
  severe_copd  = list()
)

diffexp_edgeR_func(disease_group = "general_copd", sampletype = "brush")
diffexp_edgeR_func(disease_group = "general_copd", sampletype = "biopt")

diffexp_edgeR_func(disease_group = "severe_copd", sampletype = "brush")
diffexp_edgeR_func(disease_group = "severe_copd", sampletype = "biopt")



#Save all results
saveRDS(listofresults[["general_copd"]], file = file.path(diffexp.dir, "general_copd", "results", "listofresults.rds"))
saveRDS(listofresults[["severe_copd"]], file = file.path(diffexp.dir, "severe_copd", "results","listofresults.rds"))


cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
