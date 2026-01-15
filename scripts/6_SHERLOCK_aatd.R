#SHERLOCK analysis of radiomics varaibles in relation to gene expression
# Transcriptional changes associated with emphysema severity  in COPD patients, also comparing upper and lower lobe
# Linking radiomics feature to transcriptomic signatures in COPD. Emphysema and mucus plugging scores from CT scans


options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED

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


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"


library("readxl")
# library("limma")
# library("rstatix")
# library("tibble")
# library("ggvenn")
# library("ggplot2")
# library("ggrepel")
# library("ggfortify")
library("stringr")
# library("EnsDb.Hsapiens.v79")
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


# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))

#master clinical file - output from 1a_SHERLOCK3_qc_combat.R, after removing outliers
# this file still has unecessary columns we won't use (like survey questions etc. also emphysema values are weird)
# these issues were fixed in Sherlock_database_07_25_Final.xlsx so match to that
clinical_master <- read.csv("/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/clinical_master.csv", row.names = 1)


raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #319 SEO/patient IDs
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


# saveRDS(clinical_sk_all, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds")) 
# saved again after adding z mutation data 


# Get SERPINA1 Z mutation genotypes ------------------------------------------------------------------------------------------------------------
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

row.names(serpina1_z_snp) <- serpina1_z_snp$Study.ID


matching_ids <- intersect(clinical_sk_all$Study.ID, row.names(serpina1_z_snp)) #432 patients 

clinical2 <- clinical_sk_all[clinical_sk_all$Study.ID %in% matching_ids,] #(674 samples from 432 patients)
clinical2 <- cbind(clinical2, serpina1_z_snp_GRCh38_ch14_pos94378610 = serpina1_z_snp[match(clinical2$Study.ID, serpina1_z_snp$Study.ID),
                                                          "z_mutation"])



# standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
# MM: DS < 0.2
# MZ: 0.8 ≤ DS ≤ 1.2
# ZZ: DS ≥ 1.8

#For all unique PATIENTS -------------------------------------------------------------------------------------------------------
table(clinical2[!duplicated(clinical2$Study.ID),"classification"], clinical2[!duplicated(clinical2$Study.ID),"serpina1_z_snp_GRCh38_ch14_pos94378610"])


#                      0   1   2
# Control             94   4   0
# Mild-moderate COPD 124   9   0
# Severe COPD        173  19   9



# For samples ------------------------------------------------------------------
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


sherlock1_ids <- clinical_sk_all[which(clinical_sk_all$batch == 1), "Study.ID"]
sherlock2_ids <- clinical_sk_all[which(clinical_sk_all$batch == 2), "Study.ID"]
sherlock3_ids <- clinical_sk_all[which(clinical_sk_all$batch == 3), "Study.ID"]


# SHERLOCK1 patients that don't have genotyping data
sherlock1_ids[!(sherlock1_ids %in% row.names(serpina1_z_snp))]
# [1] "A_1494" "A_922"  "A_978"  "A_1655" "A_1472" "A_1542" NA       "A_1680"
# [9] "A_1541" "A_2700" "A_2642" "A_3063" "A_2890" "A_960"


# SHERLOCK2 patients that don't have genotyping data
sherlock2_ids[!(sherlock2_ids %in% row.names(serpina1_z_snp))]
# [1] "SEO066" "SEO066" "SEO069" "SEO070" "SEO075" "SEO077" "SEO078" "SEO087"
# [9] "SEO087" "SEO185" "SEO196" "SEO196" "SEO507" "SEO507"

# SHERLOCK3 patients that don't have genotyping data
sherlock3_ids[!(sherlock3_ids %in% row.names(serpina1_z_snp))]
# [1] "SEO069" "SEO070" "SEO077" "SEO185" "SEO418" "SEO418"

clinical123_master <- clinical_sk_all
clinical123_master$serpina1_z_snp_GRCh38_ch14_pos94378610 <- serpina1_z_snp[match(clinical_sk_all$Study.ID, serpina1_z_snp$Study.ID), "z_mutation"]

##### ---------------- START EXTRA WRANGLING ------------------------ ####
#   This was done after writing SHERLOCK_sysmex_diffexp.R script but before SHERLOCK_sysmex_multivariate.R script
### 15/01/2026 SYSMEX IS UNRELATED TO THIS SCRIPT, BUT EDITING SYSMEX VARIABLE COLUMN SO THAT clinical_sherlock123_master.rds IS THE MOST UPDATED FILE THAT CAN BE USED FOR EVERYTHING ! ###
# Some rows for sysmex variable have "----" which turn it into a character. make these NA
# Sysmex first column = "Mono", last column = "EO-Z", 53 total columns
sysmex_index_start <- which(colnames(clinical123_master) == "Mono")
sysmex_index_end <- which(colnames(clinical123_master) == "EO.Z")

sapply(clinical123_master[sysmex_index_start:sysmex_index_end],class)

clinical123_master[which(clinical123_master$RELYMP.103uL == "----"),"RELYMP.103uL"] <- NA
clinical123_master$RELYMP.103uL <- as.numeric(clinical123_master$RELYMP.103uL)

clinical123_master[which(clinical123_master$RELYMP == "----"),"RELYMP"] <- NA
clinical123_master$RELYMP <- as.numeric(clinical123_master$RELYMP)

sapply(clinical123_master[sysmex_index_start:sysmex_index_end],class)

##### ---------------- END EXTRA WRANGLING ------------------------ ####
saveRDS(clinical123_master, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))


# ================================================================================== #
# 2. DIFFERENTIAL EXPRESSION =======================================================
# ================================================================================== #


































#from rachael's script
# test=cbind(getCHROM(vcf),
#            getPOS(vcf),
#            getREF(vcf),
#            getALT(vcf), 
#            getQUAL(vcf), 
#            extract.info(vcf, element = "QD", as.numeric=T), 
#            extract.info(vcf, element = "MQ", as.numeric=T),
#            extract.info(vcf, element = "FS", as.numeric=T),
#            extract.info(vcf, element = "SOR", as.numeric=T), 
#            extract.gt(vcf,element = "GT", as.numeric=FALSE),
#            extract.gt(vcf, element = "AD", as.numeric=FALSE),
#            extract.info(vcf, element = "ANN", as.numeric=FALSE), 
#            extract.info(vcf, element = "GNOMAD_AF", as.numeric=FALSE),
#            
# # Load the VCF
# vcf <- read.vcfR(file.path(data.dir, "raw", "snp_array", "20250910_RES0225_GSAv3+.imputed.vcf"), verbose = FALSE)
#                  
# # Extract genotype matrix
# gt <- extract.gt(vcf, element="GT")  # rows=SNPs, cols=samples
# 
# write.csv(gt, file = file.path(data.dir, "processed", "snp_array","genotypes.csv"), row.names = TRUE)



cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")