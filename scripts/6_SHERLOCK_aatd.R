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

library(vcfR)
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


#master - all 600+ clinical variables
raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #319 SEO/patient IDs

# INCLUDES SHERLOCK1
clinical_sk_all <- raw_clinical %>% 
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
                         clinical_brushbiopt[,c("age", "packyears", "FEV1", 
                                                "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], 
                         function(x) as.numeric(x))




# INCLUDES SHERLOCK2 and 3 only
clinical_brushbiopt_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_brushbiopt_master_named.rds"))

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))

# Get SERPINA1 Z mutation genotypes
vcf_serpina1 <- read.vcfR(
  "data/processed/snp_array/chr14_serpina1_array.vcf.gz",
  verbose = FALSE
)

head(vcf_serpina1, n= 10)


fix <- as.data.frame(getFIX(vcf_serpina1))

serpina1_snps_df <- cbind(fix[, c("CHROM","POS","REF","ALT", "QUAL", "FILTER")], 
                          extract.gt(vcf_serpina1, element = "DS", as.numeric = TRUE)
                          # extract.gt(vcf_serpina1, element = "GT", as.numeric = FALSE)
                          )
#626 samples

colnames(serpina1_snps_df)[7:ncol(serpina1_snps_df)] <- str_extract(colnames(serpina1_snps_df)[7:ncol(serpina1_snps_df)], "(?<=_)[^_]+(?=_)")

serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610),] #Z mutation
# chr14 94378610   ref=C   alt=T 

serpina1_z_snp <- as.data.frame(cbind(Study.ID = colnames(serpina1_snps_df)[-c(1:6)], 
                        z_mutation = as.numeric(serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610), -c(1:6)])))
write.csv(serpina1_z_snp, file.path(processed.data.dir, "snp_array","serpina1_z_snp_dosage.csv"), row.names = FALSE)

# standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
# MM: DS < 0.2
# MZ: 0.8 ≤ DS ≤ 1.2
# ZZ: DS ≥ 1.8

#Add mutation to clinical file
row.names(serpina1_z_snp) <- serpina1_z_snp$Study.ID
clinical_brushbiopt_master$Study.ID %in% row.names(serpina1_z_snp)
row.names(serpina1_z_snp) %in% clinical_brushbiopt_master$Study.ID
row.names(serpina1_z_snp) %in% clinical_sk_all$Study.ID

match(clinical_brushbiopt_master$Study.ID, row.names(serpina1_z_snp))
intersect(clinical_brushbiopt_master$Study.ID, row.names(serpina1_z_snp))

clinical_sk_all$Study.ID <- gsub( "A_", "A", clinical_sk_all$Study.ID )
clinical_brushbiopt_master<- cbind(clinical_brushbiopt_master, 
      serpina1_z_mutation = serpina1_z_snp[match(clinical_brushbiopt_master$Study.ID, row.names(serpina1_z_snp)),])


clinical_brush <-  clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Brush"),] #270
clinical_biopt <-  clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Biopt"),] #271


global_lung_vars <- c("RL_insp_LAA_.950HU_perc",
                      "LL_insp_LAA_.950HU_perc",
                      "RL_insp_emphysema_15percentile_HU",
                      "LL_insp_emphysema_15percentile_HU",
                      "Lungs_Pi10", 
                      "Lungs_insp_bronchial_count", 
                      "Lungs_mucus_plugs_count", 
                      "Lungs_emphysema_perc", 
                      "Lungs_airtrapping_perc")



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