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


saveRDS(clinical_sk_all, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds")) 

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))

# Get SERPINA1 Z mutation genotypes ------------------------------------------------------------------------------------------------------------
# Pull out chromosome 14
vcf_serpina1 <- read.vcfR(
  "data/processed/snp_array/chr14_serpina1_array.vcf.gz",
  verbose = FALSE
)

head(vcf_serpina1, n= 10)

fix <- as.data.frame(getFIX(vcf_serpina1))
serpina1_snps_df <- cbind(fix[, c("CHROM","POS","REF","ALT", "QUAL", "FILTER")], 
                          extract.gt(vcf_serpina1, element = "DS", as.numeric = TRUE)
                          # extract.gt(vcf_serpina1, element = "GT", as.numeric = FALSE)
                          ) #626 samples
 
# Current naming "RES0225_SEO568_GSAv3+". Convert to SEO and A number alone ("SEO568")
colnames(serpina1_snps_df)[7:ncol(serpina1_snps_df)] <- str_extract(colnames(serpina1_snps_df)[7:ncol(serpina1_snps_df)], "(?<=_)[^_]+(?=_)")

# Pull out Z mutation (ch14 pos94378610 for hg38)
serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610),] #Z mutation
# chr14 94378610   ref=C   alt=T 

#Remove the first 6 column (genotyping meta) so we are left with just dosage (ie. genotype as decimal)
serpina1_z_snp <- as.data.frame(cbind(Study.ID = colnames(serpina1_snps_df)[-c(1:6)], 
                        z_mutation = as.numeric(serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610), -c(1:6)])))
write.csv(serpina1_z_snp, file.path(processed.data.dir, "snp_array","serpina1_z_snp_dosage.csv"), row.names = FALSE)

# standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
# MM: DS < 0.2
# MZ: 0.8 ≤ DS ≤ 1.2
# ZZ: DS ≥ 1.8


#Add mutation to clinical file
serpina1_z_snp$Study.ID <- sub("^A", "A_", serpina1_z_snp$Study.ID)
row.names(serpina1_z_snp) <- serpina1_z_snp$Study.ID


matching_ids <- intersect(clinical_sk_all$Study.ID, row.names(serpina1_z_snp)) #328 patients (570 samples)

clinical2 <- clinical_sk_all[clinical_sk_all$Study.ID %in% matching_ids,] #(570 samples from 328 patients)
clinical2 <- cbind(clinical2, z_mutation = serpina1_z_snp[match(clinical2$Study.ID, serpina1_z_snp$Study.ID),
                                                          "z_mutation"])

#For all unique PATIENTS -------------------------------------------------------------------------------------------------------
table(clinical2[!duplicated(clinical2$Study.ID),"classification"], clinical2[!duplicated(clinical2$Study.ID),"z_mutation"])

#                      0   1   2
# Control             94   4   0
# Mild-moderate COPD 124   9   0
# Severe COPD         86   6   5


# For samples ------------------------------------------------------------------
clinical_brush <-  clinical2[which(clinical2$sampletype == "Brush"),] #143
table(clinical_brush$classification, clinical_brush$z_mutation)

#                      0   1   2
# Control             88   3   0
# Mild-moderate COPD 116   9   0
# Severe COPD         84   4   5



clinical_biopt <-  clinical2[which(clinical2$sampletype == "Biopt"),] #185
table(clinical_biopt$classification, clinical_biopt$z_mutation)

#                     0  1  2
# Control            69  2  0
# Mild-moderate COPD 98  6  0
# Severe COPD        75  6  5

sherlock1_ids <- clinical_sk_all[which(clinical_sk_all$batch == 1), "Study.ID"]
sherlock2_ids <- clinical_sk_all[which(clinical_sk_all$batch == 2), "Study.ID"]
sherlock3_ids <- clinical_sk_all[which(clinical_sk_all$batch == 3), "Study.ID"]

colnames(serpina1_snps_df)[-c(1:6)]

# SHERLOCK1 patients that don't have genotyping data
sherlock1_ids[!(sherlock1_ids %in% colnames(serpina1_snps_df)[-c(1:6)])]
# [1] "A_1524" "A_1446" "A_1494" "A_1024" "A_1515" "A_952"  "A_753"  "A_662"
# [9] "A_494"  "A_1189" "A_1073" "A_922"  "A_1580" "A_1537" "A_1615" "A_860"
# [17] "A_1388" "A_1402" "A_1695" "A_1694" "A_1470" "A_1584" "A_1512" "A_143"
# [25] "A_1681" "A_1209" "A_305"  "A_1735" "A_1576" "A_1790" "A_1057" "A_978"
# [33] "A_1780" "A_1477" "A_1804" "A_1772" "A_1778" "A_1122" "A_1145" "A_108"
# [41] "A_1655" "A_1002" "A_1879" "A_1119" "A_1769" "A_1896" "A_1520" "A_1657"
# [49] "A_1405" "A_53"   "A_1496" "A_1941" "A_1727" "A_1919" "A_2022" "A_1472"
# [57] "A_817"  "A_1848" "A_2028" "A_1688" "A_1911" "A_2215" "A_1542" "A_1985"
# [65] "A_1407" "A_2026" "A_727"  "A_2301" "A_1827" NA       "A_884"  "A_2304"
# [73] "A_2274" "A_1680" "A_2271" "A_2509" "A_2372" "A_1721" "A_2132" "A_2377"
# [81] "A_2321" "A_2397" "A_2135" "A_1926" "A_2453" "A_1912" "A_2315" "A_2479"
# [89] "A_2029" "A_719"  "A_2395" "A_1667" "A_1541" "A_582"  "A_1891" "A_2361"
# [97] "A_2393" "A_1954" "A_506"  "A_2448" "A_1233" "A_757"  "A_2466" "A_2557"
# [105] "A_2667" "A_2423" "A_2700" "A_2617" "A_2719" "A_2642" "A_2736" "A_3063"
# [113] "A_2709" "A_2890" "A_2625" "A_960"  "A_1664" "A_660"  "A_1840"

# SHERLOCK2 patients that don't have genotyping data
sherlock2_ids[!(sherlock2_ids %in% colnames(serpina1_snps_df)[-c(1:6)])]
# [1] "SEO066" "SEO066" "SEO069" "SEO070" "SEO075" "SEO077" "SEO078" "SEO087"
# [9] "SEO087" "SEO185" "SEO196" "SEO196" "SEO507" "SEO507"

# SHERLOCK3 patients that don't have genotyping data
sherlock3_ids[!(sherlock3_ids %in% colnames(serpina1_snps_df)[-c(1:6)])]
# [1] "SEO069" "SEO070" "SEO077" "SEO185" "SEO418" "SEO418"

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