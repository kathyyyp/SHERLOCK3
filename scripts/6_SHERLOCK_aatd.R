#SHERLOCK analysis of radiomics varaibles in relation to gene expression
# Transcriptional changes associated with emphysema severity  in COPD patients, also comparing upper and lower lobe
# Linking radiomics feature to transcriptomic signatures in COPD. Emphysema and mucus plugging scores from CT scans


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
output.dir <- file.path(output.dir, "radiomics")
if(!exists(output.dir))dir.create(output.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))

# 
# #includes sysmex data
# # Note that this file has less columns than previous version, survey columns have been removed
# raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #319 SEO/patient IDs
# 
# #master - all 600+ clinical variables
# clinical_brushbiopt_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_brushbiopt_master.rds"))
# 
# 
# hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
# 
# 
# setwd(file.path(main.dir))
# 
# 
# clinical_brush <-  clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Brush"),] #270
# clinical_biopt <-  clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Biopt"),] #271


# Load the VCF
vcf <- read.vcfR(file.path(data.dir, "raw", "snp_array", "20250910_RES0225_GSAv3+.imputed.vcf"))
                 
# Extract genotype matrix
gt <- extract.gt(vcf, element="GT")  # rows=SNPs, cols=samples

write.csv(gt, file = file.path(data.dir, "processed", "snp_array","genotypes.csv", row.names = TRUE))

# 
# # Example: convert to numeric (0,1,2)
# gt_numeric <- apply(gt, 2, function(x) {
#   sapply(x, function(g) {
#     if (g == "0/0") return(0)
#     if (g == "0/1" | g == "1/0") return(1)
#     if (g == "1/1") return(2)
#     return(NA)
#   })
# })

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")