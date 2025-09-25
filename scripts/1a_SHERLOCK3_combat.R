# SHERLOCK3 - integration of new data with SHERLOCK2
options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

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
if(!exists(postQC.data.dir)) dir.create(postQC.data.dir)

#Output directory
output.dir <- file.path(main.dir,"output")
if(!exists(output.dir)) dir.create(output.dir, recursive = TRUE)

#SHERLOCK2 dir
sherlock2.dir <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK2"

setwd(file.path(main.dir))

# ================================================================================== #
# 6. COMBAT-SEQ TO COMBINE SHERLOCK2 AND SHERLOCK3 =================================
# ================================================================================== #
cat("Starting 6. COMBAT-SEQ TO COMBINE SHERLOCK2 AND SHERLOCK3", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# NOTE: Combat needs counts in matrix class, not data frame

#Load in SHERLOCK3 data #66,138 genes
counts_sk3_brush <- readRDS(file.path(processed.data.dir, "datawrangling_qc_sk3_only", "counts_sk3_brush.rds")) #112 samples, 66138 genes
counts_sk3_biopt <- readRDS(file.path(processed.data.dir, "datawrangling_qc_sk3_only", "counts_sk3_biopt.rds")) #113 samples, 66138 genes


#Load in SHERLOCK2 data #69,972 genes
counts_sk2_brush <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts_brush.rds")) #162 samples, 69972 genes
counts_sk2_biopt <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts_biopt.rds")) #165 samples, 69972 genes

# ComBat_seq assumes main source of variation is technical. split by brush and biopsy first as combat may mistake the biological variation as batch effect = overcorrectio and loss of signal
# Combat_seq preserves raw count format (takes raw counts as input)
# Uses a negative binomial regression to model batch effects, then provide adjusted data by mapping the original data to an expected distribution if there were no batch effects.

# Brush ComBat ================================================================================
# Merge and keep matching Ensembl IDs between the two files 
counts_brush_all <- merge(counts_sk2_brush, counts_sk3_brush, by = "row.names", all = TRUE) #all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
counts_brush_merged <- merge(counts_sk2_brush, counts_sk3_brush, by = "row.names") # 274 samples,  64489 genes
row.names(counts_brush_merged) <- counts_brush_merged[,1]
counts_brush_merged <- as.matrix(counts_brush_merged[,-1])

combat.processed.data.dir <- file.path(data.dir, "processed", "datawrangling_qc", "combat_results")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir)

write.csv(counts_brush_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat.csv"))
saveRDS(counts_brush_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat.rds"))

#combat
batch_brush <- c(rep(2, ncol(counts_sk2_brush)), rep(3, ncol(counts_sk3_brush)))
counts_brush_combat <- ComBat_seq(counts_brush_merged, batch = batch_brush) #genes

write.csv(counts_brush_combat, file.path(combat.processed.data.dir, "counts_brush_combat_raw.csv"))
saveRDS(counts_brush_combat, file.path(combat.processed.data.dir, "counts_brush_combat_raw.rds"))

#Biopt Combat ================================================================================
# Merge and keep matching Ensembl IDs between the two files
counts_biopt_all <- merge(counts_sk2_biopt, counts_sk3_biopt, by = "row.names", all = TRUE) #all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
counts_biopt_merged <- merge(counts_sk2_biopt, counts_sk3_biopt, by = "row.names") # 278 samples, 64489 genes
row.names(counts_biopt_merged) <- counts_biopt_merged[,1]
counts_biopt_merged <- as.matrix(counts_biopt_merged[,-1])

write.csv(counts_biopt_merged, file.path(combat.processed.data.dir, "counts_biopt_merged_pre_combat.csv"))
saveRDS(counts_biopt_merged, file.path(combat.processed.data.dir, "counts_biopt_merged_pre_combat.rds"))

#combat
batch_biopt <- c(rep(2, ncol(counts_sk2_biopt)), rep(3, ncol(counts_sk3_biopt)))
counts_biopt_combat <- ComBat_seq(counts_biopt_merged, batch = batch_biopt) #genes

write.csv(counts_biopt_combat, file.path(combat.processed.data.dir, "counts_biopt_combat_raw.csv"))
saveRDS(counts_biopt_combat, file.path(combat.processed.data.dir, "counts_biopt_combat_raw.rds"))








# ================================================================================== #
# 7. LOAD IN DATA AGAIN==================================================================
# ================================================================================== #
cat("Starting 7. Load in data again", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

## COUNTS ---------------------------------------------------------------------------------
## Merged SHERLOCK 2 and 3 counts ##
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

##-- Pre batch correction
counts_brush_merged <- readRDS(file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat.rds"))
counts_biopt_merged <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_merged_pre_combat.rds"))

##-- Post batch correction
counts_brush_combat <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat_raw.rds"))
counts_biopt_combat <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat_raw.rds"))

## CLINICAL -------------------------------------------------------------------------------------
## MASTER CLINICAL TABLE ##
clinical_sherlock_master_preQC <- read.csv(file.path(data.dir, "processed","clinical_sherlock_full_master_preQC.csv"), row.names = 1)

# SHERLOCK2 matched GS_IDs and SEO_IDs
clinical_sk2_ids <- read_xlsx(file.path(sherlock2.dir, "data", "raw", "Samples_clinical_nameconvert.xlsx"), col_names = TRUE, .name_repair = "universal") #172 unique SEO/patient IDs



# ================================================================================== #
# 8. Create merged clinical files for brush and biopt ==============================
# ================================================================================== #
## Already have clinical file for sk3 ##

## Make sk2 clinical files

cat("Starting 8. Create merged clinical files for brush and biopt", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Include only SHERLOCK2 patients
clinical_sk2_ids <- clinical_sk2_ids[order(clinical_sk2_ids$Customer.ID),]
clinical_sk2_ids <- data.frame(clinical_sk2_ids[-which(is.na(clinical_sk2_ids$'Customer.ID')),])
clinical_sk2_ids$sampletype <- word(clinical_sk2_ids$Remarks, 1)
table(clinical_sk2_ids$sampletype)

row.names(clinical_sherlock_master_preQC) <- clinical_sherlock_master_preQC$Study.ID

# Pivot the clinical_sk2 table long to include all samples
clinical_sk2_master <- clinical_sherlock_master_preQC[clinical_sk2_ids$Customer.ID,]
row.names(clinical_sk2_master) <- clinical_sk2_ids$GS_ID
clinical_sk2_master[,1:3] 
clinical_sk2_master <- cbind(clinical_sk2_master, sampletype = clinical_sk2_ids$sampletype)

# Note: some further patients will be removed because we removed outliers for sherlock2. when we match to the combat results they'll be removed
clinical_sk3_master <- read.csv(file.path(processed.data.dir, "clinical_sk3_master_preQC.csv"), row.names = 1)

#Combine sherlock2 and sherlock3 clinical 
clinical_sherlock_master <- as.data.frame(rbind(clinical_sk2_master,clinical_sk3_master))
clinical_sherlock_master <- cbind(clinical_sherlock_master,
                                  batch = c(rep(2, nrow(clinical_sk2_master)),
                                            rep(3, nrow(clinical_sk3_master))
                                            ))

#Subset for the patients in the combat counts file for brush and biopt
clinical_brush_master <- clinical_sherlock_master[colnames(counts_brush_merged),] #274 patients (same as counts)
clinical_biopt_master <- clinical_sherlock_master[colnames(counts_biopt_merged),] #278 patients (same as counts)

#Combine brushes and biopsies (sherlock 2+3 master clinical file brush+biopt)
clinical_brushbiopt_master <- as.data.frame(rbind(clinical_brush_master, clinical_biopt_master))
counts_brushbiopt_master <- as.data.frame(cbind(counts_brush_combat, counts_biopt_combat))

# Save files

#Clinical file containing sherlock2+3 samples (brush + biopt)
write.csv(counts_brushbiopt_master, file.path(postQC.data.dir, "master","counts_brushbiopt_master.csv"))
saveRDS(counts_brushbiopt_master, file.path(postQC.data.dir,  "master","counts_brushbiopt_master.rds"))

write.csv(clinical_brushbiopt_master, file.path(postQC.data.dir, "master","clinical_brushbiopt_master.csv"))
saveRDS(clinical_brushbiopt_master, file.path(postQC.data.dir, "master", "clinical_brushbiopt_master.rds"))

#brush and biopt separated
write.csv(clinical_brush_master, file.path(postQC.data.dir, "master", "clinical_brush_master.csv"))
saveRDS(clinical_brush_master, file.path(postQC.data.dir, "master", "clinical_brush_master.rds"))

write.csv(clinical_biopt_master, file.path(postQC.data.dir,"master", "clinical_biopt_master.csv"))
saveRDS(clinical_biopt_master, file.path(postQC.data.dir, "master", "clinical_biopt_master.rds"))



# ================================================================================== #
# 8.1 Subset main clinical file ================================================
# ================================================================================== #
cat("Starting 8.1. Subset main clinical file", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

clinical_brushbiopt <- clinical_brushbiopt_master %>% 
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
    sampletype,
    batch
    
  ) %>% 
  dplyr::rename(
    sex = crf_gender,
    smoking_status = crf_smoking,
    packyears = crf_packyears,
    corticosteroid = crf_corticosteroid,
    FVC_post= postbodybox_fvc_post ,
    FEV1_FVC_post = postbodybox_fev1_fvc_post
  )


clinical_brushbiopt$classification <- make.names(clinical_brushbiopt$classification)
clinical_brushbiopt$smoking_status <- make.names(clinical_brushbiopt$smoking_status)
clinical_brushbiopt[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(clinical_brushbiopt[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], function(x) as.numeric(x))

clinical_brush <- clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),]
clinical_biopt <- clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),]

## Save clinical_brushbiopt (main patient info, smoking, ics and lung function data)
write.csv(clinical_brushbiopt, file.path(postQC.data.dir, "clinical_brushbiopt.csv"))
saveRDS(clinical_brushbiopt, file.path(postQC.data.dir, "clinical_brushbiopt.rds")) #552 samples
saveRDS(clinical_brush, file.path(postQC.data.dir, "clinical_brush.rds")) #274 samples
saveRDS(clinical_biopt, file.path(postQC.data.dir, "clinical_biopt.rds")) #278 samples

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")