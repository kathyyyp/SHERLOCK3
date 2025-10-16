# SHERLOCK - Medina - Integration of of SHERLOCK1&2&3 brush data, batch corrected with Combat_seq 
options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("sva")
library("readxl")
library("stringr")

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory #my SHERLOCK3 directory

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
# 1. CREATE SHERLOCK1 FILE FROM UDPATED MASTER CLINICAL FILE =================================
# ================================================================================== #
cat("Starting 1. CREATE SHERLOCK1 FILE FROM UDPATED MASTER CLINICAL FILE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# MASTER SHERLOCK CLINICAL FILE ===============================================================
clinical_sherlock_full <- read.csv(file.path(data.dir, "processed","clinical_sherlock_full_master_preQC.csv"), row.names = 1)
row.names(clinical_sherlock_full) <- clinical_sherlock_full$Study.ID
# NOTE: Combat needs counts in matrix class, not data frame

#Load in SHERLOCK3 data #66,138 genes
counts_sk3<- readRDS(file.path(processed.data.dir, "datawrangling_qc_sk3_only", "counts_sk3_postQC.rds")) #225 samples, 66138 genes

#Load in SHERLOCK2 data #69,972 genes #Note that 4 outliers were removed because they had low total counts, otherwise all counts available
counts_sk2 <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts.rds")) #327 samples, 69972 genes

#Load in SHERLOCK1 data #60,683 genes #Note, 4 IDs with count that didn't have clinical data. Rectified below
rawcounts_sk1 <- read.csv(file.path("data","raw","SHERLOCK1","rawcounts.txt"), sep="\t")
row.names(rawcounts_sk1)=rawcounts_sk1[,1]
rawcounts_sk1=rawcounts_sk1[,-1] # 173 samples

s_raw <- read.csv(file.path("data","raw","SHERLOCK1","mastertable_sherlock3.csv"))
sherlock1_clinical_raw <- s_raw[-which(is.na(s_raw$rna_seq.sample.id)),]
row.names(sherlock1_clinical_raw) <- sherlock1_clinical_raw$rna_seq.sample.id
# These patients didn't have counts data in original SHERLOCK1 dataset (no rnaseq sample id)
# SEO057, A1533 and A2629 ? Work out where this is !!!!!!

c(setdiff(colnames(rawcounts_sk1),sherlock1_clinical_raw$rna_seq.sample.id), setdiff(sherlock1_clinical_raw$rna_seq.sample.id, colnames(rawcounts_sk1)))
# [1] "LIB5426567_SAM24375539" "LIB5426582_SAM24375554" "LIB5426583_SAM24375555"
# [4] "LIB5426587_SAM24375559"


# Matched IDs from Alen
# LIB5426567_SAM24375539 580A A_580
# LIB5426582_SAM24375554 994A A_994
# LIB5426583_SAM24375555 984A A_984
# LIB5426587_SAM24375559 1688A A_1688

# As shown above, there were 4 patients with counts but not clinical info for sherlock1 original data 
# check this now with daan's master sherlock file

# Combine the 4 samples that had the missing clinical data originally (now available in master clinical file) with the other samples
sherlock1_ID_conversion <- data.frame(
  rnaseq_id = c("LIB5426567_SAM24375539", #580A 
                "LIB5426582_SAM24375554", #994A 
                "LIB5426583_SAM24375555", #984A
                "LIB5426587_SAM24375559"),#1688A
  patient_id = c("A_580",   #no classification ( patient not in original file )
               "A_994",     #no classification ( patient not in original file )
               "A_984",     #no classification ( patient not in original file )
               "A_1688")    
)
  
sherlock1_ID_conversion <- rbind(sherlock1_ID_conversion, 
                                 data.frame(rnaseq_id = sherlock1_clinical_raw$rna_seq.sample.id,
                                            patient_id = sherlock1_clinical_raw$sample.id))

sherlock1_ID_conversion$patient_id <- gsub("A(?=\\d)", "A_", sherlock1_ID_conversion$patient_id, perl = TRUE)


#subset for sherlock1
sherlock1_ID_conversion$patient_id %in% row.names(clinical_sherlock_full) #A_1864 was in og sherlock1 file but not the master file
clinical_sk1_master <- clinical_sherlock_full[sherlock1_ID_conversion$patient_id,]
clinical_sk1_master <- cbind(rnaseq_id = sherlock1_ID_conversion$rnaseq_id, clinical_sk1_master)


clinical_sk1_master <- cbind(clinical_sk1_master, sampletype = "brush")

row.names(clinical_sk1_master) <- clinical_sk1_master$rnaseq_id
  
row.names(clinical_sk1_master) %in% colnames(rawcounts_sk1)
colnames(rawcounts_sk1) %in% row.names(clinical_sk1_master)

clinical_sk1_master <- clinical_sk1_master[colnames(rawcounts_sk1),]

colnames(rawcounts_sk1) == row.names(clinical_sk1_master)

#Some NAs for GOLD classification
clinical_sk1_master[is.na(clinical_sk1_master$classification), "rnaseq_id"]
# [1] "LIB5426567_SAM24375539" "LIB5426582_SAM24375554" "LIB5426583_SAM24375555"
# [4] "LIB5426616_SAM24375588"

#Their sample IDs
clinical_sk1_master[is.na(clinical_sk1_master$classification), "Study.ID"]
# "A_580" 
# "A_994" 
# "A_984" 
# NA #the NA is A_1864, available in old sherlock1 clinical but not the master

na_classifications <- row.names(clinical_sk1_master)[is.na(clinical_sk1_master$classification)]
sherlock1_clinical_raw[na_classifications,]
clinical_sk1_master[na_classifications,]
sherlock1_clinical_raw[na_classifications, c("sample.id", "GOLD.stage", "FEV1.over.FVC", "FEV1.post.bronchodilator.percent.predicted")]

#In summary, only A_1864 has fev1, fev/fvc and classification info in sherlock1_clinical_raw, the other three do not have fev1 or fev/fvc data to work out classification
#Why is A_1864 missing from the master clinical file??


# Remove the 4 patients above with NA classification -----------------------------------------------------------------------
clinical_sk1_master <- clinical_sk1_master[-c(which(is.na(clinical_sk1_master$classification))),] #169 samples remaining (from 173)
counts_sk1 <- rawcounts_sk1[,row.names(clinical_sk1_master)] #60683 genes, 169 samples
clinical_sk1 <- clinical_sk1_master #incase I want to subset columns (currently have 600+ clinical_sk3 variables)

sherlock1.processed.dir <-file.path(processed.data.dir, "SHERLOCK1")
if(!exists(sherlock1.processed.dir)) dir.create(sherlock1.processed.dir)
saveRDS(counts_sk1, file.path(sherlock1.processed.dir,"counts_sk1.rds"))
write.csv(counts_sk1, file.path(sherlock1.processed.dir,"counts_sk1.csv"))

saveRDS(clinical_sk1_master, file.path(sherlock1.processed.dir,"clinical_sk1_master.rds"))
write.csv(clinical_sk1_master, file.path(sherlock1.processed.dir,"clinical_sk1_master.csv"))




# ================================================================================== #
# 1. COMBAT-SEQ TO COMBINE SHERLOCK1, SHERLOCK2 AND SHERLOCK3 =================================
# ================================================================================== #
cat("Starting 1. COMBAT-SEQ TO COMBINE SHERLOCK2 AND SHERLOCK3", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#Load in SHERLOCK3 data #66,138 genes
# counts_sk3 <- readRDS(file.path(processed.data.dir, "datawrangling_qc_sk3_only", "counts_sk3_postQC.rds")) #225 samples, 66138 genes
counts_sk3_brush <- readRDS(file.path(data.dir, "processed", "datawrangling_qc_sk3_only","counts_sk3_brush.rds")) #112 brushes

#Load in SHERLOCK2 data #69,972 genes
# counts_sk2 <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts.rds")) #327 samples, 69972 genes
counts_sk2_brush <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts_brush.rds")) #162 brushes

#Load in SHERLOCK1 data #60,683 genes #all are brushes
counts_sk1 <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds")) #169 brush samples, 60683 genes,


# ComBat ================================================================================
# Merge and keep only matching Ensembl IDs between sk1 and sk2 two files
counts_sk1_sk2_all <- merge(counts_sk1, counts_sk2_brush, by = "row.names", all = FALSE) # all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
row.names(counts_sk1_sk2_all) <- counts_sk1_sk2_all[,1]
counts_sk1_sk2_all <- as.matrix(counts_sk1_sk2_all[,-1])

# Merge the sk1_sk2 file with sk3 (keep all, including non matched Ensembl IDS - just for saving)
counts_all <- merge(counts_sk1_sk2_all,  counts_sk3_brush,  by = "row.names", all = TRUE)

# Merge the sk1_sk2 with sk3 (keep only the matched Ensembl IDS)
counts_merged <- merge(counts_sk1_sk2_all,  counts_sk3_brush,  by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1]) #48,039 genes and 443 samples

medina.data.dir <- file.path(processed.data.dir, "datawrangling_qc_allbrush_medina")
if(!exists(medina.data.dir))dir.create(medina.data.dir, recursive = TRUE)

combat.processed.data.dir <- file.path(medina.data.dir, "combat_results")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir, recursive = TRUE)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat.csv")) #48,039 genes and 721 samples
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat.rds"))

#combat
batch <- c(rep(1, ncol(counts_sk1)), rep(2, ncol(counts_sk2_brush)), rep(3, ncol(counts_sk3_brush)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat.rds"))


#Note; specifying group = sampletype in combat_seq would mean that cmbat esitmates batch effects within each sampletype and remove the within-group batch effects 
# but this would miss correcting for batch effects between sampletypes 



# ================================================================================== #
# 3. Create merged clinical files (sk1, sk2 and sk3) ===============================
# ================================================================================== #
#Combine sherlock1, sherlock2 and sherlock3 clinical 

clinical_sk1_master <- readRDS(file.path(sherlock1.processed.dir,"clinical_sk1_master.rds")) #169 samples
#remove rnaseq.id column and add batch column so colnames match sk2_sk3 clinical file
clinical_sk1_master <- clinical_sk1_master[,-1] 
clinical_sk1_master$batch <- 1
clinical_sk1_master$sampletype <- "Brush"

clinical_sk2_sk3_brush_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_brush_master.rds")) #274 brush samples (postQC on sherlock3 and 2)
clinical_sherlock_brush_master <- as.data.frame(rbind(clinical_sk1_master,clinical_sk2_sk3_brush_master)) #443 samples

#Check that patients are in same order in clinical file and combat counts file
row.names(clinical_sherlock_brush_master) == colnames(counts_combat)

# Save files

#Clinical file containing sherlock1+2+3 samples (brush only)
if(!exists(file.path(file.path(medina.data.dir, "master")))) dir.create(file.path(medina.data.dir, "master"))

saveRDS(clinical_sherlock_brush_master, file.path(medina.data.dir,  "master","clinical_sherlock_brush_master.rds"))
write.csv(clinical_sherlock_brush_master, file.path(medina.data.dir,  "master","clinical_sherlock_brush_master.csv"))



# ================================================================================== #
# 3.1 Subset master clinical file for main variables ===============================
# ================================================================================== #
cat("Starting 3.1. Subset main clinical file", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

clinical_sherlock_brush <- clinical_sherlock_brush_master %>% 
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


clinical_sherlock_brush$classification <- make.names(clinical_sherlock_brush$classification)
clinical_sherlock_brush$smoking_status <- make.names(clinical_sherlock_brush$smoking_status)
clinical_sherlock_brush[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(clinical_sherlock_brush[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], function(x) as.numeric(x))


## Save clinical_sherlock_brush (main patient info, smoking, ics and lung function data)
write.csv(clinical_sherlock_brush, file.path(medina.data.dir, "clinical_sherlock_brush_simple.csv"))
saveRDS(clinical_sherlock_brush, file.path(medina.data.dir, "clinical_sherlock_brush_simple.rds")) #443 samples





























# Files that were moved into sherlock3_studies for medina

# raw
raw_expression <- read.table(file.path("raw","20250902_Sherlock3_readcount_UMI_dedup.txt"), header = TRUE, check.names = FALSE, row.names = 1)
raw_clinical_sk3 <- read_xlsx(file.path("raw","25_08_20_SHERLOCk_1_2 mvdb.xlsx"), .name_repair = "universal") 
clinical_sk3_ids <- read_xlsx(file.path("raw", "Copy of Sample Submission Form _107165 (1).xlsx"), sheet = "sample_submission_info", col_names = TRUE, .name_repair = "universal") # unique SEO/patient IDs


# processed
write.csv(raw_clinical_sherlock_full, file.path(data.dir, "processed","clinical_sherlock_full_master_preQC.csv"))

## preQC
write.csv(counts_sk3, file.path(processed.data.dir, "counts_sk3_preQC.csv"))
write.csv(clinical_sk3_master, file.path(processed.data.dir, "clinical_sk3_master_preQC.csv"))

# postQC sk3 (removed 107165-001-146)
clinical_sk3_master <- 

saveRDS(clinical_sk3_biopt, file.path(this.processed.data.dir, "clinical_sk3_biopt.rds"))
saveRDS(counts_sk3_biopt, file.path(this.processed.data.dir,"counts_sk3_biopt.rds"))

saveRDS(clinical_sk3_brush, file.path(this.processed.data.dir,"clinical_sk3_brush.rds"))
saveRDS(counts_sk3_brush, file.path(this.processed.data.dir,"counts_sk3_brush.rds"))
  





file.path(this.processed.data.dir, "clinical_sk3_biopt.rds")
file.path(this.processed.data.dir,"counts_sk3_biopt.rds")
file.path(this.processed.data.dir,"clinical_sk3_brush.rds")
file.path(this.processed.data.dir,"counts_sk3_brush.rds")


