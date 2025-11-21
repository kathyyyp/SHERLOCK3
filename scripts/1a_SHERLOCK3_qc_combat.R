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
library("readxl")
library("stringr")
library("tidyverse")
library("edgeR")

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

#QC directory
qc.dir <- file.path(output.dir, "qc")
if(!exists(qc.dir)) dir.create(qc.dir, recursive = TRUE)

#SHERLOCK2 dir
sherlock2.dir <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK2"

setwd(file.path(main.dir))

#Saving for other UMCG users
common.dir <- "/groups/umcg-griac/tmp02/projects/SHERLOCK_2025"
if(!exists(common.dir)) dir.create(common.dir)
# ================================================================================== #
# 1. COMBAT-SEQ TO COMBINE SHERLOCK2 AND SHERLOCK3 =================================
# ================================================================================== #
cat("Starting 1. COMBAT-SEQ TO COMBINE SHERLOCK2 AND SHERLOCK3", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# NOTE: Combat needs counts in matrix class, not data frame

#Load in SHERLOCK3 data #66,138 genes
counts_sk3<- readRDS(file.path(processed.data.dir, "datawrangling_qc_sk3_only", "counts_sk3_postQC.rds")) #225 samples, 66138 genes

#Load in SHERLOCK2 data #69,972 genes
counts_sk2 <- readRDS(file.path(sherlock2.dir, "data", "processed", "datawrangling_qc", "counts.rds")) #327 samples, 69972 genes

## MASTER CLINICAL TABLE 
clinical_sherlock_master_preQC <- read.csv(file.path(data.dir, "processed","clinical_sherlock_full_master_preQC.csv"), row.names = 1)

# SHERLOCK2 matched GS_IDs and SEO_IDs
clinical_sk2_ids <- read_xlsx(file.path(sherlock2.dir, "data", "raw", "Samples_clinical_nameconvert.xlsx"), col_names = TRUE, .name_repair = "universal") #172 unique SEO/patient IDs


# ================================================================================== #
# 2. MERGE SHERLOCK2 and 3 ==========================================================
# ================================================================================== #
# 1) Merge and keep matching Ensembl IDs between the two files ----------------------------------------------

counts_all <- merge(counts_sk2, counts_sk3, by = "row.names", all = TRUE) # all = TRUE will append all non matched rows in X and Y aswell,showing NA for the unmatched
counts_merged <- merge(counts_sk2, counts_sk3, by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1])




# ================================================================================== #
# 3a. Create merged SK2 adnd SK3 clinical files for brush and biopt ==============================
# ================================================================================== #
## Already have clinical file for sk3 ##

## Make sk2 clinical files

cat("Starting 3. Create merged clinical files for brush and biopt", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

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
clinical_sk3_master <- readRDS(file.path(processed.data.dir, "clinical_sk3_master_preQC.rds"))

#Combine sherlock2 and sherlock3 clinical 
clinical_sherlock_master <- as.data.frame(rbind(clinical_sk2_master,clinical_sk3_master))
clinical_sherlock_master <- cbind(clinical_sherlock_master,
                                  batch = c(rep(2, nrow(clinical_sk2_master)),
                                            rep(3, nrow(clinical_sk3_master))
                                  ))

#Subset for the patients in the merged file
clinical_brushbiopt_master <- clinical_sherlock_master[colnames(counts_merged),]


# ================================================================================== #
# 1. CREATE SHERLOCK1 FILE FROM UDPATED MASTER CLINICAL FILE =======================
# ================================================================================== #
cat("Starting 1. CREATE SHERLOCK1 FILE FROM UDPATED MASTER CLINICAL FILE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

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
sherlock1_ID_conversion$patient_id %in% row.names(clinical_sherlock_master_preQC) #A_1864 was in og sherlock1 file but not the master file
clinical_sk1_master <- clinical_sherlock_master_preQC[sherlock1_ID_conversion$patient_id,]
clinical_sk1_master <- cbind(rnaseq_id = sherlock1_ID_conversion$rnaseq_id, clinical_sk1_master)


clinical_sk1_master <- cbind(clinical_sk1_master, sampletype = "Brush")

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
# NA #the NA is A_1864, available in old sherlock1 clinical but not the master, did not include

na_classifications <- row.names(clinical_sk1_master)[is.na(clinical_sk1_master$classification)]
sherlock1_clinical_raw[na_classifications,]
clinical_sk1_master[na_classifications,]
sherlock1_clinical_raw[na_classifications, c("sample.id", "GOLD.stage", "FEV1.over.FVC", "FEV1.post.bronchodilator.percent.predicted")]

#In summary, only A_1864 has fev1, fev/fvc and classification info in sherlock1_clinical_raw, the other three do not have fev1 or fev/fvc data to work out classification


# Remove the 4 patients above with NA classification -----------------------------------------------------------------------
clinical_sk1_master <- clinical_sk1_master[-c(which(is.na(clinical_sk1_master$classification))),] #169 samples remaining (from 173)
counts_sk1 <- rawcounts_sk1[,row.names(clinical_sk1_master)] #60683 genes, 169 samples

clinical_sk1_master <- cbind(clinical_sk1_master, batch = rep(1, nrow(clinical_sk1_master)))

# ================================================================================== #
# 4. SAMPLE EXCLUSIONS =============================================================
# ================================================================================== #
# These patients were excluded from study (Daan) and no classification available anyway - not in the clinical files here
# "A_580"
# "A_994"
# "A_984"


# A2804 (SHERLOCK1) and SEO230 (SHERLOCK2) are the same patient - keep SEO230 in most cases
clinical_sk1_master <- clinical_sk1_master[-which(clinical_sk1_master$Study.ID == "A_2804"),] #168 samples
counts_sk1 <- counts_sk1[,row.names(clinical_sk1_master)] #60683 counts and 168 samples

sherlock1.processed.dir <-file.path(processed.data.dir, "SHERLOCK1")
if(!exists(sherlock1.processed.dir)) dir.create(sherlock1.processed.dir)
saveRDS(counts_sk1, file.path(sherlock1.processed.dir,"counts_sk1.rds"))
write.csv(counts_sk1, file.path(sherlock1.processed.dir,"counts_sk1.csv"))

saveRDS(clinical_sk1_master, file.path(sherlock1.processed.dir,"clinical_sk1_master.rds"))
write.csv(clinical_sk1_master, file.path(sherlock1.processed.dir,"clinical_sk1_master.csv"))


# ================================================================================== #
# 4a. Remove 107165-001-016 and 107165-001-059 switched samples =====================
# ================================================================================== #
# 107165-001-016 and 107165-001-059 were derived from the same person (identical genotypes) but annotated as SEO267 and SEO265 respectively. Sample switch.
# - need to check when genotyping data is provided (still not provided as of 13/11/2025 )
# 107165-001-016
# 107165-001-059


clinical_brushbiopt_master <- clinical_brushbiopt_master[-which(row.names(clinical_brushbiopt_master) %in% c("107165-001-016","107165-001-059")),] #552 samples

# ================================================================================== #
# 4b. Remove patients with repeat samples - Checl Library size, gene counts ========
# ================================================================================== #
# 13/11/2025 Patient demographics table showed more counts than patients (created patient demographisc table in the next script 1b_SHERLOCK3_qc) and noticed the repeat samples
# 13/11/2025 Came back to this script to remove them. Previous revisions are on Github kathyyyp/SHERLOCK3

# sherlock1 (brush) --------------------------------------------


dup_rows <- clinical_sk1_master[clinical_sk1_master$Study.ID %in% 
                                  clinical_sk1_master$Study.ID[duplicated(clinical_sk1_master$Study.ID)], ]

write.csv(dup_rows, file.path(qc.dir, "duplicated_patients_sk1.csv"))

# view first 10 counts to check that the counts are different in these samples
dup_rows_counts <- counts_sk1[1:10,row.names(dup_rows)]
write.csv(dup_rows_counts, file.path(qc.dir, "duplicated_patients_sk1_first10counts.csv"))

# Library size and gene counts
dge <- DGEList(counts = counts_sk1)
cpm_mat <- cpm(dge, log = FALSE)

subset = counts_sk1[,row.names(dup_rows)]
cpm_mat <- cpm_mat[,colnames(subset)]

qc_df <- data.frame(
  patient = clinical_sk1_master[colnames(subset), "Study.ID"],
  batch = 1,
  sample = colnames(subset),
  lib_size = colSums(subset),
  detected_genes = colSums(cpm_mat > 1),
  prop_zero = colSums(subset == 0) / nrow(subset)
)



dup_cor <- qc_df %>%
  group_by(patient) %>%
  summarise(
    sample1 = sample[1],
    sample2 = sample[2],
    pearson = cor(log2(subset[, sample[1]] + 1),
                  log2(subset[, sample[2]] + 1),
                  method = "pearson"),
    spearman = cor(log2(subset[, sample[1]] + 1),
                   log2(subset[, sample[2]] + 1),
                   method = "spearman")
  )

qc_summary <- qc_df %>%
  left_join(dup_cor[,c("patient", "pearson", "spearman")], by = "patient") %>%
  arrange(patient)

dup_patients_qc <- data.frame()
dup_patients_qc <- rbind(dup_patients_qc, qc_summary)


# brush (sherlock 2 and 3) -------------------------------------
clinical_brush_master <- clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Brush"),] #273 samples

dup_rows <- clinical_brush_master[clinical_brush_master$Study.ID %in% 
                                    clinical_brush_master$Study.ID[duplicated(clinical_brush_master$Study.ID)], ]

write.csv(dup_rows, file.path(qc.dir, "duplicated_patients_brush.csv"))

# view first 10 counts to check that the counts are different in these samples
dup_rows_counts <- counts_merged[1:10,row.names(dup_rows)]
write.csv(dup_rows_counts, file.path(qc.dir, "duplicated_patients_brush_first10counts.csv"))


# Library size and gene counts
dge <- DGEList(counts = counts_merged)
cpm_mat <- cpm(dge, log = FALSE)

subset = counts_merged[,row.names(dup_rows)]
cpm_mat <- cpm_mat[,colnames(subset)]

qc_df <- data.frame(
  patient = clinical_brush_master[colnames(subset), "Study.ID"],
  batch = clinical_brush_master[colnames(subset), "batch"],
  sample = colnames(subset),
  lib_size = colSums(subset),
  detected_genes = colSums(cpm_mat > 1),
  prop_zero = colSums(subset == 0) / nrow(subset)
)

dup_cor <- qc_df %>%
  group_by(patient) %>%
  summarise(
    sample1 = sample[1],
    sample2 = sample[2],
    pearson = cor(log2(subset[, sample[1]] + 1),
                  log2(subset[, sample[2]] + 1),
                  method = "pearson"),
    spearman = cor(log2(subset[, sample[1]] + 1),
                   log2(subset[, sample[2]] + 1),
                   method = "spearman")
  )

qc_summary <- qc_df %>%
  left_join(dup_cor[,c("patient", "pearson", "spearman")], by = "patient") %>%
  arrange(patient)

dup_patients_qc <- rbind(dup_patients_qc, qc_summary)



#biopt (sherlock 2 and 3) -----------------------------------------
clinical_biopt_master <- clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Biopt"),] #277 samples

dup_rows <- clinical_biopt_master[clinical_biopt_master$Study.ID %in% 
                                    clinical_biopt_master$Study.ID[duplicated(clinical_biopt_master$Study.ID)], ]
write.csv(dup_rows, file.path(qc.dir, "duplicated_patients_biopt.csv"))

# view first 10 counts to check that the counts are different in these samples
dup_rows_counts <- counts_merged[1:10,row.names(dup_rows)]
write.csv(dup_rows_counts, file.path(qc.dir, "duplicated_patients_biopt_first10counts.csv"))


# Library size and gene counts
dge <- DGEList(counts = counts_merged)
cpm_mat <- cpm(dge, log = FALSE)

subset = counts_merged[,row.names(dup_rows)]
cpm_mat <- cpm_mat[,colnames(subset)]

qc_df <- data.frame(
  patient = clinical_biopt_master[colnames(subset), "Study.ID"],
  batch = clinical_biopt_master[colnames(subset), "batch"],
  sample = colnames(subset),
  lib_size = colSums(subset),
  detected_genes = colSums(cpm_mat > 1),
  prop_zero = colSums(subset == 0) / nrow(subset)
)


dup_cor <- qc_df %>%
  group_by(patient) %>%
  summarise(
    sample1 = sample[1],
    sample2 = sample[2],
    pearson = cor(log2(subset[, sample[1]] + 1),
                  log2(subset[, sample[2]] + 1),
                  method = "pearson"),
    spearman = cor(log2(subset[, sample[1]] + 1),
                   log2(subset[, sample[2]] + 1),
                   method = "spearman")
  )

qc_summary <- qc_df %>%
  left_join(dup_cor[,c("patient", "pearson", "spearman")], by = "patient") %>%
  arrange(patient)

dup_patients_qc <- rbind(dup_patients_qc, qc_summary)
write.csv(dup_patients_qc, file.path(qc.dir, "duplicated_patients_qc_sk1sk2sk3.csv"))

#to keep (higher libsize)
best_libsize_samples <- dup_patients_qc %>%
  group_by(patient) %>%
  slice_max(lib_size, n = 1) %>%
  pull(sample)

# to drop (lower lib size)
dup_samples_to_drop <- setdiff(dup_patients_qc$sample, best_libsize_samples)


# > dup_samples_to_drop
# [1] "LIB5426587_SAM24375559" "106076-002-014"         "106076-002-120"
# [4] "107165-001-041"         "106076-002-009"         "106076-002-003"
# [7] "106076-002-017"         "107165-001-079"         "106076-002-159"
# [10] "106076-002-199"

#subset to exclude
clinical_sk1_master <- clinical_sk1_master[-which(row.names(clinical_sk1_master) %in% dup_samples_to_drop),] #168 -> 167
counts_sk1

clinical_brushbiopt_master <- clinical_brushbiopt_master[-which(row.names(clinical_brushbiopt_master) %in% dup_samples_to_drop),] #541 samples
clinical_brush_master <- clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Brush"),] #270 samples
clinical_biopt_master <- clinical_brushbiopt_master[which(clinical_brushbiopt_master$sampletype == "Biopt"),] #271 samples

counts_merged <- counts_merged[,row.names(clinical_brushbiopt_master)] #541 samples


combat.processed.data.dir <- file.path(data.dir, "processed", "datawrangling_qc", "combat_results")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_merged_pre_combat.csv")) #541 samples, 64489 genes
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_merged_pre_combat.rds"))


# Combat ================================================================================

batch_name <- data.frame(id = colnames(counts_merged))
batch_name[which(batch_name$id %in% colnames(counts_sk2)), "batch"] <- 2
batch_name[which(batch_name$id %in% colnames(counts_sk3)), "batch"] <- 3

counts_combat <- ComBat_seq(counts_merged, batch = batch_name$batch) #genes

write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_combat.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_combat.rds"))


#Note; specifying group = sampletype in combat_seq would mean that cmbat esitmates batch effects within each sampletype and remove the within-group batch effects 
# but this would miss correcting for batch effects between sampletypes 


# Save files ================================================================================

# sherlock2+3 samples (brush + biopt)
if(!exists(file.path(file.path(postQC.data.dir, "master")))) dir.create(file.path(postQC.data.dir, "master"))
write.csv(clinical_brushbiopt_master, file.path(postQC.data.dir, "master","clinical_brushbiopt_master.csv"))
saveRDS(clinical_brushbiopt_master, file.path(postQC.data.dir,  "master","clinical_brushbiopt_master.rds"))

saveRDS(clinical_brush_master, file.path(postQC.data.dir,  "master","clinical_brush_master.rds"))
write.csv(clinical_brush_master, file.path(postQC.data.dir,  "master","clinical_brush_master.csv"))

saveRDS(clinical_biopt_master, file.path(postQC.data.dir,  "master","clinical_biopt_master.rds"))
write.csv(clinical_biopt_master, file.path(postQC.data.dir,  "master","clinical_biopt_master.csv"))


counts_brush_combat <- counts_combat[,row.names(clinical_brush_master)] #270 patients
saveRDS(counts_brush_combat, file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
write.csv(counts_brush_combat, file.path(combat.processed.data.dir, "counts_brush_combat.csv"))

counts_biopt_combat <- counts_combat[,row.names(clinical_biopt_master)] #271 patients 
saveRDS(counts_biopt_combat, file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
write.csv(counts_biopt_combat, file.path(combat.processed.data.dir,  "counts_biopt_combat.csv"))


# sherlock1
saveRDS(clinical_sk1_master, file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_master.rds"))
write.csv(clinical_sk1_master, file.path(processed.data.dir, "SHERLOCK1","clinical_sk1_master.csv"))

counts_sk1 <- counts_sk1[,row.names(clinical_sk1_master)] #167 patients
saveRDS(counts_sk1, file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
write.csv(counts_sk1, file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.csv"))


# Save into common folder for other researchers ================================================================================
## PRE-COMBAT
common.data.dir <- file.path(common.dir, 'data')
if(!exists(common.data.dir)) dir.create(common.data.dir)


write.csv(counts_sk1, file.path(common.data.dir, "counts_sk1.csv"), row.names = TRUE)
# I merged sk2 and sk3 first, then did my sample exclusions on the merged data. So to provide counts_sk2 and counts_sk3 WITH the excluded samples, match samples to the merged data

counts_sk2_share <- counts_sk2[,which(colnames(counts_sk2) %in% colnames(counts_merged))] #327 -> 320 
write.csv(counts_sk2_share, file.path(common.data.dir, "counts_sk2.csv"), row.names = TRUE)

counts_sk3_share <- counts_sk3[,which(colnames(counts_sk3) %in% colnames(counts_merged))] #225 -> 221
write.csv(counts_sk3_share, file.path(common.data.dir, "counts_sk3.csv"), row.names = TRUE)

clinical_master <- rbind(cbind(clinical_sk1_master[which(row.names(clinical_sk1_master) %in% colnames(counts_sk1)),-which(colnames(clinical_sk1_master) == "rnaseq_id")]),
                         cbind(clinical_sk2_master[which(row.names(clinical_sk2_master) %in% colnames(counts_sk2_share)),], batch = 2),
                         cbind(clinical_sk3_master[which(row.names(clinical_sk3_master) %in% colnames(counts_sk3_share)),], batch = 3)) #541 samples total

write.csv(clinical_master, file.path(common.data.dir, "clinical_master.csv"), row.names = TRUE)



# ================================================================================== #
# 3.1 Subset master clinical file for main variables ===============================
# ================================================================================== #
cat("Starting 3.1. Subset main clinical file", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

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
write.csv(clinical_brushbiopt, file.path(postQC.data.dir, "clinical_brushbiopt_simple.csv"))
saveRDS(clinical_brushbiopt, file.path(postQC.data.dir, "clinical_brushbiopt_simple.rds")) #541 samples

saveRDS(clinical_brush, file.path(postQC.data.dir, "clinical_brush_simple.rds")) #270 samples
saveRDS(clinical_biopt, file.path(postQC.data.dir, "clinical_biopt_simple.rds")) #271 samples



# Same for sherlock1

# ================================================================================== #
# 3.1 Subset master clinical file for main variables ===============================
# ================================================================================== #
cat("Starting 3.1. Subset main clinical file", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

clinical_sk1 <- clinical_sk1_master %>% 
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


clinical_sk1$classification <- make.names(clinical_sk1$classification)
clinical_sk1$smoking_status <- make.names(clinical_sk1$smoking_status)
clinical_sk1[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(clinical_sk1[,c("age", "packyears", "FEV1", "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], function(x) as.numeric(x))

## Save clinical_sk1 (main patient info, smoking, ics and lung function data)
write.csv(clinical_sk1, file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_simple.csv"))
saveRDS(clinical_sk1, file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_simple.rds")) 




cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")