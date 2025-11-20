# ================================================================================== #
# SCRIPT SET UP ====================================================================
# ================================================================================== #
main.dir <- "/groups/umcg-griac/tmp02/projects/SHERLOCK_2025"
# .libPaths("path/to/your/library")
setwd(main.dir)

library("sva")
library("edgeR")
library("dplyr")


# ====================================================================================================== #
# A) For SHERLOCK2 and SHERLOCK3 integration only ====================================================== #
# ====================================================================================================== #

# SUMMARY----------------------------------------------------------------------------------------------- #
# Use when SHERLOCK1 is not required for main analyses. Here, we integrate SHERLOCK2 and SHERLOCK3, while keeping SHERLOCK1 as a seperate dataset for validation
# If both biopsy and brush samples are of interest, run ComBat_seq on all samples first, then separate by sample type for downstream analyses
# 1) Subset each dataset to only include common genes
# 2) Integrate batches using ComBat_seq
# 3) Split into brushings and biopsies
# ------------------------------------------------------------------------------------------------------ #

data.dir <- file.path(main.dir, "data")
clinical_master <- read.csv(file.path(data.dir, "clinical_master.csv"), row.names = 1)
counts_sk2 <- read.csv(file.path(data.dir, "counts_sk2.csv"), row.names = 1, check.names = FALSE) #load in sherlock2 counts (320 samples)
counts_sk3 <- read.csv(file.path(data.dir, "counts_sk3.csv"), row.names = 1, check.names = FALSE) #load in sherlock3 counts (221 samples)

## 1) Subset each dataset to only include common genes ------------------------------------------------ #
counts_merged <- merge(counts_sk2, counts_sk3, by = "row.names", all = FALSE) #all = FALSE keeps only matched rows. nonmatched rows are discarded (TRUE would show NA for unmatched samples)
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1])

# Save pre-combat data, for QC or other use
combat.data.dir <- file.path(data.dir, "sherlock2_3_combat")
if(!exists(combat.data.dir))dir.create(combat.data.dir)

write.csv(counts_merged, file.path(combat.data.dir, "sherlock2_3_counts_pre_combat.csv"))


## 2) Integrate batches using ComBat_seq   -------------------------------------------------------------- #
batch <- clinical_master[colnames(counts_merged),"batch"]
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

# Save batch corrected, integrated counts
write.csv(counts_combat, file.path(combat.data.dir, "sherlock2_3_counts_integrated.csv")) 


## 3) Split into brushings and biopsies   -------------------------------------------------------------- #
brush_ids <- clinical_master %>% 
  filter(batch != 1,
         sampletype == "Brush") %>% 
  row.names()
counts_brush_combat <- counts_combat[,brush_ids] #270 samples
write.csv(counts_brush_combat, file.path(combat.data.dir, "sherlock2_3_brush_counts_integrated.csv"))

biopt_ids <- clinical_master %>% 
  filter(batch != 1,
         sampletype == "Biopt") %>% 
  row.names()
counts_biopt_combat <- counts_combat[,biopt_ids] #271 samples 
write.csv(counts_biopt_combat, file.path(combat.data.dir,  "sherlock2_3_biopt_counts_integrated.csv"))




# ===================================================================================================== #
# B) SHERLOCK1, SHERLOCK2 (Brush) and SHERLOCK3 (Brush) integration =================================== #
# ===================================================================================================== #

# SUMMARY---------------------------------------------------------------------------------------------- #
# Use when SHERLOCK1 is required for main analyses
# Integration is poor when SHERLOCK1 is included in the method above (possibly due to difference in library size) 
# Filtering out lowly expressed genes allows for better integration
# 1) Use filterByExpr() on each dataset individually
# 2) Subset for common genes
# 3) Integrate all 3 batches using ComBat_seq
# ----------------------------------------------------------------------------------------------------- #

data.dir <- file.path(main.dir, "data")
clinical_master <- read.csv(file.path(data.dir, "clinical_master.csv"), row.names = 1)
clinical_brush <- clinical_master %>% 
  filter(sampletype == "Brush")

#load in sherlock1 counts 
counts_sk1 <- read.csv(file.path(data.dir, "counts_sk1.csv"), row.names = 1, check.names = FALSE) #load in sherlock1 counts (167 samples)

#load in sherlock2 counts and subset for brushes
counts_sk2 <- read.csv(file.path(data.dir, "counts_sk2.csv"), row.names = 1, check.names = FALSE) #load in sherlock2 counts (320 samples)
sk2_brush_ids <- clinical_brush %>% 
  filter(batch == 2) %>%  
  row.names() 
counts_sk2_brush <- counts_sk2[,sk2_brush_ids] #160 brush samples

#load in sherlock3 counts and subset for brushes
counts_sk3 <- read.csv(file.path(data.dir, "counts_sk3.csv"), row.names = 1, check.names = FALSE) #load in sherlock3 counts (221 samples)
sk3_brush_ids <- clinical_brush %>% 
  filter(batch == 3) %>% 
  row.names()
counts_sk3_brush <- counts_sk3[,sk3_brush_ids] #221 brush samples


## 1) Use filterByExpr() on each dataset individually --------------------------------------------------- #

# Match clinical file with sampleIDs in counts file, use classification as grouping variable
# classification is disease severity, either "Control", "Mild.moderate.COPD" or "Severe.COPD"
  
# SHERLOCK1 
counts_sk1_keep <- filterByExpr(counts_sk1, group = clinical_brush[colnames(counts_sk1), "classification"])   
counts_sk1_filtered <- counts_sk1[counts_sk1_keep, ] 


# SHERLOCK2 Brush 
counts_sk2_brush_keep <- filterByExpr(counts_sk2_brush, group = clinical_brush[colnames(counts_sk2_brush), "classification"] )
counts_sk2_brush_filtered <- counts_sk2_brush[counts_sk2_brush_keep, ] 


# SHERLOCK3 Brush 
counts_sk3_brush_keep <- filterByExpr(counts_sk3_brush, group = clinical_brush[colnames(counts_sk3_brush), "classification"])
counts_sk3_brush_filtered <- counts_sk3_brush[counts_sk3_brush_keep, ] 

## 2) Subset for common genes  ---------------------------------------------------------------------------- #
# Get matching genes in sherlock1 and sherlock2
counts_sk1_sk2_all <- merge(counts_sk1_filtered, counts_sk2_brush_filtered, by = "row.names", all = FALSE) 
row.names(counts_sk1_sk2_all) <- counts_sk1_sk2_all[,1]
counts_sk1_sk2_all <- as.matrix(counts_sk1_sk2_all[,-1])

# Match the sk1_sk2 genes with sk3 (keep only the matched Ensembl IDS)
counts_merged <- merge(counts_sk1_sk2_all,  counts_sk3_brush_filtered,  by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1]) 

# Save pre-combat data, for QC or other use
combat.data.dir <- file.path(data.dir, "sherlock1_2_3_combat")
if(!exists(combat.data.dir))dir.create(combat.data.dir)

write.csv(counts_merged, file.path(combat.data.dir, "sherlock1_2_3_counts_brush_pre_combat.csv")) 

## 3) Integrate all 3 batches using ComBat_seq   -------------------------------------------------------------- #
batch <- clinical_brush[colnames(counts_merged),"batch"]

counts_combat <- ComBat_seq(counts_merged, batch) #genes

# Save batch corrected, integrated counts
write.csv(counts_combat, file.path(combat.data.dir, "sherlock1_2_3_counts_brush_integrated.csv")) 

