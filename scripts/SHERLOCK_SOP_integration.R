# ================================================================================== #
# SCRIPT SET UP ====================================================================
# ================================================================================== #
setwd("path/to/your/directory")
.libPaths("path/to/your/library")

library("sva")
library("edgeR")

# ====================================================================================================== #
# A) For SHERLOCK2 and SHERLOCK3 integration only ====================================================== #
# ====================================================================================================== #

# SUMMARY----------------------------------------------------------------------------------------------- #
# Use when SHERLOCK1 is not required for main analyses. Here, we integrate SHERLOCK2 and SHERLOCK3, while keeping SHERLOCK1 as a seperate dataset for validation
# If both biopsy and brush samples are of interest, run ComBat_seq on all samples first, then separate by sample type for downstream analyses
# 1) Subset each dataset to only include common genes
# 2) Integrate all 3 batches using ComBat_seq
# ------------------------------------------------------------------------------------------------------ #

counts_sk2 <- #load in sherlock2 counts 
counts_sk3 <- #load in sherlock3 counts 

## 1) Subset each dataset to only include common genes ------------------------------------------------ #
  
counts_merged <- merge(counts_sk2, counts_sk3, by = "row.names", all = FALSE) #all = FALSE keeps only matched rows. nonmatched rows are discarded (TRUE would show NA for unmatched samples)
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1])

# Save pre-combat data, for QC or other use
combat.processed.data.dir <- file.path("path/to/your/directory")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_merged_pre_combat.csv"))
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_merged_pre_combat.rds"))

## 3) Integrate all 3 batches using ComBat_seq   -------------------------------------------------------------- #
batch <- c(rep(2, ncol(counts_sk2)), rep(3, ncol(counts_sk3)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

# Save batch corrected, integrated counts
write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_combat.csv")) 
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_combat.rds"))






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

counts_sk1 <- #load in sherlock1 counts 
counts_sk2_brush <- #load in sherlock2 (brush only) 
counts_sk3_brush <- #load in sherlock3 (brush only)  

clinical_brush <- #load in clinical data

## 1) Use filterByExpr() on each dataset individually --------------------------------------------------- #

# Match clinical file with sampleIDs in counts file, use classification as grouping variable
  
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
combat.processed.data.dir <- file.path("path/to/your/directory")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir)

write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.csv")) 
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.rds"))

## 3) Integrate all 3 batches using ComBat_seq   -------------------------------------------------------------- #
batch <- c(rep(1, ncol(counts_sk1_filtered)), rep(2, ncol(counts_sk2_brush_filtered)), rep(3, ncol(counts_sk3_brush_filtered)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

# Save batch corrected, integrated counts
write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.csv")) 
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.rds"))
