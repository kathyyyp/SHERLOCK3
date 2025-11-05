


# ===================================================================================================== #
# A) For SHERLOCK2 and SHERLOCK3 integration only ====================================================== #
# ===================================================================================================== #
# 1) Subset each dataset to only include common genes
# 2) Integrate all 3 batches using ComBat_seq

# counts_sk2 <- load in sherlock2 counts #69,972 genes, 327 samples
# counts_sk3 <- load in sherlock3 counts #66,138 genes, 225 samples

# ComBat ================================================================================
# Merge and keep matching Ensembl IDs between the two files
### 1) Get common genes between all 3 datasets  ------------------------------------------------------------------------------------------------------
counts_merged <- merge(counts_sk2, counts_sk3, by = "row.names", all = FALSE) #all = FALSE keeps only matched rows. nonmatched rows are discarded (TRUE would show NA for unmatched samples)
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1])

combat.processed.data.dir <- file.path("path/to/your/directory")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir)

# Save pre-combat data, for QC or other use
write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_merged_pre_combat.csv")) #552 samples, 64489 genes
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_merged_pre_combat.rds"))

### 2) Run ComBat_seq  ------------------------------------------------------------------------------------------------------
batch <- c(rep(2, ncol(counts_sk2)), rep(3, ncol(counts_sk3)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

# Save batch corrected, integrated counts
write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_combat.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_combat.rds"))



# ===================================================================================================== #
# B) SHERLOCK1, SHERLOCK2 (Brush) and SHERLOCK3 (Brush) integration =================================== #
# ===================================================================================================== #
# Integration is poor when SHERLOCK1 is included (possibly due to difference in library size) without gene filtering
# 1) Use filterByExpr() on each dataset individually
# 2) Subset for common genes
# 3) Integrate all 3 batches using ComBat_seq

#counts_sk1 <- load in sherlock1 counts #60,683 gene, 169 samples
#counts_sk2_brush <- load in sherlock2 (brush only) counts #69,972 genes, 162 samples
#counts_sk3_brush <- load in sherlock3 (brush only) counts #66,138 genes, 112 samples


# SHERLOCK1 = 169 samples
counts_sk1_keep <- filterByExpr(counts_sk1) 
counts_sk1_filtered <- counts_sk1[counts_sk1_keep, ] 
#60,683 genes --> #18,732 genes


# SHERLOCK2 Brush = 162 samples
counts_sk2_brush_keep <- filterByExpr(counts_sk2_brush)
counts_sk2_brush_filtered <- counts_sk2_brush[counts_sk2_brush_keep, ] 
#69,972 genes -->  #16,612 genes


# SHERLOCK3 Brush = 112 samples
counts_sk3_brush_keep <- filterByExpr(counts_sk3_brush)
counts_sk3_brush_filtered <- counts_sk3_brush[counts_sk3_brush_keep, ] 
#66,138 genes  --> #19,734 genes 

### 1) Get common genes between all 3 datasets  ------------------------------------------------------------------------------------------------------
# Get matching genes in sherlock1 and sherlock2
counts_sk1_sk2_all <- merge(counts_sk1_filtered, counts_sk2_brush_filtered, by = "row.names", all = FALSE) 
row.names(counts_sk1_sk2_all) <- counts_sk1_sk2_all[,1]
counts_sk1_sk2_all <- as.matrix(counts_sk1_sk2_all[,-1])

# Match the sk1_sk2 genes with sk3 (keep only the matched Ensembl IDS)
counts_merged <- merge(counts_sk1_sk2_all,  counts_sk3_brush_filtered,  by = "row.names")
row.names(counts_merged) <- counts_merged[,1]
counts_merged <- as.matrix(counts_merged[,-1]) #15,292 genes total 

combat.processed.data.dir <- file.path("path/to/your/directory")
if(!exists(combat.processed.data.dir))dir.create(combat.processed.data.dir)

# Save pre-combat data, for QC or other use
write.csv(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.csv")) #48,039 genes and 721 samples
saveRDS(counts_merged, file.path(combat.processed.data.dir, "counts_brush_merged_pre_combat_filterbyexpr.rds"))

### 2) Run ComBat_seq  ------------------------------------------------------------------------------------------------------
batch <- c(rep(1, ncol(counts_sk1_filtered)), rep(2, ncol(counts_sk2_brush_filtered)), rep(3, ncol(counts_sk3_brush_filtered)))
counts_combat <- ComBat_seq(counts_merged, batch = batch) #genes

# Save batch corrected, integrated counts
write.csv(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.csv")) #552 samples, 64489 genes
saveRDS(counts_combat, file.path(combat.processed.data.dir, "counts_brush_combat_filterbyexpr.rds"))
