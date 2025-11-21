# SHERLOCK3 - Validation of differential expression results with SHERLOCK1
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
library("ggrepel")
library("ggfortify")
library("stringr")
library("EnsDb.Hsapiens.v79")
library("ggplot2")
library("ggpubr") #
library("edgeR")
library("DESeq2")
library("tidyverse")
library("PCAtools")
library("GSVA")

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory
setwd(main.dir)

#Data directory
data.dir <- file.path(main.dir, "data")
processed.data.dir <- file.path(data.dir,"processed")
postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir, "output")

#Diffexp directory (ensembl ids included) - running analysis on all genes (hgnc_only was just for plotting)
diffexp.dir <- file.path(output.dir, "diffexp")

diffexp.results.dir <- file.path(diffexp.dir, "results")

diffexp.figures.dir <- file.path(diffexp.dir, "figures")

#Diffexp directory (hgnc only)
diffexp.hgnconly.dir <- file.path(output.dir, "diffexp_hgnc_only")
if(!exists(diffexp.hgnconly.dir)) dir.create(diffexp.hgnconly.dir, recursive = TRUE)


sherlock1.dir <- file.path(output.dir, "sherlock1")

sherlock1.diffexp.dir <- file.path(sherlock1.dir, "diffexp")

validation.dir <- file.path(output.dir, "validation_sherlock1", "diffexp")
if(!exists(validation.dir)) dir.create(validation.dir, recursive = TRUE)

validation.figures.dir <- file.path(validation.dir, "figures")
if(!exists(validation.figures.dir)) dir.create(validation.figures.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))

# ##-- Post batch correction
counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))

clinical<- readRDS(file.path(postQC.data.dir, "clinical_brushbiopt_simple.rds")) #552 samples
clinical_brush <- readRDS(file.path(postQC.data.dir, "clinical_brush_simple.rds"))
clinical_biopt <- readRDS(file.path(postQC.data.dir, "clinical_biopt_simple.rds"))

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


#Use results with ensembl
listofresults <- readRDS(file.path(diffexp.results.dir, "listofresults.RDS"))

# ================================================================================== #
# 2. GSVA - SHERLOCK3 SIGNATURES IN SHERLOCK1 ======================================
# ================================================================================== #
cat("Starting 2. GSVA - SHERLOCK3 SIGNATURES IN SHERLOCK1", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

validation.gsva.dir <- file.path(validation.figures.dir, "gsva")
if(!exists(validation.gsva.dir)) dir.create(validation.gsva.dir)


sherlock1_counts<- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
sherlock1_clinical <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_simple.rds"))

sherlock1_expression <- voom(sherlock1_counts)

#SEVERE COPD SIGNATURE FROM SHERLOCK3
#gsva - get genes that were upregulated in brushing and biopsy severe copd-control

brush_contrast1_res <- listofresults[["brush"]][["tT2"]][["contrast1"]]
brush_contrast1_up <- row.names(brush_contrast1_res)[which(brush_contrast1_res$Legend == "Upregulated")]
brush_contrast1_down <- row.names(brush_contrast1_res)[which(brush_contrast1_res$Legend == "Downregulated")]

biopt_contrast1_res <- listofresults[["biopt"]][["tT2"]][["contrast1"]]
biopt_contrast1_up <- row.names(biopt_contrast1_res)[which(biopt_contrast1_res$Legend == "Upregulated")]
biopt_contrast1_down <- row.names(biopt_contrast1_res)[which(biopt_contrast1_res$Legend == "Downregulated")]

#MILD-MODERATE COPD SIGNATURE FROM SHERLOCK3
brush_contrast2_res <- listofresults[["brush"]][["tT2"]][["contrast2"]]
brush_contrast2_up <- row.names(brush_contrast2_res)[which(brush_contrast2_res$Legend == "Upregulated")]
brush_contrast2_down <- row.names(brush_contrast2_res)[which(brush_contrast2_res$Legend == "Downregulated")]

biopt_contrast2_res <- listofresults[["biopt"]][["tT2"]][["contrast2"]]
biopt_contrast2_up <- row.names(biopt_contrast2_res)[which(biopt_contrast2_res$Legend == "Upregulated")]
biopt_contrast2_down <- row.names(biopt_contrast2_res)[which(biopt_contrast2_res$Legend == "Downregulated")]

#UP SIGNATURE
listofgroups_up <- list(brush_contrast1_up,
                        biopt_contrast1_up,
                        brush_contrast2_up,
                        biopt_contrast2_up)

#DOWN SIGNATURE
listofgroups_down <- list(brush_contrast1_down,
                          biopt_contrast1_down,
                          brush_contrast2_down,
                          biopt_contrast2_down)


#RUN GSVA


#GSVA to see if our bronchial brushing signature (current SHERLOCK3 data) is retained in the SHERLOCK1 bronchial brushing data
gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None")


gsva_direction <- function(direction){
  
  if(direction == "up"){
    gsva_res <- gsva(as.matrix(sherlock1_expression), listofgroups_up, mx.diff = TRUE)
  }
  
  if(direction == "down"){
    gsva_res <- gsva(as.matrix(sherlock1_expression), listofgroups_down, mx.diff = TRUE)
  }
  
  
  #PLOT GSVA
  gsva_res=t(gsva_res) #the results tell us how much the set of genes was represented in each sample. ie. enrichment score of 0.9 is high- meaning the genes of interest showed up alot in sample X - now when we group the samples by copd and non copd, we can see whether certain genes are enriched in samples with or without copd
  
  
  # # FOR SHERLOCK1 VS SHERLOCK3
  # # COLNAMES FOR LISTOFGROUPS_UP
  if(direction == "up"){
    colnames(gsva_res)=c("Brush_Severe.vs.Control_Up",
                         "Biopsy_Severe.vs.Control_Up",
                         "Brush_MildModerate.vs.Control_Up",
                         "Biopsy_MildModerate.vs.Control_Up")
  }
  
  
  if(direction == "down"){
    # #COLNAMES FOR LISTOFGROUPS_DOWN
    colnames(gsva_res)=c("Brush_Severe.vs.Control_Down",
                         "Biopsy_Severe.vs.Control_Down",
                         "Brush_MildModerate.vs.Control_Down",
                         "Biopsy_MildModerate.vs.Control_Down")
  }
  
  boxplot_gsva=cbind(gsva = gsva_res,
                     disease= as.character(sherlock1_clinical$classification))
  
  boxplot_gsva <- as.data.frame(boxplot_gsva)

  
  my_comparisons <- list(c("Control", "Mild.moderate.COPD"),
                         c("Control", "Severe.COPD"),
                         c("Mild.moderate.COPD", "Severe.COPD"))
  
  # hist(as.numeric(boxplot_gsva[,1])) not all are normally distributed
  
  
  for (i in 1:4){
    
    x_order <- c('Control', 'Mild.moderate.COPD', 'Severe.COPD')
    
    
    boxplotfinal2 <- ggplot(boxplot_gsva, aes(
      x = factor(disease, level = x_order),
      y = as.numeric(boxplot_gsva[,i]),
      fill = disease)) +
      
      theme_bw()+
      
      gsva_theme +
      
      geom_boxplot(position = position_dodge(1)) +
      
      
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         paired = FALSE,
                         size = 7)+
      
      scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                 "Severe.COPD" = "#F8766D")) +
      

      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      labs(title = paste0("Signature Analysis", "(", colnames(boxplot_gsva)[i], ")")) +
      ylab (label = "Enrichment Score") +
      xlab (label = "Disease (SHERLOCK1)")
    
    
    ggsave(boxplotfinal2, file = file.path(validation.gsva.dir, paste0(colnames(boxplot_gsva)[i], ".png")), width = 3300, height = 2600, units = "px" )
    
  }
} #close function

#Run function with either upregulated or downregulated SHERLOCK3 gene set
gsva_direction(direction = "up")
gsva_direction(direction = "down")

# ================================================================================== #
# 3. BOXPLOTS (TOP SHERLOCK3 GENES ON SHERLOCK1 DATA) =======================================
# ================================================================================== #
cat("Starting 3. BOXPLOTS (TOP SHERLOCK3 GENES ON SHERLOCK1 DATA)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

validation.boxplots.dir <- file.path(validation.figures.dir, "boxplots")
if(!exists(validation.boxplots.dir)) dir.create(validation.boxplots.dir)

# WITH listoftT FROM KATHY'S SHERLOCK1 ANALYSIS
sherlock1_listofresults <- readRDS(file.path(sherlock1.diffexp.dir,"results", "listofresults.RDS")) #pvalues from sherlock1 filterbyExpr() results
listoftT_sherlock1 <- sherlock1_listofresults[["tT"]]

# SHERLOCK3 hgnc only results
listofresults_hgnc <- readRDS(file.path(diffexp.hgnconly.dir,"results","listofresults.RDS")) #top 20 from SHERLOCK3

counts_brush_voom <- voom(counts_brush)
counts_biopt_voom <- voom(counts_biopt)
listofvoomcounts <- list(brush = counts_brush_voom, biopt = counts_biopt_voom)

boxplotclinical <- sherlock1_clinical #sherlock1 clinical
boxplotcounts <- t(as.data.frame(sherlock1_expression$E)) #this is voom normalised #sherlock1 counts


if(all(sherlock1_clinical$rna_seq.sample.id == colnames(as.data.frame(sherlock1_expression$E))) == FALSE){ 
  stop("all(sherlock1_clinical$rna_seq.sample.id == colnames(as.data.frame(sherlock1_expression$E))) == FALSE") }


boxplot <- cbind(boxplotcounts,
                 classification = boxplotclinical$classification)



#for top 20
listofboxplots_contrast <- list()
for (contrast in c("contrast1", "contrast2", "contrast3")){
  top20 <- c()
  top20<- row.names(listofresults_hgnc[["brush"]][["tT2"]][[contrast]])[1:20] #get top 20 from SHERLOCK3 diffexp
  top20 <- hgnc_symbols_db[top20, "SYMBOL"]
  topgenes <- c(top20)
  
  listofboxplots_topgenes <- list()
  for (gene in c(topgenes)){
    
    geneofinterest = gene
    geneofinterestid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == geneofinterest), "GENEID"]
    
    plot <- as.data.frame(boxplot[,c(geneofinterestid,
                                     "classification")])
    
    # plot[which(plot$classification == "Not applicable"), "disease"] <- "Healthy"
    # plot[which(plot$classification == "GOLD 1" |plot$classification == "GOLD 2"), "disease"] <- "GOLD1&2"
    # plot[which(plot$classification == "GOLD 3" |plot$classification == "GOLD 4"), "disease"] <- "GOLD3&4"
    #
    
    colnames(plot)[1] <- "gene"
    
    plot$gene <- as.numeric(plot$gene)
    ### Get P-Values --------------------------------------------------------------------------------------
    # THE TABLES IN STEP 1 AND 2 ARE TO GET X POSITIONS, THE T-TEST VALUES WILL NOT BE USED #
    #### STEP 1) THIS IS  A TABLE OF THE COMPARISONS I ACTUALLY WANT --------------------------------------------------------------------------------------
    stat.table <- plot %>%
      t_test(gene ~ classification,
             comparisons = list(c("Severe.COPD", "Control"),
                                c("Mild.moderate.COPD", "Control"),
                                c("Severe.COPD", "Mild.moderate.COPD")
             ))
    stat.table<- stat.table %>%
      add_xy_position(x = "classification", dodge = 0.8)
    
    # stat.table$contrast <- paste0(stat.table$group1, "-" ,stat.table$group2)
    
    
    
    #### STEP 2) MODIFY STAT TABLE TO INCLUDE COMPARISONS OF INTEREST --------------------------------------------------------------------------------------
    stat.table3 <- stat.table
    
    stat.table3 <- cbind(stat.table3, resultsname = c("contrast1", "contrast2", "contrast3"))
    stat.table3[which(stat.table3$resultsname == "contrast1"),"p"] <- listoftT_sherlock1[["contrast1"]][geneofinterestid[1], "PValue"]
    stat.table3[which(stat.table3$resultsname == "contrast2"),"p"] <- listoftT_sherlock1[["contrast2"]][geneofinterestid[1], "PValue"]
    stat.table3[which(stat.table3$resultsname == "contrast3"),"p"] <- listoftT_sherlock1[["contrast3"]][geneofinterestid[1], "PValue"]
    stat.table3$p <- signif(as.numeric(stat.table3$p), digits = 4)
    stat.table3$y.position <- max(plot[,"gene"]) + 0.025*(max(plot[,"gene"]))
    stat.table3$y.position <- as.numeric(stat.table3$y.position)
    
    
    # x_order <- c('Healthy', 'GOLD1&2', 'GOLD3&4')
    
    boxplotimage <- ggplot(plot, aes(
      x = factor(as.factor(classification)),
      y = gene )) +
      
      #   geom_dotplot(binaxis = "y",
      #                stackdir = "center",
      #                position=position_dodge(0.5),
      #                stackratio = 0.5,
      #                dotsize = 0.7,
      #                aes(fill = smokingstatus)
      # )+
      
      
      #
    #
    geom_boxplot(position=position_dodge(),
                 aes(fill = classification)
    ) +
      
      
      scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                 "Severe.COPD" = "#F8766D")) +
      

      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      labs(title = paste0(geneofinterest, ": SHERLOCK1")) +
      ylab (label = "Expression") +
      xlab (label = "Disease") +
      
      stat_pvalue_manual(stat.table3,
                         label = "p",
                         tip.length = 0.01,
                         bracket.nudge.y = c(0, 0.5, 0.5),
                         size = 6)
    
    listofboxplots_topgenes[[geneofinterest]] <- boxplotimage
    
  }
  listofboxplots_contrast[[contrast]] <-listofboxplots_topgenes
}


### Save boxplots --------------------------------------------
cat("Saving boxplots", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

dir.create(file.path(validation.boxplots.dir, "severevscontrol"), recursive = TRUE)
dir.create(file.path(validation.boxplots.dir, "mildmoderatevscontrol"), recursive = TRUE)
dir.create(file.path(validation.boxplots.dir, "severevsmildmoderate"), recursive = TRUE)

boxplot_theme <- theme(axis.title = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       title = element_text(size = 18),
                       legend.position = "None") 

#Save brush boxplots
for (i in 1:length(listofboxplots_contrast[["contrast1"]])){
  png(file = paste0(validation.boxplots.dir,"/severevscontrol/brush_",i, ")",names(listofboxplots_contrast[["contrast1"]])[i],".png"),
      height = 14,
      width = 14,
      units = "cm",
      res = 1200)
  print(listofboxplots_contrast[["contrast1"]][[i]] + boxplot_theme) #+ theme_bw()
  dev.off()
}

for (i in 1:length(listofboxplots_contrast[["contrast2"]])){
  png(file = paste0(validation.boxplots.dir,"/mildmoderatevscontrol/brush_",i, ")",names(listofboxplots_contrast[["contrast2"]])[i],".png"),
      height = 14,
      width = 14,
      units = "cm",
      res = 1200)
  print(listofboxplots_contrast[["contrast2"]][[i]] + boxplot_theme) #+ theme_bw()
  dev.off()
}

for (i in 1:length(listofboxplots_contrast[["contrast3"]])){
  png(file = paste0(validation.boxplots.dir,"/severevsmildmoderate/brush_",i, ")",names(listofboxplots_contrast[["contrast3"]])[i],".png"),
      height = 14,
      width = 14,
      units = "cm",
      res = 1200)
  print(listofboxplots_contrast[["contrast3"]][[i]] + boxplot_theme) #+ theme_bw()
  dev.off()
}


# ================================================================================== #
# 4. SHERLOCK1 VS SHERLOCK3 VENN DIAGRAMS ==========================================
# ================================================================================== #
cat("Starting 4. SHERLOCK1 VS SHERLOCK3 VENN DIAGRAMS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

validation.venn.dir <- file.path(validation.figures.dir, "venn")
if(!exists(validation.venn.dir)) dir.create(validation.venn.dir)

# # sCOPD genes not associated with ICS or basal, cilaited and  secretory cell types
# sherlock1_sCOPD_tT2 <- read_xlsx("jos-sherlock/Supplementary Tables - Chapter2.xlsx", sheet = "Table S5", skip = 2, .name_repair = "universal")
# sherlock1_sCOPD_tT2$legend <- ifelse(
#   sherlock1_sCOPD_tT2$FDR < 0.05 & sherlock1_sCOPD_tT2$logFC > 0.25, "Upregulated",
#   ifelse(
#     sherlock1_sCOPD_tT2$FDR < 0.05 & sherlock1_sCOPD_tT2$logFC < -0.25, "Downregulated",
#     "Not Significant"))

#SHERLOCK3
brush_contrast1_up
brush_contrast1_down

#Sherlock1 tT2 results
listoftT2_sherlock1 <- sherlock1_listofresults[["tT2"]]

#Sherlock1 - SEVERE
sherlock1_sCOPD_tT2 <- as.data.frame(listoftT2_sherlock1[["contrast1"]])
sherlock1_sCOPD_up <- row.names(sherlock1_sCOPD_tT2)[which(sherlock1_sCOPD_tT2$Legend == "Upregulated")]
sherlock1_sCOPD_down <- row.names(sherlock1_sCOPD_tT2)[which(sherlock1_sCOPD_tT2$Legend == "Downregulated")]

# length(brush_contrast1_up)
# [1] 171
# > length(brush_contrast1_down)
# [1] 251
# > length(sherlock1_sCOPD_up)
# [1] 104
# > length(sherlock1_sCOPD_down)
# [1] 115

venn_up <- list(SHERLOCK1 = sherlock1_sCOPD_up,
                SHERLOCK3 = brush_contrast1_up)

venn_down <- list(SHERLOCK1 = sherlock1_sCOPD_down,
                  SHERLOCK3 = brush_contrast1_down
)

png(filename = file.path(validation.venn.dir, "severe_up.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggvenn(venn_up, stroke_size = 1.5, fill_color = c("#0073C2FF", "#EFC000FF"),text_size = 6) +
  ggtitle("Severe COPD Signature (Up)") + theme(plot.title = element_text(hjust=0.5, vjust = 3))
dev.off()

png(filename = file.path(validation.venn.dir, "severe_down.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggvenn(venn_down, stroke_size = 1.5, fill_color = c("#0073C2FF", "#EFC000FF"),text_size = 6) +
  ggtitle("Severe COPD Signature (Down)") + theme(plot.title = element_text(hjust=0.5, vjust = 3))
dev.off()

up_intersect <- intersect(sherlock1_sCOPD_up[1:30], brush_contrast1_up[1:30])
down_intersect <- intersect(sherlock1_sCOPD_down[1:50], brush_contrast1_down[1:50])

hgnc_symbols_db[match(up_intersect, hgnc_symbols_db$GENEID), "SYMBOL"]
hgnc_symbols_db[match(down_intersect, hgnc_symbols_db$GENEID), "SYMBOL"]



#Sherlock1 - MILD
sherlock1_mCOPD_tT2 <- as.data.frame(listoftT2_sherlock1[["contrast2"]])
sherlock1_mCOPD_up <- row.names(sherlock1_mCOPD_tT2)[which(sherlock1_mCOPD_tT2$Legend == "Upregulated")]
sherlock1_mCOPD_down <- row.names(sherlock1_mCOPD_tT2)[which(sherlock1_mCOPD_tT2$Legend == "Downregulated")]

# length(brush_contrast1_up)
# [1] 171
# > length(brush_contrast1_down)
# [1] 251
# > length(sherlock1_mCOPD_up)
# [1] 104
# > length(sherlock1_mCOPD_down)
# [1] 115

venn_up <- list(SHERLOCK1 = sherlock1_mCOPD_up,
                SHERLOCK3 = brush_contrast2_up)

venn_down <- list(SHERLOCK1 = sherlock1_mCOPD_down,
                  SHERLOCK3 = brush_contrast2_down
)

png(filename = paste0(validation.venn.dir, "/mild.moderate_up.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggvenn(venn_up, stroke_size = 1.5, fill_color = c("#0073C2FF", "#EFC000FF"),text_size = 6) +
  ggtitle("Mild/Moderate COPD Signature (Up)") + theme(plot.title = element_text(hjust=0.5, vjust = 3))
dev.off()

png(filename = paste0(validation.venn.dir, "/mild.moderate_down.png"), width = 15, height = 12, unit = "cm", res = 1200)
ggvenn(venn_down, stroke_size = 1.5, fill_color = c("#0073C2FF", "#EFC000FF"),text_size = 6) +
  ggtitle("Mild/Moderate COPD Signature (Down)") + theme(plot.title = element_text(hjust=0.5, vjust = 3))
dev.off()

up_intersect <- intersect(sherlock1_mCOPD_up[1:30], brush_contrast1_up[1:30])
down_intersect <- intersect(sherlock1_mCOPD_down[1:50], brush_contrast1_down[1:50])

hgnc_symbols_db[match(up_intersect, hgnc_symbols_db$GENEID), "SYMBOL"]
hgnc_symbols_db[match(down_intersect, hgnc_symbols_db$GENEID), "SYMBOL"]





# ================================================================================== #
# 5. GSVA SHERLOCK3  ===============================================================
# ================================================================================== #
cat("Starting 5. GSVA SHERLOCK3", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

gsva.dir <- file.path(diffexp.figures.dir, "gsva")
if(!exists(gsva.dir)) dir.create(gsva.dir)

#GSVA to see if our mild-moderate bronchial brushing signature is enriched in severe-copd in our same dataset
### !!! do I need to subset for only brushing data???

counts_brush_voom <- voom(counts_brush)
counts_biopt_voom <- voom(counts_biopt)


mild_brush <- list(brush_contrast2_up,
                   brush_contrast2_down)

mild_biopt <- list(biopt_contrast2_up,
                   biopt_contrast2_down)

gsva_sherlock2_func <- function(sampletype){
  
  
  if(sampletype == "brush"){
    gsva_counts_voom <- as.data.frame(counts_brush_voom$E)
    gsva_res <- gsva(as.matrix(gsva_counts_voom), mild_brush, mx.diff = TRUE)
  }
  
  if(sampletype == "biopt"){
    gsva_counts_voom <- as.data.frame(counts_biopt_voom$E)
    gsva_res <- gsva(as.matrix(gsva_counts_voom), mild_biopt, mx.diff = TRUE)
  }
  
  
  gsva_res=t(gsva_res) #the results tell us how much the set of genes was represented in each sample. ie. enrichment score of 0.9 is high- meaning the genes of interest showed up alot in sample X - now when we group the samples by copd and non copd, we can see whether certain genes are enriched in samples with or without copd
  
  
  # COLNAMES FOR MILD_BRUSH
  if(sampletype == "brush"){
    colnames(gsva_res)= c("Brush_MildModerate.vs.Control_Up",
                          "Brush_MildModerate.vs.Control_Down")
    boxplot_gsva=cbind(gsva = gsva_res,
                       disease= as.character(clinical_brush$classification))
  }
  
  if(sampletype == "biopt"){
    colnames(gsva_res)= c("Biopt_MildModerate.vs.Control_Up",
                          "Biopt_MildModerate.vs.Control_Down")
    boxplot_gsva=cbind(gsva = gsva_res,
                       disease= as.character(clinical_biopt$classification))
  }
  
  
  
  
  boxplot_gsva <- as.data.frame(boxplot_gsva)
  
  
  my_comparisons <- list(c("Control", "Severe.COPD"),
                         c("Control", "Mild.moderate.COPD"),
                         c("Severe.COPD", "Mild.moderate.COPD"))
  
  
  
  for (i in 1:2){
    
    x_order <- c('Control', 'Mild.moderate.COPD', 'Severe.COPD')
    
    
    boxplotfinal2 <- ggplot(boxplot_gsva, aes(
      x = factor(disease, level = x_order),
      y = as.numeric(boxplot_gsva[,i]),
      fill = disease)) +
      
      theme_bw()+
      
      gsva_theme +
      
      geom_boxplot(position = position_dodge(1)) +
      
      
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         paired = FALSE,
                         size = 7)+
      
      scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                 "Severe.COPD" = "#F8766D")) +
      
      # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      labs(title = paste0("Signature Analysis", "(", colnames(boxplot_gsva)[i], ")")) +
      ylab (label = "Enrichment Score") +
      xlab (label = "Disease Severity (SHERLOCK3)")
    
    
    ggsave(boxplotfinal2, file = paste0(gsva.dir,"/",colnames(boxplot_gsva)[i], ".png"), width = 3000, height = 2100, units = "px" )
    
    
  }
  
} #close function

gsva_sherlock2_func(sampletype = "brush")
gsva_sherlock2_func(sampletype = "biopt")

# ================================================================================== #
# 6. GSVA IL33 SIGNATURE =========================================================
# ====================================================================================#
cat("Starting 6. GSVA IL33 SIGNATURE", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
IL33_badi <-c("NFKB2", "RELB", "OLR1", "MAP3K8", "NFKBIZ", "PTGS2", "GBP4",
              "CSF2", "GBP2", "FEZ1", "THBS1", "STAT4", "BAMBI", "TGFBR3",
              "PCSK5", "RASL11A", "MUC20", "IRF4", "IL2RA", "MT1E", "CD74",
              "ITGA1", "IL5", "TMEM64", "PLA2G4A", "SLC27A2", "IL17RB",
              "LRIG1", "CCL3", "TNFRSF8", "IL13", "PIM2", "BATF", "C15orf48",
              "EBI3", "IRAK2", "SMAD3", "TNFRSF9", "PDLIM4", "TRAF1", "TNIP1",
              "IL7R", "MT1G", "ZFPM2", "ZC3H12A", "IL1A", "TNFAIP2", "IL32",
              "HLA-B", "LTB", "DSE", "PLA2G4C", "IL1B", "CD200", "IL4I1",
              "ST6GAL1", "NFKBIA", "GCH1", "IL15RA", "SLC7A11", "POU2F2",
              "TNFAIP3", "BIRC3", "ABCG1", "TNFAIP6", "IFIH1", "NFKB1",
              "MARCH3", "DUSP16", "NINJ1", "IL18R1", "RGS2", "TNFSF10", "GBP1",
              "TAP1", "MX1", "CXCL8", "ELL2", "TRIB1", "CSF2RB")

IL33_badi[is.na(match(IL33_badi, hgnc_symbols_db$SYMBOL))]
#IL8 = CXCL8
#SOD2 = REMOVE, couldnt find

IL33_signature <- hgnc_symbols_db[match( IL33_badi, hgnc_symbols_db$SYMBOL), "GENEID"]
IL33_signature <- list(IL33_signature)

gsva_il33_func <- function(sampletype){
  
  
  if(sampletype == "brush"){
    gsva_counts_voom <- as.data.frame(counts_brush_voom$E)
    gsva_res <- gsva(as.matrix(gsva_counts_voom), IL33_signature, mx.diff = TRUE)
    
  }
  
  if(sampletype == "biopt"){
    gsva_counts_voom <- as.data.frame(counts_biopt_voom$E)
    gsva_res <- gsva(as.matrix(gsva_counts_voom), IL33_signature, mx.diff = TRUE)
  }
  
  
  gsva_res=t(gsva_res) #the results tell us how much the set of genes was represented in each sample. ie. enrichment score of 0.9 is high- meaning the genes of interest showed up alot in sample X - now when we group the samples by copd and non copd, we can see whether certain genes are enriched in samples with or without copd
  
  if(sampletype == "brush"){
    colnames(gsva_res)= c("Brush_IL33_signature")
    boxplot_gsva=cbind(gsva = gsva_res,
                       disease= as.character(clinical_brush$classification))}
  
  if(sampletype == "biopt"){
    colnames(gsva_res)= c("Biopt_IL33_signature")
    boxplot_gsva=cbind(gsva = gsva_res,
                       disease= as.character(clinical_biopt$classification))}
  
  
  
  
  
  boxplot_gsva <- as.data.frame(boxplot_gsva)
  
  
  my_comparisons <- list(c("Control", "Severe.COPD"),
                         c("Control", "Mild.moderate.COPD"),
                         c("Severe.COPD", "Mild.moderate.COPD"))
  
  
  
  gsva.dir.il33 <- file.path(gsva.dir, "il33_geneset")
  if(!exists(gsva.dir.il33)) dir.create(gsva.dir.il33)
  
  for (i in 1){
    
    x_order <- c('Control', 'Mild.moderate.COPD', 'Severe.COPD')
    
    
    boxplotfinal2 <- ggplot(boxplot_gsva, aes(
      x = factor(disease, level = x_order),
      y = as.numeric(boxplot_gsva[,1]),
      fill = disease)) +
      
      theme_bw()+
      
      gsva_theme +
      
      geom_boxplot(position = position_dodge(1)) +
      
      
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         paired=FALSE,
                         size = 6,
                         step.increase = 0.08)+
      
      scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                 "Severe.COPD" = "#F8766D")) +
      
      # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      labs(title = paste0("Signature Analysis", "(", colnames(boxplot_gsva)[i], ")")) +
      ylab (label = "Enrichment Score") +
      xlab (label = "Disease Severity (SHERLOCK3)")
    
    
    ggsave(boxplotfinal2, file = paste0(gsva.dir.il33,"/",colnames(boxplot_gsva)[i], ".png"), width = 3000, height = 1950, units = "px" )
    
  }
}

gsva_il33_func(sampletype = "brush")
gsva_il33_func(sampletype = "biopt")

 
cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Resources consumed by job 1664106 for user umcg-kphung running on compute node nb-node-a03:
#   ================================================================================================================
#   JobID          Elapsed   AllocCPUS  AveCPU    ReqMem  MaxVMSize  MaxRSS    MaxDiskRead  MaxDiskWrite
# 1664106.batch  00:05:28  1          00:04:56          8959288K   8500612K  1040.13M     31.87M
