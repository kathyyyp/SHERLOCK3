# SHERLOCK3 - QC on merged SHERLOCK2&3 data and batch corrected with Combat_seq 
options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("rstatix")
library("ggplot2")
library("ggvenn")
library("ggrepel")
library("ggfortify")
library("edgeR")
library("tidyverse")
library("ggpubr")

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


# ================================================================================== #
# 3. CELLULAR DECONVLUTION ==========================================================
# ================================================================================== #
cat("Starting 3. CELLULAR DECONVLUTION ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Mixture/Count Matrix - CPM normalise
# raw_expression <- read.table(file.path(data.dir,"raw","20241109_Sherlock2_readcount_UMI_dedup.txt"), header = TRUE, check.names = FALSE, row.names = 1)

counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
clinical<- readRDS(file.path(postQC.data.dir, "clinical_brushbiopt_simple.rds")) #552 samples
hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))

expression_cpm <- cpm(counts)
expression_cpm <- expression_cpm[hgnc_symbols_db$GENEID,] #552 samples

hidesample <- as.data.frame(cbind(colnames(expression_cpm), seq(1:ncol(expression_cpm))))
colnames(expression_cpm) <- hidesample$V2

expression_cpm <- as.matrix(cbind(GeneSymbol = hgnc_symbols_db$SYMBOL, expression_cpm))

if(!exists(file.path(data.dir, "processed", "celldecon"))) dir.create(file.path(data.dir, "processed", "celldecon"), recursive = TRUE)
write.table(expression_cpm, file = file.path(data.dir, "processed", "celldecon", "expression_cpm.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

## RUN CIBERSORT with ref1_sigmatrix from fia (from Nov 2024, "ref1_1_1411") - check if there is a new one - this is for bronchial
## 1. Create signature matrix was done already - created Ref1_sigmatrix
## Go to 2. Impute cell fractions
## Batch correction, S-mode, 100 permutations
## refmatrix = Ref_1_1411 (single cell reference matrix)
## sigmatrix = Ref1_sigmatrix (signature matrix)
## mixture = expression_cpm


#Upload results from cibersort (results.txt) to processed data folder and run
cibersort_results <- read.delim(file.path(data.dir, "processed", "celldecon", "results_20251114.txt"))
row.names(cibersort_results) <- hidesample$V1
cibersort_results <- cibersort_results[,-1]

row.names(cibersort_results) == rownames(clinical)
cibersort_results <- cibersort_results[row.names(clinical),]

clinical_celldecon <- as.data.frame(cbind(clinical,cibersort_results))


# #split into biopt and brush
clinical_biopt_celldecon <- clinical_celldecon[which(clinical_celldecon$sampletype == "Biopt"),]

clinical_brush_celldecon <- clinical_celldecon[which(clinical_celldecon$sampletype == "Brush"),]


## Cell Decon vs Disease Plot --------------------------------------------------------------------------------------
boxplot_celldecon= cbind(id = row.names(clinical_celldecon),
                         sampletype = clinical_celldecon$sampletype,
                         classification = clinical$classification,
                         clinical_celldecon[,which(colnames(clinical_celldecon) == "Alveolar.macrophages"):ncol(clinical_celldecon)])
boxplot_celldecon <- as.data.frame(boxplot_celldecon)

boxplot_celldecon_long <- pivot_longer(boxplot_celldecon, cols = 4:21)

ciber_brush <- boxplot_celldecon_long[which(boxplot_celldecon_long$sampletype == "Brush"),]
ciber_biopt <- boxplot_celldecon_long[which(boxplot_celldecon_long$sampletype == "Biopt"),]

ciber_brush <- ciber_brush[order(ciber_brush$classification),]
colnames(ciber_brush)[7] <- "celltype"
colnames(ciber_brush)[8] <- "proportion"

ciber_biopt <- ciber_biopt[order(ciber_biopt$classification),]
colnames(ciber_biopt)[7] <- "celltype"
colnames(ciber_biopt)[8] <- "proportion"

celldecon.dir <- file.path(output.dir, "celldecon")
if(!exists(celldecon.dir)) dir.create(celldecon.dir, recursive = TRUE)

celldecon.figures.dir <- file.path(celldecon.dir, "figures")
if(!exists(celldecon.figures.dir)) dir.create(celldecon.figures.dir, recursive = TRUE)

saveRDS(ciber_biopt,file.path(celldecon.figures.dir,"celldecon_ciber_biopt.rds"))
saveRDS(ciber_brush,file.path(celldecon.figures.dir,"celldecon_ciber_brush.rds"))

#Normality test
ggqqplot(residuals(lm(proportion ~ classification, data = ciber_biopt)))

### Get P Values --------------------------------------------------------------------------------
#Step 1) Get list of comparisons
celldecon_stats <- ciber_biopt %>%
  group_by(celltype) %>%
  wilcox_test(proportion ~ classification,
              paired = FALSE,
              p.adjust.method = "BH")%>%
  add_xy_position(x = "celltype")


#Manually change y positions because right now they are unecessarily staggered
#Don't want all the p-value brackets to be adjacent in the plot because the boxplots that are far down are then very far from the brackets up top
#Get max proportion for each cell type and stagger the brackets for each cell type
maxRows <- by(data = ciber_biopt, INDICES = ciber_biopt$celltype, FUN = function(X) {X[which.max(X$proportion),]})
#by() applies a funcion to subsets of a dataframe, grouped by a factor or variable (in this case cell type)
#celltype is the grouping factor, for each unique celltype, the function will be executed
#function(X) defines a function that is applied to every susbet of data (here, X is all rows corresponding to a specific celltype), those rows become X and the function looks for the max proportion in X
max_per_celltype <- as.data.frame(do.call("rbind", maxRows))

celldecon_stats$y.position <- as.numeric(max_per_celltype[match(celldecon_stats$celltype, max_per_celltype$celltype), "proportion"])
celldecon_stats$y.position <- celldecon_stats$y.position + (celldecon_stats$y.position)*0.05

celldecon_stats$contrast <- paste0(celldecon_stats$group1, "_" ,celldecon_stats$group2)

celldecon_stats <- as.data.frame(celldecon_stats)
celldecon_stats[which(celldecon_stats$contrast == "Control_Severe.COPD"),"y.position"] <- celldecon_stats[which(celldecon_stats$contrast == "Control_Severe.COPD"),"y.position"] +0.02

#Only label the significant comparisons
celldecon_stats2 <- celldecon_stats[which(celldecon_stats$p.adj.signif != "ns"),]

png(file.path(celldecon.figures.dir,"celldecon_ciber_biopt.png"), width = 60, height = 30, units = "cm",res = 1200)

ggplot(ciber_biopt,
       aes(x=celltype,
           y=proportion))+
  geom_boxplot(aes(fill=classification))+ #have to put fill here, because if we put it in the ggplot() aes, it becomes a global aes and needs 'classification' in each subsequent step. theres no 'classification' column in celldecon_stats so it doesnt work if its global
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5,
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour =  "white"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "grey"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.margin = margin(0,0,1,2, "cm"))+
  
  stat_pvalue_manual(celldecon_stats2,
                     label = "p.adj.signif",
                     tip.length = 0.01,
                     # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                     #  bracket.shorten = 0.1,
                     size = 3)+
  scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                             "Severe.COPD" = "#F8766D")) +
  
  xlab("Cell Type") +
  ylab("Proportion")+
  labs(fill = "Classification")+
  
  # stat_compare_means(aes(group = classification),
  #                    comparisons = my_comparisons,
  #                    method = "wilcox.test")+
  
  ggtitle("Cellular Proportions (Biopt)")

dev.off()


#Normality test
ggqqplot(residuals(lm(proportion ~ classification, data = ciber_brush)))

### Get P Values --------------------------------------------------------------------------------
#Step 1) Get list of comparisons
celldecon_stats <- ciber_brush %>%
  group_by(celltype) %>%
  wilcox_test(proportion ~ classification,
              paired = FALSE,
              p.adjust.method = "BH")%>%
  add_xy_position(x = "celltype")


#Manually change y positions because right now they are unecessarily staggered
#Don't want all the p-value brackets to be adjacent in the plot because the boxplots that are far down are then very far from the brackets up top
#Get max proportion for each cell type and stagger the brackets for each cell type
maxRows <- by(data = ciber_brush, INDICES = ciber_brush$celltype, FUN = function(X) {X[which.max(X$proportion),]})
#by() applies a funcion to subsets of a dataframe, grouped by a factor or variable (in this case cell type)
#celltype is the grouping factor, for each unique celltype, the function will be executed
#function(X) defines a function that is applied to every susbet of data (here, X is all rows corresponding to a specific celltype), those rows become X and the function looks for the max proportion in X
max_per_celltype <- as.data.frame(do.call("rbind", maxRows))

celldecon_stats$y.position <- as.numeric(max_per_celltype[match(celldecon_stats$celltype, max_per_celltype$celltype), "proportion"])
celldecon_stats$y.position <- celldecon_stats$y.position + (celldecon_stats$y.position)*0.05

celldecon_stats$contrast <- paste0(celldecon_stats$group1, "_" ,celldecon_stats$group2)

celldecon_stats <- as.data.frame(celldecon_stats)
celldecon_stats[which(celldecon_stats$contrast == "Control_Severe.COPD"),"y.position"] <- celldecon_stats[which(celldecon_stats$contrast == "Control_Severe.COPD"),"y.position"] +0.02

#Only label the significant comparisons
celldecon_stats2 <- celldecon_stats[which(celldecon_stats$p.adj.signif != "ns"),]


png(file.path(celldecon.figures.dir,"celldecon_ciber_brush.png"), width = 60, height = 30, units = "cm",res = 1200)

ggplot(ciber_brush,
       aes(x=celltype,
           y=proportion))+
  geom_boxplot(aes(fill=classification))+ #have to put fill here, because if we put it in the ggplot() aes, it becomes a global aes and needs 'classification' in each subsequent step. theres no 'classification' column in celldecon_stats so it doesnt work if its global
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5,
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour =  "white"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "grey"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.margin = margin(0,0,1,2, "cm"))+
  
  stat_pvalue_manual(celldecon_stats2,
                     label = "p.adj.signif",
                     tip.length = 0.01,
                     # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                     #  bracket.shorten = 0.1,
                     size = 3)+
  scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                             "Severe.COPD" = "#F8766D")) +
  
  xlab("Cell Type") +
  ylab("Proportion")+
  labs(fill = "Classification")+
  
  # stat_compare_means(aes(group = classification),
  #                    comparisons = my_comparisons,
  #                    method = "wilcox.test")+
  
  ggtitle("Cellular Proportions (Brush)")

dev.off()


cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
#elapsed 14:08mins
#ave cpu 00:13:46
