options(error = function() { traceback(2); quit(status = 1) })
.
# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"


library("ggrepel")
library("readxl")
library(ggfortify)
library(ggplot2)
library(heatmap3)
library(gplots)
library(biomaRt)
library(edgeR)
library(limma)
library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(rstatix)
library(ggpubr)

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir, "data")
processed.data.dir <- file.path(data.dir,"processed")

#Output directory
output.dir <- file.path(main.dir, "output")

sherlock1.dir <- file.path(output.dir, "sherlock1")
if(!exists(sherlock1.dir)) dir.create(sherlock1.dir)


sherlock1.diffexp.dir <- file.path(sherlock1.dir, "diffexp")
if(!exists(sherlock1.diffexp.dir)) dir.create(sherlock1.diffexp.dir)

sherlock1.diffexp.hgnconly.dir <- file.path(sherlock1.dir, "diffexp_hgnc_only")
if(!exists(sherlock1.diffexp.hgnconly.dir)) dir.create(sherlock1.diffexp.hgnconly.dir)


# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))


sherlock1_counts<- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
sherlock1_clinical <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_simple.rds"))

#HGNC FROM ENSEMBL ID -------------------------------------------------------------------------------------------------------
my_ensembl_gene_ids <- rownames(sherlock1_counts)

#library(EnsDb.Hsapiens.v79) #
hgnc_symbols_db <- ensembldb::select(EnsDb.Hsapiens.v79, keys= my_ensembl_gene_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
row.names(hgnc_symbols_db) <- hgnc_symbols_db$GENEID


# ================================================================================== #
# 2. EdgeR DIFFERENTIAL EXPRESSION ANALYSIS ========================================
# ================================================================================== #
cat("Starting 2. EdgeR DIFFERENTIAL EXPRESSION ANALYSIS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

diffexp_edgeR <- function(this.diffexp.dir, showEnsemblID = FALSE) { 
  #defining showEnsemblID and include_ics as FALSE here just makes FALSE the default if there is no input

  
if(all(colnames(sherlock1_counts) == row.names(sherlock1_clinical)) == FALSE){ 
  stop("all(colnames(sherlock1_counts) == row.names(sherlock1_clinical) = FALSE)") }
  
DGEL<- DGEList(counts=sherlock1_counts, group = sherlock1_clinical$classification)

keep <- filterByExpr(DGEL) 


DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.

#NORMALISE
DGEL<- calcNormFactors(DGEL,method = "TMM")


design <- model.matrix(~0 + classification + age + sex + packyears, #all sherlock1 are ex smokers
                        data = sherlock1_clinical) #Design matrix

colnames(design)[1:3] <- c(levels(as.factor(sherlock1_clinical$classification)))

my.contrasts <- makeContrasts(contrast1 = Severe.COPD - Control,
              contrast2 = Mild.moderate.COPD - Control,
              contrast3 = Severe.COPD - Mild.moderate.COPD,
              levels = design)

nameconvert <- as.data.frame(cbind(
  name = c(
    "contrast1",
    "contrast2",
    "contrast3"
    ),
  contrast = c(
    "Severe.COPD - Control",
    "Mild.moderate.COPD - Control",
    "Severe.COPD - Mild.moderate.COPD"
    )
  ))

DGEL <- estimateDisp(DGEL, design) #not robust

fit <- glmQLFit(DGEL, design) #fit the GLM (design) to the DGEL(the DGEL object, which contains the counts data that has been filtered,normalised and dispersons estimtated)

listoftT <- list()
listoftT2 <- list()
listofvolcano <- list()

for (i in colnames(my.contrasts)){
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,i]) #after fitting GLM to counts data (glmQLFit), glmQLFTest perform quasi likelihood f-tests to test for differential expression. ie. hypothesis testing for differential expression. (make inferences or draw conclusions about the data)
  tT <- topTags(qlf,n=nrow(DGEL), adjust.method = "BH")$table #topTags gets top genes, here we want all of the genes, edgeR's default p.adjust method is BH
  tT$Legend <- ifelse(
    tT$FDR < 0.05 & tT$logFC > 1, "Upregulated",
    ifelse(
      tT$FDR < 0.05 & tT$logFC < -1, "Downregulated",
      "Not Significant"))
  
  tT$Legend[is.na(tT$Legend)]="Not significant"
  
  tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols

  if(showEnsemblID == TRUE){
    #for those with no hgnc symbol, label with ensembl id
    tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl
  }

  if(showEnsemblID == FALSE){
    #for those with no hgnc symbol, remove
    tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
  }
  
  
  selection <-which((tT$logFC>1|tT$logFC< -1) & 
                      tT$FDR<0.05)
  
  tT2 <- tT[selection,]
  
  listoftT[[i]] <- tT
  listoftT2[[i]] <- tT2
  
  # ================================================================================== #
  # 2.1. VOLCANO PLOT ================================================================
  # ================================================================================== #
  volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = Legend)) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"))+
    geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
    geom_vline(xintercept =-1,colour="black", linetype="dashed")+
    geom_vline(xintercept =1,colour="black", linetype="dashed")+
    geom_text_repel(data = subset(tT2[1:20,]),
                    aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines") ) +
    theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                     legend.text = element_text(size = 14),
                                     legend.title = element_text(size = 16)) 
  
  listofvolcano[[i]] <- volcano
  
} #close loop (contrast)

#Combine tT, tT2 and volcano plots
listofresults <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)


#Set results paths
diffexp.results.dir <- file.path(this.diffexp.dir, "results")
if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)


diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
if(!exists(diffexp.figures.dir)) dir.create(diffexp.figures.dir, recursive = TRUE)


#Save all results
saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))


### Save volcano plots --------------------------------------------------------------------------------------
cat("Saving VOLCANO PLOTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

volcano.dir <- file.path(diffexp.figures.dir, "volcano")
if(!exists(volcano.dir)) dir.create(volcano.dir, recursive = TRUE)

pdf(file= file.path(volcano.dir, "volcano_plots.pdf"), width = 8, height = 8)

listofresults[["volcano"]][["contrast1"]] + ggtitle(paste("Severe.COPD - Control"))
listofresults[["volcano"]][["contrast2"]] + ggtitle(paste("Mild.moderate.COPD - Control"))
listofresults[["volcano"]][["contrast3"]] + ggtitle(paste("Severe.COPD - Mild.moderate.COPD"))

dev.off()

# ================================================================================== #
# 2.2. SAVE tT AND tT2 RESULTS =====================================================
# ================================================================================== #
cat("Saving tT AND tT2 RESULTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# if(!exists(file.path(diffexp.results.dir,"tT"))) dir.create(file.path(diffexp.results.dir,"tT"), recursive = TRUE)
if(!exists(file.path(diffexp.results.dir,"tT2"))) dir.create(file.path(diffexp.results.dir,"tT2"), recursive = TRUE)

# for contrast1, contrast2 and contrast3
for (i in 1:3){
  # write.csv(listofresults[["tT"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir,"/tT/","tT_", nameconvert$contrast[i], ".csv"))
  write.csv(listofresults[["tT2"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir, "/tT2/","tT2_", nameconvert$contrast[i], ".csv"))
  
}


} #close function

#RUN THE FUNCTION ---------------------------------------------------------------#
#if show ensemblID is set to TRUE, all plots will include the genes that dont have hgnc symbols and save these to a directory folder "diffexp"
#if show ensemblID is set to FALSE, all plots will include only the genes that  have hgnc symbols and save these to a directory folder "diffexp_hgnc_only"

diffexp_edgeR(this.diffexp.dir = sherlock1.diffexp.dir, showEnsemblID  = TRUE)
diffexp_edgeR(this.diffexp.dir = sherlock1.diffexp.hgnconly.dir, showEnsemblID = FALSE)

#--------------------------------------------------------------------------------#




# ================================================================================== #
# 3. BOXPLOTS KEL AND FKBP5  =======================================================
# ================================================================================== #
cat("Starting 3. BOXPLOTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

boxplots.dir <- file.path(sherlock1.diffexp.dir, "figures", "boxplot")
if(!exists(boxplots.dir)) dir.create(boxplots.dir, recursive = TRUE)

# listoftT FROM SHERLOCK1 ANALYSIS (including ICS (so we can compare ICS and non ICS), including ensembl IDs)
sherlock1_listofresults<- readRDS(file.path(sherlock1.diffexp.dir, "results", "listofresults.RDS")) #pvalues from sherlock1 filterbyExpr() results with ICS included
listoftT_sherlock1 <- sherlock1_listofresults[["tT"]]

sherlock1_counts<- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
sherlock1_counts_voom <- voom(sherlock1_counts)

boxplotclinical <- sherlock1_clinical 
boxplotcounts <- t(as.data.frame(sherlock1_counts_voom$E)) #this is voom normalised #sherlock1 counts

#One NA
#Plotting all samples, then seperating by ICS to visualise in the plot
boxplotclinical <- boxplotclinical[!is.na(boxplotclinical$ics_use),]
boxplotcounts <- boxplotcounts[row.names(boxplotclinical),]


if(all(row.names(boxplotcounts) == row.names(boxplotclinical)) == FALSE){ 
  stop("all(row.names(boxplotcounts) == row.names(boxplotclinical)) == FALSE)") }

boxplot <- cbind(boxplotcounts,
                 classification = boxplotclinical$classification,
                 ICS = boxplotclinical$ics_use)


  for (gene in c("KEL", "FKBP5")){
    
    geneofinterest = gene
    geneofinterestid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == geneofinterest), "GENEID"]
    
    plot <- as.data.frame(boxplot[,c(geneofinterestid,
                                     "classification",
                                     "ICS")])

    colnames(plot)[1] <- "gene"
    
    plot$gene <- as.numeric(plot$gene)

    
    
    plot$classification <- factor(plot$classification, levels = c("Control", "Mild.moderate.COPD", "Severe.COPD"))
    
    boxplotimage <- ggplot(plot, aes(
      x = factor(as.factor(classification)),
      y = gene)) +
      
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
                 aes(fill = ICS)
    ) +
      
      labs(title = paste0(geneofinterest, ": SHERLOCK1")) +
      ylab (label = "Expression") +
      xlab (label = "Disease") 
  
    
png(file = file.path(boxplots.dir, paste0(gene,"_test.png")), height = 14, width = 14, units = "cm", res =1200)
print(boxplotimage)
dev.off()

  }



# ================================================================================== #
# 3.1 BOXPLOTS IL33  =======================================================
# ================================================================================== #
cat("Starting 3.1 BOXPLOTS IL33", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

boxplots.dir <- file.path(sherlock1.diffexp.dir, "figures", "boxplot")
if(!exists(boxplots.dir)) dir.create(boxplots.dir, recursive = TRUE)


# listoftT FROM SHERLOCK1 ANALYSIS (including ICS (so we can compare ICS and non ICS), including ensembl IDs)
sherlock1_listofresults<- readRDS(file.path(sherlock1.diffexp.dir, "results", "listofresults.RDS")) #pvalues from sherlock1 filterbyExpr() results with ICS included
listoftT_sherlock1 <- sherlock1_listofresults[["tT"]]

sherlock1_counts<- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
sherlock1_counts_voom <- voom(sherlock1_counts)


boxplotclinical <- sherlock1_clinical 
boxplotcounts <- t(as.data.frame(sherlock1_counts_voom$E)) #this is voom normalised #sherlock1 counts

#One NA
#Plotting all samples, then seperating by ICS to visualise in the plot
boxplotclinical <- boxplotclinical[!is.na(boxplotclinical$ics_use),]
boxplotcounts <- boxplotcounts[row.names(boxplotclinical),]

hgnc_symbols_db <- readRDS(file.path(data.dir,"processed","datawrangling_qc", "hgnc_symbols_db.rds"))


if(all(row.names(boxplotcounts) == row.names(boxplotclinical)) == FALSE){ 
  stop("all(row.names(boxplotcounts) == row.names(boxplotclinical)) == FALSE)") }

boxplot <- cbind(boxplotcounts,
                 classification = boxplotclinical$classification,
                 ICS = boxplotclinical$ics_use)


  geneofinterest = "IL33"

  geneofinterestid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == geneofinterest), "GENEID"]
  
  plot <- as.data.frame(boxplot[,c(geneofinterestid,
                                   "classification",
                                   "ICS")])
  
  colnames(plot)[1] <- "gene"
  
  plot$gene <- as.numeric(plot$gene)
  
  
  
  plot$classification <- factor(plot$classification, levels = c("Control", "Mild.moderate.COPD", "Severe.COPD"))
  
  
  
  
  ## Get P-Values --------------------------------------------------------------------------------------
  cat("Getting BOXPLOTS P-Values for annotation", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  stat.table <- plot %>%
    t_test(gene ~ classification,
           comparisons = list(c("Severe.COPD", "Control"),
                              c("Mild.moderate.COPD", "Control"),
                              c("Severe.COPD", "Mild.moderate.COPD")
           ))
  
  stat.table<- stat.table %>%
    add_xy_position(x = "classification", dodge = 0.8)
  
  
  #### STEP 2) MODIFY STAT TABLE TO INCLUDE COMPARISONS OF INTEREST --------------------------------------------------------------------------------------

  stat.table <- cbind(stat.table, resultsname = c("contrast1", "contrast2", "contrast3"))
  stat.table[which(stat.table$resultsname == "contrast1"),"p"] <- listoftT_sherlock1[["contrast1"]][geneofinterestid, "PValue"]
  stat.table[which(stat.table$resultsname == "contrast2"),"p"] <- listoftT_sherlock1[["contrast2"]][geneofinterestid, "PValue"]
  stat.table[which(stat.table$resultsname == "contrast3"),"p"] <- listoftT_sherlock1[["contrast3"]][geneofinterestid, "PValue"]
  stat.table$p <- signif(as.numeric(stat.table$p), digits = 4)
  stat.table$y.position <- max(plot[,"gene"]) + 0.025*(max(plot[,"gene"]))
  stat.table$y.position <- as.numeric(stat.table$y.position)
  
  # Grouped p value manual --------------------------------------------------------------------------------------------------------------------------
  plot_nocontrol <-plot[-which(plot$classification == "Control"),]
  
  #Compare ICS vs non ICS within disease groups
  stat.table2 <- plot_nocontrol %>%
    group_by(classification) %>%
    t_test(gene ~ ICS )
  
  # Compare disease within smoking groups
  # stat.table2 <- plot %>%
  #   group_by(smokingstatus)  %>%
  #   t_test(gene ~ classification )
  
  stat.table2<- stat.table2 %>%
    add_xy_position(x = "classification", group = "ICS", dodge = 0.8)
  
  stat.table2$p <- signif(as.numeric(stat.table2$p), digits = 4)
  stat.table2$y.position <- as.numeric(stat.table2$y.position)
  new.bracket.distance <- 0.3
  stat.table2$y.position <- max(stat.table$y.position) + new.bracket.distance*2
  
  
  
  
  
  boxplotimage <- ggplot(plot, aes(
    x = factor(as.factor(classification)),
    y = gene)) +
    
    theme_bw()+
    
    theme(plot.caption = element_text(size = 6))+
    
  geom_boxplot(position=position_dodge(),
               aes(fill = ICS)
  ) +
    
    labs(title = paste0(geneofinterest, ": SHERLOCK1"),
         caption = "Counts were voom normalised. 
       For between disease comparison, P-values from differential expression analysis shown.
       For ICS comparisons within disease groups, P-values from T-test shown.") +
    ylab (label = "Expression") +
    xlab (label = "Disease") +

    
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       bracket.nudge.y = c(0, 0.3, 0.3),
                       size = 3) +
    
    
    stat_pvalue_manual(stat.table2,
                       label = "p",
                       tip.length = 0.01,
                       size = 3) 
  
  
  png(file = file.path(boxplots.dir, paste0(geneofinterest,"_boxplot.png")), height = 14, width = 14, units = "cm", res =1200)
  print(boxplotimage)
  dev.off()
  



# ================================================================================== #
# 4. PATIENT DEMOGRAPHICS  ===================================
# ================================================================================== #
cat("Starting 4. PATIENT DEMOGRAPHICS ", Sys.time(), "\n")

qc.dir <- file.path(output.dir, "sherlock1", "qc")
if(!exists(qc.dir)) dir.create(qc.dir, recursive = TRUE)

clinical <- readRDS(file.path(data.dir, "processed", "SHERLOCK1", "clinical_sk1_simple.rds"))


# Patient demographics
#note all patient in clinical are unique
unique_patients_all <- clinical[match(unique(clinical$Study.ID), clinical$Study.ID),]


patient_demographics <- t(unique_patients_all %>% 
                            mutate(packyears = as.numeric(as.character(packyears))) %>% 
                            
                            group_by(classification) %>% 
                            
                            summarise(
                              
                              #Total patients
                              total_patients = n(),
                              
                              #sex
                              male_patients = sum(sex == "Male"), 
                              male_percentage = (male_patients / total_patients) * 100,
                              sex = paste0(male_patients, "(", round(male_percentage, digits=3), ")"),
                              
                              #age
                              mean_age = median(age),
                              range_age =paste0(min(age, na.rm = T), "-", max(age, na.rm = T)),
                              age = paste0(mean_age,"(",range_age, ")"),
                              
                              #Smoking status
                              currentsmoker = sum(smoking_status == "Current.smoker"),
                              smokerpercentage = currentsmoker/total_patients *100,
                              smoke = paste0(currentsmoker, "(", round(smokerpercentage, digits = 3), ")"),

                              #packyears
                              median_packyears = median(packyears, na.rm = TRUE),
                              range_packyears = paste0(min(packyears, na.rm = T), "-", max(packyears, na.rm = T)),
                              packyears = paste0(median_packyears, "(", range_packyears, ")"),
                              
                              
                              #FEV1
                              median_FEV1_percent_pred =  median(FEV1_percent_pred, na.rm = T),
                              range_postfev1percpred = paste0(round(min(FEV1_percent_pred, na.rm = T),digits=3), "-", round(max(FEV1_percent_pred, na.rm = T), digits = 3)),
                              fev1 = paste0(round(median_FEV1_percent_pred, digits = 3), "(", range_postfev1percpred, ")"),
                              
                              #FEV/FVC
                              median_postfev1fvcpercpred =  median(FEV1_FVC_post, na.rm = TRUE),
                              range_postfev1fvcpercpred = paste0(round(min(FEV1_FVC_post, na.rm = T)), "-", round(max(FEV1_FVC_post, na.rm = T), digits = 3 )),
                              fev_fvc = paste0(round(median_postfev1fvcpercpred, digits = 3), "(", range_postfev1fvcpercpred, ")"),
                            )
                          %>%
                            dplyr::select(
                              classification,
                              total_patients,
                              sex,
                              age,
                              smoke,
                              packyears,
                              fev1,
                              fev_fvc
                            )
)

colnames(patient_demographics) <- patient_demographics[1,]
patient_demographics <- patient_demographics[-1,]
patient_demographics <- rbind(as.numeric(table(clinical$classification)),patient_demographics)

row.names(patient_demographics) <- c(
  "Samples, n",
  "Patients, n",
  "Sex Male, n (%)",
  "Age, median (Range)",
  "Current smoker, n (%)",
  "Packyears, median (Range)",
  "FEV1 % pred. (post-bronchodilater), median (Range)",
  "FEV1/FVC % (post-bronchodilater), median (Range)"
)

write.csv(patient_demographics, file = file.path(qc.dir, "patient_demographics.csv"))
saveRDS(patient_demographics, file = file.path(qc.dir, "patient_demographics.rds"))




cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

