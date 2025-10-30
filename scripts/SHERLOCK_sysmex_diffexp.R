#SHERLOCK3 - Analysis of Sysmex data

# The Sysmex variables have been labeled in green.

# Assess 
# 1) Any association between the Sysmex parameters and clinical variables using nonlinear/multivariate analyses (age, sex, BMI, packyears, 
# lung function, emphysema scores, co-morbidities)
# 2) Any genes of which the expression is associated with 1 (or more) of the sysmex variables

# .libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary/")

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
output.dir <- file.path(output.dir, "sysmex")
if(!exists(output.dir))dir.create(output.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))


# ##-- Post batch correction
counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))


#includes sysmex data
# Note that this file has less columns than previous version, survey columns have been removed
raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #319 SEO/patient IDs

#master - all 600+ clinical variables
clinical_brushbiopt_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_brushbiopt_master.rds"))


hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


setwd(file.path(main.dir))




# ================================================================================== #
# 3) Combine clinical file with sysmex data and remove unecessary columns
# ================================================================================== #

#Match to patients in the post-QC clinical master file and counts files

# Sysmex first column = "Mono", last column = "EO-Z", 53 total columns
sysmex_index_start <- which(colnames(raw_clinical) == "Mono")
sysmex_index_end <- which(colnames(raw_clinical) == "EO-Z")

#bind the current clinical file with sysmex data, matching up the sample IDs
clinical_brushbiopt <- cbind(clinical_brushbiopt_master, raw_clinical[match(clinical_brushbiopt_master$Study.ID, 
                                                                            raw_clinical$class_incl_study_id),
                                                                      sysmex_index_start:sysmex_index_end]) #274 samples

colnames(clinical_brushbiopt) <- make.names(colnames(clinical_brushbiopt))

#rename columns
clinical_brushbiopt <- clinical_brushbiopt %>% 
  dplyr::rename(
    sex = crf_gender,
    smoking_status = crf_smoking,
    packyears = crf_packyears,
    corticosteroid = crf_corticosteroid,
    FVC_post= postbodybox_fvc_post ,
    FEV1_FVC_post = postbodybox_fev1_fvc_post
  )



clinical_brush <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),] #274
clinical_biopt <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),] #278




# ================================================================================== #
# 2) Any genes of which the expression is associated with 1 (or more) of the sysmex variables
# 2. DIFFERENTIAL EXPRESSION =======================================================
# ================================================================================== #
cat("Starting 2. DIFFERENTIAL EXPRESSION", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

diffexp.dir <- file.path(output.dir, "diffexp")
if(!exists(diffexp.dir)) dir.create(diffexp.dir, recursive = TRUE)

diffexp.hgnconly.dir <- file.path(output.dir, "diffexp_hgnc_only")
if(!exists(diffexp.hgnconly.dir)) dir.create(diffexp.hgnconly.dir, recursive = TRUE)



## MAKE LISTS -------------------------------------------------------------
#1 brush
#2 biopt
#3 brush vs biopt

#All lists must be in same order for loop 

listofclinical <- list(brush = clinical_brush, #1
                       biopt = clinical_biopt )#2


listofcounts <- list(brush = counts_brush,     #1
                     biopt = counts_biopt)     #2




## DISEASE SEVERITY EdgeR DIFFERENTIAL EXPRESSION -------------------------------------------------------------------------------------
diffexp_edgeR <- function(this.diffexp.dir, showEnsemblID = FALSE) {
  cat(paste("Starting in DIRECTORY", this.diffexp.dir), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  
  
  #Set results paths
  diffexp.results.dir <- file.path(this.diffexp.dir, "results")
  if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)
  
  
  diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir)) dir.create(diffexp.figures.dir, recursive = TRUE)
  
  
  
  listofresults <- list()
  
  for (j in names(listofclinical)){
    cat(paste("Starting SAMPLE TYPE", j), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    clinical2 <- listofclinical[[j]]
    counts2 <- listofcounts[[j]]
    
    
    
    
    listoftT <- list()
    listoftT2 <- list()
    listofvolcano <- list()
    
    
    for (this_sysmex_variable in colnames(clinical2[,619:ncol(clinical2)])){
      cat(paste("Starting SAMPLE TYPE", j, this_sysmex_variable), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      
      #DESIGN MATRIX
      # we cant correct for classification because the sysmex variable could be associated with the classification
      formula <- as.formula(paste("~ 0 +", this_sysmex_variable, " + age + sex + smoking_status"))
      
      
      design <- model.matrix(formula, #Removing the intercept by using "0 +" means “I want to estimate the expression level for each group independently, and I’ll decide later how to compare them.” instead of edgeR automatically taking 'COntrol' as the reference value since it's the first level in the factor classification (ie. control becomes the intercept)
                             data = clinical2) #Design matrix
      
      # Account for the fact that some patients are missing sysmex data
      clinical2 <- clinical2[intersect(row.names(clinical2), row.names(design)), ] 
      counts2 <- counts2[,row.names(clinical2)]
      
      # DGEList is a list-based data object. It has a matrix 'counts', a data.frame 'samples' (has info about the sample) with a column "lib.size" for library size or sequency depth, 
      DGEL<- DGEList(counts=counts2) 
      
      # FILTER
      keep <- filterByExpr(DGEL) 
      # The filterBy Expr function keeps rows that have worthwhile counts in a minimum number of samples.
      # keep <- which(rowMedians(as.matrix(DGEL))>10) #leaves 36,415 genes
      # keep <- rowSums(cpm(expression)>100) >= 2
      
      
      DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
      # expression <- DGEL$counts
      
      # NORMALISE
      DGEL<- calcNormFactors(DGEL,method = "TMM")
      
      # ESTIMATE DISPERSON
      DGEL <- estimateDisp(DGEL, design)
      
      # FIT MODEL
      fit <- glmQLFit(DGEL, design) #fit the GLM (design) to the DGEL(the DGEL object, which contains the counts data that has been filtered,normalised and dispersons estimtated)
      
      
      # coef = sysmex variable of interest
      qlf <- glmQLFTest(fit, coef=colnames(design)[1]) #after fitting GLM to counts data (glmQLFit), glmQLFTest perform quasi likelihood f-tests to test for differential expression. ie. hypothesis testing for differential expression. (make inferences or draw conclusions about the data)
      
      
      tT <- topTags(qlf,n=nrow(DGEL))$table #topTags gets top genes, here we want all of the genes, edgeR's default p.adjust method is BH
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
      
      else{
        #for those with no hgnc symbol, remove
        tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
      }
      
      selection <-which((tT$logFC>1|tT$logFC< -1)&tT$FDR<0.05)
      
      tT2 <- tT[selection,]
      
      
      listoftT[[this_sysmex_variable]] <- tT
      listoftT2[[this_sysmex_variable]] <- tT2
      
      
      # ================================================================================== #
      # 2.1. VOLCANO PLOT ================================================================
      # ================================================================================== #
      cat("Starting 2.1. VOLCANO PLOT", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      
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
      
      ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(j,"_", this_sysmex_variable,"_volcano_plot.png")), 
             width = 25, height = 25, 
             units = "cm")
      
      
      listofvolcano[[this_sysmex_variable]] <- volcano
      
      
      
    } #close loop for sysmex variable
    
    listofresults[[j]] <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)
    
    #Save all results
    saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))
    
    
    
    
    # 
    # ### Save volcano plots --------------------------------------------------------------------------------------
    # cat("Saving VOLCANO PLOTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    # 
    # volcano.dir <- file.path(diffexp.figures.dir, "volcano")
    # if(!exists(volcano.dir)) dir.create(volcano.dir, recursive = TRUE)
    # 
    # #NOTE!!!!Below code doesn't work in hpc when run as a job, but works in interactive if listofresults is loaded manually with readRDS and the code is run???
    # pdf(file= file.path(volcano.dir, "volcano_plots_hpc.pdf"), width = 8, height = 8)
    # 
    # #### Biopt Control vs Brush COntrol --------------------------------------------------------------------------------------
    # listofresults[["brushbiopt"]][["volcano"]][["contrast7"]] + ggtitle(paste("Biopt_Control - Brush_Control"))
    # 
    # #### Severe COPD - Control ---------------------------------------------------------------------------------------
    # listofresults[["brush"]][["volcano"]][["contrast1"]] + ggtitle(paste("Brush:","Severe COPD - Control"))
    # listofresults[["biopt"]][["volcano"]][["contrast1"]] + ggtitle(paste("Biopt:","Severe COPD - Control"))
    # 
    # #---Venn diagram to compare biopt and brush#
    # x <- list(Brush = row.names(listofresults[["brush"]][["tT2"]][["contrast1"]]), 
    #           Biopsy = row.names(listofresults[["biopt"]][["tT2"]][["contrast1"]])
    # )
    # 
    # (ggvenn(x, stroke_size = 1.5) + ggtitle("Severe COPD - Control"))
    # 
    # #----Delta DE to compare biopt and brush
    # listofresults[["brushbiopt"]][["volcano"]][["contrast4"]] + 
    #   ggtitle(paste("(Biopt_Severe.COPD - Brush_Severe.COPD) - \n (Biopt_Control - Brush_Control)")) + 
    #   theme(plot.title = element_text(size = 15)) 
    # 
    # #### Mild.moderate.COPD - Control --------------------------------------------------------------------------------------
    # listofresults[["brush"]][["volcano"]][["contrast2"]] + ggtitle(paste("Brush:","Mild.moderate.COPD - Control"))
    # listofresults[["biopt"]][["volcano"]][["contrast2"]] + ggtitle(paste("Biopt:","Mild.moderate.COPD - Control"))
    # 
    # #---Venn diagram to compare biopt and brush#
    # x <- list(Brush = row.names(listofresults[["brush"]][["tT2"]][["contrast2"]]), 
    #           Biopsy = row.names(listofresults[["biopt"]][["tT2"]][["contrast2"]]))
    # 
    # (ggvenn(x, stroke_size = 1.5) + ggtitle("Mild.moderate.COPD - Control"))
    # 
    # #---Delta DE to compare biopt and nasal
    # listofresults[["brushbiopt"]][["volcano"]][["contrast5"]] + 
    #   ggtitle(paste("(Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD) - \n (Biopt_Control- Brush_Control)")) + 
    #   theme(plot.title = element_text(size = 15)) 
    # 
    # #### Severe.COPD - Mild.moderate.COPD --------------------------------------------------------------------------------------
    # listofresults[["brush"]][["volcano"]][["contrast3"]] + ggtitle(paste("Brush:","Severe.COPD - Mild.moderate.COPD"))
    # listofresults[["biopt"]][["volcano"]][["contrast3"]] + ggtitle(paste("Biopt:","Severe.COPD - Mild.moderate.COPD"))
    # 
    # #---Venn diagram to compare biopt and nasal#
    # x <- list(Brush = row.names(listofresults[["brush"]][["tT2"]][["contrast3"]]), 
    #           Biopsy = row.names(listofresults[["biopt"]][["tT2"]][["contrast3"]]))
    # 
    # (ggvenn(x, stroke_size = 1.5) + ggtitle("Severe.COPD - Mild.moderate.COPD"))
    # 
    # #---Delta DE to compare biopt and nasal
    # listofresults[["brushbiopt"]][["volcano"]][["contrast6"]] + 
    #   ggtitle(paste("(Biopt_Severe.COPD - Brush_Severe.COPD) - \n (Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD)")) + 
    #   theme(plot.title = element_text(size = 15)) 
    # 
    # 
    # 
    # 
    # dev.off()
    
  } #close for loop (brush and biopt)
  
  
} #close function

diffexp_edgeR(this.diffexp.dir = diffexp.dir, showEnsemblID  = TRUE)
diffexp_edgeR(this.diffexp.dir = diffexp.hgnconly.dir, showEnsemblID = FALSE)

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

