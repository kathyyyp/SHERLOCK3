#SHERLOCK3 - Analysis of Sysmex data

# The Sysmex variables have been labeled in green in "Sherlock_database_07_25_Final.xlsx"

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

library(parallel) 
NCORES <- 8 
system.time(mclapply(1:NCORES, function(x) Sys.sleep(2), mc.cores=NCORES))

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
                                                                      sysmex_index_start:sysmex_index_end]) 

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

clinical_brushbiopt$age <- as.numeric(clinical_brushbiopt$age)

clinical_brushbiopt[which(clinical_brushbiopt$RELYMP.103uL == "----"),"RELYMP.103uL"] <- NA
clinical_brushbiopt$RELYMP.103uL <- as.numeric(clinical_brushbiopt$RELYMP.103uL)

clinical_brushbiopt[which(clinical_brushbiopt$RELYMP == "----"),"RELYMP"] <- NA
clinical_brushbiopt$RELYMP <- as.numeric(clinical_brushbiopt$RELYMP)

sapply(clinical_brushbiopt[619:ncol(clinical_brushbiopt)],class)


clinical_brush <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),] #270
clinical_biopt <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),] #271




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
    
    
    sysmex_vars <- colnames(clinical2[, 619:ncol(clinical2)])
    
    parallel_results <- mclapply(sysmex_vars, function(this_sysmex_variable) {
      
    # for (this_sysmex_variable in colnames(clinical2[,619:ncol(clinical2)])){
      cat(paste("Starting SAMPLE TYPE", j, this_sysmex_variable), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      
      
      #DESIGN MATRIX
      # we cant correct for classification because the sysmex variable could be associated with the classification
      formula <- as.formula(paste("~ 0 +", this_sysmex_variable, " + age + sex + smoking_status"))
      
      
      design <- model.matrix(formula, #Removing the intercept by using "0 +" means “I want to estimate the expression level for each group independently, and I’ll decide later how to compare them.” instead of edgeR automatically taking 'COntrol' as the reference value since it's the first level in the factor classification (ie. control becomes the intercept)
                             data = clinical2) #Design matrix
      
      
      # Account for the fact that some patients are missing sysmex data
      clinical2 <- clinical2[intersect(row.names(clinical2), row.names(design)), ] 
      counts2 <- counts2[,row.names(clinical2)]
      
      zero_proportion <- sum(clinical2[,this_sysmex_variable] == 0, na.rm = TRUE)/length(clinical2[,this_sysmex_variable])
      
      if(zero_proportion > 0.8){
        message("Skipping ", this_sysmex_variable, " (", round(zero_proportion * 100, 1), "% zeros)")
        return(list(tT = NULL, tT2 = NULL, volcano = NULL))
        
      }
      
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
        tT$FDR < 0.05 & tT$logFC > 0, "Upregulated",
        ifelse(
          tT$FDR < 0.05 & tT$logFC < 0, "Downregulated",
          "Not Significant"))
      
      tT$Legend[is.na(tT$Legend)]="Not Significant"
      
      tT$Legend <- factor(tT$Legend, levels = c("Downregulated", "Upregulated", "Not Significant"))
      
      tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols
      
      if(showEnsemblID == TRUE){
        #for those with no hgnc symbol, label with ensembl id
        tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl
      }
      
      else{
        #for those with no hgnc symbol, remove
        tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
      }
      
      selection <-which(tT$FDR<0.05)
      
      tT2 <- tT[selection,]
      
      
      # listoftT[[this_sysmex_variable]] <- tT
      # listoftT2[[this_sysmex_variable]] <- tT2
      
      
      # ================================================================================== #
      # 2.1. VOLCANO PLOT ================================================================
      # ================================================================================== #
      cat("Starting 2.1. VOLCANO PLOT", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      
      volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
        ggtitle(paste(j, "-", this_sysmex_variable)) +
        geom_point(aes(color = Legend)) +
        scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
        geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
        geom_text_repel(data = subset(tT2[1:30,]),
                        aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines") ) +
        theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                         legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 16)) 
      
      # ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(j,"_", this_sysmex_variable,"_volcano_plot.png")),
      #        width = 25, height = 25,
      #        units = "cm")
      
      

      list(
        tT = tT,
        tT2 = tT2,
        volcano = volcano
      )
      
    }, mc.cores = NCORES)
    
    
    # assign names back to results
    names(parallel_results) <- sysmex_vars
    
    # unpack into your original structure
    listoftT      <- lapply(parallel_results, `[[`, "tT")
    listoftT2     <- lapply(parallel_results, `[[`, "tT2")
    listofvolcano <- lapply(parallel_results, `[[`, "volcano")

    listofresults[[j]] <- list(
      tT      = listoftT,
      tT2     = listoftT2,
      volcano = listofvolcano
    )
    
    
    #Save all results
    saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))
    
    
    
    
    

    

    
    
    ## Grid plo t for volcano plots - 9 plots per page to save space
    
    library(patchwork)
    
    # Suppose volcano.list is your list of 40 ggplots
    total_plots <- length(listofvolcano)
    plots_per_page <- 4
    total_pages <- ceiling(total_plots / plots_per_page)
    
    
    pdf(file.path(diffexp.figures.dir,paste0(j, "_volcano_plot_panels.pdf")), width = 10, height = 10)
    
    for (i in seq_len(total_pages)) {
      start <- (i - 1) * plots_per_page + 1
      end <- min(i * plots_per_page, total_plots)
      
      page_plots <- listofvolcano[start:end]
      
      plotting_theme_func <- function(p) {
        p + theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)
        )
        
      }
      
      #map = apply a function to each element of a vector, in this case its our listofvolcanoplots
      page_plots <- Map(plotting_theme_func, page_plots)
      
      
      # Combine current 9 (or fewer) plots in a 3x3 grid
      grid_plot <- wrap_plots(page_plots, ncol = 2, nrow = 2, guides = "collect") +
        
        plot_layout(guides = "collect") +
        
        plot_annotation(
          title = paste("Sysmex variables & Gene expression - ", j),
          theme = theme(plot.title = element_text(size = 16, face = "bold"))
        ) &
        theme(legend.position = "bottom")  # or "right", "top", "left"
      
      print(grid_plot)
    }      
    
    dev.off()
    
    
  } #close for loop (brush and biopt)
  
  
} #close function

diffexp_edgeR(this.diffexp.dir = diffexp.dir, showEnsemblID  = TRUE)
# diffexp_edgeR(this.diffexp.dir = diffexp.hgnconly.dir, showEnsemblID = FALSE)




cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#Note the sysmex variables are not scaled, therefore logFC changes may not reflect biological changes
#Get all unique assoaicted genes
relymp_indexes <- which(names(listofresults[["brush"]][["tT2"]]) %in% c("RELYMP", "RELYMP.103uL"))
all_genes <- unlist(lapply(listofresults[["brush"]][["tT2"]][-c(relymp_indexes)], function(x) x$gene_symbol))

# Get unique ones
unique_genes <- unique(all_genes)

# How many unique genes? 4231 for brush
length(unique_genes)

# Create a table with genes as rows and diff exp results as columns, with the matching genes

