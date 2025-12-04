#SHERLOCK analysis of radiomics varaibles in relation to gene expression
# Transcriptional changes associated with emphysema severity  in COPD patients, also comparing upper and lower lobe
# Linking radiomics feature to transcriptomic signatures in COPD. Emphysema and mucus plugging scores from CT scans


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
library("tidyverse")


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
output.dir <- file.path(output.dir, "radiomics")
if(!exists(output.dir))dir.create(output.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(main.dir))

# ##-- Post batch correction
counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))


# Note that this file has less columns than previous version, survey columns have been removed
# IMPORTANT - radiomics data in this file is updated too. the clinical_brushbiopt_master file has weird emphysema % values 

raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #319 SEO/patient IDs

#master - all 600+ clinical variables
clinical_brushbiopt_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_brushbiopt_master.rds")) 


clinical_brushbiopt <- as.data.frame(raw_clinical[match(clinical_brushbiopt_master$Study.ID, raw_clinical$class_incl_study_id),])
clinical_brushbiopt <- cbind(clinical_brushbiopt, sampletype = clinical_brushbiopt_master$sampletype, 
                             batch = clinical_brushbiopt_master$batch,
                             classification = clinical_brushbiopt_master$classification)
row.names(clinical_brushbiopt) <- row.names(clinical_brushbiopt_master)

clinical_brushbiopt <- clinical_brushbiopt %>% 
  dplyr::rename(
    age = crf_age,
    sex = crf_gender,
    smoking_status = crf_smoking,
    packyears = crf_packyears,
    corticosteroid = crf_corticosteroid,
    FVC_post= postbodybox_fvc_post ,
    FEV1 = postbodybox_fev1_post,
    FEV1_FVC_post = postbodybox_fev1_fvc_post
  ) 


clinical_brushbiopt[,c("age", "packyears", "FEV1", 
                              "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(
                                clinical_brushbiopt[,c("age", "packyears", "FEV1", 
                                                              "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], 
                                function(x) as.numeric(x))

  


startcol <-which(colnames(clinical_brushbiopt) =="RL_insp_vol_ml")
endcol <- which(colnames(clinical_brushbiopt) =="LLL_airtrapping_emphysema_%")


#Make the names valid for R
colnames(clinical_brushbiopt)[startcol:endcol] <- gsub("%", "perc", colnames(clinical_brushbiopt)[startcol:endcol] )
colnames(clinical_brushbiopt)[startcol:endcol] <- gsub(">", "over", colnames(clinical_brushbiopt)[startcol:endcol] )
colnames(clinical_brushbiopt)[startcol:endcol] <- gsub("-", ".", colnames(clinical_brushbiopt)[startcol:endcol] )
clinical_brushbiopt[,startcol:endcol] <- sapply(clinical_brushbiopt[,startcol:endcol], 
                                                       function(x) as.numeric(x))



hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


setwd(file.path(main.dir))


clinical_brush <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),] #270
clinical_biopt <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),] #271



# Emphysema variables


# NOTE: THE BELOW ARE IDENTICAL FOR INSP AND EXP
# for emhpysema - use insp (standard bc lung is fully inflated - less variability due to volume)
# for airtrapping & small airway disease - use exp bc more relevant to small airway obstruction

#insp volume
# [480] "RL_insp_vol_ml"
# [481] "LL_insp_vol_ml"
# [482] "RUL_insp_vol_ml"
# [483] "RML_insp_vol_ml"
# [484] "RLL_insp_vol_ml"
# [485] "LUL_insp_vol_ml"
# [486] "LLL_insp_vol_ml"

# %LAAA-950HU (same as ES950?)
# LAA - (low attenuation area)
# insp = inspiratory phase
# HU - hounsfield unit
# RL_insp_LAA_.950HU - % of the right lung with density < −950 HU on inspiratory CT 
# quantifies low-density lung regions on inspiratory CT -> directly captures alveolar destruction
# [494] "RL_insp_LAA_.950HU_perc"
# [495] "LL_insp_LAA_.950HU_perc"
# [496] "RUL_insp_LAA_.950HU_perc"
# [497] "RML_insp_LAA_.950HU_perc"
# [498] "RLL_insp_LAA_.950HU_perc"
# [499] "LUL_insp_LAA_.950HU_perc"
# [500] "LLL_insp_LAA_.950HU_perc"


#mean lung density
# [501] "RL_insp_mean_density_HU"
# [502] "LL_insp_mean_density_HU"
# [503] "RUL_insp_mean_density_HU"
# [504] "RML_insp_mean_density_HU"
# [505] "RLL_insp_mean_density_HU"
# [506] "LUL_insp_mean_density_HU"
# [507] "LLL_insp_mean_density_HU"


# 15th percentile HU (PD15) - Hounsfield Unit value at the 15th percentile of inspiratory lung density; lower (more negative) values indicate more severe emphysema due to alveolar destruction.
# - sort all lung voxels from lowest to highest density -> the value at which 15% of voxels are below = PD15 - > PD15 is the density below which the most emphysematous 15% of the lung resides (more negative pd15 = more emphysema)
# more lung destroyed = pd15 is higher 
# Eg. severe copd could have 40%LAA-950 (meaning 40% of the lung has density less than -950HU = lots of low density areas = severe disease) and PD15 of -960 HU (which means the most emphysematous 15% of the lung has HU under/more negative than -960, ie. the worst regions of the lung have very low density)
# mild disease eg 5% LAA-950 and -910HU (emphysema affecting 5% of the lung and the worst portions don't have extremely low density)
# [508] "RL_insp_emphysema_15percentile_HU"
# [509] "LL_insp_emphysema_15percentile_HU"
# [510] "RUL_insp_emphysema_15percentile_HU"
# [511] "RML_insp_emphysema_15percentile_HU"
# [512] "RLL_insp_emphysema_15percentile_HU"
# [513] "LUL_insp_emphysema_15percentile_HU"
# [514] "LLL_insp_emphysema_15percentile_HU"



# PI10: Perimeter-Interpolated wall thickness at 10 mm” ie. airway wall thickness (10mM)
# - Pi10 is a way to normalise thickness bc bigger airways naturally have thicker walls & smaller airways have thinner walls
# Measures airway wall thickening (airway wall remodelling) 
# [515] "Lungs_Pi10"
# [516] "RUL_Pi10"
# [517] "RML_Pi10"
# [518] "RLL_Pi10"
# [519] "LUL_Pi10"
# [520] "LLL_Pi10"


# Bronchial count = number of visible airways (bronchi) in a lung scan
# [521] "Lungs_insp_bronchial_count"
# [522] "RUL_insp_bronchial_count"
# [523] "RML_insp_bronchial_count"
# [524] "RLL_insp_bronchial_count"
# [525] "LUL_insp_bronchial_count"
# [526] "LLL_insp_bronchial_count"


# Veins and vessels percentage
# [531] "lungs_VVP_vessels_0.1mm"
# [532] "lungs_VVP_vessels_1.2mm"
# [533] "lungs_VVP_vessels_over2mm"
# [534] "lungs_VVP_arteries_0.1mm"
# [535] "lungs_VVP_arteries_1.2mm"
# [536] "lungs_VVP_arteries_over2mm"
# [537] "lungs_VVP_veins_0.1mm"
# [538] "lungs_VVP_veins_1.2mm"
# [539] "lungs_VVP_veins_over2mm"



# mucus plugs
# [540] "Lungs_mucus_vol_ml"
# [541] "Lungs_mucus_plugs_count"
# [542] "RUL_mucus_plugs_count"
# [543] "RML_mucus_plugs_count"
# [544] "RLL_mucus_plugs_count"
# [545] "LUL_mucus_plugs_count"
# [546] "LLL_mucus_plugs_count"

# % lung normal
# [575] "Lungs_normal_perc"
# [576] "RL_normal_perc"
# [577] "LL_normal_perc"
# [578] "RUL_normal_perc"
# [579] "RML_normal_perc"
# [580] "RLL_normal_perc"
# [581] "LUL_normal_perc"
# [582] "LLL_normal_perc"

# % emphysema
# [583] "Lungs_emphysema_perc"
# [584] "RL_emphysema_perc"
# [585] "LL_emphysema_perc"
# [586] "RUL_emphysema_perc"
# [587] "RML_emphysema_perc"
# [588] "RLL_emphysema_perc"
# [589] "LUL_emphysema_perc"
# [590] "LLL_emphysema_perc"
# 
# [591] "Lungs_airtrapping_perc"
# [592] "RL_airtrapping_perc"
# [593] "LL_airtrapping_perc"
# [594] "RUL_airtrapping_perc"
# [595] "RML_airtrapping_perc"
# [596] "RLL_airtrapping_perc"
# [597] "LUL_airtrapping_perc"
# [598] "LLL_airtrapping_perc"

# [599] "Lungs_airtrapping_emphysema_perc"
# [600] "RL_airtrapping_emphysema_perc"
# [601] "LL_airtrapping_emphysema_perc"
# [602] "RUL_airtrapping_emphysema_perc"
# [603] "RML_airtrapping_emphysema_perc"
# [604] "RLL_airtrapping_emphysema_perc"
# [605] "LUL_airtrapping_emphysema_perc"
# [606] "LLL_airtrapping_emphysema_perc"


# Comparisons to do

# Associate emphysema severity with gene expression -----------------------------------------------
# A) Whole-lung emphysema 
    # - "Lungs_Pi10", "Lungs_insp_bronchial_count", "Lungs_mucus_plugs_count", "Lungs_emphysema_perc", "Lungs_airtrapping_perc", "Lungs_airtrapping_emphysema_perc"
    # - "RL_insp_LAA_.950HU_perc", "LL_insp_emphysema_15percentile_HU", "RL_insp_mean_density_HU", "LL_insp_mean_density_HU", "RL_emphysema_perc", "LL_emphysema_perc",  "RL_airtrapping_perc",  "LL_airtrapping_perc", "RL_airtrapping_emphysema_perc", "LL_airtrapping_emphysema_perc"
  

# Compare upper vs lower lobes???? .-------------------------------------------------------------------------------------------
# B) Right: Upper vs middle vs lower - block for patients?
    # - "RUL_insp_LAA_.950HU_perc", "RML_insp_LAA_.950HU_perc", "RLL_insp_LAA_.950HU_perc"
    # - "RUL_insp_emphysema_15percentile_HU", "RML_insp_emphysema_15percentile_HU", "RLL_insp_emphysema_15percentile_HU"
    # - RUL_Pi10", "RML_Pi10", "RLL_Pi10"
    # mean density,  insp vol, bronchial count, mucus plugging, airtrapping, %emphysema, VVP, % airtrapping

# C) Left: Upper vs lower - block for patients id
    # - "LUL_insp_LAA_.950HU_perc", "LLL_insp_LAA_.950HU_perc"
    # - "LUL_insp_emphysema_15percentile_HU", "LLL_insp_emphysema_15percentile_HU"
    # - "LUL_Pi10", "LLL_Pi10"


# ================================================================================== #
# 2. DIFFERENTIAL EXPRESSION =======================================================
# ================================================================================== #
cat("Starting 2. DIFFERENTIAL EXPRESSION", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

diffexp.dir <- file.path(output.dir, "diffexp")
if(!exists(diffexp.dir)) dir.create(diffexp.dir, recursive = TRUE)

# diffexp.hgnconly.dir <- file.path(output.dir, "diffexp_hgnc_only")
# if(!exists(diffexp.hgnconly.dir)) dir.create(diffexp.hgnconly.dir, recursive = TRUE)



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
    
    

    #set up variables of interest (for emphysema, insp) - do air trapping exp seperately?
    global_lung_vars <- c("RL_insp_LAA_.950HU_perc",
                          "LL_insp_LAA_.950HU_perc",
                          "RL_insp_emphysema_15percentile_HU",
                          "LL_insp_emphysema_15percentile_HU",
                          "Lungs_Pi10", 
                          "Lungs_insp_bronchial_count", 
                          "Lungs_mucus_plugs_count", 
                          "Lungs_emphysema_perc", 
                          "Lungs_airtrapping_perc")

    
    # right_lung_vars <- c("RUL_insp_LAA_.950HU_perc", "RML_insp_LAA_.950HU_perc", "RLL_insp_LAA_.950HU_perc",
    #                      "RUL_insp_emphysema_15percentile_HU", "RML_insp_emphysema_15percentile_HU", "RLL_insp_emphysema_15percentile_HU",
    #                      "RUL_Pi10", "RML_Pi10", "RLL_Pi10")
    # 
    # 
    # left_lung_vars <-c("LUL_insp_LAA_.950HU_perc", "LLL_insp_LAA_.950HU_perc",
    #                    "LUL_insp_emphysema_15percentile_HU", "LLL_insp_emphysema_15percentile_HU",
    #                    "LUL_Pi10", "LLL_Pi10")

    
    # =============================
    
    
    for (this_variable in global_lung_vars){
      cat(paste("Starting SAMPLE TYPE", j, this_variable), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      
      
      # #Remove the NAs -  outliers
      # clinical2 <- clinical2[!is.na(clinical2[,this_variable]),]
      # counts2 <- counts2[,row.names(clinical2)]
      # 
      # 
      # #Filter out variables where there are too many zeros (cutoff at >80% zeros)
      # zero_proportion <- sum(clinical2[,this_variable] == 0, na.rm = TRUE)/length(clinical2[,this_variable])
      # 
      # if(zero_proportion > 0.8){
      #   cat(paste0("Skipping ", this_variable, " (", round(zero_proportion * 100, 1), "% zeros)"))
      #   next        
      # }
      # 
      # 
      # # Filter out outliers (3sd from the mean) - these are falsely driving the diffexp results
      # stdev  <- sd(clinical2[,this_variable], na.rm = TRUE)
      # 
      # #if value is 3 standard deviations away from the mean, make it NA
      # if (!is.na(stdev) && stdev > 0) {
      #   mean_val <- mean(clinical2[,this_variable], na.rm = TRUE)
      #   outliers <- abs(clinical2[,this_variable] - mean_val) > 3 * stdev
      #   clinical2[outliers, this_variable] <- NA
      # }
      # 
      
       
      #Remove the NAs - some are missing and some are outliers
      clinical2 <- clinical2[!is.na(clinical2[,this_variable]),]
      counts2 <- counts2[,row.names(clinical2)]
      
      
      
      #DESIGN MATRIX
      # we cant correct for classification because the emphysema variable could be associated with the classification
      formula <- as.formula(paste("~ 0 +", this_variable, " + age + sex + smoking_status"))
      
      
      design <- model.matrix(formula, #Removing the intercept by using "0 +" means “I want to estimate the expression level for each group independently, and I’ll decide later how to compare them.” instead of edgeR automatically taking 'COntrol' as the reference value since it's the first level in the factor classification (ie. control becomes the intercept)
                             data = clinical2) #Design matrix
      
      
      
      # DGEList is a list-based data object. It has a matrix 'counts', a data.frame 'samples' (has info about the sample) with a column "lib.size" for library size or sequency depth, 
      DGEL<- DGEList(counts=counts2, group = clinical2$classification) 
      
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
      
      
      # coef = variable of interest
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
      
      
      listoftT[[this_variable]] <- tT
      listoftT2[[this_variable]] <- tT2
      
      
      # ================================================================================== #
      # 2.1. VOLCANO PLOT ================================================================
      # ================================================================================== #
      cat("Starting 2.1. VOLCANO PLOT", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      
      volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
        ggtitle(paste(j, "-", this_variable)) +
        geom_point(aes(color = Legend)) +
        scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
        geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
        geom_text_repel(data = subset(tT2[1:30,]),
                        aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines") ) +
        theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                         legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 16)) 
      
      ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(j,"_", this_variable,"_volcano_plot.png")),
             width = 25, height = 25,
             units = "cm")
      
      
      listofvolcano[[this_variable]] <- volcano
      
      
    } #close loop for variable
    
    
    
    listofresults[[j]] <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)
    
    
    
    
    
    
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
          title = paste("Emphysema variables & Gene expression - ", j),
          theme = theme(plot.title = element_text(size = 16, face = "bold"))
        ) &
        theme(legend.position = "bottom")  # or "right", "top", "left"
      
      print(grid_plot)
    }      
    
    dev.off()
    
    
  } #close for loop (brush and biopt)
  #Save all results
  
  
  
  saveRDS(listofresults, file = file.path(diffexp.results.dir, "listofresults.RDS"))
  
  
} #close function

diffexp_edgeR(this.diffexp.dir = diffexp.dir, showEnsemblID  = TRUE)
# diffexp_edgeR(this.diffexp.dir = diffexp.hgnconly.dir, showEnsemblID = FALSE)




cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


# [480] "RL_insp_vol_ml"
# [481] "LL_insp_vol_ml"
# [482] "RUL_insp_vol_ml"
# [483] "RML_insp_vol_ml"
# [484] "RLL_insp_vol_ml"
# [485] "LUL_insp_vol_ml"
# [486] "LLL_insp_vol_ml"
# [487] "RL_insp_LAA_.910HU_."
# [488] "LL_insp_LAA_.910HU_."
# [489] "RUL_insp_LAA_.910HU_."
# [490] "RML_insp_LAA_.910HU_."
# [491] "RLL_insp_LAA_.910HU_."
# [492] "LUL_insp_LAA_.910HU_."
# [493] "LLL_insp_LAA_.910HU_."
# [494] "RL_insp_LAA_.950HU_."
# [495] "LL_insp_LAA_.950HU_."
# [496] "RUL_insp_LAA_.950HU_."
# [497] "RML_insp_LAA_.950HU_."
# [498] "RLL_insp_LAA_.950HU_."
# [499] "LUL_insp_LAA_.950HU_."
# [500] "LLL_insp_LAA_.950HU_."
# [501] "RL_insp_mean_density_HU"
# [502] "LL_insp_mean_density_HU"
# [503] "RUL_insp_mean_density_HU"
# [504] "RML_insp_mean_density_HU"
# [505] "RLL_insp_mean_density_HU"
# [506] "LUL_insp_mean_density_HU"
# [507] "LLL_insp_mean_density_HU"
# [508] "RL_insp_emphysema_15percentile_HU"
# [509] "LL_insp_emphysema_15percentile_HU"
# [510] "RUL_insp_emphysema_15percentile_HU"
# [511] "RML_insp_emphysema_15percentile_HU"
# [512] "RLL_insp_emphysema_15percentile_HU"
# [513] "LUL_insp_emphysema_15percentile_HU"
# [514] "LLL_insp_emphysema_15percentile_HU"
# [515] "Lungs_Pi10"
# [516] "RUL_Pi10"
# [517] "RML_Pi10"
# [518] "RLL_Pi10"
# [519] "LUL_Pi10"
# [520] "LLL_Pi10"
# [521] "Lungs_insp_bronchial_count"
# [522] "RUL_insp_bronchial_count"
# [523] "RML_insp_bronchial_count"
# [524] "RLL_insp_bronchial_count"
# [525] "LUL_insp_bronchial_count"
# [526] "LLL_insp_bronchial_count"
# [527] "Lungs_G1.G6_BwaBoa"
# [528] "Lungs_G1.G6_BoutA"
# [529] "Lungs_G1.G6_BinA"
# [530] "Lungs_G1.G6_BwtA"
# [531] "lungs_VVP_vessels_0.1mm"
# [532] "lungs_VVP_vessels_1.2mm"
# [533] "lungs_VVP_vessels_.2mm"
# [534] "lungs_VVP_arteries_0.1mm"
# [535] "lungs_VVP_arteries_1.2mm"
# [536] "lungs_VVP_arteries_.2mm"
# [537] "lungs_VVP_veins_0.1mm"
# [538] "lungs_VVP_veins_1.2mm"
# [539] "lungs_VVP_veins_.2mm"
# [540] "Lungs_mucus_vol_ml"
# [541] "Lungs_mucus_plugs_count"
# [542] "RUL_mucus_plugs_count"
# [543] "RML_mucus_plugs_count"
# [544] "RLL_mucus_plugs_count"
# [545] "LUL_mucus_plugs_count"
# [546] "LLL_mucus_plugs_count"
# [547] "RL_exp_vol_ml"
# [548] "LL_exp_vol_ml"
# [549] "RUL_exp_vol_ml"
# [550] "RML_exp_vol_ml"
# [551] "RLL_exp_vol_ml"
# [552] "LUL_exp_vol_ml"
# [553] "LLL_exp_vol_ml"
# [554] "RL_exp_mean_density_HU"
# [555] "LL_exp_mean_density_HU"
# [556] "RUL_exp_mean_density_HU"
# [557] "RML_exp_mean_density_HU"
# [558] "RLL_exp_mean_density_HU"
# [559] "LUL_exp_mean_density_HU"
# [560] "LLL_exp_mean_density_HU"
# [561] "RL_exp_emphysema_15percentile_HU"
# [562] "LL_exp_emphysema_15percentile_HU"
# [563] "RUL_exp_emphysema_15percentile_HU"
# [564] "RML_exp_emphysema_15percentile_HU"
# [565] "RLL_exp_emphysema_15percentile_HU"
# [566] "LUL_exp_emphysema_15percentile_HU"
# [567] "LLL_exp_emphysema_15percentile_HU"
# [568] "RL_exp_LAA_.856HU_."
# [569] "LL_exp_LAA_.856HU_."
# [570] "RUL_exp_LAA_.856HU_."
# [571] "RML_exp_LAA_.856HU_."
# [572] "RLL_exp_LAA_.856HU_."
# [573] "LUL_exp_LAA_.856HU_."
# [574] "LLL_exp_LAA_.856HU_."
# [575] "Lungs_normal_."
# [576] "RL_normal_."
# [577] "LL_normal_."
# [578] "RUL_normal_."
# [579] "RML_normal_."
# [580] "RLL_normal_."
# [581] "LUL_normal_."
# [582] "LLL_normal_."
# [583] "Lungs_emphysema_."
# [584] "RL_emphysema_."
# [585] "LL_emphysema_."
# [586] "RUL_emphysema_."
# [587] "RML_emphysema_."
# [588] "RLL_emphysema_."
# [589] "LUL_emphysema_."
# [590] "LLL_emphysema_."
# [591] "Lungs_airtrapping_."
# [592] "RL_airtrapping_."
# [593] "LL_airtrapping_."
# [594] "RUL_airtrapping_."
# [595] "RML_airtrapping_."
# [596] "RLL_airtrapping_."
# [597] "LUL_airtrapping_."
# [598] "LLL_airtrapping_."
# [599] "Lungs_airtrapping_emphysema_."
# [600] "RL_airtrapping_emphysema_."
# [601] "LL_airtrapping_emphysema_."
# [602] "RUL_airtrapping_emphysema_."
# [603] "RML_airtrapping_emphysema_."
# [604] "RLL_airtrapping_emphysema_."
# [605] "LUL_airtrapping_emphysema_."
# [606] "LLL_airtrapping_emphysema_."
