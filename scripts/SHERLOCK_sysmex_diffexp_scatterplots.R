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
library("patchwork")

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
# 2) Combine clinical file with sysmex data and remove unecessary columns
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

#Some rows for sysmex varaible have "----" which turn it into a character. make these NA
clinical_brushbiopt[which(clinical_brushbiopt$RELYMP.103uL == "----"),"RELYMP.103uL"] <- NA
clinical_brushbiopt$RELYMP.103uL <- as.numeric(clinical_brushbiopt$RELYMP.103uL)

clinical_brushbiopt[which(clinical_brushbiopt$RELYMP == "----"),"RELYMP"] <- NA
clinical_brushbiopt$RELYMP <- as.numeric(clinical_brushbiopt$RELYMP)

sapply(clinical_brushbiopt[619:ncol(clinical_brushbiopt)],class)

clinical_brush <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),] #270
clinical_biopt <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),] #271




# # ================================================================================== #
# # 3. SCATTERPLOTS (individual) =====================================================
# # ================================================================================== #
# diffexp.dir <- file.path(output.dir, "diffexp_withgroup")
# diffexp.results.dir <- file.path(diffexp.dir, "results")
# diffexp.figures.dir <- file.path(diffexp.dir, "figures")
# 
# 
# listofresults <- readRDS(file.path(diffexp.results.dir, "listofresults3.RDS"))
# 
# 
# clinical <- clinical_brush
# counts <- counts_brush
# counts_voom <- voom(counts_brush)
# 
# gene <- "CAV1"
# gene_symbol=hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == gene), "GENEID"]
# sysmex_variable <- "IG"
# 
# clinical2 <- clinical[!is.na(clinical[,sysmex_variable]),]
# counts2 <- counts_voom$E[,row.names(clinical2)]
# 
# #Filter out variables where there are too many zeros (cutoff at >80% zeros)
# zero_proportion <- sum(clinical2[,this_sysmex_variable] == 0, na.rm = TRUE)/length(clinical2[,this_sysmex_variable])
# 
# 
# 
# # Filter out outliers (3sd from the mean) - these are falsely driving the diffexp results
# stdev  <- sd(clinical2[,this_sysmex_variable], na.rm = TRUE)
# 
# #if value is 3 standard deviations away from the mean, make it NA
# if (!is.na(stdev) && stdev > 0) {
#   mean_val <- mean(clinical2[,this_sysmex_variable], na.rm = TRUE)
#   outliers <- abs(clinical2[,this_sysmex_variable] - mean_val) > 3 * stdev
#   clinical2[outliers, this_sysmex_variable] <- NA
#   if(sum(outliers > 0)){
#     caption <- paste("Removed outliers",paste0(clinical2[outliers, "Study.ID"], collapse = ","))}
#   else{
#     caption <- NULL}
# }
# 
# 
# 
# 
# scatterplot_data <- as.data.frame(cbind(sysmex_variable = clinical2[,sysmex_variable],
#                                         gene = counts2[gene_symbol,],
#                                         classification = clinical2$classification,
#                                         sample = clinical2$Study.ID))
# 
# 
# scatterplot_theme <- theme(axis.title = element_text(size = 24),
#                            axis.text = element_text(size = 24),
#                            title = element_text(size = 20),
#                            legend.text = element_text(size = 18))
# 
# 
# 
# 
# 
# 
# #geom_point, split by disease
# boxplotfinal2 <- ggplot(scatterplot_data, aes(
#   x = as.numeric(sysmex_variable),
#   y = as.numeric(gene))) +
#   
#   theme_bw()+
#   scatterplot_theme +
#   geom_point(aes(colour=classification)) +
#   geom_text(aes(label = sample)) +
#   
#   # stat_pvalue_manual(stat.table.gsva,
#   #                    label = "p",
#   #                    tip.length = 0.01,
#   #                    size = 7)+
#   #
#   # # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
#   # stat_summary(fun.y = mean, fill = "red",
#   #              geom = "point", shape = 21, size =4,
#   #              show.legend = TRUE) +
#   #
# 
# 
# theme(axis.text.x = element_text(size = 18))+
#   labs(title = paste0(gene,"_vs_", sysmex_variable),
#        color = "Disease Severity" #legend title
#   ) +
#   scale_color_manual(values = c("Control" = "#00BA38",
#                                 "Mild-moderate COPD" = "#619CFF",
#                                 "Severe COPD" = "#F8766D"))+
#   ylab (label = gene) +
#   xlab (label = sysmex_variable) +
#   labs(caption = caption)
# 
# 
# ggsave(boxplotfinal2, filename = file.path(diffexp.figures.dir, paste0(gene,"_vs_", sysmex_variable,"_labelled_3sd.png")),
#        width = 30, height = 25,
#        units = "cm")
# 



# ================================================================================== #
# 4. SCATTERPLOTS LOOP ==================================================================
# ================================================================================== #
diffexp.dir <- file.path(output.dir, "diffexp_withgroup")
diffexp.results.dir <- file.path(diffexp.dir, "results")
diffexp.figures.dir <- file.path(diffexp.dir, "figures")


scatterplot_theme <- theme(axis.title = element_text(size = 24),
                           axis.text = element_text(size = 24),
                           title = element_text(size = 20),
                           legend.text = element_text(size = 18))


listofresults <- readRDS(file.path(diffexp.results.dir, "listofresults3.RDS"))


#for brush or biopt
for (i in c("brush", "biopt")){
  
# i = "brush"
  if(i == "brush"){
    listoftT2 = listofresults[["brush"]][["tT2"]]
    clinical = clinical_brush
    counts = counts_brush
  }else if (i == "biopt"){
    listoftT2 = listofresults[["biopt"]][["tT2"]]
    clinical = clinical_brush
    counts = counts_brush
  }
  
  message("Running voom once for: ", i)
  counts_voom <- voom(counts)
  
  
  
  # open pdf for this sample type
  pdf(file.path(diffexp.figures.dir,
                paste0(i,"_scatterplot_panels.pdf")),
      width = 25, height = 20)
  
  # loop through Sysmex varaibles
  for (sysmex_variable in names(listoftT2)) {
    
    message("Processing Sysmex variable: ", sysmex_variable)
    
    #get top 12 genes from tT
    tT2 <- listoftT2[[sysmex_variable]]
    if (nrow(tT2) >= 12) {
      top12 <- row.names(tT2)[1:12]
    } else {
      top12 <- row.names(tT2)
    }
    
    if(length(top12) == 0){next}
    
  
  #filter out the NAs  
    clinical2 <- clinical[!is.na(clinical[,sysmex_variable]),]
    counts2 <- counts_voom$E[,row.names(clinical2)]
  
    
    plotting_func <- function(gene){
      
      geneid <- gene
      gene_hgnc <- ifelse(gene %in% hgnc_symbols_db$GENEID,
                          hgnc_symbols_db[which(hgnc_symbols_db$GENEID == gene), "SYMBOL"],
                          gene)
      
      
      # Filter out outliers (3sd from the mean) - these are falsely driving the diffexp results
      stdev  <- sd(clinical2[,sysmex_variable], na.rm = TRUE)
      
      #if value is 3 standard deviations away from the mean, make it NA
      if (!is.na(stdev) && stdev > 0) {
        mean_val <- mean(clinical2[,sysmex_variable], na.rm = TRUE)
        outliers <- abs(clinical2[,sysmex_variable] - mean_val) > 3 * stdev
        clinical2[outliers, sysmex_variable] <- NA
      }
      

      
      scatterplot_data <- as.data.frame(cbind(sysmex_variable = clinical2[,sysmex_variable],
                                              gene = counts2[geneid,],
                                              classification = clinical2$classification,
                                              sample = clinical2$Study.ID))
      
      
      logFC_res <- tT2[geneid,"logFC"]
      pval_res <- tT2[geneid,"PValue"]
      
      #geom_point, split by disease
      scatterplotfinal <- ggplot(scatterplot_data, aes(
        x = as.numeric(sysmex_variable),
        y = as.numeric(gene))) +
        
        theme_bw()+
        scatterplot_theme +
        geom_point(aes(colour=classification)) +
        # geom_text(aes(label = sample)) +
        
        
        
        theme(axis.text.x = element_text(size = 18))+
        labs(title = paste0(gene_hgnc,"_vs_", sysmex_variable),
             color = "Disease Severity" #legend title
        ) +
        scale_color_manual(values = c("Control" = "#00BA38",
                                      "Mild-moderate COPD" = "#619CFF",
                                      "Severe COPD" = "#F8766D"))+
        ylab (label = gene_hgnc) +
        xlab (label = sysmex_variable) +
        
        
        annotate(
          "text",
          x = Inf,# adjust horizontally
          y = Inf,  # adjust vertically
          hjust = 1.05,
          vjust = 1.15,
          label = paste0("logFC = ", signif(logFC_res, 3), "\n", 
                         "p = ", signif(pval_res, 3)),
          size = 5
        )
      
    } #close plotting func
    
    
    # generate the page immediately
    plots <- map(top12, plotting_func)
    
    page <- wrap_plots(plots, ncol = 4, nrow = 3) +
    plot_layout(guides = "collect") +
      plot_annotation(
        title = paste("Sysmex variables & Gene expression - ", i),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      ) &
      theme(legend.position = "bottom")  # or "right", "top", "left"
    
    print(page)  # PDF page written immediately
    
    rm(plots, page, tT2, clinical2, counts2, top12)
    gc()
  } #close sysmex variable loop
  
  dev.off()
  
  rm(counts_voom)
  gc()
} #close sample type loop



  # 
  
  # 
  # # this works when testing one
  # library(patchwork)
  # 
  # pdf(file.path(diffexp.figures.dir,paste0(i,"_", sysmex_variable, "_scatterplot_plot_panels.pdf")), width = 10, height = 10)
  # wrap_plots(listofboxplots, ncol = 4, nrow = 3, guides = "collect") +
  #   
  #   # plot_layout(guides = "collect") +
  #   
  #   plot_annotation(
  #     # title = paste(sysmex_variable, "&", gene,"-", i),
  #     theme = theme(plot.title = element_text(size = 16, face = "bold"))
  #   ) &
  #   theme(legend.position = "bottom")  # or "right", "top", "left"
  # dev.off()
  # 
  
  
cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")












