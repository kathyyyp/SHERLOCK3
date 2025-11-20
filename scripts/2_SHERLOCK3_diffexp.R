# SHERLOCK3 - Differential Expression
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
library("sva")


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



# ================================================================================== #
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
clinical$treat <- paste0(clinical$sampletype,"_",clinical$classification)

listofclinical <- list(brush = clinical_brush, #1
                       biopt = clinical_biopt, #2
                       brushbiopt = clinical)  #3

listofcounts <- list(brush = counts_brush,     #1
                     biopt = counts_biopt,     #2
                     brushbiopt = counts)      #3


listofvariables <- list(brush ="classification",  #1
                        biopt ="classification",  #2
                        brushbiopt="treat")       #3
# design 1 and 2 - compare disease severity, within brush and biopsy
design1 <- model.matrix(~0 + classification + age + sex + smoking_status, #Removing the intercept by using "0 +" means “I want to estimate the expression level for each group independently, and I’ll decide later how to compare them.” instead of edgeR automatically taking 'COntrol' as the reference value since it's the first level in the factor classification (ie. control becomes the intercept)
                        data = clinical_brush) #Design matrix

colnames(design1)[1:3] <- c(levels(as.factor(clinical_brush$classification)))

design2 <- model.matrix(~0 + classification + age + sex + smoking_status,
                        data = clinical_biopt) #Design matrix

colnames(design2)[1:3] <- c(levels(as.factor(clinical_biopt$classification)))


# design 3 - compare between brush and biopsy
design3 <- model.matrix(~0 + treat + age + sex + smoking_status,
                        data = clinical) #Design matrix


print(colnames(design3))
colnames(design3)[1:6] <- c(levels(as.factor(clinical$treat)))


listofdesigns <- list(brush = design1, #1 brush
                      biopt = design2, #2 biopt
                      brushbiopt = design3 #3 brush vs biopt
)

listofcontrasts <- list(
  
  #1
  brush = makeContrasts(contrast1 = Severe.COPD - Control,
                        contrast2 = Mild.moderate.COPD - Control,
                        contrast3 = Severe.COPD - Mild.moderate.COPD,
                        levels = listofdesigns[[1]]),
  
  #2
  biopt = makeContrasts(contrast1 = Severe.COPD - Control,
                        contrast2 = Mild.moderate.COPD - Control,
                        contrast3 = Severe.COPD - Mild.moderate.COPD,
                        levels = listofdesigns[[2]]),
  
  #3
  brushbiopt = makeContrasts(contrast4 = (Biopt_Severe.COPD - Brush_Severe.COPD) - (Biopt_Control - Brush_Control), #tells us what genes are more biopsy-specific than brushes in Severe COPD compared to controls. DOES not say that they are unique to biopsies. To do that we have to biopsy_severe - biopsy_control and then remove any genes that were also in brush_severe - brush_control and subtract
                             contrast5 = (Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD) - (Biopt_Control- Brush_Control),
                             contrast6 = (Biopt_Severe.COPD - Brush_Severe.COPD) - (Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD),
                             contrast7 = Biopt_Control - Brush_Control,
                             levels = listofdesigns[[3]])
  
)




nameconvert <- cbind(name = 
                       c("contrast1", 
                         "contrast2", 
                         "contrast3",
                         "contrast4",
                         "contrast5",
                         "contrast6",
                         "contrast7"),
                     
                     contrast=
                       c("Severe.COPD - Control",
                         "Mild.moderate.COPD - Control",
                         "Severe.COPD - Mild.moderate.COPD",
                         "(Biopt_Severe.COPD - Brush_Severe.COPD) - (Biopt_Control - Brush_Control)",
                         "(Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD) - (Biopt_Control- Brush_Control)",
                         "(Biopt_Severe.COPD - Brush_Severe.COPD) - (Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD)",
                         "Biopt_Control - Brush_Control")
)
nameconvert <- as.data.frame(nameconvert)

listofresults <- list()


## DISEASE SEVERITY EdgeR DIFFERENTIAL EXPRESSION -------------------------------------------------------------------------------------
diffexp_edgeR <- function(this.diffexp.dir, showEnsemblID = FALSE) {
  
  for (j in names(listofclinical)){
    
    clinical2 <- listofclinical[[j]]
    counts2 <- listofcounts[[j]]
    
    
    # FILTER - Remove genes that are unexpressed or lowly expressed
    DGEL<- DGEList(counts=counts2, group = clinical2[,listofvariables[[j]]])
    keep <- filterByExpr(DGEL) #leaves 22967 genes #filterbyExpr needs grouping variable otherwise it will filter toomuch (we grouped it in the DGEL list)
    # DGEList is a list-based data object. It has a matrix 'counts', a data.frame 'samples' (has info about the sample) with a column "lib.size" for library size or sequency depth, 
    # The filterBy Expr function keeps rows that have worthwhile counts in a minimum number of samples.
    # keep <- which(rowMedians(as.matrix(DGEL))>10) #leaves 36,415 genes
    # keep <- rowSums(cpm(expression)>100) >= 2
    
    
    DGEL<-DGEL[keep, , keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
    # expression <- DGEL$counts
    
    #NORMALISE
    # normalise the library sizes (If a small proportion of highly expressed genes consume a substantial proportion of the total library size for a particular sample, this will cause the remaining genes to be under-sampled for that sample and appear falsely as downregulated)
    # the default method is trimmed mean of M-values (TMM) between each pair of samples
    DGEL<- calcNormFactors(DGEL,method = "TMM")
    
    #DESIGN MATRIX
    design <- listofdesigns[[j]] 
    
    # ESTIMATE DISPERSON
    # estimate common dispersion, trended dispersions and tagwise dispersions (dispersion is the variation in gene expression across samples
    # estimatedisp() calculates
    # common dispersion (global measure of disperson across dataset)
    # tagwise disperson (disperson as a function of mean expression level of genes. genes with higher expression tend to have higher variability/larger dispersion so this accounts for that)
    # trended disperson (calculated seperately for each gene)
    DGEL <- estimateDisp(DGEL, design)
    
    #FIT MODEL
    # Fit a negative binomial generalized log-linear model (GLM) to the read counts for each gene.
    fit <- glmQLFit(DGEL, design) #fit the GLM (design) to the DGEL(the DGEL object, which contains the counts data that has been filtered,normalised and dispersons estimtated)
    #fits a GLM with quasi-likelihood approach(QL approach allows for the overdispersion present in RNA-Seq count data by modeling it appropriately.)
    # A GLM is a flexible class of models that generalizes linear regression to accommodate data that may not necessarily follow a normal distribution.
    # Set up contrasts between our variables of interest and recalculate model coefficients.
    # Experimental vs Control
    # makeContrasts(contrasts=variable being tested , levels=model).
    
    my.contrasts <- listofcontrasts[[j]]
    
    
    listoftT <- list()
    listoftT2 <- list()
    listofvolcano <- list()
    
    for (i in colnames(my.contrasts)){
      qlf <- glmQLFTest(fit, contrast=my.contrasts[,i]) #after fitting GLM to counts data (glmQLFit), glmQLFTest perform quasi likelihood f-tests to test for differential expression. ie. hypothesis testing for differential expression. (make inferences or draw conclusions about the data)
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
      
      
      listoftT[[i]] <- tT
      listoftT2[[i]] <- tT2
      
      
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
      listofvolcano[[i]] <- volcano
      
    } #close loop (contrast)
    
    listofresults[[j]] <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)
    
  } #close loop (brush or biopt)
  
  
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
  
  #NOTE!!!!Below code doesn't work in hpc when run as a job, but works in interactive if listofresults is loaded manually with readRDS and the code is run???
  pdf(file= file.path(volcano.dir, "volcano_plots_hpc.pdf"), width = 8, height = 8)
  
  #### Biopt Control vs Brush COntrol --------------------------------------------------------------------------------------
  listofresults[["brushbiopt"]][["volcano"]][["contrast7"]] + ggtitle(paste("Biopt_Control - Brush_Control"))
  
  #### Severe COPD - Control ---------------------------------------------------------------------------------------
  listofresults[["brush"]][["volcano"]][["contrast1"]] + ggtitle(paste("Brush:","Severe COPD - Control"))
  listofresults[["biopt"]][["volcano"]][["contrast1"]] + ggtitle(paste("Biopt:","Severe COPD - Control"))
  
  #---Venn diagram to compare biopt and brush#
  x <- list(Brush = row.names(listofresults[["brush"]][["tT2"]][["contrast1"]]), 
            Biopsy = row.names(listofresults[["biopt"]][["tT2"]][["contrast1"]])
  )
  
  (ggvenn(x, stroke_size = 1.5) + ggtitle("Severe COPD - Control"))
  
  #----Delta DE to compare biopt and brush
  listofresults[["brushbiopt"]][["volcano"]][["contrast4"]] + 
    ggtitle(paste("(Biopt_Severe.COPD - Brush_Severe.COPD) - \n (Biopt_Control - Brush_Control)")) + 
    theme(plot.title = element_text(size = 15)) 
  
  #### Mild.moderate.COPD - Control --------------------------------------------------------------------------------------
  listofresults[["brush"]][["volcano"]][["contrast2"]] + ggtitle(paste("Brush:","Mild.moderate.COPD - Control"))
  listofresults[["biopt"]][["volcano"]][["contrast2"]] + ggtitle(paste("Biopt:","Mild.moderate.COPD - Control"))
  
  #---Venn diagram to compare biopt and brush#
  x <- list(Brush = row.names(listofresults[["brush"]][["tT2"]][["contrast2"]]), 
            Biopsy = row.names(listofresults[["biopt"]][["tT2"]][["contrast2"]]))
  
  (ggvenn(x, stroke_size = 1.5) + ggtitle("Mild.moderate.COPD - Control"))
  
  #---Delta DE to compare biopt and nasal
  listofresults[["brushbiopt"]][["volcano"]][["contrast5"]] + 
    ggtitle(paste("(Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD) - \n (Biopt_Control- Brush_Control)")) + 
    theme(plot.title = element_text(size = 15)) 
  
  #### Severe.COPD - Mild.moderate.COPD --------------------------------------------------------------------------------------
  listofresults[["brush"]][["volcano"]][["contrast3"]] + ggtitle(paste("Brush:","Severe.COPD - Mild.moderate.COPD"))
  listofresults[["biopt"]][["volcano"]][["contrast3"]] + ggtitle(paste("Biopt:","Severe.COPD - Mild.moderate.COPD"))
  
  #---Venn diagram to compare biopt and nasal#
  x <- list(Brush = row.names(listofresults[["brush"]][["tT2"]][["contrast3"]]), 
            Biopsy = row.names(listofresults[["biopt"]][["tT2"]][["contrast3"]]))
  
  (ggvenn(x, stroke_size = 1.5) + ggtitle("Severe.COPD - Mild.moderate.COPD"))
  
  #---Delta DE to compare biopt and nasal
  listofresults[["brushbiopt"]][["volcano"]][["contrast6"]] + 
    ggtitle(paste("(Biopt_Severe.COPD - Brush_Severe.COPD) - \n (Biopt_Mild.moderate.COPD - Brush_Mild.moderate.COPD)")) + 
    theme(plot.title = element_text(size = 15)) 
  
  
  
  
  dev.off()
  
  
  # ================================================================================== #
  # 2.2. SAVE tT AND tT2 RESULTS =====================================================
  # ================================================================================== #
  cat("Saving tT AND tT2 RESULTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # diffexp.results.dir <- file.path(data.dir, "results")
  # if(!exists(diffexp.results.dir)) dir.create(diffexp.results.dir, recursive = TRUE)
  
  # diffexp.results.dir <- file.path(data.dir, "results", "hgnc_only")
  
  if(!exists(file.path(diffexp.results.dir,"tT"))) dir.create(file.path(diffexp.results.dir,"tT"), recursive = TRUE)
  if(!exists(file.path(diffexp.results.dir,"tT2"))) dir.create(file.path(diffexp.results.dir,"tT2"), recursive = TRUE)
  
  # for contrast1, contrast2 and contrast3
  for (i in 1:3){
    write.csv(listofresults[["brush"]][["tT"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir,"/tT/","brush", "_tT_", nameconvert$contrast[i], ".csv"))
    write.csv(listofresults[["brush"]][["tT2"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir, "/tT2/","brush", "_tT2_", nameconvert$contrast[i], ".csv"))
    
    write.csv(listofresults[["biopt"]][["tT"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir,"/tT/","biopt", "_tT_", nameconvert$contrast[i], ".csv"))
    write.csv(listofresults[["biopt"]][["tT2"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir, "/tT2/","biopt", "_tT2_", nameconvert$contrast[i], ".csv"))
  }
  
  #for contrast4, contrast5, contrast6 and contrast7
  for (i in 4:7){
    write.csv(listofresults[["brushbiopt"]][["tT"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir, "/tT/","brushbiopt", "_tT_", nameconvert$contrast[i], ".csv"))
    write.csv(listofresults[["brushbiopt"]][["tT2"]][[nameconvert$name[i]]], file = paste0(diffexp.results.dir, "/tT2/","brushbiopt", "_tT2_", nameconvert$contrast[i], ".csv"))
  }
  
  
  # ================================================================================== #
  # 2.3. BOXPLOTS ====================================================================
  # ================================================================================== #
  cat("Starting 2.3. BOXPLOTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  counts_brush_voom <- voom(counts_brush)
  counts_biopt_voom <- voom(counts_biopt)
  listofvoomcounts <- list(brush = counts_brush_voom, biopt = counts_biopt_voom)
  
  listofboxplots <- list()
  for (i in names(listofvoomcounts)){
    
    # # for specific gene
    # geneofinterest = "SERPINA3"
    # 
    # geneofinterestid = hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == geneofinterest), "GENEID"]
    # 
    
    boxplotclinical <- listofclinical[[i]]
    boxplotcounts <- listofvoomcounts[[i]]
    
    boxplotdata <- as.data.frame(t(boxplotcounts$E))
    
    boxplot <- cbind(boxplotdata,
                     classification = boxplotclinical$classification,
                     smokingstatus = boxplotclinical$smoking_status)
    #for top 10
    
    listofboxplots_contrast <- list()
    for (contrast in c("contrast1", "contrast2", "contrast3")){
      top10 <- c()
      
      tT2_genes<- row.names(listofresults[[i]][["tT2"]][[contrast]])
      
      top10 <- if (length(tT2_genes) >= 10) {
        tT2_genes[1:10]
      } else {
        tT2_genes
      }
      
      
      if(showEnsemblID == TRUE){#with ensembl 
        top10 <- ifelse(is.na(hgnc_symbols_db[top10, "SYMBOL"]),
                        top10,
                        hgnc_symbols_db[top10, "SYMBOL"])
        
        topgenes <- c(top10, "IL33","IL1RL1")
        
      }
      
      else{ #hgnc_only
        top10 <- hgnc_symbols_db[top10, "SYMBOL"]
        topgenes <- c(top10, "IL33","IL1RL1")
      }
      
      print(i)
      print(contrast)
      print(topgenes)  
      
      listofboxplots_topgenes <- list()
      for (gene in c(topgenes)){
        
        geneofinterest <- gene
        
        if(showEnsemblID == TRUE){
          #if some genes are ensembl and some are hgnc, just get ensembl ids (still need hgnc symbols saved to topgenes object so we can label the plot later)
          geneofinterestid <- ifelse(geneofinterest %in% hgnc_symbols_db$SYMBOL,
                                     hgnc_symbols_db[which(hgnc_symbols_db == geneofinterest), "GENEID"],
                                     geneofinterest)
          
          
        }
        
        else{
          #if all genes are hgnc symbols, match to hgnc_symbols_db and get ensembl IDs
          geneofinterestid <- hgnc_symbols_db[which(hgnc_symbols_db == geneofinterest), "GENEID"]
        }
        
        
        plot <- boxplot[,c(geneofinterestid,
                           "classification",
                           "smokingstatus")]
        
        colnames(plot)[1] <- "gene"
        

        plot <- as.data.frame(plot)
        
        ## Get P-Values --------------------------------------------------------------------------------------
        cat("Getting BOXPLOTS P-Values for annotation", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        
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
        stat.table3[which(stat.table3$resultsname == "contrast1"),"p"] <- listofresults[[i]][["tT"]][["contrast1"]][geneofinterestid[1], "PValue"]
        stat.table3[which(stat.table3$resultsname == "contrast2"),"p"] <- listofresults[[i]][["tT"]][["contrast2"]][geneofinterestid[1], "PValue"]
        stat.table3[which(stat.table3$resultsname == "contrast3"),"p"] <- listofresults[[i]][["tT"]][["contrast3"]][geneofinterestid[1], "PValue"]
        stat.table3$p <- signif(as.numeric(stat.table3$p), digits = 4)
        stat.table3$y.position <- max(plot[,"gene"]) + 0.025*(max(plot[,"gene"]))
        stat.table3$y.position <- as.numeric(stat.table3$y.position)
        
        # plot2 <- plot[-which(row.names(plot) == "106076-002-311"),] #remove outlier for il33
        
        cat("Making BOXPLOTS ggplot2", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        
        boxplotimage <- ggplot(plot, aes(
          x = as.factor(classification),
          y = gene
          # fill = smokingstatus
        )) +
          
          theme_bw()+
          
          
          
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
                     aes(fill = classification)) +
          
          
          scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                     "Severe.COPD" = "#F8766D")) +
          
          scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
          scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
          
          labs(title = paste0(geneofinterest,": SHERLOCK2&3 - ", i)) +
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
    listofboxplots[[i]] <- listofboxplots_contrast
    
  }
  
  ## Save boxplots --------------------------------------------
  cat("Saving BOXPLOTS", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  boxplot.dir <- file.path(diffexp.figures.dir, "boxplots")
  if (!exists(boxplot.dir)) dir.create(boxplot.dir, recursive = TRUE)
  
  dir.create(file.path(boxplot.dir, "severevscontrol"), recursive = TRUE)
  dir.create(file.path(boxplot.dir, "mildmoderatevscontrol"), recursive = TRUE)
  dir.create(file.path(boxplot.dir, "severevsmildmoderate"), recursive = TRUE)
  
  boxplot_theme <- theme(axis.title = element_text(size = 22),
                         axis.text = element_text(size = 22),
                         title = element_text(size = 18),
                         legend.position = "None") 
  
  #Save brush boxplots
  for (i in 1:length(listofboxplots[["brush"]][["contrast1"]])){
    png(file = paste0(boxplot.dir, "/severevscontrol/brush_",i, ")",names(listofboxplots[["brush"]][["contrast1"]])[i],".png"),
        height = 20,
        width= 20,
        units = "cm",
        res = 800)
    print(listofboxplots[["brush"]][["contrast1"]][[i]] + boxplot_theme)
    dev.off()
  }
  
  for (i in 1:length(listofboxplots[["brush"]][["contrast2"]])){
    png(file = paste0(boxplot.dir, "/mildmoderatevscontrol/brush_",i, ")",names(listofboxplots[["brush"]][["contrast2"]])[i],".png"),
        height = 20,
        width= 20,
        units = "cm",
        res = 800)
    print(listofboxplots[["brush"]][["contrast2"]][[i]] +  boxplot_theme)
    dev.off()
  }
  
  for (i in 1:length(listofboxplots[["brush"]][["contrast3"]])){
    png(file = paste0(boxplot.dir, "/severevsmildmoderate/brush_",i, ")",names(listofboxplots[["brush"]][["contrast3"]])[i],".png"),
        height = 20,
        width= 20,
        units = "cm",
        res = 800)
    print(listofboxplots[["brush"]][["contrast3"]][[i]] +  boxplot_theme)
    dev.off()
  }
  
  #Save biopsy boxplots
  for (i in 1:length(listofboxplots[["biopt"]][["contrast1"]])){
    png(file = paste0(boxplot.dir, "/severevscontrol/biopt_",i, ")",names(listofboxplots[["biopt"]][["contrast1"]])[i],".png"),
        height = 20,
        width= 20,
        units = "cm",
        res = 800)
    print(listofboxplots[["biopt"]][["contrast1"]][[i]] + boxplot_theme)
    dev.off()
  }
  
  for (i in 1:length(listofboxplots[["biopt"]][["contrast2"]])){
    png(file = paste0(boxplot.dir, "/mildmoderatevscontrol/biopt_",i, ")",names(listofboxplots[["biopt"]][["contrast2"]])[i],".png"),
        height = 20,
        width= 20,
        units = "cm",
        res = 800)
    print(listofboxplots[["biopt"]][["contrast2"]][[i]] + boxplot_theme)
    dev.off()
  }
  
  for (i in 1:length(listofboxplots[["biopt"]][["contrast3"]])){
    png(file = paste0(boxplot.dir, "/severevsmildmoderate/biopt_",i, ")",names(listofboxplots[["biopt"]][["contrast3"]])[i],".png"),
        height = 20,
        width= 20,
        units = "cm",
        res = 800)
    print(listofboxplots[["biopt"]][["contrast3"]][[i]] +  boxplot_theme)
    dev.off()
  }
  
  
} #close function

#RUN THE FUNCTION ---------------------------------------------------------------#
#if show ensemblID is set to TRUE, all plots will include the genes that dont have hgnc symbols and save these to a directory folder "diffexp"
#if show ensemblID is set to FALSE, all plots will include only the genes that  have hgnc symbols and save these to a directory folder "diffexp_hgnc_only"

diffexp_edgeR(this.diffexp.dir = diffexp.dir, showEnsemblID  = TRUE)
diffexp_edgeR(this.diffexp.dir = diffexp.hgnconly.dir, showEnsemblID = FALSE)
#--------------------------------------------------------------------------------#

cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

