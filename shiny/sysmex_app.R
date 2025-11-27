#Shiny app for Sherlock RNA-Seq - Sysmex - Analysis Nov 2025
# ml RPlus
# R
# library(shiny)
# runApp("shiny/sysmex_app.R", host = "0.0.0.0", port = 3838,launch.browser = FALSE)
# 
# wsl
# cd ${HOME} #to redirect to the linux directory
# 
# 
# #for login node (head node)
# ssh  -i ~/.ssh/private_key_ssh \
# -N -J umcg-kphung@tunnel.hpc.rug.nl \
# -L 3838:localhost:3838 \
# umcg-kphung@nibbler
# 
# # for compute node (interactive job) - match the node ID
# ssh  -i ~/.ssh/private_key_ssh \
# -N -J umcg-kphung@tunnel.hpc.rug.nl \
# -L 3838:localhost:3838 \
# umcg-kphung@nb-node-a02


# In Chrome: http://localhost:3838/

# ================================================================================================================================== #
# ================================== ******** INSTRUCTIONS FOR USER :) ********* ===================================================
# ==================================================================================================================================#

# STEP 1: Change the line below to your directory (path to where the 'shiny' folder is located on your computer)
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"


# STEP 2: Click "Run App" in the top right of the source pane OR highlight all of the script after this line to the end of the page and press Ctrl + Enter to run

setwd(my_directory)
# .libPaths(file.path(my_directory,"sherlock_app_library")) 

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

library("shiny")
library("DT")
library("shinydashboard")
library("tidyverse")
library("ggplot2")
library("ggrepel")
library("rstatix")
library("limma")
library("ggpubr")
library("readxl")
library("edgeR")

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir, "data")

# #Output directory
output.dir <- file.path(main.dir, "output")

processed.data.dir <- file.path(data.dir, "processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")


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

patient_demographics <- readRDS(file.path(output.dir,"qc", "patient_demographics_postcombat.rds"))
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

#Some rows for sysmex varaible have "----" which turn it into a character. make these NA
clinical_brushbiopt[which(clinical_brushbiopt$RELYMP.103uL == "----"),"RELYMP.103uL"] <- NA
clinical_brushbiopt$RELYMP.103uL <- as.numeric(clinical_brushbiopt$RELYMP.103uL)

clinical_brushbiopt[which(clinical_brushbiopt$RELYMP == "----"),"RELYMP"] <- NA
clinical_brushbiopt$RELYMP <- as.numeric(clinical_brushbiopt$RELYMP)


clinical_brush <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),] #270
clinical_biopt <-  clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),] #271

cat("Finished loading data")


# ================================================================================== #
# A. SHINY (UI) ==================================================================
# ================================================================================== #

ui <- dashboardPage(
  
  dashboardHeader(title = "SHERLOCK Sysmex"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("1) View Data (post QC)", tabName = "viewdata", icon = icon("search")),
      menuItem("2) Differential Expression", tabName = "comparison", icon = icon("bar-chart")),
      menuItem("3) Gene Expression", tabName = "gene_expression", icon = icon("bar-chart"))
    )
    
  ), #closedashboardSidebar
  
  
  dashboardBody(
    title = "SHERLOCK - Sysmex RNA-Seq",
    tabItems(
      
      
      ## UI 1) View Data (post QC)  ---------------------------------------------------------------------------------
      tabItem("viewdata",
              
              tabsetPanel(
                
                tabPanel("Clinical",
                         DTOutput("viewclinical")), #show table
                
                tabPanel("Counts",
                         DTOutput("viewcounts")), #show table
                
                tabPanel("Patient Demographics",
                         DTOutput("viewdemographics")) #show table
                
              ) #close tabset Panel
              
              
              
      ), #close tabItem
      
      
      ## UI 2) Differential Expression -------------------------------------------------------------------------------
      tabItem("comparison",
              selectInput("study2", 
                          label = "1) Select study & sample type of interest",
                          choices = c("SHERLOCK - Brush", "SHERLOCK - Biopsy"),
                          width = "50%"),
              
              selectInput("sysmex_variable",
                          label = "2) Select Sysmex variable of interest",
                          choices = colnames(clinical_brushbiopt)[619:ncol(clinical_brushbiopt)],
                          width = "50%"),
              
              h4("Press backspace to search, leave blank to include all genes"),
              
              tabsetPanel(
                tabPanel("Results",
                         DT::DTOutput("deg2") #show table
                         
                ),
                tabPanel("Volcano",
                         plotOutput("volcano2",
                                    width = "800px",
                                    height = "800px")
                         
                         
                )
                
              ) #close TabsetPanel
      ), #close Tab Item
      
      ## UI 3) Gene Expression -------------------------------------------------------------------------------
      tabItem("gene_expression",
              selectInput("study3", 
                          label = "1) Select study & sample type of interest",
                          choices = c("SHERLOCK - Brush", "SHERLOCK - Biopsy"),
                          width = "50%"),
              
              selectizeInput("gene3",
                             label = "2) Enter gene (HGNC symbol) of interest (case sensitive)",
                             width = "50%",
                             choices = NULL),
              
              selectInput("sysmex_variable3",
                          label = "3) Select Sysmex variable of interest",
                          choices = colnames(clinical_brushbiopt)[619:ncol(clinical_brushbiopt)],
                          width = "50%"),
              
              downloadButton("downloadscatterplot3",
                             "Download expression data for this gene"),
              
              plotOutput("scatterplot3",
                         width = "650px",
                         height = "650px"),
              
              p("Nominal P values from differential expression shown. Counts were voom normalised (limma version 3.54.0)") 
              
              
      )
      
      
      
    ) #close tabItems
  )#close Dashboard body
) #close UI




server <- function(input, output, session){
  
  
  # ================================================================================== #
  ## SERVER 1) VIEW DATA  =============================================================
  # ================================================================================== #
  
  
  ### output$viewcounts View counts file -------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  output$viewcounts <-  renderDT({
    
    datatable(counts,
              extensions = c('Buttons'),
              options = list(pageLength = 25,
                             scrollX = TRUE,
                             scrollY = "1000px",
                             dom = 'Bfrtip', 
                             buttons = c('copy', 'csv', 'excel'))
    ) #close datatable
  }) #close render DT
  
  
  ### output$viewclinical  View clinical file -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$viewclinical <-  renderDT({
    
    datatable(clinical_brushbiopt,
              extensions = c('Buttons'),
              filter = list(position = "top", 
                            clear = TRUE),
              options = list(pageLength = 25,
                             scrollX = TRUE,
                             scrollY = "1000px",
                             dom = 'Bfrtip', 
                             buttons = c('copy', 'csv', 'excel'))
    ) #close data table
  })#close render DT
  
  
  ### output$viewdemographics View patient demographics file ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$viewdemographics<-  renderDT({
    patient_demographics <- readRDS(file.path(output.dir,"qc", "patient_demographics_postcombat.rds"))
    
    datatable(patient_demographics,
              extensions = c('Buttons'),
              options = list(
                dom = 'Bfrtip', 
                buttons = c('copy', 'csv', 'excel', 'print'))
    ) #close data table
  }) #close render DT
  
  
  # ================================================================================== #
  ## SERVER 2) DIFFERENTIAL EXPRESSION ===============================================
  # ================================================================================== #
  
  
  ### Reactive expression to load the correct data ---------------------------------------------------------------------------------------------------------------------------------------------------------
  study_data2 <- reactive({
    req(input$study2)  # ensure input is not NULL
    
    listofresults <- readRDS(file.path(output.dir, "sysmex", "diffexp_withgroup", "results", "listofresults.RDS"))
    if(input$study2 == "SHERLOCK - Brush"){
      
      
      list(
        clinical2 = clinical_brush,
        counts2 = counts_brush,
        tT = listofresults[["brush"]][["tT"]],
        tT2 = listofresults[["brush"]][["tT2"]]
      )
    }
    
    else if(input$study2 == "SHERLOCK - Biopsy"){
      
      
      list(
        clinical2 = clinical_biopt,
        counts2 = counts_biopt,
        tT = listofresults[["biopt"]][["tT"]],
        tT2 = listofresults[["biopt"]][["tT2"]]
      )
    }
    
    
  }) #close reactive study_data2
  
  
  
  
  ### Reactive differential expression script -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  ### output$deg2 Render results tT table ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$deg2 <- DT::renderDT({
    
    this_sysmex_variable <- input$sysmex_variable
    
    withProgress(message = "Loading edgeR differential expression results...", value = 0.5, {
      
      
      tT <- as.data.frame(study_data2()$tT[[this_sysmex_variable]])
      
      incProgress(0.8)
      
      datatable(tT,
                extensions = c('Buttons'),
                filter = list(position = "top", 
                              clear = TRUE),
                options = list(pageLength = 25,
                               scrollX = TRUE,
                               scrollY = "1000px",
                               dom = 'Bfrtip', 
                               buttons = c('copy', 'csv', 'excel')))
      
    }) #close withProgress
  }) #close renderDT
  
  
  ### output$volcano2 Render volcano plot  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$volcano2 <-renderPlot({
    
    withProgress(message = "Loading plot ...", value = 0.5, {
      
      this_sysmex_variable <- input$sysmex_variable
      
      clinical2 = study_data2()$clinical2
      counts2 = study_data2()$counts2 
      
      #Remove the NAs -  outliers
      clinical2 <- clinical2[!is.na(clinical2[,this_sysmex_variable]),]
      counts2 <- counts2[,row.names(clinical2)]
      
      
      # Filter out outliers (3sd from the mean) - these are falsely driving the diffexp results
      stdev  <- sd(clinical2[,this_sysmex_variable], na.rm = TRUE)
      
      #if value is 3 standard deviations away from the mean, make it NA
      if (!is.na(stdev) && stdev > 0) {
        mean_val <- mean(clinical2[,this_sysmex_variable], na.rm = TRUE)
        outliers <- abs(clinical2[,this_sysmex_variable] - mean_val) > 3 * stdev
        clinical2[outliers, this_sysmex_variable] <- NA
        if(sum(outliers > 0)){
          caption <- paste("Removed outliers",paste0(clinical2[outliers, "Study.ID"], collapse = ","))}
        else{
          caption <- NULL}
      }
      
      
      
      tT <- as.data.frame(study_data2()$tT[[this_sysmex_variable]])
      tT2 <- as.data.frame(study_data2()$tT2[[this_sysmex_variable]])
      
      
      
      incProgress(0.8)
      
      volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
        geom_point(aes(color = Legend)) +
        scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
        geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
        geom_text_repel(data = subset(tT2[1:30,]),
                        aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines") ) +
        theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                         legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 16)) +
        labs(title = paste0(input$study2,"\n",input$sysmex_variable)) +
        labs(caption = caption)
    })
    print(volcano)
    
  }) #close renderPlot volcano2
  
  
  
  
  # ================================================================================== #
  ## SERVER 3) GENE EXPRESSION (Boxplots) ===============================================
  # ================================================================================== #
  
  ### Reactive expression to load the correct data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  study_data3 <- reactive({
    req(input$study3)  # ensure input is not NULL
    
    withProgress(message = "Loading required data ...", value = 0.2, {
      
      listofresults <- readRDS(file.path(output.dir, "sysmex", "diffexp_withgroup", "results", "listofresults.RDS"))
      
      
      if(input$study3 == "SHERLOCK - Brush"){
        
        listoftT <- listofresults[["brush"]][["tT"]]
        clinical <- clinical_brush
        counts <- counts_brush
        counts_voom <- voom(counts_brush)
        
        this_sysmex_variable <- input$sysmex_variable3
        
        #Remove the NAs
        clinical <- clinical[!is.na(clinical[,this_sysmex_variable]),]
        counts <- counts_voom$E[,row.names(clinical)]
        incProgress(0.7)
        
        print("test1")
        print(head(clinical))
        
        list(
          clinical2 = clinical,
          counts2 = counts,
          listoftT = listoftT)
      }
      
      else if(input$study3 == "SHERLOCK - Biopsy"){
        
        listoftT <- listofresults[["biopt"]][["tT"]]
        clinical <- clinical_biopt
        counts <- counts_biopt
        counts_voom <- voom(counts_biopt)
        
        this_sysmex_variable <- input$sysmex_variable3
        
        clinical <- clinical[!is.na(clinical[,this_sysmex_variable]),]
        counts <- counts_voom$E[,row.names(clinical)]
        incProgress(0.7)
        
        
        list(
          clinical2 = clinical,
          counts2 = counts,
          listoftT = listoftT)
      }
      
      
      
    }) #close progress message
    
  }) #close reactive study_data3
  
  ### Update selectizeInput for gene3  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### Include the ones with hgnc symbols and without
  observe({
    updateSelectizeInput(session, "gene3", choices = c(hgnc_symbols_db$SYMBOL, setdiff(row.names(counts), hgnc_symbols_db$GENEID)), server = TRUE)
  }) 
  
  
  
  ### output$scatterplot3 Render scatterplots ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$scatterplot3 <- renderPlot({
    
    
    req(input$gene3)
    req(input$sysmex_variable3)
    

      listoftT = study_data3()$listoftT
      clinical2 = study_data3()$clinical2
      counts2 = study_data3()$counts2 #these are voom counts!
      this_sysmex_variable <- input$sysmex_variable3
      gene <- input$gene3
      
      
      withProgress(message = "Creating boxplot...", value = 0.5, {
        
      print("test2")
      print(head(clinical2))
      
      gene_symbol <- ifelse(gene %in% hgnc_symbols_db$SYMBOL,
                            hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == gene), "GENEID"],
                            gene)
      
      
      
      
      # Filter out outliers (3sd from the mean) - these are falsely driving the diffexp results
      stdev  <- sd(clinical2[,this_sysmex_variable], na.rm = TRUE)
      
      #if value is 3 standard deviations away from the mean, make it NA
      if (!is.na(stdev) && stdev > 0) {
        mean_val <- mean(clinical2[,this_sysmex_variable], na.rm = TRUE)
        outliers <- abs(clinical2[,this_sysmex_variable] - mean_val) > 3 * stdev
        clinical2[outliers, this_sysmex_variable] <- NA
        if(sum(outliers > 0)){
          caption <- paste("Removed outliers",paste0(clinical2[outliers, "Study.ID"], collapse = ","))}
        else{
          caption <- NULL}
      }
      
      
      print(this_sysmex_variable)
      print(clinical2[,this_sysmex_variable])
      
      scatterplot_data <- as.data.frame(cbind(sysmex_variable = clinical2[,this_sysmex_variable],
                                              gene = counts2[gene_symbol,],
                                              classification = clinical2$classification,
                                              sample = clinical2$Study.ID))
      
      
      scatterplot_theme <- theme(axis.title = element_text(size = 24),
                                 axis.text = element_text(size = 24),
                                 title = element_text(size = 20),
                                 legend.text = element_text(size = 16),
                                 legend.position = "bottom")
      
      
      
      logFC_res <- listoftT[[this_sysmex_variable]][gene_symbol,"logFC"]
      pval_res <- listoftT[[this_sysmex_variable]][gene_symbol,"PValue"]
      
      
      #geom_point, split by disease
      boxplotimage <- ggplot(scatterplot_data, aes(
        x = as.numeric(sysmex_variable),
        y = as.numeric(gene))) +
        
        theme_bw()+
        scatterplot_theme +
        geom_point(aes(colour=classification)) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(colour = "black")) +
        
        
        theme(axis.text.x = element_text(size = 18))+
        labs(title = paste0(gene,"_vs_", this_sysmex_variable),
             color = "Disease Severity" #legend title
        ) +
        scale_color_manual(values = c("Control" = "#00BA38",
                                      "Mild-moderate COPD" = "#619CFF",
                                      "Severe COPD" = "#F8766D"))+
        ylab (label = gene) +
        xlab (label = this_sysmex_variable) +
        
        
        
        annotate(
        "text",
        x = -Inf,# adjust horizontally
        y = -Inf,  # adjust vertically
        hjust = -0,
        vjust = -0.5,
        label = paste0("logFC = ", signif(logFC_res, 3), "\n", 
                       "p = ", signif(pval_res, 3)),
        size = 5
      )
      
      # labs(caption = caption)
      
    }) #close withProgress
    
    print(boxplotimage)
    
  }) #Close render plot
  
  
  
  ### output$downloadboxplot3  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$downloadscatterplot3<-downloadHandler(
    
    filename = function() {
      paste0(input$gene3,"_data.csv")
    },
    content = function(file) { 
      
      
      
      clinical2 = study_data3()$clinical2
      counts2 = study_data3()$counts2 #these are voom counts!
      this_sysmex_variable <- input$sysmex_variable3
      gene <- input$gene3
      
      gene_symbol <- ifelse(gene %in% hgnc_symbols_db$SYMBOL,
                            hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == gene), "GENEID"],
                            gene)
      
      
      scatterplot_data <- as.data.frame(cbind(sysmex_variable = clinical2[,this_sysmex_variable],
                                              gene = counts2[gene_symbol,],
                                              classification = clinical2$classification,
                                              sample = clinical2$Study.ID))
      
      
      write.csv(scatterplot_data, file, row.names = T)
    } #close content
  ) #close downloadHandler
  
  
}#close servere


shinyApp(ui = ui, server = server)
