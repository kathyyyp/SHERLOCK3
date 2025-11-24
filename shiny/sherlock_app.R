#Shiny app for Sherlock RNA-Seq Analysis Nov 2025

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


# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir, "data")

# #Output directory
output.dir <- file.path(main.dir, "output")

processed.data.dir <- file.path(data.dir,"processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


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


# ================================================================================== #
# A. SHINY (UI) ==================================================================
# ================================================================================== #

ui <- dashboardPage(
  
  dashboardHeader(title = "SHERLOCK Transcriptomics Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("1) View Data (post QC)", tabName = "viewdata", icon = icon("search")),
      menuItem("2) Differential Expression", tabName = "comparison", icon = icon("bar-chart")),
      menuItem("3) Gene Expression", tabName = "gene_expression", icon = icon("bar-chart")),
      menuItem("4) TF Analysis", tabName = "tf_analysis", icon = icon("bar-chart")),
      menuItem("5) Cellular Deconvolution", tabName = "celldecon", icon = icon("th"))
      
    ) #close sidebarMenu
  ), #closedashboardSidebar
  
  
  dashboardBody(
    title = "SHERLOCK",
    tabItems(
      
      ## UI 1) View Data (post QC)  ---------------------------------------------------------------------------------
      tabItem("viewdata",
              selectInput("study1", label = "Select study of interest",
                          choices = c("SHERLOCK1", "SHERLOCK2&3")),
              tabsetPanel(
                tabPanel("Clinical",
                         DTOutput("viewclinical")), #show table
                
                tabPanel("Counts",
                         DTOutput("viewcounts")), #show table
                
                tabPanel("Patient Demographics",
                         DTOutput("viewdemographics"),
                         column(width = 8, 
                                fluidRow("Note - SHERLOCK2&3: 4 outliers were removed following QC. Sherlock1: All samples from ex-smokers") 
                         )
                ), #show table
                
                tabPanel("QC details",
                         uiOutput("QCinfo")
                         
                )
                                
              ) #close tabset Panel
              
              
              
      ), #close tabItem
      
      
      ## UI 2) Differential Expression -------------------------------------------------------------------------------
      tabItem("comparison",
              selectInput("study2", 
                          label = "1) Select study & sample type of interest",
                          choices = c("SHERLOCK2&3 - Brush", "SHERLOCK1 - Brush", "SHERLOCK2&3 - Biopsy", "SHERLOCK2&3 - Brush & Biopsy"),
                          width = "50%"),
              
              selectInput("comparison2", 
                          label = "2) Select comparison of interest",
                          choices = nameconvert$contrast,
                          width = "50%"),
              
              h4("Press backspace to search, leave blank to include all genes"),
              
              tabsetPanel(
                tabPanel("Results",
                         DT::DTOutput("deg2"), #show table
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
                          choices = c("SHERLOCK2&3 - Brush", "SHERLOCK1 - Brush", "SHERLOCK2&3 - Biopsy"),
                          width = "50%"),
              
              selectizeInput("gene3",
                             label = "2) Enter gene (HGNC symbol) of interest (case sensitive)",
                             width = "50%",
                             choices = NULL),
              
              
              downloadButton("downloadboxplot3",
                             "Download expression data for this gene"),
              
              plotOutput("boxplot3",
                         width = "650px",
                         height = "650px"),
              
              p("Nominal P values from unpaired T-test shown. Counts were voom normalised (limma version 3.54.0)") 
              
              
      ),
      
      ## UI 4) TF Analysis Expression -------------------------------------------------------------------------------
      tabItem("tf_analysis",
              selectInput("study4", 
                          label = "1) Select study & sample type of interest",
                          choices = c("SHERLOCK2&3 - Brush", "SHERLOCK1 - Brush", "SHERLOCK2&3 - Biopsy", "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs"),
                          width = "50%"),
              
              selectInput("comparison4", 
                          label = "2) Select comparison of interest",
                          choices = nameconvert$contrast[1:3],
                          width = "50%"),
              
              selectizeInput("gene4",
                             label = "For Boxplot tab: Enter TF of interest (case sensitive)",
                             width = "50%",
                             choices = NULL),
              
              
              tabsetPanel(
                
                tabPanel("Volcano",
                         plotOutput("volcano4",
                                    width = "800px",
                                    height = "800px")),
                
                
                tabPanel("Results",
                         # downloadButton("downloadtable2","Download Results Table") #download table=
                         # div(style = 'overflow-x: scroll',=
                         DT::DTOutput("table4") #show table
                ),
                
                tabPanel("Boxplot",
                         
                         downloadButton("downloadboxplot4",
                                        "Download TF activity data for this TF"),
                         
                         plotOutput("boxplot4",
                                    width = "800px",
                                    height = "800px"),
                         
                         p("Nominal P values from unpaired T-test shown. TF activity inferred from RNA-data using",
                           br(),
                           "VIPER pipeline (version 1.32.0), with regulons inferred from ARACNe (Califano Lab .jar release).")
                         
                ),
                
              ) #close TabsetPanel
              
      ), #close tabItem
      
      ## UI 5) Cell Decon tab  ---------------------------------------------------------------------------------------------
      tabItem("celldecon",
              selectInput("study5", 
                          label = "Select study & sample type of interest",
                          choices = c("SHERLOCK2&3 - Brush", "SHERLOCK2&3 - Biopsy")),
              column(width = 8, 
                     plotOutput("celldeconplot5",
                                width = "1300px",
                                height = "700px"),
                     
                     fluidRow( "BH-corrected P values from unpaired T-test shown. Counts were CPM normalised for CIBERSORT") 
              )
              
      ) #close tabItem
      
    ) #close tabItems
  )#close Dashboard body
) #close UI





# ================================================================================== #
# B. SERVER ========================================================================
# ================================================================================== #

server <- function(input, output, session){
  
  # ================================================================================== #
  ## SERVER 1) VIEW DATA  =============================================================
  # ================================================================================== #
  
  ### Get files for selected study -----------------------------------------------------------------
  study_data1 <- reactive({
    req(input$study1)  # ensure input is not NULL

    
    if(input$study1 == "SHERLOCK1"){
      clinical_filepath <- file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_master.rds")
      counts_filepath <- file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds")
      
      list(
        clinical = readRDS(clinical_filepath),
        counts = readRDS(counts_filepath),
        patientdemographics = readRDS(file.path(output.dir, "sherlock1", "qc", "patient_demographics.rds"))
      )
    }
    
    else if(input$study1 == "SHERLOCK2&3"){
      
      clinical_filepath <- file.path(postQC.data.dir, "master","clinical_brushbiopt_master.rds")
      counts_filepath <- file.path(combat.processed.data.dir, "counts_combat.rds")
      
      list(
        clinical = readRDS(clinical_filepath),
        counts = readRDS(counts_filepath),
        patientdemographics = readRDS(file.path(qc.dir, "patient_demographics_postcombat.rds"))
        
      )
    }
  }) #close reactive study_data1
  
  ### output$viewcounts View counts file -------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$viewcounts <-  renderDT({
    
    datatable(study_data1()$counts,
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
    
    datatable(study_data1()$clinical, 
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
    datatable(study_data1()$patientdemographics,
              extensions = c('Buttons'),
              options = list(
                dom = 'Bfrtip', 
                buttons = c('copy', 'csv', 'excel', 'print'))
    ) #close data table
  }) #close render DT
  
  ### output$viewdemographics View patient demographics file ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$QCinfo <- renderUI({
    if (input$study1 == "SHERLOCK1") {
      HTML("Samples excluded:<br>
      - A_2804: A2804 (SHERLOCK1) and SEO230 (SHERLOCK2) are the same patient - advised to keep SEO230 in most cases <br>
      - LIB5426587_SAM24375559: Duplicated sample from patient A_1688. Kept LIB5426624_SAM24375596 due to greater library size <br>
      - A_580, A_994, A_984: Patients were excluded from SHERLOCK study")
    }
      
    else if (input$study1 == "SHERLOCK2&3") {
        HTML("Samples excluded:<br>
      - 106076-002-232, 106076-002-311, 106076-002-200, 106076-002-015: Excluded from SHERLOCK2 as <500,000 total reads <br>
      - 107165-001-146: Excluded from SHERLOCK3 as sample is labelled as female but has high number of counts (3.129 CPM) for DDX3Y (male-specific gene on Y-chromosome) <br>
      - 107165-001-016 and 107165-001-059 were derived from the same person (identical genotypes) but annotated as SEO267 and SEO265 respectively. Removed as possible sample switch. <br>
      <br>
      The below are duplicate samples from the same patient. Kept the sample with greater library size <br>
      - 106076-002-014: Kept duplicate sample 107165-001-010, Patient SEO077 <br>
      - 106076-002-120: Kept duplicate sample 107165-001-032, Patient SEO185 <br>
      - 107165-001-041: Kept duplicate sample 106076-002-168, Patient SEO310 <br>
      - 106076-002-009: Kept duplicate sample 107165-001-161, Patient SEO069 <br>
      - 106076-002-003: Kept duplicate sample 107165-001-160, Patient SEO070 <br>
      - 106076-002-017: Kept duplicate sample 107165-001-163, Patient SEO084 <br>
      - 107165-001-079: Kept duplicate sample 106076-002-023, Patient SEO096 <br>
      - 106076-002-159: Kept duplicate sample 107165-001-125, Patient SEO308 <br>
      - 106076-002-199: Kept duplicate sample 107165-001-071, Patient SEO317")
    } 
  })
  
  # ================================================================================== #
  ## SERVER 2) DIFFERENTIAL EXPRESSION ===============================================
  # ================================================================================== #
  
  ### Observe selectInput and change selection choices ---------------------------------------------------------------------------------------------------------------------------------------------------------
  
  observeEvent(input$study2, { #observe() automatically observes any reactive expressions in its body changes whereas observeEvent runs only when the specified event changes
    if (input$study2 == "SHERLOCK2&3 - Brush & Biopsy") {
      choices <- nameconvert$contrast[4:nrow(nameconvert)]
    } else {
      choices <- nameconvert$contrast[1:3]
    }
    
    updateSelectInput(session, "comparison2", choices = choices)
  }) #close observeEvent
  
  
  
  ### Reactive expression to load the correct data ---------------------------------------------------------------------------------------------------------------------------------------------------------
  study_data2 <- reactive({
    req(input$study2)  # ensure input is not NULL
    
    if(input$study2 == "SHERLOCK1 - Brush"){
      
      results_filepath <- file.path(output.dir, "sherlock1", "diffexp", "results", "listofresults.RDS")
      listofresults <- readRDS(results_filepath)
      
      list(
        listoftT = listofresults$tT,
        listoftT2 = listofresults$tT2
      )
    }
    
    else if(input$study2 == "SHERLOCK2&3 - Brush"){
      
      results_filepath <- file.path(output.dir, "diffexp", "results", "listofresults.RDS")
      listofresults <- readRDS(results_filepath)
      
      list(
        listoftT =listofresults$brush$tT,
        listoftT2 = listofresults$brush$tT2
      )
      
    }
    else if(input$study2 == "SHERLOCK2&3 - Biopsy"){
      
      results_filepath <- file.path(output.dir, "diffexp", "results", "listofresults.RDS")
      listofresults <- readRDS(results_filepath)
      
      list(
        listoftT = listofresults$biopt$tT,
        listoftT2 = listofresults$biopt$tT2
      )
      
    }
    else if(input$study2 == "SHERLOCK2&3 - Brush & Biopsy"){
      
      results_filepath <- file.path(output.dir, "diffexp", "results", "listofresults.RDS")
      listofresults <- readRDS(results_filepath)
      
      list(
        listoftT = listofresults$brushbiopt$tT,
        listoftT2 = listofresults$brushbiopt$tT2
      )
    }
    
  }) #close reactive study_data2
  
  
  
  ### output$deg2 Render results tT table ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$deg2 <- DT::renderDT({
    
    withProgress(message = "Loading data ...", value = 0.5, {
      c <- nameconvert[which(nameconvert$contrast == input$comparison2),"name"]
      listoftT <- study_data2()$listoftT
      listoftT2 <- study_data2()$listoftT2
      
      incProgress(0.8)
      
      tT <- as.data.frame(listoftT[[c]])
      tT2 <- as.data.frame(listoftT2[[c]])
      
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
      
      c <- nameconvert[which(nameconvert$contrast == input$comparison2),"name"]
      listoftT <- study_data2()$listoftT
      listoftT2 <- study_data2()$listoftT2
      
      tT <- as.data.frame(listoftT[[c]])
      tT2 <- as.data.frame(listoftT2[[c]])
      
      incProgress(0.8)
      
      volcano <- ggplot(tT, aes(x = logFC, y = -log10(PValue))) +
        geom_point(aes(color = Legend)) +
        scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"))+
        geom_hline(yintercept =-log10(max(tT2$PValue)),colour="black", linetype="dashed")+
        geom_vline(xintercept =-1,colour="black", linetype="dashed")+
        geom_vline(xintercept =1,colour="black", linetype="dashed")+
        geom_text_repel(data = tT2[1:20,],
                        aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines") ) +
        theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                         legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 16)) +
        labs(title = paste0(input$study2,"\n",input$comparison2))
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
      
      if(input$study3 == "SHERLOCK1 - Brush"){
        
        clinical <- readRDS(file.path(data.dir, "processed", "SHERLOCK1", "clinical_sk1_simple.rds"))
        counts <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
        
        
        results_filepath <- file.path(output.dir, "sherlock1", "diffexp", "results", "listofresults.RDS")
        listofresults <- readRDS(results_filepath)
        
        list(
          clinical = clinical,
          counts = counts,
          listoftT = listofresults$tT,
          listoftT2 = listofresults$tT2
        )
      }
      
      
      
    
      
      else if(input$study3 == "SHERLOCK2&3 - Brush"){
        
        clinical <- readRDS(file.path(postQC.data.dir, "clinical_brush_simple.rds"))
        counts <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
        
        results_filepath <- file.path(output.dir, "diffexp", "results", "listofresults.RDS")
        listofresults <- readRDS(results_filepath)
        
        list(
          clinical = clinical,
          counts = counts,
          listoftT =listofresults$brush$tT,
          listoftT2 = listofresults$brush$tT2
        )
        
      }
      else if(input$study3 == "SHERLOCK2&3 - Biopsy"){
        
        clinical <- readRDS(file.path(postQC.data.dir, "clinical_biopt_simple.rds"))
        counts <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
        results_filepath <- file.path(output.dir, "diffexp", "results", "listofresults.RDS")
        listofresults <- readRDS(results_filepath)
        
        list(
          clinical = clinical,
          counts = counts,
          listoftT = listofresults$biopt$tT,
          listoftT2 = listofresults$biopt$tT2
        )
        
      }
    }) #close progress message
    
  }) #close reactive study_data3
  
  
  ### Update selectizeInput for gene3  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  gene_name_conversion <- reactive ({
    genes_table <- as.data.frame(row.names(study_data3()$counts))
    row.names(genes_table) <- row.names(study_data3()$counts)
    genes_table$gene_symbol <- hgnc_symbols_db[row.names(study_data3()$counts), "SYMBOL"] #add hgnc symbols
    genes_table[which(is.na(genes_table$gene_symbol)), "gene_symbol"] <- row.names(genes_table)[(which(is.na(genes_table$gene_symbol)))] 
    colnames(genes_table)[1] <- "ensembl_id"
    genes_table
    
  }) #close gene_name_conversion reactive
  
  observe({
    updateSelectizeInput(session, "gene3", choices = gene_name_conversion()$gene_symbol, server = TRUE)
  }) 
  
  
  ### output$boxplot3 Render boxplots ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$boxplot3 <- renderPlot({
    
    req(input$gene3)
    
    withProgress(message = "Creating boxplot...", value = 0.5, {
      
      boxplotcounts <- study_data3()$counts
      boxplotclinical <- study_data3()$clinical
      
      boxplotcounts_voom <- voom(boxplotcounts)
      incProgress(0.7)
      
      boxplotdata <- as.data.frame(t(boxplotcounts_voom$E))
      
      
      boxplot <- cbind(boxplotdata,
                       classification = boxplotclinical$classification
      )
      
      geneofinterest <- input$gene3
      
      gene_name_conversion_table <- gene_name_conversion()
      
      geneofinterestid <- gene_name_conversion_table[which(gene_name_conversion_table$gene_symbol == input$gene3), "ensembl_id"]
      
      plot <- boxplot[,c(geneofinterestid,
                         "classification")]
      
      colnames(plot)[1] <- "gene"
      
      plot <- as.data.frame(plot)
      
      
    }) #close progressmessage
    
    stat.table <- plot %>%
      t_test(gene ~ classification,
             comparisons = list(c("Severe.COPD", "Control"),
                                c("Mild.moderate.COPD", "Control"),
                                c("Severe.COPD", "Mild.moderate.COPD")
             ))
    
    
    stat.table<- stat.table %>%
      add_xy_position(x = "classification", dodge = 0.8)
    
    stat.table3 <- stat.table
    
    # Replace stat table p values with diff exp p values
    listoftT <- study_data3()$listoftT
    stat.table3 <- cbind(stat.table3, resultsname = c("contrast1", "contrast2", "contrast3"))
    stat.table3[which(stat.table3$resultsname == "contrast1"),"p"] <- listoftT[["contrast1"]][geneofinterestid, "PValue"]
    stat.table3[which(stat.table3$resultsname == "contrast2"),"p"] <- listoftT[["contrast2"]][geneofinterestid, "PValue"]
    stat.table3[which(stat.table3$resultsname == "contrast3"),"p"] <- listoftT[["contrast3"]][geneofinterestid, "PValue"]
    stat.table3$p <- signif(as.numeric(stat.table3$p), digits = 4)
    stat.table3$y.position <- max(plot[,"gene"]) + 0.025*(max(plot[,"gene"]))
    stat.table3$y.position <- as.numeric(stat.table3$y.position)
    
    # Boxplot 
    boxplotimage <- ggplot(plot, aes(
      x = as.factor(classification),
      y = gene,
      group = classification
    )) +
      
      theme_bw(base_size =  18)+
      
      theme(axis.text = element_text(size = 18),
            title = element_text(size = 20),
            legend.text = element_text(size = 16),
            legend.position = "None") +
      
      
      geom_boxplot(position=position_dodge(1)) +
      
      geom_jitter(aes(color = classification),
                  alpha = 0.5,
                  size = 2.5, 
                  width = 0.2) +
      
      scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                 "Severe.COPD" = "#F8766D")) +
      
      # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      labs(title = paste0(geneofinterest, ": ",input$study3)) +
      ylab (label = "Expression") +
      xlab (label = "Disease") +
      
      
      stat_pvalue_manual(stat.table3,
                         label = "p",
                         tip.length = 0.01,
                         bracket.nudge.y = c(0, 0.5, 0.5),
                         size = 6) +
      
      stat_summary(fun = mean, fill = "red",
                   geom = "point", shape = 21, size =4,
                   show.legend = TRUE) 
    
    
    print(boxplotimage)
    
  }) #Close render plot
  
  
  
  ### output$downloadboxplot3  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$downloadboxplot3<-downloadHandler(
    
    filename = function() {
      paste0(input$gene3,"_data.csv")
    },
    content = function(file) { 
      
      
      boxplotcounts <- study_data3()$counts
      boxplotclinical <- study_data3()$clinical
      
      boxplotcounts_voom <- voom(boxplotcounts)
      incProgress(0.7)
      
      boxplotdata <- as.data.frame(t(boxplotcounts_voom$E))
      
      
      boxplot <- cbind(boxplotdata,
                       classification = boxplotclinical$classification
      )
      
      geneofinterest <- input$gene3
      
      gene_name_conversion_table <- gene_name_conversion()
      
      geneofinterestid <- gene_name_conversion_table[which(gene_name_conversion_table$gene_symbol == input$gene3), "ensembl_id"]
      
      plot <- boxplot[,c(geneofinterestid,
                         "classification")]
      
      colnames(plot)[1] <- "gene"
      
      plot <- as.data.frame(plot)
      
      write.csv(plot, file, row.names = T)
    } #close content
  ) #close downloadHandler
  



  # ================================================================================== #
  ## SERVER 4) TF ANALYSIS ===========================================================
  # ================================================================================== #

  ### Reactive expression to load the correct data ---------------------------------------------------------------------------------
  study_data4 <- reactive({
    req(input$study4)  # ensure input is not NULL

    if(input$study4 == "SHERLOCK1 - Brush"){

      clinical <- readRDS(file.path(data.dir, "processed", "SHERLOCK1", "clinical_sk1_simple.rds"))
      TF_results <- readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "TF_analysis_results.rds"))
      listoftT = readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "listoftT.rds"))
      listoftT2 = readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "listoftT2.rds"))
      
      names(listoftT2)[names(listoftT2) == "Mild.moderate.COPDvsControl"] <- "Mild.moderate.COPD - Control"
      names(listoftT2)[names(listoftT2) == "Severe.COPDvsControl"] <- "Severe.COPD - Control"
      names(listoftT2)[names(listoftT2) == "Severe.COPDvsMild.moderate.COPD"] <- "Severe.COPD - Mild.moderate.COPD"
      
      names(listoftT)[names(listoftT) == "Mild.moderate.COPDvsControl"] <- "Mild.moderate.COPD - Control"
      names(listoftT)[names(listoftT) == "Severe.COPDvsControl"] <- "Severe.COPD - Control"
      names(listoftT)[names(listoftT) == "Severe.COPDvsMild.moderate.COPD"] <- "Severe.COPD - Mild.moderate.COPD"

      list(
        TF_results = TF_results,
        listoftT = listoftT,
        listoftT2 = listoftT2,
        clinical = clinical
      )
    }

    else if(input$study4 == "SHERLOCK2&3 - Brush"){
      all_clinical <- readRDS(file.path(postQC.data.dir, "master","clinical_brushbiopt_master.rds"))
      sherlock2_brush_clinical <- all_clinical[which(all_clinical$sampletype == "Brush"),]
      TF_results <- readRDS(file.path(output.dir, "tf_analysis", "brush", "results", "TF_analysis_results.rds"))
      listoftT = readRDS(file.path(output.dir, "tf_analysis", "brush", "results", "listoftT.rds"))
      listoftT2 = readRDS(file.path(output.dir, "tf_analysis", "brush", "results", "listoftT2.rds"))

      names(listoftT2)[names(listoftT2) == "Mild.moderate.COPDvsControl"] <- "Mild.moderate.COPD - Control"
      names(listoftT2)[names(listoftT2) == "Severe.COPDvsControl"] <- "Severe.COPD - Control"
      names(listoftT2)[names(listoftT2) == "Severe.COPDvsMild.moderate.COPD"] <- "Severe.COPD - Mild.moderate.COPD"

      names(listoftT)[names(listoftT) == "Mild.moderate.COPDvsControl"] <- "Mild.moderate.COPD - Control"
      names(listoftT)[names(listoftT) == "Severe.COPDvsControl"] <- "Severe.COPD - Control"
      names(listoftT)[names(listoftT) == "Severe.COPDvsMild.moderate.COPD"] <- "Severe.COPD - Mild.moderate.COPD"

      list(
        TF_results = TF_results,
        listoftT = listoftT,
        listoftT2 = listoftT2,
        clinical = sherlock2_brush_clinical

      )
    }
    else if(input$study4 == "SHERLOCK2&3 - Biopsy"){
      all_clinical <- readRDS(file.path(postQC.data.dir, "master","clinical_brushbiopt_master.rds"))
      sherlock3_biopsy_clinical <- all_clinical[which(all_clinical$sampletype == "Biopt"),]
      TF_results <- readRDS(file.path(output.dir, "tf_analysis", "biopt", "results", "TF_analysis_results.rds"))
      listoftT = readRDS(file.path(output.dir, "tf_analysis", "biopt", "results", "listoftT.rds"))
      listoftT2 = readRDS(file.path(output.dir, "tf_analysis", "biopt", "results", "listoftT2.rds"))

      names(listoftT2)[names(listoftT2) == "Mild.moderate.COPDvsControl"] <- "Mild.moderate.COPD - Control"
      names(listoftT2)[names(listoftT2) == "Severe.COPDvsControl"] <- "Severe.COPD - Control"
      names(listoftT2)[names(listoftT2) == "Severe.COPDvsMild.moderate.COPD"] <- "Severe.COPD - Mild.moderate.COPD"

      names(listoftT)[names(listoftT) == "Mild.moderate.COPDvsControl"] <- "Mild.moderate.COPD - Control"
      names(listoftT)[names(listoftT) == "Severe.COPDvsControl"] <- "Severe.COPD - Control"
      names(listoftT)[names(listoftT) == "Severe.COPDvsMild.moderate.COPD"] <- "Severe.COPD - Mild.moderate.COPD"

      list(
        TF_results = TF_results,
        listoftT = listoftT,
        listoftT2 = listoftT2,
        clinical = sherlock3_biopsy_clinical

      )
    }

    else if(input$study4 == "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs"){
      Sherlock1_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "sherlock1", "results", "listoftT2.rds"))
      Sherlock3_tT2 <- readRDS(file.path(output.dir, "tf_analysis", "brush","results", "listoftT2.rds"))

      list(
        Sherlock1_tT2 = Sherlock1_tT2,
        Sherlock3_tT2 = Sherlock3_tT2
      )

    }



  }) #close reactive study_data4

  ### output$volcano4 Render volcano ---------------------------------------------------------------------------------
  output$volcano4 <- renderPlot({

    validate(
      need(input$study4 != "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs",
           "Go to Results tab to see overlapping TFs between SHERLOCK2&3 (Brush) and SHERLOCK1")
    )

    withProgress(message = "Loading plot ...", value = 0.5, {

      c <-input$comparison4 #these files are saved as listoftT$mildmoderate - control instead of listoftT$contrast1

      listoftT <- study_data4()$listoftT
      listoftT2 <- study_data4()$listoftT2

      tT <- as.data.frame(listoftT[[c]])
      tT2 <- as.data.frame(listoftT2[[c]])

      tT$Legend <- ifelse(
        tT$adj.P.Val < 0.05 & tT$logFC > 0, "Upregulated",
        ifelse(
          tT$adj.P.Val < 0.05 & tT$logFC < -0, "Downregulated",
          "Not Significant"))

      tT2$gene_symbol <- row.names(tT2)

      incProgress(0.8)


      volcano <- ggplot(tT, aes(x = logFC, y = -log10(P.Value))) +
        geom_point(aes(color = Legend)) +
        scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"))+
        geom_hline(yintercept =-log10(max(tT2$P.Value)),colour="black", linetype="dashed")+
        geom_text_repel(data = tT2[1:20,],
                        aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines") ) +
        theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                         legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 16)) +
        labs(title = paste0(input$study4,"\n",input$comparison4))



    }) #close Progress message

    print(volcano)
  }) #close renderPlot

  ### output$results 4 Render list of tT2 TFs ---------------------------------------------------------------------------------
  output$table4 <- renderDT({

    withProgress(message = "Loading results ...", value = 0.5, {

      if(input$study4 == "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs"){

        Sherlock1_tT2 <- study_data4()$Sherlock1_tT2
        Sherlock3_tT2 <- study_data4()$Sherlock3_tT2

        #add _sherlock1 or 2 suffix after colnames so can differentiate the studies in the cbind table after
        Sherlock1_tT2 <- lapply(Sherlock1_tT2, function(df) {  #apply function to each element of the list
          colnames(df) <- paste0(colnames(df), "_Sherlock1")
          df #return the modified df
        })

        Sherlock3_tT2 <- lapply(Sherlock3_tT2, function(df) {  #apply function to each element of the list
          colnames(df) <- paste0(colnames(df), "_Sherlock3")
          df
        })

        #get the overlapping TFs between hserlock1 and SHERLOCK2&3
        if(input$comparison4 == "Mild.moderate.COPD - Control"){
          common1 <- intersect(row.names(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]]), row.names(Sherlock3_tT2[["Mild.moderate.COPDvsControl"]]))
          results_table <- cbind(Sherlock1_tT2[["Mild.moderate.COPDvsControl"]][common1,],Sherlock3_tT2[["Mild.moderate.COPDvsControl"]][common1,])
        }

        else if(input$comparison4 == "Severe.COPD - Control"){
          common2 <- intersect(row.names(Sherlock1_tT2[["Severe.COPDvsControl"]]), row.names(Sherlock3_tT2[["Severe.COPDvsControl"]]))
          results_table <- cbind(Sherlock1_tT2[["Severe.COPDvsControl"]][common2,],Sherlock3_tT2[["Severe.COPDvsControl"]][common2,])
        }

        else if(input$comparison4 == "Severe.COPD - Mild.moderate.COPD"){
          common3 <- intersect(row.names(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]]), row.names(Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]]))
          results_table <- cbind(Sherlock1_tT2[["Severe.COPDvsMild.moderate.COPD"]][common3,],Sherlock3_tT2[["Severe.COPDvsMild.moderate.COPD"]][common3,])
        }

        validate(
          need(nrow(results_table) >= 1, "No overlapping TFs found")
        )

      } #close main if statement (if overlap)


      else{
        results_list <- study_data4()$listoftT
        results_table <- results_list[[input$comparison4]]
      }
    })

    #datatable
    datatable(results_table,
              extensions = c('Buttons'),
              filter = list(position = "top",
                            clear = TRUE),
              options = list(pageLength = 25,
                             scrollX = TRUE,
                             scrollY = "1000px",
                             dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel')))

  }) #close render DT


  ### Update selectizeInput for gene4  ---------------------------------------------------------------------------------

  observe({

    if(input$study4 == "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs"){
      updateSelectizeInput(session, "gene4", choices = "not applicable",
                           selected = "not applicable",
                           server = FALSE)
    }

    else{
      listoftT <- study_data4()$listoftT[["Severe.COPD - Control"]] #listoftT and TF_results have same number of TFs, but sherlock1, SHERLOCK2&3 brush and sherlock3 biopsy have different numbers
      updateSelectizeInput(session, "gene4", choices = row.names(listoftT), server = TRUE)
    }

  }) #close observe




  ### output$boxplot4 Render boxplot   ---------------------------------------------------------------------------------
  output$boxplot4 <- renderPlot({

    validate(
      need(input$study4 != "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs",
           "Go to Results tab to see overlapping TFs between SHERLOCK2&3 (Brush) and SHERLOCK1")
    )

    withProgress(message = "Loading plot ...", value = 0.5, {


      TF_results = study_data4()$TF_results
      clinical = study_data4()$clinical


      boxplot_data <- as.data.frame(cbind(tf_activity = t(TF_results[input$gene4,]),
                                          group = clinical$classification))


      colnames(boxplot_data)[1] <- "tf_activity"


      my_comparisons <- combn(unique(clinical$classification), 2, simplify = FALSE)


      x_order <- c("Control","Mild.moderate.COPD","Severe.COPD")


      boxplot_data$tf_activity <- as.numeric(boxplot_data$tf_activity)
      boxplot_data$group <- factor(boxplot_data$group, levels = x_order)


      stat.table <- boxplot_data %>%
        t_test(tf_activity ~ group,
               comparisons = list(c("Severe.COPD", "Control"),
                                  c("Mild.moderate.COPD", "Control"),
                                  c("Severe.COPD", "Mild.moderate.COPD")
               ))


      stat.table<- stat.table %>%
        add_xy_position(x = "group", dodge = 0.8)

      stat.table3 <- stat.table

      #replace stat.table p values with pvalues from diffexp
      listoftT <- study_data4()$listoftT
      stat.table3 <- cbind(stat.table3, resultsname = c("Severe.COPD - Control", "Mild.moderate.COPD - Control", "Severe.COPD - Mild.moderate.COPD"))
      stat.table3[which(stat.table3$resultsname == "Severe.COPD - Control"),"p"] <- listoftT[["Severe.COPD - Control"]][input$gene4, "P.Value"]
      stat.table3[which(stat.table3$resultsname == "Mild.moderate.COPD - Control"),"p"] <- listoftT[["Mild.moderate.COPD - Control"]][input$gene4, "P.Value"]
      stat.table3[which(stat.table3$resultsname == "Severe.COPD - Mild.moderate.COPD"),"p"] <- listoftT[["Severe.COPD - Mild.moderate.COPD"]][input$gene4, "P.Value"]
      stat.table3$p <- signif(as.numeric(stat.table3$p), digits = 4)
      stat.table3$y.position <- max(boxplot_data[,"tf_activity"]) + 0.025*(max(boxplot_data[,"tf_activity"]))
      stat.table3$y.position <- as.numeric(stat.table3$y.position)

      boxplot_theme <- theme(axis.title = element_text(size = 24),
                             axis.text = element_text(size = 24),
                             title = element_text(size = 20),
                             legend.position = "None")


      #boxplot
      boxplotfinal <- ggplot(boxplot_data, aes(
        x = factor(group, level = x_order),
        y = tf_activity,
        group = group)) +

        theme_bw()+

        boxplot_theme +

        geom_boxplot(position = position_dodge(1)) +

        geom_jitter(aes(color = group),
                    alpha = 0.5,
                    size = 2.5,
                    width = 0.2) +

        {if(nrow(stat.table) >0 )
          stat_pvalue_manual(stat.table3,
                             label = "p",
                             tip.length = 0.01,
                             bracket.nudge.y = c(0, 1, 1),
                             size = 6)
        } +
        stat_summary(fun = mean, fill = "red",
                     geom = "point", shape = 21, size =4,
                     show.legend = TRUE) +



        theme(axis.text.x = element_text(size = 18))+
        labs(title = paste(input$study4, "TF Activity:", input$gene4)

        ) +
        ylab (label = "TF Activity Score") +
        xlab (label = "Disease")

    }) #close message


    print(boxplotfinal)

  })



  ### output$downloadboxplot 4  ---------------------------------------------------------------------------------
  output$downloadboxplot4<-downloadHandler(

    filename = function() {
      paste0(input$gene4,"_data.csv")
    },
    content = function(file) {

      validate(
        need(input$study4 != "SHERLOCK2&3 (Brush) & SHERLOCK1 - Commonly Activated TFs",
             "Go to Results tab to see overlapping TFs between SHERLOCK2&3 (Brush) and SHERLOCK1")
      )


      TF_results = study_data4()$TF_results
      clinical = study_data4()$clinical


      boxplot_data <- as.data.frame(cbind(tf_activity = t(TF_results[input$gene4,]),
                                          group = clinical$classification))


      colnames(boxplot_data)[1] <- "tf_activity"



      my_comparisons <- combn(unique(clinical$classification), 2, simplify = FALSE)


      x_order <- c("Control","Mild.moderate.COPD","Severe.COPD")


      boxplot_data$tf_activity <- as.numeric(boxplot_data$tf_activity)
      boxplot_data$group <- factor(boxplot_data$group, levels = x_order)


      write.csv(boxplot_data,file, row.names = T)
    } #close content
  ) #close download handler

  
  # ================================================================================== #
  ## SERVER 5) Cell Decon ===========================================================
  # ================================================================================== #
  
  ### Reactive expression to load the correct data ---------------------------------------------------------------------------------
  study_data5 <- reactive({
    withProgress(message = "Loading data ...", value = 0.5, {
      
      if(input$study5 == "SHERLOCK2&3 - Brush"){
        readRDS(file.path(output.dir, "celldecon", "figures", "celldecon_ciber_brush.rds"))
      }
      
      else if(input$study5 == "SHERLOCK2&3 - Biopsy"){
        readRDS(file.path(output.dir, "celldecon", "figures", "celldecon_ciber_biopt.rds"))
      }
      
    }) #close message
  }) #close reactive
  
  ### output$boxplot5 Render cell decon plot   ---------------------------------------------------------------------------------
  
  output$celldeconplot5 <- renderPlot({
    
    ciber_data <- study_data5()
    
    celldecon_stats <- ciber_data %>%
      group_by(celltype) %>%
      wilcox_test(proportion ~ classification,
                  paired = FALSE,
                  p.adjust.method = "BH")%>%
      add_xy_position(x = "celltype")
    
    maxRows <- by(data = ciber_data, INDICES = ciber_data$celltype, FUN = function(X) {X[which.max(X$proportion),]})
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
    
    
    ggplot(ciber_data,
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
                         size = 3)+
      scale_fill_manual(values=c("Control" = "#00BA38" , "Mild.moderate.COPD" = "#619CFF",
                                 "Severe.COPD" = "#F8766D")) +
      
      xlab("Cell Type") +
      ylab("Proportion")+
      labs(fill = "Classification")+
      
      ggtitle(paste0("Cellular Proportions: ", input$study5))
    
  }) #close renderplot
  
  
} #close shiny server


shinyApp(ui = ui, server = server)



# runApp("shiny/sherlock_app.R", host = "0.0.0.0", port = 3838,launch.browser = FALSE)
# In Chrome: http://localhost:3838/

