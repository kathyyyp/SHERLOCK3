#Shiny app for Sherlock RNA-Seq - Sysmex - Analysis Nov 2025
# runApp("sysmex_app.R", host = "0.0.0.0", port = 3838,launch.browser = FALSE)
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



# ================================================================================== #
# A. SHINY (UI) ==================================================================
# ================================================================================== #

ui <- dashboardPage(
  
  dashboardHeader(title = "SHERLOCK2"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("2) Differential Expression", tabName = "comparison", icon = icon("bar-chart")),
      menuItem("3) Gene Expression", tabName = "gene_expression", icon = icon("bar-chart"))
      )

  ), #closedashboardSidebar
  
  
  dashboardBody(
    title = "SHERLOCK - Sysmex RNA-Seq",
    tabItems(
      
      
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
                          choices = c("SHERLOCK - Brush", "SHERLOCK - Biopsy"),
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
              
              
      )
      
     
              
    ) #close tabItems
  )#close Dashboard body
) #close UI




server <- function(input, output, session){
  
  
  # ================================================================================== #
  ## SERVER 2) DIFFERENTIAL EXPRESSION ===============================================
  # ================================================================================== #
  
  
  ### Reactive expression to load the correct data ---------------------------------------------------------------------------------------------------------------------------------------------------------
  study_data2 <- reactive({
    req(input$study2)  # ensure input is not NULL
    
    results_filepath <- file.path(output.dir, "diffexp_withgroup", "results", "listofresults.RDS")
    listofresults <- readRDS(results_filepath)
    
    if(input$study2 == "SHERLOCK - Brush"){
      
      
      list(
        listoftT = listofresults$tT,
        listoftT2 = listofresults$tT2
      )
    }
    
    else if(input$study2 == "SHERLOCK - Biopsy"){
      
      
      list(
        listoftT =listofresults$brush$tT,
        listoftT2 = listofresults$brush$tT2
      )
    }
 
    
  }) #close reactive study_data2
  
  
  
  ### output$deg2 Render results tT table ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  output$deg2 <- DT::renderDT({
    
    withProgress(message = "Loading data ...", value = 0.5, {
      c <- input$sysmex_variable
      
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
      
      c <- nameconvert[which(nameconvert$contrast == input$sysmex_variable),"name"]
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
        labs(title = paste0(input$study2,"\n",input$sysmex_variable))
    })
    print(volcano)
    
  }) #close renderPlot volcano2
  
}


shinyApp(ui = ui, server = server)
