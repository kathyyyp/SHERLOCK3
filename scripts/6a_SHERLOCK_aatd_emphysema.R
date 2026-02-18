# SHERLOCK analysis of genotyping data (specifically SERPINA1 mutations) in relation to emphysema and gene expression
# Aim: comparing emphysema in SERPINA1 different genotypes

options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED

# ================================================================================== #
## A. SCRIPT SET UP ===============================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("readxl")
library("ggplot2")
library("DESeq2")
library("ggrepel")
library("ggfortify")
library("stringr")
library("tidyverse")
library("ggpubr")
library("rstatix")
library("GSVA")

# ================================================================================== #
## B. SET UP DIRECTORY & OUTPUT PATHS =============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir,"data")

processed.data.dir <- file.path(data.dir,"processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir, "output","aatd", "aatd_emphysema")
if(!exists(output.dir))dir.create(output.dir, recursive = TRUE)

setwd(file.path(main.dir))

# ================================================================================== #
## 1. LOAD IN DATA ===============================================================
# ================================================================================== #

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
clinical123_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
colnames(clinical123_master)[which(colnames(clinical123_master) == "serpina1_z_snp_GRCh38_ch14_pos94378610")] <- "serpina1_z_snp"

clinical123_master[which(clinical123_master$serpina1_genotype_group == "SS"), "serpina1_genotype_group"] <- NA


to_new_col <- which(clinical123_master$serpina1_expanded_genotypes %in% c("MM", "MZ", "MS", "ZZ"))

#new column with cleaner genotypes (excluding MI, MNull, MP393S, ZI, SS)
clinical123_master[to_new_col, "serpina1_genotype_group_clean"] <- clinical123_master[to_new_col,"serpina1_expanded_genotypes"]


# Brush SHERLOCK1, 2 and 3 counts
counts123_brush <- read.csv(file.path(combat.processed.data.dir, "sherlock1_2_3_counts_brush_integrated.csv"), check.names = FALSE, row.names = 1) #from SHERLOCK_SOP_integration.R script (copied from  "/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/sherlock1_2_3_combat)
clinical_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush"),] 

# Biopsy SHERLOCK 2 and 3 counts
counts23_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
clinical_biopt <- clinical123_master[which(clinical123_master$sampletype == "Biopt"),]
counts23_biopt <- counts23_biopt[,row.names(clinical_biopt)]


# ================================================================================== #
# 2. SERPINA1 GENOTYPE VS EMPHYSEMA ================================================
# ================================================================================== #

# #Loaded up clinical_subset from antonia_sherlock1.R (SevereCOPD data)
# clinical_severeCOPD <- clinical_subset
# row.names(clinical_severeCOPD) <- paste0("A_", row.names(clinical_severeCOPD))
# matched_ids <- intersect(row.names(clinical_severeCOPD), clinical_brush$Study.ID)
# 
# emphysema_df <- cbind(clinical_severeCOPD[matched_ids, c("age","sex", "fev1fvc", "es950", "upperlobes_es950", "lowerlobes_es950")], 
#                       clinical_brush[match( matched_ids, clinical_brush$Study.ID),c(3:4,158, 277:291,374:381)]) 
# write.csv(emphysema_df, file = file.path(output.dir, "matching_emphysema_variables_attempt.csv"))          
# 17/02/26 Cannot match up SevereCOPD emphysema variables with SHERLOCK database. Eg. "es950" for patient A_1002 is "40" but there is nothing matching this number in SHERLOCK database    
# 18/02 Daan said to just use the variables we have. The severeCOPD sohort was measured with different method than the SHERLOCK dataset so can't figure out how they got their values "es950", "lowerlobe_es950" and "upperlobe_es950"

#Sanity check, mean emphysema perc in each group
t(clinical_brush %>% group_by(classification)%>% 
  summarise(
    
    Lungs_emphysema_perc = mean(Lungs_emphysema_perc, na.rm = TRUE),
    RUL_emphysema_perc = mean(RUL_emphysema_perc, na.rm = TRUE),
    RLL_emphysema_perc = mean(RLL_emphysema_perc, na.rm = TRUE),
  
    LUL_emphysema_perc = mean(LUL_emphysema_perc, na.rm = TRUE),
    LLL_emphysema_perc = mean(LLL_emphysema_perc, na.rm = TRUE),
    
    RUL_insp_LAA_.910HU_perc = mean(RUL_insp_LAA_.910HU_perc, na.rm = TRUE),
    RLL_insp_LAA_.910HU_perc = mean(RLL_insp_LAA_.910HU_perc, na.rm = TRUE),
    
    LUL_insp_LAA_.910HU_perc = mean(LUL_insp_LAA_.910HU_perc, na.rm = TRUE),
    LLL_insp_LAA_.910HU_perc = mean(LLL_insp_LAA_.910HU_perc, na.rm = TRUE)
    
  ) )

# classification           "Control"   "Mild-moderate COPD" "Severe COPD"
# Lungs_emphysema_perc     " 4.331688" "10.816088"          "26.963802"
# RUL_emphysema_perc       " 3.772633" "11.153599"          "27.496931"
# RLL_emphysema_perc       " 3.511897" " 9.011973"          "26.295175"
# LUL_emphysema_perc       " 5.182521" "11.797256"          "25.907441"
# LLL_emphysema_perc       " 3.727967" " 9.390359"          "25.143759"
# RUL_insp_LAA_.910HU_perc "31.73705"  "37.85742"           "66.03790"
# RLL_insp_LAA_.910HU_perc "29.39791"  "33.60580"           "62.95580"
# LUL_insp_LAA_.910HU_perc "34.74721"  "39.50909"           "63.15542"
# LLL_insp_LAA_.910HU_perc "28.65119"  "33.21299"           "61.47269"



#clinical123_master contains all samples (so same patient is in there twice, brush& biopsy)
clinical123_unique <- clinical123_master %>% distinct(Study.ID, .keep_all = TRUE)

#Remove NAs and calculate upper, lower and upper/lower (cant take average if only upper or lower is available)
clinical123_unique$upperlobe_emphysema_perc <-  ifelse(is.na(clinical123_unique$RUL_emphysema_perc) | is.na(clinical123_unique$LUL_emphysema_perc),
                                                      NA,
                                                      (clinical123_unique$RUL_emphysema_perc+clinical123_unique$LUL_emphysema_perc)/2)

clinical123_unique$lowerlobe_emphysema_perc <-  ifelse(is.na(clinical123_unique$RLL_emphysema_perc) |is.na(clinical123_unique$LLL_emphysema_perc),
                                                      NA,
                                                      (clinical123_unique$RLL_emphysema_perc+clinical123_unique$LLL_emphysema_perc)/2)


clinical123_unique$upperoverlowerlobe_emphysema_perc <-ifelse(is.na(clinical123_unique$upperlobe_emphysema_perc) | is.na(clinical123_unique$lowerlobe_emphysema_perc),
                                                             NA,
                                                             clinical123_unique$upperlobe_emphysema_perc/clinical123_unique$lowerlobe_emphysema_perc)





#Remove NAs and calculate upper, lower and upper/lower (cant take average if only upper or lower is available)
clinical123_unique$upperlobe_insp_LAA950_perc <-  ifelse(is.na(clinical123_unique$RUL_insp_LAA_.910HU_perc) | is.na(clinical123_unique$LUL_insp_LAA_.910HU_perc),
                                                        NA,
                                                        (clinical123_unique$RUL_insp_LAA_.910HU_perc+clinical123_unique$LUL_insp_LAA_.910HU_perc)/2)

clinical123_unique$lowerlobe_insp_LAA950_perc <-  ifelse(is.na(clinical123_unique$RLL_insp_LAA_.910HU_perc) | is.na(clinical123_unique$LLL_insp_LAA_.910HU_perc),
                                                        NA,
                                                        (clinical123_unique$RLL_insp_LAA_.910HU_perc+clinical123_unique$LLL_insp_LAA_.910HU_perc)/2)


clinical123_unique$upperoverlowerlobe_insp_LAA950_perc <-ifelse(is.na(clinical123_unique$upperlobe_insp_LAA950_perc) | is.na(clinical123_unique$lowerlobe_insp_LAA950_perc),
                                                               NA,
                                                               clinical123_unique$upperlobe_insp_LAA950_perc/clinical123_unique$lowerlobe_insp_LAA950_perc)


#Sanity check, how many samples for each genotype have emphysema values
sample_num <- t(clinical123_unique[!is.na(clinical123_unique$serpina1_expanded_genotypes),] %>% 
    group_by(serpina1_expanded_genotypes)%>% 

    summarise(
      
      upperlobe_emphysema_perc = sum(!is.na(upperlobe_emphysema_perc)),
      lowerlobe_emphysema_perc = sum(!is.na(lowerlobe_emphysema_perc)),
      upperoverlowerlobe_emphysema_perc = sum(!is.na(upperoverlowerlobe_emphysema_perc)),
      
      upperlobe_insp_LAA950_perc = sum(!is.na(upperlobe_insp_LAA950_perc)),
      lowerlobe_insp_LAA950_perc = sum(!is.na(lowerlobe_insp_LAA950_perc)),
      upperoverlowerlobe_insp_LAA950_perc = sum(!is.na(upperoverlowerlobe_insp_LAA950_perc))
      
      
    ) )

write.csv(sample_num, file.path(output.dir, "samples_with_emphysema_data.csv"))

clinical123_unique$serpina1_LOF <- ifelse(clinical123_unique$serpina1_expanded_genotypes %in% c("ZZ", "Znull", "ZI"),
                                          "SERPINA1_LOF",
                                          clinical123_unique$serpina1_expanded_genotypes)
  

plot <- clinical123_unique[!is.na(clinical123_unique$serpina1_expanded_genotypes),]
plot <- plot[,c("upperlobe_emphysema_perc",
                "lowerlobe_emphysema_perc",
                "upperoverlowerlobe_emphysema_perc",
                "upperlobe_insp_LAA950_perc",
                "lowerlobe_insp_LAA950_perc",
                "upperoverlowerlobe_insp_LAA950_perc",
                "serpina1_expanded_genotypes",
                "classification",
                "Study.ID")]

plot$Study.ID <- row.names(plot)
colnames(plot)[which(colnames(plot) %in% "serpina1_expanded_genotypes")] <- "genotype"
plot <- as.data.frame(plot)

plot$genotype <- as.factor(plot$genotype)


# Duplicate the LOF rows

LOF_df <- plot[which(plot$genotype %in%  c("ZZ", "Znull", "ZI")),]
LOF_df$genotype <- "LOF"
plot <- rbind(
  plot,
  LOF_df)


#boxplot
boxplot_theme <- theme(axis.title = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       title = element_text(size = 18),
                       legend.position = "bottom") 

figures.dir <- file.path(output.dir, "figures", "serpina1_genotype_vs_emphysema")
if(!exists(figures.dir)) dir.create(figures.dir)


for (i in colnames(plot)[1:6]){
  
  plot2 <- plot[,c(i,
                    "genotype",
                    "classification",
                    "Study.ID")]
  
  
  plot2 <- plot2[!is.na(plot[,i]),] 
  plot2 <- plot2[!is.infinite(plot2[,i]),] 
  
  colnames(plot2)[1] <- "emphysema_variable"
  
  
  if(i == "lowerlobe_emphysema_perc" | i == "upperoverlowerlobe_emphysema_perc"){ #The Znull is missing lowerlobe_emphysema_perc
    plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "ZI","LOF"))
  }else{
  plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "Znull", "ZI","LOF"))}

  
  plot_data <- plot2 %>%
  select(emphysema_variable, genotype) %>%
  na.omit() %>% 
  group_by(genotype) %>%
  mutate(row_id = row_number()) %>%
  
  ungroup() %>%
  pivot_wider(
    names_from = genotype,
    values_from = emphysema_variable
  )
  
  plot_data <- plot_data[,c(levels(plot2$genotype))]
write.csv(plot_data, file = file.path(prism.files.dir, paste0(i,"_prism_plot_data.csv")))
  
  ### Get P-Values --------------------------------------------------------------------------------------
  # Run dunn_test (should give warning if kruskal-wallis test returns false)
  stat.table <- plot2 %>%
    dunn_test(emphysema_variable ~ genotype,
              p.adjust.method = "bonferroni") #do i need this?? this isnt multiple testing
  
  # Filter to only the comparisons you want (optional)
  # e.g. remove redundant ones involving Non-wildtype vs MZ/ZZ
  stat.table <- stat.table %>%
    filter(p < 0.05) %>%          # only significant ones
    add_xy_position(x = "genotype")   # auto-calculates bracket positions
  
  stat.table$p <- signif(as.numeric(stat.table$p), digits = 4)
  y_max = max(plot2[,"emphysema_variable"])
  stat.table$y.position <- y_max + 0.025*y_max
  y_step = y_max* 0.08
  stat.table$y.position <- y_max + (seq_len(nrow(stat.table)) * y_step)

  
  
boxplotimage <- ggplot(plot2, aes(
  x = as.factor(genotype),
  y = emphysema_variable)) +
  
  theme_bw()+
  
  boxplot_theme +
  
  geom_boxplot(position = position_dodge(1),
               aes(fill = genotype,
                   alpha = 0.5)) +
  
  
  geom_jitter(aes(color = genotype),
              alpha = 0.5,
              size = 2.5,
              width = 0.3) +
  
  # 
  # scale_fill_manual(values=c("MM" = "#00BA38" , "MZ" = "#619CFF",
  #                            "ZZ" = "#F8766D")) +
  # 
  
  # scale_color_manual(values=c("Mild-moderate COPD" = "#E68613" , 
  #                             "Severe COPD" = "#C77CFF")) +
  
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
    theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))+
  
    labs(title = paste(i),
         caption = "p = <0.05 shown. Kruskal-Wallis test with Dunn’s Multiple comparisons test.") +
    ylab (label =  paste(i)) +
    xlab (label = paste("SERPINA1 Genotype")) +
    
    # stat_compare_means(
    #   method = "anova",
    #   label = "p.signif"
    # )+
    
  stat_pvalue_manual(
    stat.table,
    label = "p",
    # label = "p.adj",   #p.aj.signif gives the stars **
    tip.length = 0.01 ) +
  
    guides(
      fill  = "none",
      alpha = "none"
    ) 
  

ggsave(boxplotimage,
       file = file.path(figures.dir, paste0(i,"_boxplot.png")),
       height = 20, width = 30, units = "cm", dpi = 800)

  print(i)
}




# ================================================================================== #
# 2. SERPINA1 GENOTYPE VS SERPINA1 EXPRESSION ===================================
# ================================================================================== #

clinical123_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush"),]

#Remove NAs and calculate upper, lower and upper/lower (cant take average if only upper or lower is available)
clinical123_brush$upperlobe_emphysema_perc <-  ifelse(is.na(clinical123_brush$RUL_emphysema_perc) | is.na(clinical123_brush$LUL_emphysema_perc),
                                                      NA,
                                                      (clinical123_brush$RUL_emphysema_perc+clinical123_brush$LUL_emphysema_perc)/2)

clinical123_brush$lowerlobe_emphysema_perc <-  ifelse(is.na(clinical123_brush$RLL_emphysema_perc) |is.na(clinical123_brush$LLL_emphysema_perc),
                                                      NA,
                                                      (clinical123_brush$RLL_emphysema_perc+clinical123_brush$LLL_emphysema_perc)/2)


clinical123_brush$upperoverlowerlobe_emphysema_perc <-ifelse(is.na(clinical123_brush$upperlobe_emphysema_perc) | is.na(clinical123_brush$lowerlobe_emphysema_perc),
                                                             NA,
                                                             clinical123_brush$upperlobe_emphysema_perc/clinical123_brush$lowerlobe_emphysema_perc)





#Remove NAs and calculate upper, lower and upper/lower (cant take average if only upper or lower is available)
clinical123_brush$upperlobe_insp_LAA950_perc <-  ifelse(is.na(clinical123_brush$RUL_insp_LAA_.910HU_perc) | is.na(clinical123_brush$LUL_insp_LAA_.910HU_perc),
                                                      NA,
                                                      (clinical123_brush$RUL_insp_LAA_.910HU_perc+clinical123_brush$LUL_insp_LAA_.910HU_perc)/2)

clinical123_brush$lowerlobe_insp_LAA950_perc <-  ifelse(is.na(clinical123_brush$RLL_insp_LAA_.910HU_perc) | is.na(clinical123_brush$LLL_insp_LAA_.910HU_perc),
                                                        NA,
                                                        (clinical123_brush$RLL_insp_LAA_.910HU_perc+clinical123_brush$LLL_insp_LAA_.910HU_perc)/2)


clinical123_brush$upperoverlowerlobe_insp_LAA950_perc <-ifelse(is.na(clinical123_brush$upperlobe_insp_LAA950_perc) | is.na(clinical123_brush$lowerlobe_insp_LAA950_perc),
                                                             NA,
                                                             clinical123_brush$upperlobe_insp_LAA950_perc/clinical123_brush$lowerlobe_insp_LAA950_perc)




clinical123_brush$serpina1_LOF <- ifelse(clinical123_brush$serpina1_expanded_genotypes %in% c("ZZ", "Znull", "ZI"),
                                          "SERPINA1_LOF",
                                          clinical123_brush$serpina1_expanded_genotypes)
  

plot <- clinical123_brush[!is.na(clinical123_brush$serpina1_expanded_genotypes),]
plot <- plot[,c("upperlobe_emphysema_perc",
                "lowerlobe_emphysema_perc",
                "upperoverlowerlobe_emphysema_perc",
                "upperlobe_insp_LAA950_perc",
                "lowerlobe_insp_LAA950_perc",
                "upperoverlowerlobe_insp_LAA950_perc",
                "serpina1_expanded_genotypes",
                "classification",
                "Study.ID")]



row.names(clinical123_brush) == colnames(counts123_brush)

# normalise gene expression values !!!!!
dds <- DESeqDataSetFromMatrix(countData = counts123_brush,
                              colData = clinical123_brush,
                              design = ~ 1) #no design needed, just need to make dds object so i can vst normalise

counts123_brush_vst <-  assay(vst(dds))

#oull out serpina1
serpina1_ensemblid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == "SERPINA1"),"GENEID"]

serpina1_expression <- as.data.frame(counts123_brush_vst[serpina1_ensemblid, ])
row.names(serpina1_expression)
serpina1_expression <- serpina1_expression[row.names(plot),]


plot <- cbind(plot,
              serpina1_expression_vst = serpina1_expression)


plot$Study.ID <- row.names(plot)
colnames(plot)[which(colnames(plot) %in% "serpina1_expanded_genotypes")] <- "genotype"
plot <- as.data.frame(plot)

plot$genotype <- as.factor(plot$genotype)


# Duplicate the LOF rows

LOF_df <- plot[which(plot$genotype %in%  c("ZZ", "Znull", "ZI")),]
LOF_df$genotype <- "LOF"
plot <- rbind(
  plot,
  LOF_df)


#boxplot
boxplot_theme <- theme(axis.title = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       title = element_text(size = 18),
                       legend.position = "bottom") 

figures.dir <- file.path(output.dir, "figures", "serpina1_genotype_vs_expression")
if(!exists(figures.dir)) dir.create(figures.dir)


  plot2 <- plot[,c("serpina1_expression_vst",
                    "genotype",
                    "classification",
                    "Study.ID")]
  
  
  plot2 <- plot2[!is.na(plot[,"serpina1_expression_vst"]),] 

  
  if(i == "lowerlobe_emphysema_perc" | i == "upperoverlowerlobe_emphysema_perc"){ #The Znull is missing lowerlobe_emphysema_perc
    plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "ZI","LOF"))
  }else{
  plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "Znull", "ZI","LOF"))}

  
#Save prism files for anotnia
  
prism.files.dir <- file.path(output.dir, "data_for_prism")
if(!exists(prism.files.dir))dir.create(prism.files.dir)


plot_data <- plot2 %>%
  select(serpina1_expression_vst, genotype) %>%
  na.omit() %>% 
  group_by(genotype) %>%
  mutate(row_id = row_number()) %>%
  
  ungroup() %>%
  pivot_wider(
    names_from = genotype,
    values_from = serpina1_expression_vst
  )

plot_data <- plot_data[,c(levels(plot2$genotype))]
write.csv(plot_data, file = file.path(prism.files.dir, "serpina1_expression_prism_plot_data.csv"))


  ### Get P-Values --------------------------------------------------------------------------------------
  # Run dunn_test (should give warning if kruskal-wallis test returns false)
  stat.table <- plot2 %>%
    dunn_test(serpina1_expression_vst ~ genotype,
              p.adjust.method = "bonferroni") #do i need this?? this isnt multiple testing
  
  # Filter to only the comparisons you want (optional)
  # e.g. remove redundant ones involving Non-wildtype vs MZ/ZZ
  stat.table <- stat.table %>%
    filter(p < 0.05) %>%          # only significant ones
    add_xy_position(x = "genotype")   # auto-calculates bracket positions
  
  stat.table$p <- signif(as.numeric(stat.table$p), digits = 4)
  y_max = max(plot2[,"serpina1_expression_vst"])
  stat.table$y.position <- y_max + 0.025*y_max
  y_step = y_max* 0.08
  stat.table$y.position <- y_max + (seq_len(nrow(stat.table)) * y_step)

  
  
boxplotimage <- ggplot(plot2, aes(
  x = as.factor(genotype),
  y = serpina1_expression_vst)) +
  
  theme_bw()+
  
  boxplot_theme +
  
  geom_boxplot(position = position_dodge(1),
               aes(fill = genotype,
                   alpha = 0.5)) +
  
  
  geom_jitter(aes(color = genotype),
              alpha = 0.5,
              size = 2.5,
              width = 0.3) +
  
  # 
  # scale_fill_manual(values=c("MM" = "#00BA38" , "MZ" = "#619CFF",
  #                            "ZZ" = "#F8766D")) +
  # 
  
  # scale_color_manual(values=c("Mild-moderate COPD" = "#E68613" , 
  #                             "Severe COPD" = "#C77CFF")) +
  
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
    theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))+
  
    labs(title = paste("SERPINA1 Expression"),
         caption = "p = <0.05 shown. Kruskal-Wallis test with Dunn’s Multiple comparisons test.",
         "Gene expression counts were vst normalised") +
    ylab (label =  paste("SERPINA1 Expression")) +
    xlab (label = paste("SERPINA1 Genotype")) +
    
    # stat_compare_means(
    #   method = "anova",
    #   label = "p.signif"
    # )+
    
  stat_pvalue_manual(
    stat.table,
    label = "p",
    # label = "p.adj",   #p.aj.signif gives the stars **
    tip.length = 0.01 ) +
  
    guides(
      fill  = "none",
      alpha = "none"
    ) 
  

ggsave(boxplotimage,
       file = file.path(figures.dir, paste0("serpina1_expression","_boxplot.png")),
       height = 20, width = 30, units = "cm", dpi = 800)








# ================================================================================== #
# 3. SERPINA1 EXPRESSION VS EMPHYSEMA ================================================
# ================================================================================== #

clinical123_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush"),]

#Remove NAs and calculate upper, lower and upper/lower (cant take average if only upper or lower is available)
clinical123_brush$upperlobe_emphysema_perc <-  ifelse(is.na(clinical123_brush$RUL_emphysema_perc) | is.na(clinical123_brush$LUL_emphysema_perc),
                                                      NA,
                                                      (clinical123_brush$RUL_emphysema_perc+clinical123_brush$LUL_emphysema_perc)/2)

clinical123_brush$lowerlobe_emphysema_perc <-  ifelse(is.na(clinical123_brush$RLL_emphysema_perc) |is.na(clinical123_brush$LLL_emphysema_perc),
                                                      NA,
                                                      (clinical123_brush$RLL_emphysema_perc+clinical123_brush$LLL_emphysema_perc)/2)


clinical123_brush$upperoverlowerlobe_emphysema_perc <-ifelse(is.na(clinical123_brush$upperlobe_emphysema_perc) | is.na(clinical123_brush$lowerlobe_emphysema_perc),
                                                             NA,
                                                             clinical123_brush$upperlobe_emphysema_perc/clinical123_brush$lowerlobe_emphysema_perc)





#Remove NAs and calculate upper, lower and upper/lower (cant take average if only upper or lower is available)
clinical123_brush$upperlobe_insp_LAA950_perc <-  ifelse(is.na(clinical123_brush$RUL_insp_LAA_.910HU_perc) | is.na(clinical123_brush$LUL_insp_LAA_.910HU_perc),
                                                      NA,
                                                      (clinical123_brush$RUL_insp_LAA_.910HU_perc+clinical123_brush$LUL_insp_LAA_.910HU_perc)/2)

clinical123_brush$lowerlobe_insp_LAA950_perc <-  ifelse(is.na(clinical123_brush$RLL_insp_LAA_.910HU_perc) | is.na(clinical123_brush$LLL_insp_LAA_.910HU_perc),
                                                        NA,
                                                        (clinical123_brush$RLL_insp_LAA_.910HU_perc+clinical123_brush$LLL_insp_LAA_.910HU_perc)/2)


clinical123_brush$upperoverlowerlobe_insp_LAA950_perc <-ifelse(is.na(clinical123_brush$upperlobe_insp_LAA950_perc) | is.na(clinical123_brush$lowerlobe_insp_LAA950_perc),
                                                             NA,
                                                             clinical123_brush$upperlobe_insp_LAA950_perc/clinical123_brush$lowerlobe_insp_LAA950_perc)


plot <- clinical123_brush[,c("upperlobe_emphysema_perc",
                "lowerlobe_emphysema_perc",
                "upperoverlowerlobe_emphysema_perc",
                "upperlobe_insp_LAA950_perc",
                "lowerlobe_insp_LAA950_perc",
                "upperoverlowerlobe_insp_LAA950_perc",
                "classification",
                "Study.ID")]

row.names(clinical123_brush) == colnames(counts123_brush)


# normalise gene expression values !!!!!
dds <- DESeqDataSetFromMatrix(countData = counts123_brush,
                              colData = clinical123_brush,
                              design = ~ 1) #no design needed, just need to make dds object so i can vst normalise

counts123_brush_vst <-  assay(vst(dds))

#oull out serpina1
serpina1_ensemblid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == "SERPINA1"),"GENEID"]

serpina1_expression <- as.data.frame(counts123_brush_vst[serpina1_ensemblid, ])
row.names(serpina1_expression)
serpina1_expression <- serpina1_expression[row.names(plot),]


plot <- cbind(plot,
              serpina1_expression_vst = serpina1_expression)

plot$Study.ID <- row.names(plot)
plot <- as.data.frame(plot)


boxplot_theme <- theme(axis.title = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       title = element_text(size = 18),
                       legend.position = "bottom") 

figures.dir <- file.path(output.dir, "figures", "serpina1_expression_vs_emphysema")
if(!exists(figures.dir)) dir.create(figures.dir)


for (i in colnames(plot)[1:6]){
  
  plot2 <- plot[,c(i,
                   "serpina1_expression_vst",
                   "classification",
                   "Study.ID")]
  
  
  plot2 <- plot2[!is.na(plot[,i]),] 
  plot2 <- plot2[!is.infinite(plot2[,i]),] 
  
  
  
  colnames(plot2)[1] <- "emphysema_variable"
  

  
  
  scatterplotimage <- ggplot(plot2, aes(
    x = serpina1_expression_vst,
    y = emphysema_variable)) +
    
    theme_bw()+
    
    boxplot_theme +
    
    geom_point(aes(colour = classification)) +
    
    
    stat_cor(method = "spearman")  +
    
    scale_colour_manual(values=c("Control" = "#00BA38" , "Mild-moderate COPD" = "#619CFF",
                               "Severe COPD" = "#F8766D")) +
    
    
    labs(title = paste(i),
         caption = "Gene expression counts were vst normalised. Spearman Correlation.") +
    ylab (label =  paste(i)) +
    xlab (label = paste("SERPINA1 Expression")) #+
    
    # stat_compare_means(
    #   method = "anova",
    #   label = "p.signif"
    # )+
    
    # scale_x_continuous(labels = scales::label_number(accuracy = 0.01))
    
  
  ggsave(scatterplotimage,
         file = file.path(figures.dir, paste0(i,"_scatterplot.png")),
         height = 20, width = 20, units = "cm", dpi = 800)
  
  print(i)
}




# ================================================================================== #
# 4. PATIENT DEMOGRAPHICS for paper ================================================
# ================================================================================== #

patient_demographics <- t(clinical123_unique %>% 
                            group_by(classification) %>% 
                            summarise(
                              
                              #Total patients
                              total_patients = n(),
                              
                              #sex
                              male_patients = sum(sex == "Male"), 
                              male_percentage = (male_patients / total_patients) * 100,
                              sex = paste0(male_patients, "(", round(male_percentage, digits=3), ")"),
                              
                              #age
                              mean_age = mean(age,  na.rm = TRUE),
                              sd_age = sd(age, na.rm = TRUE),
                              age = paste0(round(mean_age,digits =3),"(",round(sd_age, digits =3),")"),
                              
                              
                              #packyears
                              mean_packyears = mean(packyears, na.rm = TRUE),
                              sd_packyears = sd(packyears, na.rm = T),
                              packyears = paste0(round(mean_packyears,digits =3), "(",round(sd_packyears,digits=3), ")"),
                              
                              
                              #FEV1
                              mean_FEV1_pred =  mean(FEV1_pred, na.rm = T),
                              sd_FEV1_pred = sd(FEV1_pred, na.rm = T),
                              fev1 = paste0(round(mean_FEV1_pred, digits = 3), "(", round(sd_FEV1_pred,digits=3),  ")"),
                              
                              
                              #FEV1 perc pred
                              mean_FEV1_percent_pred =  mean(FEV1_percent_pred, na.rm = T),
                              sd_FEV1_percent_pred = sd(FEV1_percent_pred, na.rm = T),
                              fev1_percent_pred = paste0(round(mean_FEV1_percent_pred, digits = 3), "(", round(sd_FEV1_percent_pred,digits=3),  ")"),
                              
                              
                              #FEV1FVC
                              mean_FEV1FVC =  mean(FEV1_FVC_post, na.rm = T),
                              sd_FEV1FVC = sd(FEV1_FVC_post, na.rm = T),
                              fev1fvc = paste0(round(mean_FEV1FVC, digits = 3), "(", round(sd_FEV1FVC,digits=3),  ")"),
                              
                    
                              #upperoverlowerlobe_insp_LAA950_perc
                              mean_upperoverlowerlobe_insp_LAA950_perc =  mean(upperoverlowerlobe_insp_LAA950_perc, na.rm = T),
                              sd_upperoverlowerlobe_insp_LAA950_perc = sd(upperoverlowerlobe_insp_LAA950_perc, na.rm = T),
                              upperoverlowerlobe_insp_LAA950_perc = paste0(round(mean_upperoverlowerlobe_insp_LAA950_perc, digits = 3), "(", round(sd_upperoverlowerlobe_insp_LAA950_perc,digits=3),  ")"),
                              
                              #upperlobe_insp_LAA950_perc
                              mean_upperlobe_insp_LAA950_perc =  mean(upperlobe_insp_LAA950_perc, na.rm = T),
                              sd_upperlobe_insp_LAA950_perc = sd(upperlobe_insp_LAA950_perc, na.rm = T),
                              upperlobe_insp_LAA950_perc = paste0(round(mean_upperlobe_insp_LAA950_perc, digits = 3), "(", round(sd_upperlobe_insp_LAA950_perc,digits=3),  ")"),
                              
                              #lowerlobe_insp_LAA950_perc
                              mean_lowerlobe_insp_LAA950_perc =  mean(lowerlobe_insp_LAA950_perc, na.rm = T),
                              sd_lowerlobe_insp_LAA950_perc = sd(lowerlobe_insp_LAA950_perc, na.rm = T),
                              lowerlobe_insp_LAA950_perc = paste0(round(mean_lowerlobe_insp_LAA950_perc, digits = 3), "(", round(sd_lowerlobe_insp_LAA950_perc,digits=3),  ")"),
                              
                              
                            )
                          %>%
                            dplyr::select(
                              classification,
                              total_patients,
                              sex,
                              age,
                              packyears,
                              fev1,
                              fev1_percent_pred,
                              fev1fvc,
                              upperoverlowerlobe_insp_LAA950_perc,
                              upperlobe_insp_LAA950_perc,
                              lowerlobe_insp_LAA950_perc#,
                              # aatd  #check sherlock table for this
                            )
)




colnames(patient_demographics) <- patient_demographics[1,]
patient_demographics <- patient_demographics[-1,]

row.names(patient_demographics) <- c(
  "Samples, n",
  "Gender, male n(%)",
  "Age, mean (SD)",
  "Pack years, mean (SD)",
  "FEV1, mean (SD)",
  "FEV1 (% predicted), mean (SD)",
  "FEV1/FVC, mean (SD)",
  "Upper/Lower Lobes insp %LAA-950, % mean (SD)",
  "Upper Lobes %LAA-950, % mean (SD)",
  "Lower Lobes %LAA-950, % mean (SD)")



# Save 
write.csv(patient_demographics, file = file.path(file.path(output.dir ,"patient_demographics.csv")))


#Use chi squared test for categorical variables (how to do this for 3 variables??????)

table_to_insert <-table(clinical123_unique$sex, clinical123_unique$classification)
sex_chi <- chisq.test(table_to_insert)$p.value


#Use t.test or wilcoxon test for continuous variables (shapiro.test shows all p<0.05, therefore not normally distributed?? but paper used T-test)
continuous_variables <- c("age", 
                          "packyears", 
                          "FEV1_pred",
                          "FEV1_percent_pred", 
                          "FEV1_FVC_post", 
                          "upperlobe_insp_LAA950_perc",
                          "lowerlobe_insp_LAA950_perc",
                          "upperoverlowerlobe_insp_LAA950_perc")




stats_kruskal<- data.frame()
for (i in continuous_variables){
  formula <- as.formula(paste(i, "~ classification"))
  test_res <- kruskal.test(formula = formula, data = clinical123_unique)
  pval <- test_res$p.value
  stats_kruskal <- rbind(stats_kruskal, data.frame(variable = i, pval))
}


library(broom)
# 
# stats_aov <- data.frame()
# for (i in continuous_variables){
#   formula <- as.formula(paste(i, "~ classification"))
#   clinical123_unique_clean <- clinical123_unique[!is.na(clinical123_unique[,i]), ]
#   test_res <- aov(formula = formula, data = clinical123_unique_clean)
#   pval <- tidy(test_res)$p.value[1]
#   stats_aov <- rbind(stats_aov, data.frame(variable = i, pval))
# }

stats_col_kruskal <- rbind("classification", data.frame(variable = "sex", pval = sex_chi), stats_kruskal)
# stats_col_aov <- rbind("classification", data.frame(variable = "sex", pval = sex_chi), stats_aov)



patient_demographics <-cbind(patient_demographics,
                             # pval_aov = stat_col_wilcox$pval,
                             pval_kruskal = stats_col_kruskal$pval
                             )

write.csv(patient_demographics, file = file.path(file.path(output.dir, "patient_demographics_withstats.csv")))

#manual tests
# wilcox_res <- c(
# wilcox.test(formula = age ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = packyears ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = fev1_pred ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = fev1fvc ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = es950 ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = upperlobes_es950 ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = lowerlobes_es950 ~ group, data = clinical_subset, paired = FALSE)$p.value,
# wilcox.test(formula = aatd ~ group, data = clinical_subset, paired = FALSE)$p.value
# )




# 
# # ================================================================================== #
# # 3. PRISM: SERPINA1 GENOTYPE VS EMPHYSEMA (FOR ANTONIA) ===========================
# # ================================================================================== #
# plot <- clinical123_unique[!is.na(clinical123_unique$serpina1_expanded_genotypes),]
# plot <- plot[,c("upperlobe_emphysema_perc",
#                 "lowerlobe_emphysema_perc",
#                 "upperoverlowerlobe_emphysema_perc",
#                 "upperlobe_insp_LAA950_perc",
#                 "lowerlobe_insp_LAA950_perc",
#                 "upperoverlowerlobe_insp_LAA950_perc",
#                 "serpina1_expanded_genotypes")]
# 
# plot$Study.ID <- row.names(plot)
# colnames(plot)[which(colnames(plot) %in% "serpina1_expanded_genotypes")] <- "genotype"
# plot <- as.data.frame(plot)
# 
# plot$genotype <- as.factor(plot$genotype)
# 
# 
# # Duplicate the LOF rows
# 
# LOF_df <- plot[which(plot$genotype %in%  c("ZZ", "Znull", "ZI")),]
# LOF_df$genotype <- "LOF"
# plot <- rbind(
#   plot,
#   LOF_df)
# 
# 
# prism.files.dir <- file.path(output.dir, "data_for_prism")
# if(!exists(prism.files.dir))dir.create(prism.files.dir)
# for (i in colnames(plot)[1:6]){
#  
# plot_data <- plot %>%
#   select(i, genotype) %>%
#   na.omit() %>% 
#   group_by(genotype) %>%
#   mutate(row_id = row_number()) %>%
#   
#   ungroup() %>%
#   pivot_wider(
#     names_from = genotype,
#     values_from = i
#   )
# write.csv(plot_data, file = file.path(prism.files.dir, paste0(i,"_prism_plot_data.csv")))
# }
# 
# # ================================================================================== #
# # 4. PRISM: SERPINA1 EXPRESSION VS EMPHYSEMA (FOR ANTONIA) =========================
# # ================================================================================== #
# 
# ## serpina1 expression 
# plot <- clinical123_brush[,c("upperlobe_emphysema_perc",
#                 "lowerlobe_emphysema_perc",
#                 "upperoverlowerlobe_emphysema_perc",
#                 "upperlobe_insp_LAA950_perc",
#                 "lowerlobe_insp_LAA950_perc",
#                 "upperoverlowerlobe_insp_LAA950_perc",
#                 "classification",
#                 "Study.ID")]
# 
# row.names(clinical123_brush) == colnames(counts123_brush)
# 
# #normalise gene epression
# dds <- DESeqDataSetFromMatrix(countData = counts123_brush,
#                               colData = clinical123_brush,
#                               design = ~ 1) #no design needed, just need to make dds object so i can vst normalise
# 
# counts123_brush_vst <-  assay(vst(dds))
# 
# #pull out serpina1
# serpina1_ensemblid <- hgnc_symbols_db[which(hgnc_symbols_db$SYMBOL == "SERPINA1"),"GENEID"]
# 
# serpina1_expression <- as.data.frame(counts123_brush_vst[serpina1_ensemblid, ])
# row.names(serpina1_expression)
# serpina1_expression <- serpina1_expression[row.names(plot),]
# 
# 
# plot <- cbind(plot,
#               serpina1_expression_vst = serpina1_expression)
# 
# plot$Study.ID <- row.names(plot)
# plot <- as.data.frame(plot)
# 
# serpina1_plot_data <- plot %>%
#   select(serpina1_expression, genotype) %>%
#   na.omit() %>% 
#   group_by(genotype) %>%
#   mutate(row_id = row_number()) %>%
#   
#   ungroup() %>%
#   pivot_wider(
#     names_from = genotype,
#     values_from = serpina1_expression
#   )
# write.csv(serpina1_plot_data, file = file.path(output.dir, "serpina1_plot_data.csv"))
# 
# 
