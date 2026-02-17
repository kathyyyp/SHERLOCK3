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

#Loaded up clinical_subset from antonia_sherlock1.R (SevereCOPD data)
clinical_severeCOPD <- clinical_subset
row.names(clinical_severeCOPD) <- paste0("A_", row.names(clinical_severeCOPD))
matched_ids <- intersect(row.names(clinical_severeCOPD), clinical_brush$Study.ID)

emphysema_df <- cbind(clinical_severeCOPD[matched_ids, c("age","sex", "fev1fvc", "es950", "upperlobes_es950", "lowerlobes_es950")], 
                      clinical_brush[match( matched_ids, clinical_brush$Study.ID),c(3:4,158, 277:291,374:381)]) 
write.csv(emphysema_df, file = file.path(output.dir, "matching_emphysema_variables_attempt.csv"))          
#17/02/25
# Cannot match up SevereCOPD emphysema variables with SHERLOCK database. Eg. "es950" for patient A_1002 is "40" but there is nothing matching this number in SHERLOCK database    


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

clinical123_unique$upperlobe_emphysema_perc <-  (clinical123_unique$RUL_emphysema_perc+clinical123_unique$LUL_emphysema_perc)/2
clinical123_unique$lowerlobe_emphysema_perc <-  (clinical123_unique$RLL_emphysema_perc+clinical123_unique$LLL_emphysema_perc)/2
clinical123_unique$upperoverlowerlobe_emphysema_perc <-  clinical123_unique$upperlobe_emphysema_perc/clinical123_unique$lowerlobe_emphysema_perc

clinical123_unique$upperlobe_insp_LAA950_perc <- (clinical123_unique$RUL_insp_LAA_.910HU_perc+clinical123_unique$LUL_insp_LAA_.910HU_perc)/2
clinical123_unique$lowerlobe_insp_LAA950_perc <- (clinical123_unique$RLL_insp_LAA_.910HU_perc+clinical123_unique$LLL_insp_LAA_.910HU_perc)/2
clinical123_unique$upperoverlowerlobe_insp_LAA950_perc <- clinical123_unique$upperlobe_insp_LAA950_perc/clinical123_unique$lowerlobe_insp_LAA950_perc 



#Sanity check, mean emphysema perc in each group
t(clinical123_unique[!is.na(clinical123_unique$serpina1_expanded_genotypes),] %>% 
    group_by(serpina1_expanded_genotypes)%>% 

    summarise(
      
      upperlobe_emphysema_perc = sum(!is.na(upperlobe_emphysema_perc)),
      lowerlobe_emphysema_perc = sum(!is.na(lowerlobe_emphysema_perc)),

      upperlobe_insp_LAA950_perc = sum(!is.na(upperlobe_insp_LAA950_perc)),
      lowerlobe_insp_LAA950_perc = sum(!is.na(lowerlobe_insp_LAA950_perc))
      
    ) )



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


boxplot_theme <- theme(axis.title = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       title = element_text(size = 18),
                       legend.position = "bottom") 

figures.dir <- file.path(output.dir, "figures")
if(!exists(figures.dir)) dir.create(figures.dir)


for (i in colnames(plot)[1:6]){
  
  plot2 <- plot[,c(i,
                    "genotype",
                    "classification",
                    "Study.ID")]
  
  
  plot2 <- plot2[!is.na(plot[,i]),] 
  plot2 <- plot2[!is.infinite(plot2[,i]),] 
  
  
  
  colnames(plot2)[1] <- "emphysema_variable"
  
  if(i == "lowerlobe_emphysema_perc"){ #The Znull is missing lowerlobe_emphysema_perc
    plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "ZI","LOF"))
  }else{
  plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "Znull", "ZI","LOF"))}

  
  ## Get P-Values --------------------------------------------------------------------------------------
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
# 2. SERPINA1 GENE EXPRESSION VS EMPHYSEMA ================================================
# ================================================================================== #


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


boxplot_theme <- theme(axis.title = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       title = element_text(size = 18),
                       legend.position = "bottom") 

figures.dir <- file.path(output.dir, "figures")
if(!exists(figures.dir)) dir.create(figures.dir)


for (i in colnames(plot)[1:6]){
  
  plot2 <- plot[,c(i,
                   "genotype",
                   "classification",
                   "Study.ID")]
  
  
  plot2 <- plot2[!is.na(plot[,i]),] 
  plot2 <- plot2[!is.infinite(plot2[,i]),] 
  
  
  
  colnames(plot2)[1] <- "emphysema_variable"
  
  if(i == "lowerlobe_emphysema_perc"){ #The Znull is missing lowerlobe_emphysema_perc
    plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "ZI","LOF"))
  }else{
    plot2$genotype <- factor(plot2$genotype, levels = c("MM",  "MZ", "MS", "MI", "Mnull", "MP393S", "SS", "ZZ", "Znull", "ZI","LOF"))}
  
  
  ## Get P-Values --------------------------------------------------------------------------------------
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