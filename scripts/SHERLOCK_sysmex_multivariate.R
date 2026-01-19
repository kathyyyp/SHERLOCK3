#SHERLOCK3 - Analysis of Sysmex data
#This script was run after SHERLOCK_sysmex_diffexp, 5_SHERLOCK_radiomics and 6_SHERLOCK_aatd. Refer to those for clinical master table generation details

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

#clinical_sherlock123_master.rds - created in 6_SHERLOCK_aat.R script (same as clinical_brushbiopt_master.rds but extra unecessary columns were removed and emphysema values were fixed up (some were missing decimal places before))
clinical123_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))

clinical_brushbiopt <- clinical123_master[which(clinical123_master$batch != 1),]

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))


setwd(file.path(main.dir))


# ================================================================================== #
# 3) Combine clinical file with sysmex data and remove unecessary columns
# ================================================================================== #
# All radiomics varialbes should be numeric
radiomics_index_start <-which(colnames(clinical123_master) =="RL_insp_vol_ml")
radiomics_index_end <- which(colnames(clinical123_master) =="LLL_airtrapping_emphysema_perc")


clinical_brushbiopt[, radiomics_index_start:radiomics_index_end] <- lapply(
  clinical_brushbiopt[, radiomics_index_start:radiomics_index_end],
  function(x) {
    if (is.list(x)) {
      as.numeric(unlist(x))
    } else {
      as.numeric(x)
    }
  }
)

sapply(clinical_brushbiopt[radiomics_index_start:radiomics_index_end],class)


# # All sysmex varialbes should be numeric
sysmex_index_start <- which(colnames(clinical_brushbiopt) == "Mono")
sysmex_index_end <- which(colnames(clinical_brushbiopt) == "EO.Z")
sapply(clinical_brushbiopt[sysmex_index_start:sysmex_index_end],class)

saveRDS(clinical_brushbiopt, file.path(postQC.data.dir,  "master","clinical_brushbiopt_master_202601.rds"))
write.csv(clinical_brushbiopt, file.path(postQC.data.dir,  "master","clinical_brushbiopt_master_202601.csv"))

# ================================================================================== #
# 2) Any association between the Sysmex parameters and clinical variables using 
# nonlinear/multivariate analyses (age, sex, BMI, packyears, lung function, 
# emphysema scores, co-morbidities)
# ================================================================================== #

# Aim #
# lm(predictor ~ outcome)

# Linear models for each sysmex variable and clinical variable
# mono ~ age
# mono ~ sex
# mono ~ smoking_status 
# etc
# rbc ~ age
# rbc ~ sex
# etc


# Make sure data types are coorrect
clinical_brushbiopt$sex <- factor(clinical_brushbiopt$sex)
clinical_brushbiopt$smoking_status <- factor(clinical_brushbiopt$smoking_status)
clinical_brushbiopt$ics_use <- factor(clinical_brushbiopt$ics_use)
clinical_brushbiopt$classification <- factor(clinical_brushbiopt$classification)

numeric_vars <- c("age", 
                  "packyears", 
                  "FEV1_percent_pred", 
                  "FEV1_FVC_post", 
                  # "crf_years_of_cessation_num", #need to convert this to years but dont know when this info was collected (year collected - year of cessation)
                  "RL_insp_LAA_.950HU_perc",
                  "LL_insp_LAA_.950HU_perc",
                  "RL_insp_emphysema_15percentile_HU",
                  "LL_insp_emphysema_15percentile_HU",
                  "Lungs_Pi10", 
                  "Lungs_insp_bronchial_count", 
                  "Lungs_mucus_plugs_count", 
                  "Lungs_emphysema_perc", 
                  "Lungs_airtrapping_perc")


categorical_vars  <- c("sex", "smoking_status", "ics_use", "classification")

sysmex_vars <- colnames(clinical_brushbiopt)[sysmex_index_start:sysmex_index_end]




clinical_brush <- clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Brush"),]
clinical_biopt <- clinical_brushbiopt[which(clinical_brushbiopt$sampletype == "Biopt"),]

#How to interpret lm() results
# Example code
fit <- lm(Mono ~ sex, data = clinical_brush)
coef(summary(fit))

#Results of above code        
# Coefficients:
#               Estimate  Std. Error  t value  Pr(>|t|)
# (Intercept)     0.62       0.05      12.4     <0.001
# sexMale         0.18       0.07       2.6      0.01


# (Intercept) = mean monocyte count in females
# sexMale = different in males compared to females
# So, males have 0.18 higher monocyte count on average compared to females and difference is significant (p = 0.01)
# Note the p < 0.001 is from testing whether the mean monocyte count in the reference group (Females) different from 0? - this is not biologically relevant or useful for us (bc it should be different to 0)



results <- list()

clinical <- clinical_brush

#For every sysmex variable, run lm() against a clinical variable
results <- list()

for (y in sysmex_vars) {
  for (x in c(numeric_vars, categorical_vars)) {
    
    fit <- lm(as.formula(paste(y, "~", x)), data = clinical)
    coefs <- as.data.frame(coef(summary(fit)))
    coefs$term <- rownames(coefs)
    
    # drop intercept
    coefs <- subset(coefs, term != "(Intercept)")
    
    # reference level (only for factors)
    ref <- if (x %in% categorical_vars) {
      fit$xlevels[[x]][1]
    } else {
      NA
    }
    
    results[[paste(y, x, sep = "_")]] <- data.frame(
      outcome   = y,
      predictor = x,
      term      = coefs$term,
      estimate  = coefs$Estimate,
      p_value   = coefs$`Pr(>|t|)`,
      reference = ref
    )
  }
}

results_df <- do.call(rbind, results)
row.names(results_df) <- NULL

results_df$FDR <- p.adjust(results_df$p_value, method = "BH")
sig <- subset(results_df, FDR < 0.05)




for (i in this_sysmex_variable)
#Filter out variables where there are too many zeros (cutoff at >80% zeros)
zero_proportion <- sum(clinical2[,this_sysmex_variable] == 0, na.rm = TRUE)/length(clinical2[,this_sysmex_variable])

if(zero_proportion > 0.8){
  cat(paste0("Skipping ", this_sysmex_variable, " (", round(zero_proportion * 100, 1), "% zeros)"))
  next        
}


# Filter out outliers (3sd from the mean) - these are falsely driving the diffexp results
stdev  <- sd(clinical2[,this_sysmex_variable], na.rm = TRUE)

#if value is 3 standard deviations away from the mean, make it NA
if (!is.na(stdev) && stdev > 0) {
  mean_val <- mean(clinical2[,this_sysmex_variable], na.rm = TRUE)
  outliers <- abs(clinical2[,this_sysmex_variable] - mean_val) > 3 * stdev
  clinical2[outliers, this_sysmex_variable] <- NA
}





sysmex_vars <- clinical_brushbiopt
clinical_vars <- c("age", "sex", "smoking", "packyears", "BMI", "FEV1_percent_pred", "emphysema")

df <- clinical_data   

# a) Visualise correlations/scatterplots
ggplot(data, aes(y = mean_SD, x = Var1))+ 
  geom_point(size = 3)+
  facet_grid(Var3 ~ Var2, scales = "free_x", space = "free")+
  labs(title =" SD across all parameter sets",
       x = "Var1 ", 
       y= "Mean SD")+ 
  ylim(0, 500)







fit <- lm(Mono ~ sex, data = clinical_brush)








# 
# # b) Seperate linear models
# 
# lm_screen <- expand.grid(
# sysmex = sysmex_vars,
# clinical = clinical_vars,
# stringsAsFactors = FALSE
# ) %>%
#   mutate(
#     formula = paste(sysmex, "~", clinical),
#     model = map(formula, ~ lm(as.formula(.x), data = df)),
#     tidy = map(model, tidy)
#   ) %>%
#   unnest(tidy) %>%
#   filter(term != "(Intercept)") %>%
#   select(sysmex, clinical, term, estimate, p.value) %>%
#   group_by(sysmex) %>%
#   mutate(FDR = p.adjust(p.value, method = "BH")) %>%
#   ungroup()




# 
# 
# # c) Multivariate model
# # this is adjustd for other variables?
# multi_models <- map(
#   sysmex_vars,
#   ~ lm(
#     as.formula(
#       paste(.x, "~", paste(clinical_vars, collapse = " + "))
#     ),
#     data = df
#   )
# )
# 
# 
# 
# 
# #extract betas and fdr (per-sysmex fdr)
# multi_results <- map_df(
#   names(multi_models),
#   ~ tidy(multi_models[[.x]]) %>%
#     filter(term != "(Intercept)") %>%
#     mutate(sysmex = .x)
# )
# 
# multi_results <- multi_results %>%
#   group_by(sysmex) %>%
#   mutate(FDR = p.adjust(p.value, method = "BH")) %>%
#   ungroup()
# names(multi_models) <- sysmex_vars
# 
# #global fdr (across all tests)?
# multi_results$FDR_global <- p.adjust(multi_results$p.value, method = "BH")

