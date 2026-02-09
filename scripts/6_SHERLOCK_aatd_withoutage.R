# SHERLOCK analysis of genotyping data (sepcifically SERPINA1 mutations) in relation to gene expression
# Aim: comparing gene expression for MZ and ZZ patients (in COPD)

# SETTING FOR HPC BATCH JOB (run script from 2.DIFFERENTIAL EXPRESSION as a batch job, everything before that was run in interactive)
hpc_batch_job = TRUE


# ### BASH ####
# # note the .vcf.gz file needs the accompanying .tbi file
# # View header
# module load BCFtools
# bcftools view -h data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz
# 
# #see one variant
# bcftools view -H data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz | head -n 1
# 
# # Subset chr1
# bcftools view -r chr1:1-5e6 data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz -Oz -o data/processed/snp_array/chr1_subset.vcf.gz
# 
# # -Oz means Output format is bgZipped vcf
# 
# Subset chr14
# bcftools view -r chr14 data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz -Oz -o data/processed/snp_array/chr14_array.vcf.gz
# bcftools view -r chr14:94376747-94388602 data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gz -Oz -o data/processed/snp_array/chr14_serpina1_array.vcf.gz


#Check size
# bcftools stats data/processed/snp_array/chr14_serpina1_array.vcf.gz | grep "^SN"

# should be hg38 ! according to alen
# Subset samples of interest?
# bcftools view -S samples.txt data/raw/snp_array/20250910_RES0225_GSAv3+.imputed.vcf.gzz -Oz -o data/processed/snp_array/samples_subset.vcf.gz

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT 
# chr1    
# 10390   chr1:10390_CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA_C       CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA    
# C       
# .       
# LOWIMPQUAL   
# DR2=0.16;AF=0.0009;INFO=0.163865;MAF=0.0009;AF_EUR=0.00085644;AF_ADMIXED=0.001;AF_EAS=0;AF_AFR=0;IMPUTED;RSID=rs1557426845      
# GT:DS:GP    
# 0|0:0:1,0,0 

# Z rs28929474		
# S rs17580		
# I rs28931570		
# M(Herleen) rs199422209	
# Q0amersfoort rs199422210		

if(hpc_batch_job == FALSE){ #don't need this part for hpc
  
  
  # ================================================================================== #
  # A. SCRIPT SET UP =================================================================
  # ================================================================================== #
  my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"
  
  
  library("readxl")
  library("stringr")
  library("vcfR")
  
  
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
  output.dir <- file.path(output.dir, "aatd")
  if(!exists(output.dir))dir.create(output.dir)
  
  # ================================================================= #
  # 1. LOAD IN DATA ==================================================================
  # ================================================================================== #
  setwd(file.path(main.dir))
  
  hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
  
  #master clinical file - output from 1a_SHERLOCK3_qc_combat.R, after removing outliers
  # this file still has unecessary columns we won't use (like survey questions etc. also emphysema values are weird)
  # these issues were fixed in Sherlock_database_07_25_Final.xlsx so match to that
  clinical_master <- read.csv("/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/clinical_master.csv", row.names = 1)
  
  
  raw_clinical <- read_xlsx(file.path(data.dir,"raw","Sherlock_database_07_25_Final.xlsx")) #598 SEO/patient IDs
  # Note that this file has less columns than previous version, survey columns have been removed
  # Also, this file has one row per patient, clinical_sk_all has one row per sample (same patients can have more than one sample bc brush and biopt)
  # IMPORTANT - radiomics data in this file is updated too. the clinical_brushbiopt_master file has weird emphysema % values 
  
  
  clinical_sk_all <- as.data.frame(raw_clinical[match(clinical_master$Study.ID, raw_clinical$class_incl_study_id),])
  # clinical_sk_all <- clinical_sk_all[match(clinical_master$Study.ID, clinical_sk_all$Study.ID),]
  
  clinical_sk_all <- cbind(clinical_sk_all, sampletype = clinical_master$sampletype, 
                           batch = clinical_master$batch,
                           classification = clinical_master$classification)
  row.names(clinical_sk_all) <- row.names(clinical_master)
  
  clinical_sk_all <- clinical_sk_all %>% 
    dplyr::rename(
      Study.ID = class_incl_study_id,
      age = crf_age,
      sex = crf_gender,
      smoking_status = crf_smoking,
      packyears = crf_packyears,
      corticosteroid = crf_corticosteroid,
      FVC_post= postbodybox_fvc_post ,
      FEV1 = postbodybox_fev1_post,
      FEV1_FVC_post = postbodybox_fev1_fvc_post
    ) 
  
  
  clinical_sk_all[,c("age", "packyears", "FEV1", 
                     "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")] <- sapply(
                       clinical_sk_all[,c("age", "packyears", "FEV1", 
                                          "FEV1_percent_pred", "FEV1_FVC_post", "FVC_post")], 
                       function(x) as.numeric(x))
  
  
  # LIB5426634_SAM24375606 
  #this patient had no data in raw_clinical but is in clinical_sk1_master.rds(generated from mastertable_sherlock3.csv ie. the SHERLOCK1 data alen gave me in 2024 
  #(mastertable_sherlock3.csv does not mean sherlock3 study, just the name of alen's sherlock1 file version 3))
  clinical_sk1_master <- readRDS(file.path(processed.data.dir, "SHERLOCK1", "clinical_sk1_master.rds"))
  clinical_sk1_master["LIB5426634_SAM24375606",]
  clinical_sk1_master <- clinical_sk1_master %>% 
    dplyr::rename(
      sex = crf_gender,
      smoking_status = crf_smoking,
      packyears = crf_packyears,
      corticosteroid = crf_corticosteroid,
      FVC_post= postbodybox_fvc_post ,
      FEV1_FVC_post = postbodybox_fev1_fvc_post
    ) 
  clinical_sk1_master$age <- as.numeric(clinical_sk1_master$age)
  common_cols <- intersect(colnames(clinical_sk1_master), colnames(clinical_sk_all))
  clinical_sk_all["LIB5426634_SAM24375606", common_cols] <- clinical_sk1_master["LIB5426634_SAM24375606", common_cols]
  
  
  ## Get SERPINA1 Z mutation genotypes ------------------------------------------------------------------------------------------------------------
  # Pull out chromosome 14
  vcf_serpina1 <- read.vcfR(
    "data/processed/snp_array/chr14_serpina1_array.vcf.gz",
    verbose = FALSE
  )
  
  head(vcf_serpina1, n= 10)
  
  fix <- as.data.frame(getFIX(vcf_serpina1))
  serpina1_snps_df <- cbind(fix[, c("CHROM","POS", "ID","REF","ALT", "QUAL", "FILTER")], 
                            extract.gt(vcf_serpina1, element = "DS", as.numeric = TRUE)
                            # extract.gt(vcf_serpina1, element = "GT", as.numeric = FALSE)
  ) #626 samples
  
  # Current naming "RES0225_SEO568_GSAv3+". Convert to SEO and A number alone ("SEO568")
  colnames(serpina1_snps_df)[8:ncol(serpina1_snps_df)] <- str_extract(colnames(serpina1_snps_df)[8:ncol(serpina1_snps_df)], "(?<=_)[^_]+(?=_)")
  
  # Pull out Z mutation (ch14 pos94378610 for hg38)
  # This should be rs28929474 (SERPINA1 Pi*Z mutation) but RS IDs aren't labelled
  serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610),] #Z mutation
  # chr14 94378610   ref=C   alt=T 
  
  #Remove the first 7 column (genotyping meta) so we are left with just dosage (ie. genotype as decimal)
  serpina1_z_snp <- as.data.frame(cbind(Study.ID = colnames(serpina1_snps_df)[-c(1:7)], 
                                        z_mutation = as.numeric(serpina1_snps_df[which(serpina1_snps_df$POS ==  94378610), -c(1:7)])))
  
  # standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
  # MM: DS < 0.2
  # MZ: 0.8 ≤ DS ≤ 1.2
  # ZZ: DS ≥ 1.8
  
  #Convert names in clinical file to new names provided by Daan (the A numbers in genotyping file are different to clinical file)
  A_number_conversion_file <- as.data.frame(read_excel(file.path(data.dir, "raw", "A_number_conversion_file.xlsx")))
  A_number_conversion_file$genotyping_ID <- paste0("A", A_number_conversion_file$'A-numbers Genotyping')
  A_number_conversion_file$clinical_ID <- paste0("A", A_number_conversion_file$'A-numbers Sherlock database')
  
  
  
  #Add mutation to clinical file
  
  A_IDs <- serpina1_z_snp[grep("^A", serpina1_z_snp$Study.ID),]
  not_A_IDs <- serpina1_z_snp[!serpina1_z_snp$Study.ID %in% A_IDs$Study.ID,]
  
  
  A_IDs$Study.ID %in% A_number_conversion_file$genotyping_ID
  A_IDs$Study.ID.corrected <- A_number_conversion_file[match(A_IDs$Study.ID, A_number_conversion_file$genotyping_ID), "clinical_ID"]
  
  
  A_IDs$Study.ID.corrected <- sub("^A", "A_", A_IDs$Study.ID.corrected )
  
  
  serpina1_z_snp <- rbind(cbind(Study.ID = A_IDs$Study.ID.corrected, 
                                z_mutation = A_IDs$z_mutation),
                          not_A_IDs)
  
  row.names(serpina1_z_snp) <- serpina1_z_snp$Study.ID #There are 626 patients with genotyping data
  length(unique(serpina1_z_snp$Study.ID))
  
  write.csv(serpina1_z_snp, file.path(processed.data.dir, "snp_array","serpina1_z_snp_dosage.csv"), row.names = FALSE)
  
  
  matching_ids <- intersect(clinical_sk_all$Study.ID, row.names(serpina1_z_snp)) #432 patients that have genotyping data AND gene expresion data. 194 have genotyping but not gene expression
  
  clinical2 <- clinical_sk_all[clinical_sk_all$Study.ID %in% matching_ids,] #(674 samples from 432 patients)
  clinical2 <- cbind(clinical2, serpina1_z_snp_GRCh38_ch14_pos94378610 = serpina1_z_snp[match(clinical2$Study.ID, serpina1_z_snp$Study.ID),
                                                                                        "z_mutation"])
  
  
  
  # standard for definitions based on dosage DS, but the dosages for this snp are all whole numbers anyway
  # MM: DS < 0.2
  # MZ: 0.8 ≤ DS ≤ 1.2
  # ZZ: DS ≥ 1.8
  
  #For all unique PATIENTS -------------------------------------------------------------------------------------------------------#
  table(clinical2[!duplicated(clinical2$Study.ID),"classification"], clinical2[!duplicated(clinical2$Study.ID),"serpina1_z_snp_GRCh38_ch14_pos94378610"])
  
  
  #                      0   1   2
  # Control             94   4   0
  # Mild-moderate COPD 124   9   0
  # Severe COPD        173  19   9
  
  
  
  # For samples ------------------------------------------------------------------#
  clinical_brush <-  clinical2[which(clinical2$sampletype == "Brush"),] #143 (sherlock1, 2 and 3)
  table(clinical_brush$classification, clinical_brush$serpina1_z_snp_GRCh38_ch14_pos94378610)
  
  
  #                      0   1   2
  # Control             88   3   0
  # Mild-moderate COPD 116   9   0
  # Severe COPD        171  17   9
  
  
  
  clinical_biopt <-  clinical2[which(clinical2$sampletype == "Biopt"),] #185 (sherlock 2 and 3)
  table(clinical_biopt$classification, clinical_biopt$serpina1_z_snp_GRCh38_ch14_pos94378610)
  
  
  #                     0  1  2
  # Control            69  2  0
  # Mild-moderate COPD 98  6  0
  # Severe COPD        75  6  5
  
  
  
  # Make master table!
  clinical123_master <- clinical_sk_all
  clinical123_master$serpina1_z_snp_GRCh38_ch14_pos94378610 <- serpina1_z_snp[match(clinical_sk_all$Study.ID, serpina1_z_snp$Study.ID), "z_mutation"]
  
  
  
  ### START EXTRA WRANGLING ------------------------------------------------------------------------------------------------------
  # 15/01/26: Came back to add this section after running 5_SHERLOCK_radiomics and SHERLOCK_sysmex_diffexp.R scripts,  before starting SHERLOCK_sysmex_multivariate.R script
  # The below wrangling is unrelated to this script. editing radiomics and sysmex variables columns so that clinical_sherlock123_master.rds IS THE MOST UPDATED FILE THAT CAN BE USED FOR EVERYTHING ! ###
  
  
  ## FOR RADIOMICS COLUMNS
  #Make the names valid for R
  radiomics_index_start <-which(colnames(clinical123_master) =="RL_insp_vol_ml")
  radiomics_index_end <- which(colnames(clinical123_master) =="LLL_airtrapping_emphysema_%")
  
  colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] <- gsub("%", "perc", colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] )
  colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] <- gsub(">", "over", colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] )
  colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] <- gsub("-", ".", colnames(clinical123_master)[radiomics_index_start:radiomics_index_end] )
  clinical123_master[, radiomics_index_start:radiomics_index_end] <- lapply(
    clinical123_master[, radiomics_index_start:radiomics_index_end],
    function(x) {
      if (is.list(x)) {
        as.numeric(unlist(x))
      } else {
        as.numeric(x)
      }
    }
  )
  sapply(clinical123_master[radiomics_index_start:radiomics_index_end],class)
  
  
  # FOR SYSMEX COLUMNS
  # Some rows for sysmex variable have "----" which turn it into a character. make these NA
  # Sysmex first column = "Mono", last column = "EO-Z", 53 total columns
  sysmex_index_start <- which(colnames(clinical123_master) == "Mono")
  sysmex_index_end <- which(colnames(clinical123_master) == "EO-Z")
  
  colnames(clinical123_master)[sysmex_index_start:sysmex_index_end] <- make.names(colnames(clinical123_master)[sysmex_index_start:sysmex_index_end])
  
  sapply(clinical123_master[sysmex_index_start:sysmex_index_end],class)
  
  clinical123_master[which(clinical123_master$RELYMP.103uL == "----"),"RELYMP.103uL"] <- NA
  clinical123_master$RELYMP.103uL <- as.numeric(clinical123_master$RELYMP.103uL)
  
  clinical123_master[which(clinical123_master$RELYMP == "----"),"RELYMP"] <- NA
  clinical123_master$RELYMP <- as.numeric(clinical123_master$RELYMP)
  
  sapply(clinical123_master[sysmex_index_start:sysmex_index_end],class)
  
  ### END EXTRA WRANGLING ---------------------------------------------------------------------------
  
  saveRDS(clinical123_master, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
  write.csv(clinical123_master, file.path(postQC.data.dir,  "master","clinical_sherlock123_master.csv"))
  
  ### SUMMARISING NUMBERS ----------------------------------------------------------------------------------
  
  
  sherlock1_ids <- unique(clinical123_master[which(clinical123_master$batch == 1), "Study.ID"])
  sherlock2_ids <- unique(clinical123_master[which(clinical123_master$batch == 2), "Study.ID"])
  sherlock3_ids <- unique(clinical123_master[which(clinical123_master$batch == 3), "Study.ID"])
  
  
  ### 626 Patients have genotyping data
  ### 598 patients have clinical data in master file
  ### 541 patients have expression data
  ### 24 PATIENTS THAT HAVE EXPRESSION DATA AND CLINICAL DATA BUT NO GENOTYPING data  ------------------------------------------------- #
  setdiff(unique(clinical123_master$Study.ID), row.names(serpina1_z_snp))
  cat(setdiff(unique(clinical123_master$Study.ID), row.names(serpina1_z_snp)), sep = "\n")
  
  # [1] "A_1494" "A_922"  "A_978"  "A_1655" "A_1472" "A_1542" "A_1680" "A_1541"
  # [9] "A_2700" "A_2642" "A_3063" "A_2890" "A_960"  "SEO066" "SEO069" "SEO070"
  # [17] "SEO075" "SEO077" "SEO078" "SEO087" "SEO185" "SEO196" "SEO507" "SEO418"
  
  # Breakdown------------------------------------------#
  # SHERLOCK1
  sherlock1_ids[!(sherlock1_ids %in% row.names(serpina1_z_snp))]
  # [1] "A_1494" "A_922"  "A_978"  "A_1655" "A_1472" "A_1542" "A_1680" "A_1541"
  # [9] "A_2700" "A_2642" "A_3063" "A_2890" "A_960"
  
  
  # SHERLOCK2 
  sherlock2_ids[!(sherlock2_ids %in% row.names(serpina1_z_snp))]
  # [1] "SEO066" "SEO069" "SEO070" "SEO075" "SEO077" "SEO078" "SEO087" "SEO185"
  # [9] "SEO196" "SEO507"
  
  
  # SHERLOCK3 
  sherlock3_ids[!(sherlock3_ids %in% row.names(serpina1_z_snp))]
  # [1] "SEO069" "SEO070" "SEO077" "SEO185" "SEO418"
  
  ### 432 PATIENTS THAT HAVE GENOTYPING DATA BUT NO EXPRESSION data ### -------------------------------------------------#
  cat(intersect(clinical_sk_all$Study.ID, row.names(serpina1_z_snp)), sep = "\n")
  
  ### 194 PATIENTS THAT HAVE GENOTYPING DATA BUT NO EXPRESSION data ### -------------------------------------------------#
  setdiff(row.names(serpina1_z_snp), clinical123_master$Study.ID)
  cat(setdiff(row.names(serpina1_z_snp), clinical123_master$Study.ID), sep = "\n")
  
  
  ### 179 PATIENTS THAT HAVE GENOTYPING DATA BUT NOT CLINICAL DATA  ### -------------------------------------------------#
  setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id))
  cat(setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id)), sep = "\n")
  
  
  ### 74 PATIENTS THAT HAVE GENOTYPING DATA AND CLINICAL DATA BUT NO EXPRESSION DATA ### ---------------------------------#
  setdiff(setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id)), matching_ids)
  cat(setdiff(setdiff(row.names(serpina1_z_snp), sub("^A", "A_", raw_clinical$class_incl_study_id)), matching_ids), sep = "\n")
  
  
  # ##-- Post batch correction
  counts <- readRDS(file.path(combat.processed.data.dir, "counts_combat.rds"))
  counts_brush <- readRDS(file.path(combat.processed.data.dir, "counts_brush_combat.rds"))
  counts_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
  sherlock1_counts<- readRDS(file.path(processed.data.dir, "SHERLOCK1", "counts_sk1.rds"))
  
  unique(clinical123_master[intersect(colnames(counts_brush), row.names(clinical123_master)), "Study.ID"]) #270 patients have brushes (sherlock 2,3)
  unique(clinical123_master[intersect(colnames(counts_biopt), row.names(clinical123_master)), "Study.ID"]) #271 patients have biopsies (sherlock 2,3)
  unique(clinical123_master[intersect(colnames(sherlock1_counts), row.names(clinical123_master)), "Study.ID"]) #167 patients have brushes (sherlock 1)
  
  
  intersect(row.names(serpina1_z_snp), raw_clinical$class_incl_study_id)
  ### --------------------------------------------------------------------------------------------------------------------------------------#
  




# ================================================================================== #
# 1A. SAMPLE/PATIENT DEMOGRAPHICS TABLE ============================================
# ================================================================================== #


demographics_func <- function(sampletype, disease_group){
  
  if(sampletype == "brush"){
    clinical <- clinical_brush
  }
  
  if(sampletype == "biopt"){
    clinical <- clinical_biopt
  }
  
  #Remove samples with NA serpina1 snp data
  clinical <- clinical[-which(is.na(clinical$serpina1_snp)),]
  
  #make names valid
  clinical$smoking_status<- make.names(clinical$smoking_status)
  
  if(disease_group == "mildmoderate_copd"){
    clinical <- clinical[which(clinical$classification == "Mild-moderate COPD"),]}
  
  if(disease_group == "severe_copd"){
    clinical <- clinical[which(clinical$classification == "Severe COPD"),]}
  
  
  demographics <- t(
    
    clinical %>%   
      
      mutate(
        serpina1_snp = recode(
          serpina1_snp,
          "0" = "MM",
          "1" = "MZ",
          "2" = "ZZ"
        )
      ) %>% 
      
      group_by(serpina1_snp) %>% 
      
      summarise(
        total_patients = n(), #total patents per genotype group
        
        #sex
        male_patients = sum(sex == "Male"), 
        male_percentage = (male_patients / total_patients) * 100,
        sex = paste0(male_patients, "(", round(male_percentage, digits=3), ")"),
        
        #age
        mean_age = median(age),
        range_age =paste0(min(age, na.rm = T), "-", max(age, na.rm = T)),
        age = paste0(mean_age,"(",range_age, ")"),
        
        #Smoking status
        currentsmoker = sum(smoking_status == "Current.smoker"),
        smokerpercentage = currentsmoker/total_patients *100,
        smoke = paste0(currentsmoker, "(", round(smokerpercentage, digits = 3), ")"),
        
        #packyears
        median_packyears = median(packyears, na.rm = TRUE),
        range_packyears = paste0(min(packyears, na.rm = T), "-", max(packyears, na.rm = T)),
        packyears = paste0(median_packyears, "(", range_packyears, ")"),
        
        
        #FEV1
        median_FEV1_percent_pred =  median(FEV1_percent_pred, na.rm = T),
        range_postfev1percpred = paste0(round(min(FEV1_percent_pred, na.rm = T),digits=3), "-", round(max(FEV1_percent_pred, na.rm = T), digits = 3)),
        fev1 = paste0(round(median_FEV1_percent_pred, digits = 3), "(", range_postfev1percpred, ")"),
        
        #FEV/FVC
        median_postfev1fvcpercpred =  median(FEV1_FVC_post, na.rm = TRUE),
        range_postfev1fvcpercpred = paste0(round(min(FEV1_FVC_post, na.rm = T)), "-", round(max(FEV1_FVC_post, na.rm = T), digits = 3 )),
        fev_fvc = paste0(round(median_postfev1fvcpercpred, digits = 3), "(", range_postfev1fvcpercpred, ")")
        
      ) %>%
      
      
      dplyr::select(
        serpina1_snp,
        total_patients,
        serpina1_snp,
        sex,
        age,
        smoke,
        packyears,
        fev1,
        fev_fvc
      )
  ) 
  
  colnames(demographics) <- demographics[1,]
  
  row.names(demographics) <- c(
    "AATD Pi Genotype",
    "Patients, n",
    "Sex Male, n (%)",
    "Age, median (Range)",
    "Current smoker, n (%)",
    "Packyears, median (Range)",
    "FEV1 % pred. (post-bronchodilater), median (Range)",
    "FEV1/FVC % (post-bronchodilater), median (Range)"
  )
  return(demographics)
}

#Bind together mild and severe and save
brush_demographics <- cbind(demographics_func(sampletype = "brush", disease_group = "mildmoderate_copd"), 
                            demographics_func(sampletype = "brush", disease_group = "severe_copd"))
write.csv(brush_demographics, file.path(output.dir, "brush_aatd_demographics.csv"))

biopt_demographics <- cbind(demographics_func(sampletype = "biopt", disease_group = "mildmoderate_copd"), 
                            demographics_func(sampletype = "biopt", disease_group = "severe_copd"))
write.csv(biopt_demographics, file.path(output.dir, "biopt_aatd_demographics.csv"))

} #close the hpc_batch_job. script from here on will be run


# ================================================================================== #
# 2. DIFFERENTIAL EXPRESSION =======================================================
# ================================================================================== #

options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED

# ================================================================================== #
## 2A. SCRIPT SET UP ===============================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("readxl")
library("ggplot2")
library("DESeq2")
library("ggrepel")
library("ggfortify")
library("stringr")
library("tidyverse")

# ================================================================================== #
## 2B. SET UP DIRECTORY & OUTPUT PATHS =============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir,"data")

processed.data.dir <- file.path(data.dir,"processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir, "output","aatd")

#DESEQ DIRECTORY
diffexp.dir <- file.path(output.dir, "diffexp_serpina1_deseq_withoutage", "rm_LIB5426585_SAM24375557")
if(!exists(diffexp.dir))dir.create(diffexp.dir)

setwd(file.path(main.dir))

# ================================================================================== #
## 2.1. LOAD IN DATA ===============================================================
# ================================================================================== #

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
clinical123_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
colnames(clinical123_master)[which(colnames(clinical123_master) == "serpina1_z_snp_GRCh38_ch14_pos94378610")] <- "serpina1_snp"

# Brush SHERLOCK1, 2 and 3 counts
counts123_brush <- read.csv(file.path(combat.processed.data.dir, "sherlock1_2_3_counts_brush_integrated.csv"), check.names = FALSE, row.names = 1) #from SHERLOCK_SOP_integration.R script (copied from  "/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/sherlock1_2_3_combat)
clinical_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush"),] 

# Biopsy SHERLOCK 2 and 3 counts
counts23_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
clinical_biopt <- clinical123_master[which(clinical123_master$sampletype == "Biopt"),]
counts23_biopt <- counts23_biopt[,row.names(clinical_biopt)]

#General COPD vs Severe COPD
# MM = 0
# MZ = 1
# ZZ = 2

table(clinical_brush$classification, clinical_brush$serpina1_snp)
#                      0   1   2
# Control             88   3   0
# Mild-moderate COPD 116   9   0
# Severe COPD        171  17   9

table(clinical_biopt$classification, clinical_biopt$serpina1_snp)
#                     0  1  2
# Control            69  2  0
# Mild-moderate COPD 98  6  0
# Severe COPD        75  6  5

#make valid names




# ================================================================================== #
# 3. DIFFERENTIAL EXPRESSION (DESeq2) =============================================
# General COPD (Subset only for Mild/Mod and Severe COPD. Compare the genotypes)
# Severe COPD (Subset for only Severe COPD. Compare the genotypes)
# ================================================================================== #

# Create empty results lists to save tT, tT2 and volcano into
listofresults <- list(
  general_copd = list(),
  severe_copd  = list()
)

# Identify outliers to be removed - identified after running diffexp and plotting top genes (boxplots)
brush_outliers <- c("LIB5426585_SAM24375557")
biopt_outliers <- NULL

# Create function to run differential expression with DESEq2
diffexp_deseq_func <- function(sampletype, disease_group){
  cat(paste("Starting 3. DIFFERENTIAL EXPRESSION (DESeq)", sampletype), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  if(disease_group == "general_copd"){
    this.diffexp.dir <- file.path(diffexp.dir, "general_copd")}

  if(disease_group == "severe_copd"){
    this.diffexp.dir <- file.path(diffexp.dir, "severe_copd")}

  if(!exists(this.diffexp.dir))dir.create(this.diffexp.dir)

  #Results and figures directory
  diffexp.results.dir <- file.path(this.diffexp.dir, "results")
  if(!exists(diffexp.results.dir))dir.create(diffexp.results.dir)

  diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
  if(!exists(diffexp.figures.dir))dir.create(diffexp.figures.dir)


  if(sampletype == "brush"){
    clinical <- clinical_brush
    counts <- counts123_brush
    outliers <- brush_outliers
  }

  if(sampletype == "biopt"){
    clinical <- clinical_biopt
    counts <- counts23_biopt
    outliers <- biopt_outliers
  }


  if(all(colnames(counts) == row.names(clinical)) == FALSE){
    stop("all(colnames(counts) == row.names(clinical) = FALSE)") }

  if(length(outliers >= 1)){
  clinical <- clinical[-which(row.names(clinical) %in% outliers),]
  }
  
  #Remove samples with NA serpina1 snp data
  clinical <- clinical[-which(is.na(clinical$serpina1_snp)),]
  counts <- counts[,row.names(clinical)]

  #make names valid
  clinical$smoking_status<- make.names(clinical$smoking_status)

  if(disease_group == "general_copd"){
    #Only include COPD samples
    clinical <- clinical[-which(clinical$classification == "Control"),]
    counts <- counts[,row.names(clinical)]

    # DGE
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = clinical,
                                  design = ~ 0 + serpina1_snp + sex + smoking_status)
  }



  if(disease_group == "severe_copd"){
    #Only include COPD samples
    clinical <- clinical[which(clinical$classification == "Severe COPD"),]
    counts <- counts[,row.names(clinical)]

    # DGE
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = clinical,
                                  design = ~ 0 + serpina1_snp + sex ) #no smokers in severe group

  }



  # Filter the genes that are lowly expressed and normalize
  # Low exp genes affects statistics - can make p value very significant even if only a few samples are lowly expressed amongst other samples with no expression.  also affects multiple testing.
  # One method = keep row medians that are greater than 10 (ie. half of the samples for a gene must have a minimum number of 10 counts)
  keep <- rowMedians(counts(dds)) >= 10

  dds <- dds[keep,]
  dds <- DESeq(dds)

  # results extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values;
  resultsNames(dds)

  listofcontrasts <- list(
    MZ_MM = c("serpina1_snp", "1", "0"),#This means we have set snp1 (MZ) as control and we define fold change based on snp1 as baseline (log fold change = MZ - MMM)
    ZZ_MM = c("serpina1_snp", "2", "0"),
    ZZ_MZ = c("serpina1_snp", "2", "1")

  )

  listoftT <- list()
  listoftT2 <- list()
  listofvolcano <- list()

  for (contrast in names(listofcontrasts)){
    cat(paste(sampletype, contrast), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

    results <- results(dds, contrast = c(listofcontrasts[[contrast]][1], #variable name
                                         listofcontrasts[[contrast]][2], #reference
                                         listofcontrasts[[contrast]][3])) #other comparison

    tT <- as.data.frame(results)

    # is this done auto
    # tT$padj=p.adjust(tT$pvalue,method="BH")

    tT <- as.data.frame(results)
    tT <-tT[order(tT$pvalue),]
    # baseMean is the average of the normalized count values, dividing by size factors, taken over all samples.


    tT$Legend <- ifelse(
      tT$padj < 0.05 & tT$log2FoldChange > 0, "Upregulated",
      ifelse(
        tT$padj < 0.05 & tT$log2FoldChange < 0, "Downregulated",
        "Not Significant"))

    tT$Legend[is.na(tT$Legend)]="Not Significant"

    tT$Legend <- factor(tT$Legend, levels = c("Downregulated", "Upregulated", "Not Significant"))

    tT$gene_symbol=hgnc_symbols_db[row.names(tT), "SYMBOL"] #add hgnc symbols

    # if(showEnsemblID == TRUE){
    #for those with no hgnc symbol, label with ensembl id
    tT[which(is.na(tT$gene_symbol)), "gene_symbol"] <- row.names(tT)[(which(is.na(tT$gene_symbol)))] #listofresults_withensembl
    # }
    #
    # else{
    #   #for those with no hgnc symbol, remove
    #   tT <- tT[-which(is.na(tT$gene_symbol)), ] #listofresults_hgnconly
    # }

    selection <-which(tT$padj<0.05)

    tT2 <- tT[selection,]

    listoftT[[contrast]] <- tT
    listoftT2[[contrast]] <- tT2
    
    write.csv(tT, file = file.path(diffexp.results.dir, paste0(sampletype, "_", contrast, "_tT.csv")))
    write.csv(tT2, file = file.path(diffexp.results.dir, paste0(sampletype, "_", contrast, "_tT2.csv")))


    # ================================================================================== #
    # 3.1. VOLCANO PLOT ================================================================
    # ================================================================================== #
    cat(paste("Starting 3.1. VOLCANO PLOT", contrast), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

    volcano <- ggplot(tT, aes(x = log2FoldChange, y = -log10(pvalue))) +
      ggtitle(paste0(sampletype,": ",contrast)) +
      geom_point(aes(color = Legend)) +
      scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red"), drop = FALSE)+
      geom_hline(yintercept =-log10(max(tT2$pvalue)),colour="black", linetype="dashed")+
      geom_text_repel(data = subset(tT2[1:30,]),
                      aes(label= gene_symbol),size = 4, box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines") ) +
      theme_bw(base_size = 18) + theme(legend.position = "bottom",
                                       legend.text = element_text(size = 14),
                                       legend.title = element_text(size = 16))

    listofvolcano[[contrast]] <- volcano

    ggsave(volcano, filename = file.path(diffexp.figures.dir, paste0(sampletype, "_", contrast, "_volcano.png")),
           width = 25, height = 25,
           units = "cm")

  } #close listofcontrasts loop

  listofresults <- list(tT = listoftT, tT2 = listoftT2, volcano = listofvolcano)

} #end diffexp deseq function


# Run function on general COPD data
listofresults_brush_general <- diffexp_deseq_func(disease_group = "general_copd", sampletype = "brush")
listofresults_biopt_general <- diffexp_deseq_func(disease_group = "general_copd", sampletype = "biopt")

# Run function on severe COPD data
listofresults_brush_severe <- diffexp_deseq_func(disease_group = "severe_copd", sampletype = "brush")
listofresults_biopt_severe <- diffexp_deseq_func(disease_group = "severe_copd", sampletype = "biopt")


#Save all results (Made seperate directoriess for general COPD and severe COPD.)
#listofresults contains the lists "brush" and "biopt", which then each contain lists "listoftT", "listoftT2" and "listofvolcano")
saveRDS(list(brush = listofresults_brush_general,
             biopt = listofresults_biopt_general),
        file = file.path(diffexp.dir, "general_copd", "results", "listofresults.rds"))

saveRDS(list(brush = listofresults_brush_severe,
             biopt = listofresults_biopt_severe),
        file = file.path(diffexp.dir, "severe_copd", "results","listofresults.rds"))


# } #close the hpc_batch_job. script from here on will be run


# options(error = function() { traceback(); quit(status = 1) })
#options(error = ...) tells r to run the function
#traceback()	prints the call stack (what functions were running in what order at the time of failure. traceback(2) means skip the top frame (the error handler itself).
#quit(status = 1) tells R to exit and that the script has failed (1=FAIL and 0 = SUCCESS) - sacct command will show job status as FAILED

# ================================================================================== #
## 6A. SCRIPT SET UP ===============================================================
# ================================================================================== #
my_directory <- "/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3"

library("readxl")
library("ggplot2")
library("DESeq2")
library("ggrepel")
library("ggfortify")
library("stringr")
library("tidyverse")
library("rstatix")
library("ggpubr")

# ================================================================================== #
## 6B. SET UP DIRECTORY & OUTPUT PATHS =============================================
# ================================================================================== #
main.dir <- my_directory

#Data directory
data.dir <- file.path(main.dir,"data")

processed.data.dir <- file.path(data.dir,"processed")

postQC.data.dir <- file.path(processed.data.dir, "datawrangling_qc")
combat.processed.data.dir <- file.path(postQC.data.dir, "combat_results")

#Output directory
output.dir <- file.path(main.dir, "output","aatd")

#DESEQ DIRECTORY
diffexp.dir <- file.path(output.dir, "diffexp_serpina1_deseq_withoutage", "rm_LIB5426585_SAM24375557")
if(!exists(diffexp.dir))dir.create(diffexp.dir)

#EDGER DIRECTORY
# diffexp.dir <- file.path(output.dir, "diffexp_serpina1_edgeR_withoutage")
# if(!exists(diffexp.dir))dir.create(diffexp.dir)

setwd(file.path(main.dir))

# ================================================================================== #
## 6.1. LOAD IN DATA ===============================================================
# ================================================================================== #

hgnc_symbols_db <- readRDS(file.path(postQC.data.dir,"hgnc_symbols_db.rds"))
clinical123_master <- readRDS(file.path(postQC.data.dir,  "master","clinical_sherlock123_master.rds"))
colnames(clinical123_master)[which(colnames(clinical123_master) == "serpina1_z_snp_GRCh38_ch14_pos94378610")] <- "serpina1_snp"

# Brush SHERLOCK1, 2 and 3 counts
clinical_brush <- clinical123_master[which(clinical123_master$sampletype == "Brush"),] 
counts123_brush <- read.csv(file.path(combat.processed.data.dir, "sherlock1_2_3_counts_brush_integrated.csv"), check.names = FALSE, row.names = 1) #from SHERLOCK_SOP_integration.R script (copied from  "/groups/umcg-griac/tmp02/projects/SHERLOCK_2025/data/sherlock1_2_3_combat)

# Biopsy SHERLOCK 2 and 3 counts
clinical_biopt <- clinical123_master[which(clinical123_master$sampletype == "Biopt"),]
counts23_biopt <- readRDS(file.path(combat.processed.data.dir, "counts_biopt_combat.rds"))
counts23_biopt <- counts23_biopt[,row.names(clinical_biopt)]

# ================================================================================== #
# 6. GSVA ==========================================================================
# GSVA of ZZ signature on rest of data (to see if those genes go up in MZ at all)
# ================================================================================== #
library(GSVA)
library(ggpubr)

#Get tT2 for severe copd
listofresults <- readRDS(file.path(diffexp.dir, "severe_copd", "results", "listofresults.rds"))

#Make dds objects and vst normalise the countss

#Run GSVA of gene signature that was upregulated in ZZ (compared to MM) in Severe COPD on
#  i) Severe COPD patients
# ii) General COPD (Severe COPD and Mild COPD combined)

gsva_func <- function(sampletype, disease_group){
  
  
  if(sampletype == "brush"){
    counts <- counts123_brush
    clinical <- clinical_brush
    outliers <- brush_outliers}
  
  if(sampletype == "biopt"){
    counts <- counts23_biopt
    clinical <- clinical_biopt
    outliers <- biopt_outliers}
  
  
  if(disease_group == "general_copd"){
    #Only include COPD samples
    clinical <- clinical[-which(clinical$classification == "Control"),]
    this.diffexp.dir <- file.path(diffexp.dir, "general_copd")}
  
  
  if(disease_group == "severe_copd"){
    #Only include COPD samples
    clinical <- clinical[which(clinical$classification == "Severe COPD"),]
    this.diffexp.dir <- file.path(diffexp.dir, "severe_copd")}
  
  
  diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
  gsva.dir <- file.path(diffexp.figures.dir, "gsva")
  if(!exists(gsva.dir))dir.create(gsva.dir)
  
  
  
  if(length(outliers >= 1)){
  clinical <- clinical[-which(row.names(clinical) %in% outliers),]
  }

  #Remove samples with NA serpina1 snp data and change numbers to letters
  clinical <- clinical[-which(is.na(clinical$serpina1_snp)),] %>% 
    mutate(serpina1_snp = recode(
      serpina1_snp,
        "0" = "MM",
        "1" = "MZ",
        "2" = "ZZ" ))
  
  #make names valid
  clinical$smoking_status<- make.names(clinical$smoking_status)
  
  #match counts to clinical
  counts <- counts[,row.names(clinical)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = clinical,
                                design = ~ 1) #no design needed, just need to make dds aobject so i can vst normalise
  
  counts_vst <- assay(vst(dds)) #assay extracts the counts
  
  # Get gene signature set
  tT2 <- listofresults[[sampletype]][["tT2"]] #severe results
  severe_ZZ_MM_tT2 <- tT2[["ZZ_MM"]]
  
  # When including outlier (for ZZ-MM, genes are upregulated only) !!!!!!!!!!!!!!!!! ----------------
  # severe_ZZ_MM_up <- row.names(severe_ZZ_MM_tT2)[which(severe_ZZ_MM_tT2$Legend == "Upregulated")]
  # gsva_res <- gsva(as.matrix(counts_vst), 
  #                  list(severe_ZZ_MM_up)
  #                  mx.diff = TRUE)
  # gsva_res <- t(gsva_res) #the results tell us how much the set of genes was represented in each sample. ie. enrichment score of 0.9 is high- meaning the genes of interest showed up alot in sample X - now when we group the samples by copd and non copd, we can see whether certain genes are enriched in samples with or without copd
  # colnames(gsva_res) <- paste0(sampletype,"_ZZ_MM_Up")
  # --------------------------------------------------------------------------------- #
  
  #When not including outlier (for ZZ-MM, genes are downregulated only) !!!!!!!!!!!!!!!!! -----------
  severe_ZZ_MM_down <- row.names(severe_ZZ_MM_tT2)[which(severe_ZZ_MM_tT2$Legend == "Downregulated")]
  gsva_res <- gsva(as.matrix(counts_vst), 
                   list(severe_ZZ_MM_down), 
                   mx.diff = TRUE)
  
  gsva_res <- t(gsva_res)
  colnames(gsva_res) <- paste0(sampletype,"_ZZ_MM_Down")
  # --------------------------------------------------------------------------------- #
  boxplot_gsva=cbind(gsva = gsva_res,
                       genotype= as.character(clinical$serpina1_snp),
                       disease = as.character(clinical$classification))
  
  boxplot_gsva <- as.data.frame(boxplot_gsva)
  
  
  my_comparisons <- list(c("MM", "ZZ"),
                         c("MM", "MZ"),
                         c("ZZ", "MZ"))
  
  
  
    x_order <- c("MM", "MZ", "ZZ")
  #make this a loop if there are multiple gene set lists         y = as.numeric(boxplot_gsva[,i]),
    
    gsva_theme <- theme(axis.title = element_text(size = 24),
                        axis.text = element_text(size = 24),
                        title = element_text(size = 20),
                        legend.position = "bottom")
    
    
    boxplot <- ggplot(boxplot_gsva, aes(
      x = factor(genotype, levels = x_order),
      y = as.numeric(boxplot_gsva[,1]),
      fill = genotype)) +
      
      theme_bw()+
      
      gsva_theme +
      
      geom_boxplot(position = position_dodge(1),
                   aes(alpha = 0.5)) +
      
      
      geom_jitter(aes(color = disease),
                  alpha = 0.5,
                  size = 2.5,
                  width = 0.3) +
      

      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         paired = FALSE,
                         size = 7)+

      scale_fill_manual(values=c("MM" = "#00BA38" , "MZ" = "#619CFF",
                                 "ZZ" = "#F8766D")) +
      
      
      scale_color_manual(values=c("Mild-moderate COPD" = "#E68613" , 
                                  "Severe COPD" = "#C77CFF")) +
      
      # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
     
       guides(
        fill  = "none",
        alpha = "none"
      ) +
      
      labs(title = paste0("Signature Analysis", "(", paste0("SevereCOPD_",colnames(boxplot_gsva)), ")")) +
      ylab (label = "Enrichment Score") +
      xlab (label = "AATD Pi Genotype")
    
    
    ggsave(boxplot, file = file.path(gsva.dir,paste0("SevereSig_",colnames(boxplot_gsva)[1], "_", disease_group,"data", ".png")), width = 3000, height = 2100, units = "px" )
    
} #close function


gsva_func(sampletype = "brush", disease_group = "severe_copd")
# gsva_func(sampletype = "biopt", disease_group = "severe_copd")

gsva_func(sampletype = "brush", disease_group =  "general_copd")
# gsva_func(sampletype = "biopt", disease_group =  "general_copd")



# ================================================================================== #
# 7. BOXPLOTS ==========================================================================
# ================================================================================== #

boxplot_func <- function(sampletype, disease_group){
  
  listofresults <- readRDS(file.path(diffexp.dir, disease_group, "results", "listofresults.rds"))
  
  
  if(sampletype == "brush"){
    counts <- counts123_brush
    clinical <- clinical_brush
    outliers <- brush_outliers}
  
  if(sampletype == "biopt"){
    counts <- counts23_biopt
    clinical <- clinical_biopt
    outliers <- biopt_outliers}
  
  
  if(disease_group == "general_copd"){
    #Only include COPD samples
    clinical <- clinical[-which(clinical$classification == "Control"),]
    this.diffexp.dir <- file.path(diffexp.dir, "general_copd")}
  
  
  if(disease_group == "severe_copd"){
    #Only include COPD samples
    clinical <- clinical[which(clinical$classification == "Severe COPD"),]
    this.diffexp.dir <- file.path(diffexp.dir, "severe_copd")}
  
  
  diffexp.figures.dir <- file.path(this.diffexp.dir, "figures")
  boxplot.dir <- file.path(diffexp.figures.dir, "boxplot")
  if(!exists(boxplot.dir))dir.create(boxplot.dir)
  
  
  if(length(outliers >= 1)){
  clinical <- clinical[-which(row.names(clinical) %in% outliers),]
  }
  
  #Remove samples with NA serpina1 snp data and change numbers to letters
  clinical <- clinical[-which(is.na(clinical$serpina1_snp)),] %>% 
    mutate(serpina1_snp = recode(
      serpina1_snp,
      "0" = "MM",
      "1" = "MZ",
      "2" = "ZZ" ))
  
  #make names valid
  clinical$smoking_status<- make.names(clinical$smoking_status)
  
  #match counts to clinical
  counts <- counts[,row.names(clinical)]
  
  #Make dds objects and vst normalise the countss
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = clinical,
                                design = ~ 1) #no design needed, just need to make dds object so i can vst normalise
  
  counts_vst <- assay(vst(dds)) #assay extracts the counts
  
  # Get top genes 
  tT2 <- listofresults[[sampletype]][["tT2"]][["ZZ_MM"]] #severe or general results depending on function selection
  
  
  boxplotdata <- as.data.frame(t(counts_vst))
  
  boxplot <- cbind(boxplotdata,
                   genotype = as.character(clinical$serpina1_snp),
                   disease = as.character(clinical$classification),
                   age = clinical$age,
                   smoking_status = clinical$smoking_status,
                   packyears = clinical$packyears,
                   ics = clinical$ics_use)
  
  tT2_genes<- row.names(tT2)
  
  top10 <- if (length(tT2_genes) >= 10) {
    tT2_genes[1:10]
  } else {
    tT2_genes
  }
  
  
  
  pdf(file = file.path(boxplot.dir, paste0(disease_group,"_",sampletype, "_boxplot",".pdf")),
      height = 8,
      width= 6)
  
  cat("Starting plots", disease_group, sampletype, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  
  
  for (gene in c(top10)){
    
    gene_hgnc <- ifelse(is.na(hgnc_symbols_db[gene, "SYMBOL"]),
                        gene,
                        hgnc_symbols_db[gene, "SYMBOL"])
    
    
    plot <- boxplot[,c(gene,
                       "genotype",
                       "disease",
                       "age",
                       "packyears",
                       "ics")]
    
    
    plot$Study.ID <- row.names(plot)
    
    
    colnames(plot)[1] <- "gene"
    
    plot <- as.data.frame(plot)
    
    
    ## Get P-Values --------------------------------------------------------------------------------------
    
    # THE TABLES IN STEP 1 AND 2 ARE TO GET X POSITIONS, THE T-TEST VALUES WILL NOT BE USED #
    #### STEP 1) THIS IS  A TABLE OF THE COMPARISONS I ACTUALLY WANT --------------------------------------------------------------------------------------
    stat.table <- plot %>%
      t_test(gene ~ genotype)
    
    stat.table<- stat.table %>%
      add_xy_position(x = "genotype", dodge = 0.8)
    
    # stat.table$contrast <- paste0(stat.table$group1, "-" ,stat.table$group2)
    
    
    
    #### STEP 2) MODIFY STAT TABLE TO INCLUDE COMPARISONS OF INTEREST --------------------------------------------------------------------------------------
    stat.table3 <- stat.table
    
    
    stat.table3 <- cbind(stat.table3, resultsname = c("MZ_MM", "ZZ_MM", "ZZ_MZ"))
    stat.table3[which(stat.table3$resultsname == "MZ_MM"),"p"] <- listofresults[[sampletype]][["tT"]][["MZ_MM"]][gene, "pvalue"]
    stat.table3[which(stat.table3$resultsname == "ZZ_MM"),"p"] <- listofresults[[sampletype]][["tT"]][["ZZ_MM"]][gene, "pvalue"]
    stat.table3[which(stat.table3$resultsname == "ZZ_MZ"),"p"] <- listofresults[[sampletype]][["tT"]][["ZZ_MZ"]][gene, "pvalue"]
    stat.table3$p <- signif(as.numeric(stat.table3$p), digits = 4)
    stat.table3$y.position <- max(plot[,"gene"]) + 0.025*(max(plot[,"gene"]))
    stat.table3$y.position <- as.numeric(stat.table3$y.position)
    
    
    boxplot_theme <- theme(axis.title = element_text(size = 22),
                           axis.text = element_text(size = 22),
                           title = element_text(size = 18),
                           legend.position = "bottom") 
    
    boxplotimage <- ggplot(plot, aes(
      x = as.factor(genotype),
      y = gene)) +
      
      theme_bw()+
      
      boxplot_theme +
      
      geom_boxplot(position = position_dodge(1),
                   aes(fill = genotype,
                       alpha = 0.5)) +
      
      
      geom_jitter(aes(color = disease),
                  alpha = 0.5,
                  size = 2.5,
                  width = 0.3) +
      
      
      scale_fill_manual(values=c("MM" = "#00BA38" , "MZ" = "#619CFF",
                                 "ZZ" = "#F8766D")) +
      
      
      scale_color_manual(values=c("Mild-moderate COPD" = "#E68613" , 
                                  "Severe COPD" = "#C77CFF")) +
      
      scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      
      labs(title = paste(gene_hgnc, "expression in", disease_group, "patients")) +
      ylab (label =  paste(gene_hgnc, "expression")) +
      xlab (label = paste("AATD Pi Genotype")) +
    
      
      stat_pvalue_manual(stat.table3,
                         label = "p",
                         tip.length = 0.01,
                         bracket.nudge.y = c(0, 0.5, 0),
                         size = 6)   +
    
      
      guides(
        fill  = "none",
        alpha = "none"
      ) 
    
    
    ggsave(boxplotimage,
           file = file.path(boxplot.dir, paste0(sampletype, "_", gene_hgnc, "_boxplot.png")),
           height = 20, width = 20, units = "cm", dpi = 800)
    
    
    
    boxplotimage_labelled <- boxplotimage +
      geom_text_repel(
        aes(label = Study.ID, color = disease),
        size = 4,
        show.legend = FALSE) 
    
    
    
    image <- ggarrange(plotlist = list(boxplotimage,
                                       boxplotimage_labelled),
                       nrow = 2,
                       ncol = 1)
    
    print(image)
    
    
  } #close loop of top10genes
  
  dev.off()
} # close function

  


boxplot_func(sampletype = "brush", disease_group = "severe_copd")
# boxplot_func(sampletype = "biopt", disease_group = "severe_copd")

boxplot_func(sampletype = "brush", disease_group =  "general_copd")
# boxplot_func(sampletype = "biopt", disease_group =  "general_copd")


cat("END OF THIS JOB", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#check these genes in the tT from diffexp results where I still included the outlier 9are they still in same direction)

tT2 <- listofresults[[sampletype]][["tT2"]][["ZZ_MM"]] #outlier removed, ZZ_MM, severe_copd brush
# 4 genes all downregulated, check these genes in the tT from diffexp still including the outlier
# > tT2
#                  baseMean   log2FoldChange lfcSE      stat       pvalue
# ENSG00000136689 372.89814      -3.428353 0.7124784 -4.811869 1.495255e-06
# ENSG00000166920  53.57579      -2.961956 0.6402094 -4.626542 3.718218e-06
# ENSG00000205420 276.10489      -2.801264 0.6184296 -4.529642 5.908380e-06
# ENSG00000171401 173.24922      -3.684852 0.8197240 -4.495235 6.949306e-06
#                    padj          Legend       gene_symbol
# ENSG00000136689 0.02552700 Downregulated       IL1RN
# ENSG00000166920 0.02965964 Downregulated    C15orf48
# ENSG00000205420 0.02965964 Downregulated       KRT6A
# ENSG00000171401 0.02965964 Downregulated       KRT13

old_tT <- readRDS("/groups/umcg-griac/tmp02/projects/KathyPhung/SHERLOCK3/output/aatd/diffexp_serpina1_deseq_withoutage/severe_copd/results/listofresults.rds")
old_tT <- old_tT[["brush"]][["tT"]][["ZZ_MM"]]

# > old_tT[row.names(tT2),]
#                  baseMean   log2FoldChange lfcSE      stat       pvalue
# ENSG00000136689 374.10901      -2.271335 0.6723549 -3.378179 7.296764e-04 
# ENSG00000166920  53.68061      -2.322292 0.5880406 -3.949203 7.841195e-05
# ENSG00000205420 277.16105      -2.008045 0.5833357 -3.442349 5.766859e-04
# ENSG00000171401 174.13736      -2.252625 0.7664434 -2.939063 3.292066e-03
#                    padj          Legend       gene_symbol
# ENSG00000136689 0.10156835 Not Significant       IL1RN
# ENSG00000166920 0.03933512   Downregulated    C15orf48
# ENSG00000205420 0.09196181 Not Significant       KRT6A
# ENSG00000171401 0.21029766 Not Significant       KRT13

## Match up with exome data
clinical_exome <- fread(file.path(data.dir, "raw", "exome", "clinfile.txt"))
clinical_exome$Study.ID <- paste0("A_", clinical_exome$aanmnr)

serpina1_z_snp <- read.csv(file.path(processed.data.dir, "snp_array","serpina1_z_snp_dosage.csv"))
intersect(clinical_exome$Study.ID, serpina1_z_snp$Study.ID)

#179 genotyping patients have exome seq data (does clinfile have all exome seq data though????)
# > intersect(clinical_exome$Study.ID, serpina1_z_snp$Study.ID)
# [1] "A_1520" "A_1470" "A_1515" "A_1496" "A_1477" "A_1507" "A_1512" "A_1537"
# [9] "A_1524" "A_1446" "A_1584" "A_1024" "A_1580" "A_1615" "A_952"  "A_1657"
# [17] "A_753"  "A_1167" "A_1687" "A_580"  "A_662"  "A_1664" "A_1353" "A_1688"
# [25] "A_1695" "A_494"  "A_1189" "A_1073" "A_1716" "A_1694" "A_1681" "A_860"
# [33] "A_1388" "A_1435" "A_1576" "A_1721" "A_1402" "A_1727" "A_994"  "A_984"
# [41] "A_1119" "A_1735" "A_1754" "A_1756" "A_143"  "A_1769" "A_1209" "A_1778"
# [49] "A_1780" "A_305"  "A_1790" "A_1772" "A_1804" "A_1057" "A_817"  "A_1826"
# [57] "A_1856" "A_1840" "A_1848" "A_1849" "A_1860" "A_1879" "A_1868" "A_1877"
# [65] "A_1122" "A_108"  "A_1145" "A_1002" "A_1896" "A_1902" "A_1912" "A_884"
# [73] "A_1919" "A_1827" "A_1911" "A_1926" "A_1941" "A_1405" "A_1954" "A_1407"
# [81] "A_1961" "A_1950" "A_1985" "A_2028" "A_53"   "A_2036" "A_2026" "A_2022"
# [89] "A_1391" "A_660"  "A_727"  "A_1983" "A_2154" "A_2029" "A_2139" "A_2132"
# [97] "A_2257" "A_597"  "A_2135" "A_2095" "A_2140" "A_2215" "A_2274" "A_2271"
# [105] "A_1233" "A_2304" "A_2301" "A_2315" "A_2307" "A_2423" "A_2326" "A_2321"
# [113] "A_1891" "A_2328" "A_2352" "A_2358" "A_2365" "A_2360" "A_2361" "A_1904"
# [121] "A_2372" "A_2377" "A_2393" "A_2292" "A_2395" "A_582"  "A_2397" "A_2552"
# [129] "A_1560" "A_506"  "A_2479" "A_1277" "A_2448" "A_2453" "A_2509" "A_2489"
# [137] "A_2507" "A_2554" "A_757"  "A_2557" "A_719"  "A_1667" "A_2625" "A_2629"
# [145] "A_2594" "A_2466" "A_2607" "A_2617" "A_2638" "A_2666" "A_2719" "A_2667"
# [153] "A_2677" "A_2684" "A_2709" "A_2803" "A_2736" "A_2804" "A_2785" "A_2872"
# [161] "A_2981" "A_3172" "A_3119" "A_3062" "A_3130" "A_3145" "A_3167" "A_3241"
# [169] "A_86"   "A_3335" "A_3394" "A_3566" "A_3542" "A_3618" "A_3469" "A_3733"
# [177] "A_3823" "A_3830" "A_3695"


###   venn diagram --------------------------------------------------------------------------------------------------------------------------------------#

library(ggvenn)
#---Venn diagram to compare biopt and brush#
x <- list(expression = clinical_sk_all$Study.ID, 
          genotyping = serpina1_z_snp$Study.ID,
          clinical = raw_clinical$class_incl_study_id)



vennplot = ggvenn(x, stroke_size = 1.5) + ggtitle("SHERLOCk Samples")
ggsave(
  
  filename = "SHERLOCK_samples_venn.png",
  plot = vennplot,
  width = 6,
  height = 6,
  dpi = 300
)

clinical_only <- setdiff(x$clinical, union(x$genotyping, x$expression))
serp_only     <- setdiff(x$genotyping, union(x$clinical, x$expression))
expression_only      <- setdiff(x$expression, union(x$clinical, x$genotyping))
