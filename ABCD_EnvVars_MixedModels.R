## load required packages ## 

package_list <- c("cvTools", "summarytools", "psych", "lme4", "gvlma", "extrafont", "lattice", 
                  "ade4", "fastDummies", "tidyverse", "FactoMineR", "factoextra", "magrittr", 
                  "data.table", "corrplot", "effectsize", "dplyr", "Hmisc", "foreach", "tableone", 
                  "gridExtra", "RColorBrewer", "lmerTest", "optimx", "broom", "MuMIn", "effects", 
                  "sjPlot", "JWileymisc", "multilevelTools")

lapply(package_list, require, character.only=TRUE)

## set working directory ##

setwd("working_directory_path")

## Create paths to files ## 

rds_path <- c("path_to_rds_file") #note: DEAP Study release 2.0.1 update (DOI 10.15154/1506087) 
first_release_txt_path <- c("path_to_release1.1_acspsw02_txt_file") #note: we use this acspsw02 txt file to determine subjects available for Release 1.1
dcan_lab_matched_split_tsv_path <- c("path_to_DCAN_labs_matched_case_splits_participant_list_tsv_file") #note:  (https://nda.nih.gov/edit_collection.html?id=3165). 

## Import Rds file ##

full_data <- read_rds(rds_path)

## Missingness for other culture vars ##

acculturation_data <- full_data %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  dplyr::select(via_accult_ss_hc_p, via_accult_ss_amer_p, accult_phenx_q1, accult_phenx_q2, accult_phenx_q3_dropdwn,
                accult_phenx_q4, accult_phenx_q5, accult_phenx_q1_p, accult_phenx_q2_p, accult_phenx_q3_dropdwn_p,
                accult_phenx_q4_p, accult_phenx_q5_p)

acculturation_table <- CreateTableOne(data = acculturation_data)
summary(acculturation_table)

## 1. Import DEAP Rds file and set up working directories, get subject IDs for Release 1 ##

abcd_data <- full_data %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  dplyr::select(subjectid, abcd_site, rel_family_id,age, demo_gender_id_p, demo_sex_p, female, race_ethnicity, #general participant info/ sociodemo 
                demo_comb_income_p,demo_prnt_ed_p, demo_prtnr_ed_p, 
                meim_ss_exp_p, meim_ss_com_p, macvs_ss_fo_p, macvs_ss_fr_p, macvs_ss_fs_p, macvs_ss_isr_p, macvs_ss_r_p, #culture environment
                fes_ss_fc_p_pr, parental_monitoring_ss_mean, fes_ss_fc_pr, crpbi_acceptance_ss_studycaregiver, # home environment
                school_risk_phenx_ss_ses, school_risk_phenx_ss_iiss, school_risk_phenx_ss_dfs, #school environment
                neighb_phenx_ss_mean_p, neighb_phenx, reshist_addr1_walkindex, reshist_addr1_no2, reshist_addr1_pm25, #neighborhood environment
                reshist_addr1_proxrd,reshist_addr1_d1a, reshist_addr1_p1vlnt, reshist_addr1_drugtot, reshist_addr1_drgsale, 
                reshist_addr1_drgposs, reshist_addr1_dui, reshist_addr1_mjsale, reshist_addr1_adi_edu_l, reshist_addr1_adi_edu_h, 
                reshist_addr1_adi_work_c, reshist_addr1_adi_income, reshist_addr1_adi_in_dis, reshist_addr1_adi_home_v,
                reshist_addr1_adi_rent, reshist_addr1_adi_mortg, reshist_addr1_adi_home_o, reshist_addr1_adi_unemp,
                reshist_addr1_adi_pov, reshist_addr1_adi_b138, reshist_addr1_adi_sp, reshist_addr1_adi_ncar,
                reshist_addr1_adi_ntel, reshist_addr1_adi_nplumb,  reshist_addr1_adi_crowd,
                nihtbx_picvocab_uncorrected, nihtbx_flanker_uncorrected, nihtbx_list_uncorrected, nihtbx_cardsort_uncorrected, # cog behavioral scores
                nihtbx_pattern_uncorrected, nihtbx_picture_uncorrected, nihtbx_reading_uncorrected, pea_ravlt_sd_trial_i_tc,
                pea_ravlt_sd_trial_ii_tc, pea_ravlt_sd_trial_iii_tc, pea_ravlt_sd_trial_iv_tc, pea_ravlt_sd_trial_v_tc, lmt_scr_perc_correct,
                neurocog_pc1.bl, neurocog_pc2.bl, neurocog_pc3.bl # cog BPPCs
  )

ind_pea_ravlt = c(which(names(abcd_data)=="pea_ravlt_sd_trial_i_tc"),which(names(abcd_data)=="pea_ravlt_sd_trial_ii_tc"),
                  which(names(abcd_data)=="pea_ravlt_sd_trial_iii_tc"),which(names(abcd_data)=="pea_ravlt_sd_trial_iv_tc"),
                  which(names(abcd_data)=="pea_ravlt_sd_trial_v_tc"))

abcd_data$pea_ravlt_ld = apply(abcd_data[,ind_pea_ravlt],1,sum)

first_release_subjects <- read.delim(first_release_txt_path, na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>%
  slice(-1) %>%
  dplyr::select(src_subject_id)

## dummy variables for sex and race/ethnicity ##
## convert age to years; convert environment factor variables to continuous variables ## 

abcd_data <- abcd_data %>% 
  dummy_cols(ignore_na = TRUE, select_columns = c("race_ethnicity")) %>%
  mutate(rel_family_id = factor(rel_family_id),
         age = age/12,
         male = ifelse(female == "yes", 0, 1),
         neighb_phenx = as.numeric(factor(neighb_phenx, 
                                          levels = c("Strongly Disagree", "Disagree", 
                                                     "Neutral (neither agree nor disagree)",
                                                     "Agree", "Strongly Agree"),
                                          labels = c(1:5))))

## log transform neighborhood environment variables ##

abcd_data <- abcd_data %>%
  mutate(Res_density_log = log1p(reshist_addr1_d1a), #residential density
         UCR_ViolCrime_log = log1p(reshist_addr1_p1vlnt), # total adult violent crime
         UCR_DrugAbuse_log = log1p(reshist_addr1_drugtot), # drug abuse violations
         UCR_DrugSale_log = log1p(reshist_addr1_drgsale), #drug sale total
         UCR_DrugPoss_log = log1p(reshist_addr1_drgposs), # drug possession
         UCR_DUI_log = log1p(reshist_addr1_dui), # dui reports
         UCR_MJSale_log = log1p(reshist_addr1_mjsale), # marijuana sales
         ProxRoads_log = log(reshist_addr1_proxrd), # proximity to roads
         reshist_addr1_adi_ncar_log = log1p(reshist_addr1_adi_ncar),
         reshist_addr1_adi_ntel_log = log1p(reshist_addr1_adi_ntel),
         reshist_addr1_adi_nplumb_log = log1p(reshist_addr1_adi_nplumb),
         reshist_addr1_adi_crowd_log = log1p(reshist_addr1_adi_crowd))

## transform income and education variables into continuous integer values ## 

household.income = as.character(abcd_data$demo_comb_income_p)  
household.income[household.income == "Less than $5,000"] = 1 #  Less than $5,000
household.income[household.income == "$5,000 through $11,999"] = 2 # $5,000 through $11,999
household.income[household.income == "$12,000 through $15,999"] = 3 # $12,000 through $15,999
household.income[household.income == "$16,000 through $24,999"] = 4 # $16,000 through $24,999
household.income[household.income == "$25,000 through $34,999"] = 5 # $25,000 through $34,999
household.income[household.income == "$35,000 through $49,999"] = 6 # $35,000 through $49,999
household.income[household.income == "$50,000 through $74,999"] = 7 # $50,000 through $74,999
household.income[household.income == "$75,000 through $99,999"] = 8 # $75,000 through $99,999
household.income[household.income == "$100,000 through $199,999"] = 9 # $100,000 through $199,999
household.income[household.income == "$200,000 and greater"] = 10 # $200,000 and greater
household.income[household.income == "Don't know"] = NA
household.income[household.income == "Refuse to answer"] = NA
household.income[household.income %in% c(NA, "Don't know", "Refuse to answer")] = NA
abcd_data$household.income = factor( household.income, levels= 1:10, 
                                     labels = c("LT5K", "8.5K", "14K", "20.5K", "30K",
                                                "42.5K", "62.5K", "87.5K", "150K", "200KP") )

highest.education = rep("999", length(abcd_data$demo_prnt_ed_p))
highest.education[abcd_data$demo_prnt_ed_p == "Never attended/Kindergarten only"] = 0
highest.education[abcd_data$demo_prnt_ed_p == "1st grade"] = 1
highest.education[abcd_data$demo_prnt_ed_p == "2nd grade"] = 2
highest.education[abcd_data$demo_prnt_ed_p == "3rd grade"] = 3
highest.education[abcd_data$demo_prnt_ed_p == "4th grade"] = 4
highest.education[abcd_data$demo_prnt_ed_p == "5th grade"] = 5
highest.education[abcd_data$demo_prnt_ed_p == "6th grade"] = 6
highest.education[abcd_data$demo_prnt_ed_p == "7th grade"] = 7
highest.education[abcd_data$demo_prnt_ed_p == "8th grade"] = 8
highest.education[abcd_data$demo_prnt_ed_p == "9th grade"] = 9
highest.education[abcd_data$demo_prnt_ed_p == "10th grade"] = 10
highest.education[abcd_data$demo_prnt_ed_p == "11th grade"] = 11
highest.education[(abcd_data$demo_prnt_ed_p == "12th grade, no diploma")] = 12
highest.education[(abcd_data$demo_prnt_ed_p == "High school graduate")] = 13
highest.education[abcd_data$demo_prnt_ed_p == "GED or equivalent"] = 14
highest.education[abcd_data$demo_prnt_ed_p == "Some college, no degree"] = 15
highest.education[(abcd_data$demo_prnt_ed_p == "Associate degree: Occupational, Technical, or Vocational")] = 16
highest.education[(abcd_data$demo_prnt_ed_p == "Associate degree: Academic Program")] = 17
highest.education[abcd_data$demo_prnt_ed_p == "Bachelor's degree (ex. BA, AB, BS, BBS)"] = 18
highest.education[abcd_data$demo_prnt_ed_p == "Master's degree (ex. MA, MS, MEng, MEd, MBA)"] = 19
highest.education[abcd_data$demo_prnt_ed_p == "Professional School degree (ex. MD, DDS, DVN, JD)"] = 20
highest.education[abcd_data$demo_prnt_ed_p == "Doctoral degree (ex. PhD, EdD)"] = 21
highest.education[abcd_data$demo_prnt_ed_p == "Refused to answer"] = 999
highest.education[highest.education == 999] = NA

highest.education2 = rep("999", length(abcd_data$demo_prtnr_ed_p))
highest.education2[abcd_data$demo_prnt_ed_p == "Never attended/Kindergarten only"] = 0
highest.education2[abcd_data$demo_prnt_ed_p == "1st grade"] = 1
highest.education2[abcd_data$demo_prnt_ed_p == "2nd grade"] = 2
highest.education2[abcd_data$demo_prnt_ed_p == "3rd grade"] = 3
highest.education2[abcd_data$demo_prnt_ed_p == "4th grade"] = 4
highest.education2[abcd_data$demo_prnt_ed_p == "5th grade"] = 5
highest.education2[abcd_data$demo_prnt_ed_p == "6th grade"] = 6
highest.education2[abcd_data$demo_prnt_ed_p == "7th grade"] = 7
highest.education2[abcd_data$demo_prnt_ed_p == "8th grade"] = 8
highest.education2[abcd_data$demo_prnt_ed_p == "9th grade"] = 9
highest.education2[abcd_data$demo_prnt_ed_p == "10th grade"] = 10
highest.education2[abcd_data$demo_prnt_ed_p == "11th grade"] = 11
highest.education2[abcd_data$demo_prnt_ed_p == "12th grade, no diploma"] = 12 
highest.education2[abcd_data$demo_prnt_ed_p == "High school graduate"] = 13
highest.education2[abcd_data$demo_prnt_ed_p == "GED or equivalent"] = 14
highest.education2[abcd_data$demo_prnt_ed_p == "Some college, no degree"] = 15
highest.education2[abcd_data$demo_prnt_ed_p == "Associate degree: Occupational, Technical, or Vocational"] = 16 
highest.education2[abcd_data$demo_prnt_ed_p == "Associate degree: Academic Program"] = 17
highest.education2[abcd_data$demo_prnt_ed_p == "Bachelor's degree (ex. BA, AB, BS, BBS)"] = 18
highest.education2[abcd_data$demo_prnt_ed_p == "Master's degree (ex. MA, MS, MEng, MEd, MBA)"] = 19
highest.education2[abcd_data$demo_prnt_ed_p == "Professional School degree (ex. MD, DDS, DVN, JD)"] = 20
highest.education2[abcd_data$demo_prnt_ed_p == "Doctoral degree (ex. PhD, EdD)"] = 21
highest.education2[abcd_data$demo_prnt_ed_p == "Refused to answer"] = 999
highest.education2[highest.education2 == 999] = NA

abcd_data$highest.education = (pmax(as.numeric(highest.education), as.numeric(highest.education2),na.rm=T)) 

high.educ1 = highest.education
high.educ2 = highest.education2
high.educ1[which(high.educ1 == "999")] = NA
high.educ2[which(high.educ2 == "999")] = NA
high.educ1[which(high.educ1 == "777")] = NA
high.educ2[which(high.educ2 == "777")] = NA
high.educ = pmax(as.numeric(as.character(high.educ1)), as.numeric(as.character(high.educ2)), na.rm=T)
idx <- which(high.educ %in% 0:12, arr.ind = TRUE)
high.educ[idx] = 1 # "< HS Diploma"
idx <- which(high.educ %in% 13:14, arr.ind = TRUE)
high.educ[idx] = 2 # "HS Diploma/GED"
idx <- which(high.educ %in% 15:17, arr.ind = TRUE)
high.educ[idx] = 3 # "Some College"
idx <- which(high.educ == 18, arr.ind = TRUE)
high.educ[idx] = 4 # "Bachelor"
idx <- which(high.educ %in% 19:21, arr.ind = TRUE)
high.educ[idx] = 5 # "Post Graduate Degree"
high.educ[which(high.educ == "999")]=NA
high.educ[which(high.educ == "777")]=NA
abcd_data$high.educ = factor( high.educ, levels= 1:5, labels = c("LTHSDiploma","HSDiplomaGED","SomeCollege","Bachelor","PostGradDegree") )

household.income_num <- paste(abcd_data$household.income)
household.income_num[household.income_num == "LT5K"] = 5000 #  Less than $5,000
household.income_num[household.income_num == "8.5K"] = 8500 # $5,000 through $11,999
household.income_num[household.income_num == "14K"] = 14000 # $12,000 through $15,999
household.income_num[household.income_num == "20.5K"] = 20500 # $16,000 through $24,999
household.income_num[household.income_num == "30K"] = 30000 # $25,000 through $34,999
household.income_num[household.income_num == "42.5K"] = 42500 # $35,000 through $49,999
household.income_num[household.income_num == "62.5K"] = 62500 # $50,000 through $74,999
household.income_num[household.income_num == "87.5K"] = 87500 # $75,000 through $99,999
household.income_num[household.income_num == "150K"] = 150000 # $100,000 through $199,999
household.income_num[household.income_num == "200KP"] = 200000 # $200,000 and greater
abcd_data$household.income_num <- as.numeric(household.income_num)

## Gather and rename relevant variables ## 

abcd_data_final <- abcd_data %>%
  dplyr::select(subjectid, abcd_site, rel_family_id, age, male, race_ethnicity_White, race_ethnicity_Black, race_ethnicity_Hispanic, race_ethnicity_Asian,
                race_ethnicity_Other, household.income_num, highest.education, meim_ss_exp_p, meim_ss_com_p, macvs_ss_fo_p, macvs_ss_fr_p,
                macvs_ss_fs_p, macvs_ss_isr_p, macvs_ss_r_p, fes_ss_fc_p_pr, parental_monitoring_ss_mean, fes_ss_fc_pr, crpbi_acceptance_ss_studycaregiver,
                school_risk_phenx_ss_ses, school_risk_phenx_ss_iiss, school_risk_phenx_ss_dfs, neighb_phenx_ss_mean_p, neighb_phenx, reshist_addr1_walkindex,
                reshist_addr1_no2, reshist_addr1_pm25, ProxRoads_log, Res_density_log, UCR_ViolCrime_log, UCR_DrugAbuse_log, UCR_DrugSale_log, UCR_DrugPoss_log,
                UCR_DUI_log, UCR_MJSale_log, reshist_addr1_adi_edu_l, reshist_addr1_adi_edu_h, reshist_addr1_adi_work_c, reshist_addr1_adi_income, 
                reshist_addr1_adi_in_dis, reshist_addr1_adi_home_v, reshist_addr1_adi_rent, reshist_addr1_adi_mortg, reshist_addr1_adi_home_o, reshist_addr1_adi_unemp,
                reshist_addr1_adi_pov, reshist_addr1_adi_b138, reshist_addr1_adi_sp, reshist_addr1_adi_ncar, reshist_addr1_adi_ntel, reshist_addr1_adi_nplumb,
                reshist_addr1_adi_crowd, nihtbx_picvocab_uncorrected, nihtbx_flanker_uncorrected, nihtbx_list_uncorrected, nihtbx_cardsort_uncorrected, 
                nihtbx_pattern_uncorrected, nihtbx_picture_uncorrected, nihtbx_reading_uncorrected, pea_ravlt_ld, lmt_scr_perc_correct, neurocog_pc1.bl, neurocog_pc2.bl, neurocog_pc3.bl)

colnames(abcd_data_final) <- c("subject_id",  "abcd_site", "family_id",  "Age", "Male", "White", "Black", "Hispanic", "Asian", "Other", "Income",
                               "Education", "EthIdentExploration", "EthIdentCommitment",  "FamObligation", "FamReferent", "FamSupport", "IndepSelfReliance",
                               "Religion", "FamConflict_P", "ParentMonitoring", "FamConflict_Y",  "ParentAcceptance", "SchoolEnvironment", "SchoolInvolvement",
                               "SchoolDisengagement", "NeighSafety_P",  "NeighSafety_Y",  "Walkability", "NO2Exposure", "PM25Exposure", "ProxRoads_log",
                               "ResDensity_log", "UCR_ViolCrime_log", "UCR_DrugAbuse_log",  "UCR_DrugSale_log", "UCR_DrugPoss_log", "UCR_DUI_log",
                               "UCR_MJSale_log", "ADI_Edu_LTHS",  "ADI_Edu_HSDip", "ADI_Occ_WhiteCollar", "ADI_MedianFamInc", "ADI_IncDisparityIdx",
                               "ADI_MedianHomeValue", "ADI_MedianGrossRent",  "ADI_MedianMonthMortg", "ADI_PercHOwner", "ADI_PercUnemp",
                               "ADI_PercFamPoverty", "ADI_PercPovBelow138",  "ADI_PercSingle", "ADI_PercLogNCar", "ADI_PercLogNTel",
                               "ADI_PercLogNPlumb", "ADI_PercCrowding", "NIHtb_PicVocab", "NIHtb_Flanker", "NIHtb_List", "NIHtb_CardSort",
                               "NIHtb_Pattern", "NIHtb_Picture", "NIHtb_Reading", "RAVLT", "LMT", "Neurocog_BPPC1_GenAbility",
                               "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")

nmissing_abcd_data_final <- colSums(is.na(abcd_data_final))

abcd_data_final <- abcd_data_final[complete.cases(abcd_data_final),]
write_csv(abcd_data_final, "abcd_data_final.csv")

abcd_data_final_loginc <- abcd_data_final %>%
  mutate(Income = log(Income))

## Create descriptive statistics table with missingness ## 

abcd_stats_total <- abcd_data %>%
  dplyr::select(subjectid, abcd_site, age, male, race_ethnicity, household.income_num, highest.education, 
                #home 
                fes_ss_fc_p_pr, parental_monitoring_ss_mean, fes_ss_fc_pr, crpbi_acceptance_ss_studycaregiver,
                #school
                school_risk_phenx_ss_ses, school_risk_phenx_ss_iiss, school_risk_phenx_ss_dfs, 
                #neighborhood
                neighb_phenx_ss_mean_p, neighb_phenx, reshist_addr1_walkindex,
                reshist_addr1_no2, reshist_addr1_pm25, ProxRoads_log, Res_density_log, UCR_ViolCrime_log, UCR_DrugAbuse_log, UCR_DrugSale_log, UCR_DrugPoss_log,
                UCR_DUI_log, UCR_MJSale_log, reshist_addr1_adi_edu_l, reshist_addr1_adi_edu_h, reshist_addr1_adi_work_c, reshist_addr1_adi_income, 
                reshist_addr1_adi_in_dis, reshist_addr1_adi_home_v, reshist_addr1_adi_rent, reshist_addr1_adi_mortg, reshist_addr1_adi_home_o, reshist_addr1_adi_unemp,
                reshist_addr1_adi_pov, reshist_addr1_adi_b138, reshist_addr1_adi_sp, reshist_addr1_adi_ncar, reshist_addr1_adi_ntel, reshist_addr1_adi_nplumb,
                reshist_addr1_adi_crowd,
                
                #culture
                meim_ss_exp_p, meim_ss_com_p, macvs_ss_fo_p, macvs_ss_fr_p, macvs_ss_fs_p, macvs_ss_isr_p, macvs_ss_r_p,
                # neurocog 
                nihtbx_picvocab_uncorrected, nihtbx_flanker_uncorrected, nihtbx_list_uncorrected, nihtbx_cardsort_uncorrected, 
                nihtbx_pattern_uncorrected, nihtbx_picture_uncorrected, nihtbx_reading_uncorrected, pea_ravlt_ld, lmt_scr_perc_correct, neurocog_pc1.bl, neurocog_pc2.bl, neurocog_pc3.bl) %>%
  mutate(male = factor(male, levels = c(0,1), labels = c("female", "male")))


colnames(abcd_stats_total) <- c("Subject_ID", "ABCD Site", "Age", "Sex", "Race/Ethnicity", "Household Income", "Caregiver Education", 
                                
                                "Family Conflict (Parent)","Parent Monitoring", "Family Conflict (Youth)", "Parent Acceptance", 
                                
                                "School Environment", "School Involvement", "School Disengagement",
                                
                                "Neighborhood Safety (Parent)", "Neighborhood Safety (Youth)", "Walkability", "NO2 Exposure", "PM25 Exposure","Proximity to Roads (log)",
                                "Residential Density (log)", "Violent Crime Rate (log)", "Drug Abuse Rate (log)", "Drug Sales Rate (log)", "Drug Possession Rate (log)",
                                "DUI Rate (log)", "Marijuana Sales Rate (log)", "Education < HS", "Education - HS Diploma", "Occupation- White Collar",
                                "Median Family Income",  "Income Disparity Index", "Median Home Value", "Median Gross Rent", "Median Monthly Mortgage",  
                                "Home Ownership (%)","Unemployed (%)","Families in Poverty (%)", "Families Below 138% Povery Line (%)","Single/Not Married (%)",
                                "Car Ownership (%, log)", "Telephone in Home (%, log)", "Plumbing in Home (%, log)", "Overcrowding (%)", 
                                "Ethnic Identy Exploration", "Ethnic Identity Commitment", "Family Obligation", "Family as Referent", "Family Support",
                                "Indepence/SelfReliance", "Religion", 
                                "NIHtb Pic Vocab", "NIHtb Flanker", "NIHtb List Sorting", 
                                "NIHtb Card Sort", "NIHtb Pattern Recognition", "NIHtb Picture", "NIHtb Reading", "RAVLT", "Little Man", "General Cog Ability", 
                                "Executive Functioning", "Learning/Memory")

## releqse split stats ## 

abcd_stats_total_rel1 <- abcd_stats_total %>%
  filter(Subject_ID %in% first_release_subjects[[1]]) %>%
  dplyr::select(2:ncol(abcd_stats_total))

abcd_stats_total_rel2 <- abcd_stats_total %>%
  filter(Subject_ID %nin% first_release_subjects[[1]]) %>%
  dplyr::select(2:ncol(abcd_stats_total))

abcd_stats_totaltab_rel1 <- CreateTableOne(data = abcd_stats_total_rel1)
abcd_stats_summarytab_rel1 <- summary(abcd_stats_totaltab_rel1)
sink("abcd_missingtable_rel1.csv")
print(summary(abcd_stats_totaltab_rel1), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
sink()

abcd_stats_totaltab_rel2 <- CreateTableOne(data = abcd_stats_total_rel2)
sink("abcd_missingtable_rel2.csv")
print(summary(abcd_stats_totaltab_rel2), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
sink()

### Complete case release split stats ## 

abcd_stats_rel1_complete <- abcd_stats_total_rel1[complete.cases(abcd_stats_total_rel1),]
abcd_stats_rel2_complete <- abcd_stats_total_rel2[complete.cases(abcd_stats_total_rel2),]

(abcd_stats_completetab_rel1 <- CreateTableOne(data = abcd_stats_rel1_complete))
stats_table_rel1 <- print(abcd_stats_completetab_rel1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(stats_table_rel1, file = "abcd_stats_table_rel1.csv")

(abcd_stats_completetab_rel2 <- CreateTableOne(data = abcd_stats_rel2_complete))
stats_table_rel2 <- print(abcd_stats_completetab_rel2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(stats_table_rel2, file = "abcd_stats_table_rel2.csv")

## full sample stats ##

abcd_stats_full <- abcd_stats_total %>%
  dplyr::select(2:ncol(abcd_stats_total))

(abcd_stats_completetab_full <- CreateTableOne(data = abcd_stats_full))
stats_table_full <- print(abcd_stats_completetab_full, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(stats_table_full, file = "abcd_stats_table_full.csv")

### full sample complete case stats ##

abcd_stats_fullcomplete <- abcd_stats_full[complete.cases(abcd_stats_full),]

(abcd_stats_completetab_fullcomplete <- CreateTableOne(data = abcd_stats_fullcomplete))
stats_table_fullcomplete <- print(abcd_stats_completetab_fullcomplete, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(stats_table_fullcomplete, file = "abcd_stats_table_fullcomplete.csv")

### excluded subsample stats ##

abcd_stats_full_subsample <- abcd_stats_total %>%
  filter(Subject_ID %nin% abcd_data_final$subject_id) %>%
  dplyr::select(-c("Subject_ID"))

(abcd_stats_completetab_excluded <- CreateTableOne(data = abcd_stats_full_subsample))
stats_table_excluded <- print(abcd_stats_completetab_excluded, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(stats_table_excluded, file = "abcd_stats_table_excluded.csv")


## compare mean scores across full and complete case samples ##

abcd_stats_ttest <- data.table(Variable = colnames(abcd_stats_full),
                               tstat = NA,
                               df = NA,
                               pval = NA) 

for (i in 1:ncol(abcd_stats_full)){
  if (i %in% c(1,3,4)){
    next
  }
  test <- t.test(abcd_stats_full[,i], abcd_stats_fullcomplete[,i], paired = FALSE, var.equal = FALSE) 
  temp_tstat <- test$statistic
  temp_pval <- test$p.value
  temp_df <- test$parameter
  abcd_stats_ttest$tstat[i] <- temp_tstat
  abcd_stats_ttest$df[i] <- temp_df
  abcd_stats_ttest$pval[i] <- temp_pval
}

rm(temp_tstat, temp_pval, temp_df)
write.csv(abcd_stats_ttest, file = "abcd_stats_ttable.csv")

## 2. Gather necessary variables for 4 environment categories, sociodemographic, and neurocognitive BPPCA weights ##

CultureEnvVar_list <- data.table(c("EthIdentExploration", "EthIdentCommitment", "FamObligation", "FamReferent", 
                                   "FamSupport", "IndepSelfReliance", "Religion"), "Culture")

HomeEnvVar_list <- data.table(c("FamConflict_P", "ParentMonitoring", "FamConflict_Y", "ParentAcceptance"), "Home")

SchoolEnvVar_list <- data.table(c("SchoolEnvironment", "SchoolInvolvement", "SchoolDisengagement"), "School")

NeighborhoodEnvVar_list <- data.table(c("NeighSafety_P", "NeighSafety_Y", "Walkability", "NO2Exposure", "PM25Exposure",
                                        "ProxRoads_log", "ResDensity_log", "UCR_ViolCrime_log",  "UCR_DrugAbuse_log",  "UCR_DrugSale_log", 
                                        "UCR_DrugPoss_log", "UCR_DUI_log", "UCR_MJSale_log", "ADI_Edu_LTHS", "ADI_Edu_HSDip", 
                                        "ADI_Occ_WhiteCollar",  "ADI_MedianFamInc", "ADI_IncDisparityIdx",  "ADI_MedianHomeValue", 
                                        "ADI_MedianGrossRent", "ADI_MedianMonthMortg", "ADI_PercHOwner", "ADI_PercUnemp", "ADI_PercFamPoverty", 
                                        "ADI_PercPovBelow138", "ADI_PercSingle", "ADI_PercLogNCar", "ADI_PercLogNTel", "ADI_PercLogNPlumb",
                                        "ADI_PercCrowding"), "Neighborhood")

NeurocogVar_list <- data.table(c("NIHtb_PicVocab", "NIHtb_Flanker","NIHtb_List", "NIHtb_CardSort", "NIHtb_Pattern",
                                 "NIHtb_Picture", "NIHtb_Reading", "RAVLT", "LMT"), "Neurocog")

SociodemoVar_list <- data.table(c("White", "Black", "Hispanic", "Asian", "Other", "Income","Education"), "Sociodemo")

## Create subject splits based on release (1.1 and 2.0.1), randomized site groups, and DCAN matched splits ##

first_release_data <- abcd_data_final %>%
  filter(subject_id %in% first_release_subjects[[1]]) %>%
  standardize()

first_release_bppca <- abcd_data_final %>%
  filter(subject_id %in% first_release_subjects[[1]]) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

second_release_data <- abcd_data_final %>%
  filter(subject_id %nin% first_release_subjects[[1]]) %>%
  standardize()

second_release_bppca <- abcd_data_final %>%
  filter(subject_id %nin% first_release_subjects[[1]]) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

## create site splits ##

set.seed(1992)

nrand <- 100
n_sites <- length(unique(abcd_data_final$abcd_site))
n_site_split <- n_sites/2
site_list <- as.character(unique(abcd_data_final$abcd_site))
site_splits_rand <- data.frame(matrix(nrow = n_site_split, ncol = nrand))

for (i in 1:nrand){
  tmp_sitesplit <- sample(site_list, n_site_split, replace = FALSE, prob = NULL)
  site_splits_rand[,i] <- tmp_sitesplit
  
  rm(tmp_sitesplit)
}

## even splits using DCAN labs' matched splits ##

matched_splits <- read_tsv(dcan_lab_matched_split_tsv_path) %>%
  dplyr::select(participant_id, matched_group) %>%
  mutate(participant_id = substring(participant_id, 9)) %>%
  distinct()

matched_data <- abcd_data_final %>%
  mutate(subject_id = as.character(subject_id),
         subject_id = substring(subject_id,6)) %>%
  left_join(matched_splits, by = c("subject_id" = "participant_id"))

matched1_data <- matched_data %>%
  filter(matched_group == 1) %>%
  dplyr::select(-c("matched_group")) %>%
  standardize()

matched1_bppca <- matched_data %>%
  filter(matched_group == 1) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

matched2_data <- matched_data %>%
  filter(matched_group == 2) %>%
  dplyr::select(-c("matched_group")) %>%
  standardize()

matched2_bppca <- matched_data %>%
  filter(matched_group == 2) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

matched3_data <- matched_data %>%
  filter(matched_group == 3) %>%
  dplyr::select(-c("matched_group")) %>%
  standardize()

matched3_bppca <- matched_data %>%
  filter(matched_group == 3) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

## 3. Run PCA for environmental factors, sociodemographic factors, and neurocognitive scores in each subsample ##

## Release 1 and Release 2 group PCAs ##

## Release 1 ##

rel1_culture_pca     <- PCA(first_release_data[CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_home_pca        <- PCA(first_release_data[HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_school_pca      <- PCA(first_release_data[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_neighborhood_pca <- PCA(first_release_data[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_sociodemo_pca    <- PCA(first_release_data[SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_neurocog_pca    <- PCA(first_release_data[NeurocogVar_list$V1], scale.unit = TRUE, graph = TRUE)

first_release_data <- first_release_data %>%
  mutate(culture_PC1 = -1 * rel1_culture_pca$ind$coord[,1],
         home_PC1 = rel1_home_pca$ind$coord[,1],
         school_PC1 = rel1_school_pca$ind$coord[,1],
         neighborhood_PC1 = -1 * rel1_neighborhood_pca$ind$coord[,1],
         sociodemo_PC1 = rel1_sociodemo_pca$ind$coord[,1],
         neurocog_PC1 = rel1_neurocog_pca$ind$coord[,1],
         neurocog_PC2 = rel1_neurocog_pca$ind$coord[,2],
         neurocog_PC3 = rel1_neurocog_pca$ind$coord[,3], 
         Neurocog_BPPC1_GenAbility = first_release_bppca$Neurocog_BPPC1_GenAbility,
         Neurocog_BPPC2_ExecFunc = first_release_bppca$Neurocog_BPPC2_ExecFunc,
         Neurocog_BPPC3_LearnMem = first_release_bppca$Neurocog_BPPC3_LearnMem)

# write_csv(first_release_data, "abcd_enviro_firstrelease.csv")

## Release 2 ##

rel2_culture_pca     <- PCA(second_release_data[CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_home_pca        <- PCA(second_release_data[HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_school_pca      <- PCA(second_release_data[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_neighborhood_pca <- PCA(second_release_data[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_sociodemo_pca    <- PCA(second_release_data[SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_neurocog_pca         <- PCA(second_release_data[NeurocogVar_list$V1], scale.unit = TRUE, graph = TRUE)

second_release_data <- second_release_data %>%
  mutate(culture_PC1 = -1 * rel2_culture_pca$ind$coord[,1],
         home_PC1 = rel2_home_pca$ind$coord[,1],
         school_PC1 = rel2_school_pca$ind$coord[,1],
         neighborhood_PC1 = -1 * rel2_neighborhood_pca$ind$coord[,1],
         sociodemo_PC1 = rel2_sociodemo_pca$ind$coord[,1],
         neurocog_PC1 = rel2_neurocog_pca$ind$coord[,1],
         neurocog_PC2 = rel2_neurocog_pca$ind$coord[,2],
         neurocog_PC3 = rel2_neurocog_pca$ind$coord[,3],
         Neurocog_BPPC1_GenAbility = second_release_bppca$Neurocog_BPPC1_GenAbility,
         Neurocog_BPPC2_ExecFunc = second_release_bppca$Neurocog_BPPC2_ExecFunc,
         Neurocog_BPPC3_LearnMem = second_release_bppca$Neurocog_BPPC3_LearnMem)

# write_csv(second_release_data, "abcd_enviro_secondrelease.csv")

## Random site groups 1 and 2 PCAs ## 

site1_culture_PC_list <- vector("list", nrand)
site1_home_PC_list <- vector("list", nrand)
site1_school_PC_list <- vector("list", nrand)
site1_neighborhood_PC_list <- vector("list", nrand)
site1_sociodemo_PC_list <- vector("list", nrand)

site2_culture_PC_list <- vector("list", nrand)
site2_home_PC_list <- vector("list", nrand)
site2_school_PC_list <- vector("list", nrand)
site2_neighborhood_PC_list <- vector("list", nrand)
site2_sociodemo_PC_list <- vector("list", nrand)

for (i in 1:nrand){
  
  site1_data <- abcd_data_final %>%
    filter(abcd_site %in% site_splits_rand[,i]) %>%
    standardize()
  
  site2_data <- abcd_data_final %>%
    filter(abcd_site %nin% site_splits_rand[,i]) %>%
    standardize()
  
  site1_culture_PC_list[[i]]    <- PCA(site1_data[CultureEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site1_home_PC_list[[i]]    <- PCA(site1_data[HomeEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site1_school_PC_list[[i]]    <- PCA(site1_data[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site1_neighborhood_PC_list[[i]]    <- PCA(site1_data[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site1_sociodemo_PC_list[[i]]    <- PCA(site1_data[SociodemoVar_list$V1], scale.unit = TRUE, graph = FALSE)
  
  site2_culture_PC_list[[i]]    <- PCA(site2_data[CultureEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site2_home_PC_list[[i]]    <- PCA(site2_data[HomeEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site2_school_PC_list[[i]]    <- PCA(site2_data[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site2_neighborhood_PC_list[[i]]    <- PCA(site2_data[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = FALSE)
  site2_sociodemo_PC_list[[i]]    <- PCA(site2_data[SociodemoVar_list$V1], scale.unit = TRUE, graph = FALSE)
  
}

#Gather variable loadings for each site split ## 

site1_culture_varload <- data.frame(matrix(NA, nrow = 7, ncol = 100), row.names = CultureEnvVar_list$V1)
site1_home_varload <- data.frame(matrix(NA, nrow = 4, ncol = 100), row.names = HomeEnvVar_list$V1)
site1_school_varload <- data.frame(matrix(NA, nrow = 3, ncol = 100), row.names = SchoolEnvVar_list$V1)
site1_neighborhood_varload <- data.frame(matrix(NA, nrow = 30, ncol = 100), row.names = NeighborhoodEnvVar_list$V1)
site1_sociodemo_varload <- data.frame(matrix(NA, nrow = 7, ncol = 100), row.names = SociodemoVar_list$V1)

site2_culture_varload <- data.frame(matrix(NA, nrow = 7, ncol = 100), row.names = CultureEnvVar_list$V1)
site2_home_varload <- data.frame(matrix(NA, nrow = 4, ncol = 100), row.names = HomeEnvVar_list$V1)
site2_school_varload <- data.frame(matrix(NA, nrow = 3, ncol = 100), row.names = SchoolEnvVar_list$V1)
site2_neighborhood_varload <- data.frame(matrix(NA, nrow = 30, ncol = 100), row.names = NeighborhoodEnvVar_list$V1)
site2_sociodemo_varload <- data.frame(matrix(NA, nrow = 7, ncol = 100), row.names = SociodemoVar_list$V1)

for (i in 1:nrand){
  site1_culture_varload[,i] <- site1_culture_PC_list[[i]]$var$contrib[,1]
  site1_home_varload[,i] <- site1_home_PC_list[[i]]$var$contrib[,1]
  site1_school_varload[,i] <- site1_school_PC_list[[i]]$var$contrib[,1]
  site1_neighborhood_varload[,i] <- site1_neighborhood_PC_list[[i]]$var$contrib[,1]
  site1_sociodemo_varload[,i] <- site1_sociodemo_PC_list[[i]]$var$contrib[,1]
  
  site2_culture_varload[,i] <- site2_culture_PC_list[[i]]$var$contrib[,1]
  site2_home_varload[,i] <- site2_home_PC_list[[i]]$var$contrib[,1]
  site2_school_varload[,i] <- site2_school_PC_list[[i]]$var$contrib[,1]
  site2_neighborhood_varload[,i] <- site2_neighborhood_PC_list[[i]]$var$contrib[,1]
  site2_sociodemo_varload[,i] <- site2_sociodemo_PC_list[[i]]$var$contrib[,1]
  
}

## grab the variance explained for site split PCS ##
site1_home_variance <- c(rep(NA, 100))
site1_school_variance <- c(rep(NA, 100))
site1_neighborhood_variance <- c(rep(NA, 100))
site1_culture_variance <- c(rep(NA, 100))
site1_sociodemo_variance <- c(rep(NA, 100))

site2_home_variance <- c(rep(NA, 100))
site2_school_variance <- c(rep(NA, 100))
site2_neighborhood_variance <- c(rep(NA, 100))
site2_culture_variance <- c(rep(NA, 100))
site2_sociodemo_variance <- c(rep(NA, 100))

for (i in 1:nrand){
  site1_home_variance[i] <- site1_home_PC_list[[i]]$eig[1,2]
  site1_school_variance[i] <- site1_school_PC_list[[i]]$eig[1,2]
  site1_neighborhood_variance[i] <- site1_neighborhood_PC_list[[i]]$eig[1,2]
  site1_culture_variance[i] <- site1_culture_PC_list[[i]]$eig[1,2]
  site1_sociodemo_variance[i] <- site1_sociodemo_PC_list[[i]]$eig[1,2]
  
  site2_home_variance[i] <- site2_home_PC_list[[i]]$eig[1,2]
  site2_school_variance[i] <- site2_school_PC_list[[i]]$eig[1,2]
  site2_neighborhood_variance[i] <- site2_neighborhood_PC_list[[i]]$eig[1,2]
  site2_culture_variance[i] <- site2_culture_PC_list[[i]]$eig[1,2]
  site2_sociodemo_variance[i] <- site2_sociodemo_PC_list[[i]]$eig[1,2]
}

site1_home_avgvariance <- mean(site1_home_variance)
site1_school_avgvariance <- mean(site1_school_variance)
site1_neighborhood_avgvariance <- mean(site1_neighborhood_variance)
site1_culture_avgvariance <- mean(site1_culture_variance)
site1_sociodemo_avgvariance <- mean(site1_sociodemo_variance)

site2_home_avgvariance <- mean(site2_home_variance)
site2_school_avgvariance <- mean(site2_school_variance)
site2_neighborhood_avgvariance <- mean(site2_neighborhood_variance)
site2_culture_avgvariance <- mean(site2_culture_variance)
site2_sociodemo_avgvariance <- mean(site2_sociodemo_variance)

(sitesplit_home_avgvariance <- mean(site1_home_avgvariance, site2_home_avgvariance))
(sitesplit_school_avgvariance <- mean(site1_school_avgvariance, site2_school_avgvariance))
(sitesplit_neighborhood_avgvariance <- mean(site1_neighborhood_avgvariance, site2_neighborhood_avgvariance))
(sitesplit_culture_avgvariance <- mean(site1_culture_avgvariance, site2_culture_avgvariance))
(sitesplit_sociodemo_avgvariance <- mean(site1_sociodemo_avgvariance, site2_sociodemo_avgvariance))

## Matched-group PCA ##

## Matched group 1 ## 

matched1_culture_pca     <- PCA(matched1_data[CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched1_home_pca        <- PCA(matched1_data[HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched1_school_pca      <- PCA(matched1_data[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched1_neighborhood_pca <- PCA(matched1_data[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched1_sociodemo_pca    <- PCA(matched1_data[SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched1_neurocog_pca    <- PCA(matched1_data[NeurocogVar_list$V1], scale.unit = TRUE, graph = TRUE)

matched1_data <- matched1_data %>%
  mutate(culture_PC1 = -1 * matched1_culture_pca$ind$coord[,1],
         home_PC1 = matched1_home_pca$ind$coord[,1],
         school_PC1 = matched1_school_pca$ind$coord[,1],
         neighborhood_PC1 = -1 * matched1_neighborhood_pca$ind$coord[,1],
         sociodemo_PC1 = matched1_sociodemo_pca$ind$coord[,1],
         neurocog_PC1 = matched1_neurocog_pca$ind$coord[,1],
         neurocog_PC2 = matched1_neurocog_pca$ind$coord[,2],
         neurocog_PC3 = matched1_neurocog_pca$ind$coord[,3], 
         Neurocog_BPPC1_GenAbility = matched1_bppca$Neurocog_BPPC1_GenAbility,
         Neurocog_BPPC2_ExecFunc = matched1_bppca$Neurocog_BPPC2_ExecFunc,
         Neurocog_BPPC3_LearnMem = matched1_bppca$Neurocog_BPPC3_LearnMem)

## Matched group 2 ##

matched2_culture_pca     <- PCA(matched2_data[CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched2_home_pca        <- PCA(matched2_data[HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched2_school_pca      <- PCA(matched2_data[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched2_neighborhood_pca <- PCA(matched2_data[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched2_sociodemo_pca    <- PCA(matched2_data[SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)
matched2_neurocog_pca    <- PCA(matched2_data[NeurocogVar_list$V1], scale.unit = TRUE, graph = TRUE)

matched2_data <- matched2_data %>%
  mutate(culture_PC1 = -1 * matched2_culture_pca$ind$coord[,1],
         home_PC1 = matched2_home_pca$ind$coord[,1],
         school_PC1 = matched2_school_pca$ind$coord[,1],
         neighborhood_PC1 = -1 * matched2_neighborhood_pca$ind$coord[,1],
         sociodemo_PC1 = matched2_sociodemo_pca$ind$coord[,1],
         neurocog_PC1 = matched2_neurocog_pca$ind$coord[,1],
         neurocog_PC2 = matched2_neurocog_pca$ind$coord[,2],
         neurocog_PC3 = matched2_neurocog_pca$ind$coord[,3], 
         Neurocog_BPPC1_GenAbility = matched2_bppca$Neurocog_BPPC1_GenAbility,
         Neurocog_BPPC2_ExecFunc = matched2_bppca$Neurocog_BPPC2_ExecFunc,
         Neurocog_BPPC3_LearnMem = matched2_bppca$Neurocog_BPPC3_LearnMem)

## Scree plots for Release 1 and Release 2 groups ##

## Release 1.1 ## 

fviz_screeplot(rel1_culture_pca, addlabels =TRUE) +
  ggtitle("Release 1 Culture")
fviz_screeplot(rel1_home_pca, addlabels =TRUE) +
  ggtitle("Release 1 Home")
fviz_screeplot(rel1_school_pca, addlabels =TRUE) +
  ggtitle("Release 1 School")
fviz_screeplot(rel1_neighborhood_pca, addlabels =TRUE) +
  ggtitle("Release 1 Neighborhood")
fviz_screeplot(rel1_sociodemo_pca, addlabels =TRUE) +
  ggtitle("Release 1 Sociodemo")
fviz_screeplot(rel1_neurocog_pca, addlabels =TRUE) +
  ggtitle("Release 1 Neurocog")

## Release 2.0.1 ## 

fviz_screeplot(rel2_culture_pca, addlabels =TRUE) +
  ggtitle("Release 2 Culture")
fviz_screeplot(rel2_home_pca, addlabels =TRUE) +
  ggtitle("Release 2 Home")
fviz_screeplot(rel2_school_pca, addlabels =TRUE) +
  ggtitle("Release 2 School")
fviz_screeplot(rel2_neighborhood_pca, addlabels =TRUE) +
  ggtitle("Release 2 Neighborhood")
fviz_screeplot(rel2_sociodemo_pca, addlabels =TRUE) +
  ggtitle("Release 2 Sociodemo")
fviz_screeplot(rel2_neurocog_pca, addlabels =TRUE) +
  ggtitle("Release 2 Neurocog")

## Scree plots for matched groups ##

## matched group 1 ##
fviz_screeplot(matched1_culture_pca, addlabels =TRUE) +
  ggtitle("Matched Group 1 Culture")
fviz_screeplot(matched1_home_pca, addlabels =TRUE) +
  ggtitle("Matched Group 1 Home")
fviz_screeplot(matched1_school_pca, addlabels =TRUE) +
  ggtitle("Matched Group 1 School")
fviz_screeplot(matched1_neighborhood_pca, addlabels =TRUE) +
  ggtitle("Matched Group 1 Neighborhood")
fviz_screeplot(matched1_sociodemo_pca, addlabels =TRUE) +
  ggtitle("Matched Group 1 Sociodemo")
fviz_screeplot(matched1_neurocog_pca, addlabels =TRUE) +
  ggtitle("Matched Group 1 Neurocog")

## matched group 2 ## 

fviz_screeplot(matched2_culture_pca, addlabels =TRUE) +
  ggtitle("Matched Group 2 Culture")
fviz_screeplot(matched2_home_pca, addlabels =TRUE) +
  ggtitle("Matched Group 2 Home")
fviz_screeplot(matched2_school_pca, addlabels =TRUE) +
  ggtitle("Matched Group 2 School")
fviz_screeplot(matched2_neighborhood_pca, addlabels =TRUE) +
  ggtitle("Matched Group 2 Neighborhood")
fviz_screeplot(matched2_sociodemo_pca, addlabels =TRUE) +
  ggtitle("Matched Group 2 Sociodemo")
fviz_screeplot(matched2_neurocog_pca, addlabels =TRUE) +
  ggtitle("Matched Group 2 Neurocog")

## Regress out age and sex from PC participant loadings ##

## first release resids ## 

X <- as.matrix(first_release_data[c("Age",
                                    "Male")])

Y <- as.matrix(first_release_data[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                                    "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = first_release_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

first_release_data <- first_release_data %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"])

first_release_resid <- first_release_data %>%
  dplyr::select(subject_id, abcd_site, culture_PC1_resid, home_PC1_resid, school_PC1_resid, neighborhood_PC1_resid, 
                sociodemo_PC1_resid, Neurocog_BPPC1_GenAbility_resid, Neurocog_BPPC2_ExecFunc_resid, Neurocog_BPPC3_LearnMem_resid)

write.csv(first_release_resid, 'ABCD_release1_PCresid.csv', row.names = FALSE)

## second release resids ## 

X <- as.matrix(second_release_data[c("Age",
                                     "Male")])

Y <- as.matrix(second_release_data[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                                     "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])


lm.ControlDem <- lm(Y ~ X, data = second_release_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

second_release_data <- second_release_data %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"])

second_release_resid <- second_release_data %>%
  dplyr::select(subject_id, abcd_site, culture_PC1_resid, home_PC1_resid, school_PC1_resid, neighborhood_PC1_resid, 
                sociodemo_PC1_resid, Neurocog_BPPC1_GenAbility_resid, Neurocog_BPPC2_ExecFunc_resid, Neurocog_BPPC3_LearnMem_resid)

write.csv(second_release_resid, 'ABCD_release2_PCresid.csv', row.names = FALSE)

## Match group 1 ## 

X <- as.matrix(matched1_data[c("Age",
                               "Male")])

Y <- as.matrix(matched1_data[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                               "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = matched1_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

matched1_data <- matched1_data %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"])

## Match group 2 ##

X <- as.matrix(matched2_data[c("Age",
                               "Male")])

Y <- as.matrix(matched2_data[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                               "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = matched2_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

matched2_data <- matched2_data %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"])

## Match group 3 ## 

X <- as.matrix(matched3_data[c("Age",
                               "Male")])

Y <- as.matrix(matched3_data[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                               "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = matched3_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

matched3_data <- matched3_data %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"])

## Regress out age and sex from enviro/sociodem vars first and compare participant loadings with the regressed PC loadings above ##

## first release ##

first_release_varresid_data <- first_release_data %>%
  dplyr::select(c(Age, Male, CultureEnvVar_list$V1, HomeEnvVar_list$V1, SchoolEnvVar_list$V1, NeighborhoodEnvVar_list$V1, SociodemoVar_list$V1))

X <- as.matrix(first_release_varresid_data[c("Age", "Male")])

Y <- as.matrix(first_release_varresid_data[c(CultureEnvVar_list$V1, HomeEnvVar_list$V1, SchoolEnvVar_list$V1, NeighborhoodEnvVar_list$V1, SociodemoVar_list$V1)])

lm.ControlDem <- lm(Y ~ X, data = first_release_varresid_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

rel1_culturevarres_pca     <- PCA(CombMeasures_allcog_resid[,CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_homevarres_pca        <- PCA(CombMeasures_allcog_resid[,HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_schoolvarres_pca      <- PCA(CombMeasures_allcog_resid[,SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_neighborhoodvarres_pca <- PCA(CombMeasures_allcog_resid[,NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_sociodemovarres_pca    <- PCA(CombMeasures_allcog_resid[,SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)

## second release ## 

second_release_varresid_data <- second_release_data %>%
  dplyr::select(c(Age, Male, CultureEnvVar_list$V1, HomeEnvVar_list$V1, SchoolEnvVar_list$V1, NeighborhoodEnvVar_list$V1, SociodemoVar_list$V1))

X <- as.matrix(second_release_varresid_data[c("Age", "Male")])

Y <- as.matrix(second_release_varresid_data[c(CultureEnvVar_list$V1, HomeEnvVar_list$V1, SchoolEnvVar_list$V1, NeighborhoodEnvVar_list$V1, SociodemoVar_list$V1)])

lm.ControlDem <- lm(Y ~ X, data = second_release_varresid_data)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

rel2_culturevarres_pca     <- PCA(CombMeasures_allcog_resid[,CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_homevarres_pca        <- PCA(CombMeasures_allcog_resid[,HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_schoolvarres_pca      <- PCA(CombMeasures_allcog_resid[,SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_neighborhoodvarres_pca <- PCA(CombMeasures_allcog_resid[,NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_sociodemovarres_pca    <- PCA(CombMeasures_allcog_resid[,SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)

## Compare participant loadings from the two PC sets with regressed age and sex ## 

## release 1 ## 

(rel1_home_rescor <- cor.test(first_release_data$home_PC1_resid, rel1_homevarres_pca$ind$coord[,1], method = "spearman"))
(rel1_school_rescor <- cor.test(first_release_data$school_PC1_resid, rel1_schoolvarres_pca$ind$coord[,1], method = "spearman"))
(rel1_neighborhood_rescor <- cor.test(first_release_data$neighborhood_PC1_resid, rel1_neighborhoodvarres_pca$ind$coord[,1], method = "spearman"))
(rel1_culture_rescor <- cor.test(first_release_data$culture_PC1_resid, rel1_culturevarres_pca$ind$coord[,1], method = "spearman"))
(rel1_sociodemo_rescor <- cor.test(first_release_data$sociodemo_PC1_resid, rel1_sociodemovarres_pca$ind$coord[,1], method = "spearman"))

## release 2 ##

(rel2_home_rescor <- cor.test(second_release_data$home_PC1_resid, rel2_homevarres_pca$ind$coord[,1], method = "spearman"))
(rel2_school_rescor <- cor.test(second_release_data$school_PC1_resid, rel2_schoolvarres_pca$ind$coord[,1], method = "spearman"))
(rel2_neighborhood_rescor <- cor.test(second_release_data$neighborhood_PC1_resid, rel2_neighborhoodvarres_pca$ind$coord[,1], method = "spearman"))
(rel2_culture_rescor <- cor.test(second_release_data$culture_PC1_resid, rel2_culturevarres_pca$ind$coord[,1], method = "spearman"))
(rel2_sociodemo_rescor <- cor.test(second_release_data$sociodemo_PC1_resid, rel2_sociodemovarres_pca$ind$coord[,1], method = "spearman"))

## Plot variable contributions for environment and sociodemographic categories ## 

## set colors for figs and plot theme ## 

color_sociodem <- "#EBF094"
color_home <- "#B1D8B7"
color_school <-  "#7CC644"
color_culture <- "#116530"
color_neighborhood <- "#16A637"
color_genabil <- "#CB6CE6"
color_execfunc <- "#FF66C4"
color_learnmem <- "#8C52FF"
color_cog <- "#77AAAD"

pc_theme <- theme_minimal() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(), 
                                    axis.title.x = element_blank(), 
                                    axis.title.y = element_blank(),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_blank())

pc_theme_toptitles <- theme_minimal() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                              panel.grid.minor = element_blank(), 
                                              axis.title.x = element_blank(), 
                                              axis.title.y = element_blank(),
                                              axis.text.x = element_text(size = 14),
                                              axis.text.y = element_blank())

pc_theme2 <- theme_minimal() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(), 
                                     axis.title.x = element_text(size = 14), 
                                     axis.title.y = element_blank(),
                                     axis.text.x = element_text(size = 14),
                                     axis.text.y = element_blank())

## Culture ##

## rel 1.1 ##

rel1_culture_pccontrib <- data.table(Contribution = rel1_culture_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_culture_pca$var$contrib))

rel1_culture_pcdirect <- data.table(Coord = rel1_culture_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_culture_pca$var$coord))

rel1_negative_dir_culture <- which(rel1_culture_pcdirect$Coord < 0)
rel1_culture_pccontrib$Contribution[rel1_negative_dir_culture] <- -1 * rel1_culture_pccontrib$Contribution[rel1_negative_dir_culture]

(rel1_culture_pccontrib_plot <- rel1_culture_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_culture) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("rel1_culture_pccontrib_plot.png")

## rel 2.0.1 ##

rel2_culture_pccontrib <- data.table(Contribution = rel2_culture_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_culture_pca$var$contrib),
         rel1_order = rel1_culture_pccontrib$Contribution)

rel2_culture_pcdirect <- data.table(Coord =  rel2_culture_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_culture_pca$var$coord))

rel2_negative_dir_culture <- which(rel2_culture_pcdirect$Coord < 0)
rel2_culture_pccontrib$Contribution[rel2_negative_dir_culture] <- -1 * rel2_culture_pccontrib$Contribution[rel2_negative_dir_culture]

(rel2_culture_pccontrib_plot <- rel2_culture_pccontrib %>%
    ggplot(aes(x = reorder(Variable, rel1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_culture) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("rel2_culture_pccontrib_plot.png")

## Home ## 

## rel 1.1 ## 

rel1_home_pccontrib <- data.table(Contribution = rel1_home_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_home_pca$var$contrib))

rel1_home_pcdirect <- data.table(Coord = rel1_home_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_home_pca$var$coord))

rel1_negative_dir_home <- which(rel1_home_pcdirect$Coord < 0)
rel1_home_pccontrib$Contribution[rel1_negative_dir_home] <- -1 * rel1_home_pccontrib$Contribution[rel1_negative_dir_home]

(rel1_home_pccontrib_plot <- rel1_home_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_home) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("rel1_home_pccontrib_plot.png")

## rel 2.0.1 ##

rel2_home_pccontrib <- data.table(Contribution = rel2_home_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_home_pca$var$contrib),
         rel1_order = rel1_home_pccontrib$Contribution)

rel2_home_pcdirect <- data.table(Coord = rel2_home_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_home_pca$var$coord))

rel2_negative_dir_home <- which(rel2_home_pcdirect$Coord < 0)
rel2_home_pccontrib$Contribution[rel2_negative_dir_home] <- -1 * rel2_home_pccontrib$Contribution[rel2_negative_dir_home]

(rel2_home_pccontrib_plot <- rel2_home_pccontrib %>%
    ggplot(aes(x = reorder(Variable, rel1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_home) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("rel2_home_pccontrib_plot.png")

## school ## 

## rel 1.1 ## 

rel1_school_pccontrib <- data.table(Contribution = rel1_school_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_school_pca$var$contrib))

rel1_school_pcdirect <- data.table(Coord = rel1_school_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_school_pca$var$coord))

rel1_negative_dir_school <- which(rel1_school_pcdirect$Coord < 0)
rel1_school_pccontrib$Contribution[rel1_negative_dir_school] <- -1 * rel1_school_pccontrib$Contribution[rel1_negative_dir_school]

(rel1_school_pccontrib_plot <- rel1_school_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_school) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("rel1_school_pccontrib_plot.png")

## rel 2.0.1 ## 

rel2_school_pccontrib <- data.table(Contribution = rel2_school_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_school_pca$var$contrib),
         rel1_order = rel1_school_pccontrib$Contribution)

rel2_school_pcdirect <- data.table(Coord = rel2_school_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_school_pca$var$coord))

rel2_negative_dir_school <- which(rel2_school_pcdirect$Coord < 0)
rel2_school_pccontrib$Contribution[rel2_negative_dir_school] <- -1 * rel2_school_pccontrib$Contribution[rel2_negative_dir_school]

(rel2_school_pccontrib_plot <- rel2_school_pccontrib %>%
    ggplot(aes(x = reorder(Variable, rel1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_school) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("rel2_school_pccontrib_plot.png")

## neighborhood ## 

## rel 1.1 ## 

rel1_neighborhood_pccontrib <- data.table(Contribution = rel1_neighborhood_pca$var$contrib[,1])%>%
  mutate(Variable = rownames(rel1_neighborhood_pca$var$contrib))

rel1_neighborhood_pcdirect <- data.table(Coord = rel1_neighborhood_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_neighborhood_pca$var$coord))

rel1_negative_dir_neighborhood <- which(rel1_neighborhood_pcdirect$Coord < 0)
rel1_neighborhood_pccontrib$Contribution[rel1_negative_dir_neighborhood] <- -1 * rel1_neighborhood_pccontrib$Contribution[rel1_negative_dir_neighborhood]

(rel1_neighborhood_pccontrib_plot <- rel1_neighborhood_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_neighborhood) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("rel1_neighborhood_pccontrib_plot.png")

## rel 2.0.1 ## 

rel2_neighborhood_pccontrib <- data.table(Contribution = rel2_neighborhood_pca$var$contrib[,1])%>%
  mutate(Variable = rownames(rel2_neighborhood_pca$var$contrib),
         rel1_order = rel1_neighborhood_pccontrib$Contribution)

rel2_neighborhood_pcdirect <- data.table(Coord =  rel2_neighborhood_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_neighborhood_pca$var$coord))

rel2_negative_dir_neighborhood <- which(rel2_neighborhood_pcdirect$Coord < 0)
rel2_neighborhood_pccontrib$Contribution[rel2_negative_dir_neighborhood] <- -1 * rel2_neighborhood_pccontrib$Contribution[rel2_negative_dir_neighborhood]

(rel2_neighborhood_pccontrib_plot <- rel2_neighborhood_pccontrib %>%
    ggplot(aes(x = reorder(Variable, rel1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_neighborhood) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("rel2_neighborhood_pccontrib_plot.png")

## sociodemographics ## 

## rel 1.1 ## 

rel1_sociodemo_pccontrib <- data.table(Contribution = rel1_sociodemo_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_sociodemo_pca$var$contrib))

rel1_sociodemo_pcdirect <- data.table(Coord = rel1_sociodemo_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_sociodemo_pca$var$coord))

rel1_negative_dir_sociodemo <- which(rel1_sociodemo_pcdirect$Coord < 0)
rel1_sociodemo_pccontrib$Contribution[rel1_negative_dir_sociodemo] <- -1 * rel1_sociodemo_pccontrib$Contribution[rel1_negative_dir_sociodemo]

(rel1_sociodemo_pccontrib_plot <- rel1_sociodemo_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_sociodem) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("rel1_sociodemo_pccontrib_plot.png")

## rel 2.0.1 ## 

rel2_sociodemo_pccontrib <- data.table(Contribution = rel2_sociodemo_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_sociodemo_pca$var$contrib),
         rel1_order = rel1_sociodemo_pccontrib$Contribution)

rel2_sociodemo_pcdirect <- data.table(Coord = rel2_sociodemo_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_sociodemo_pca$var$coord))

rel2_negative_dir_sociodemo <- which(rel2_sociodemo_pcdirect$Coord < 0)
rel2_sociodemo_pccontrib$Contribution[rel2_negative_dir_sociodemo] <- -1 * rel2_sociodemo_pccontrib$Contribution[rel2_negative_dir_sociodemo]

(rel2_sociodemo_pccontrib_plot <- rel2_sociodemo_pccontrib %>%
    ggplot(aes(x = reorder(Variable, rel1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_sociodem) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("rel2_sociodemo_pccontrib_plot.png")

## Combine plots for Release 1 and 2 PCs ## 

(relplot_pcfig <- grid.arrange(rel1_sociodemo_pccontrib_plot,
                               rel2_sociodemo_pccontrib_plot,
                               rel1_home_pccontrib_plot,
                               rel2_home_pccontrib_plot,
                               rel1_school_pccontrib_plot,
                               rel2_school_pccontrib_plot,
                               rel1_neighborhood_pccontrib_plot,
                               rel2_neighborhood_pccontrib_plot,
                               rel1_culture_pccontrib_plot,
                               rel2_culture_pccontrib_plot,
                               heights = c(1, 1, 1, 3, 1),
                               nrow = 5))
ggsave('relplot_pcfig3.png', relplot_pcfig, height = 17, width = 14)

## Plot PC contributions for matched groups ## 

## Culture ## 

## Matched group 1 ## 

matched1_culture_pccontrib <- data.table(Contribution = matched1_culture_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched1_culture_pca$var$contrib))

matched1_culture_pcdirect <- data.table(Coord = matched1_culture_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched1_culture_pca$var$coord))

matched1_negative_dir_culture <- which(matched1_culture_pcdirect$Coord < 0)
matched1_culture_pccontrib$Contribution[matched1_negative_dir_culture] <- -1 * matched1_culture_pccontrib$Contribution[matched1_negative_dir_culture]

(matched1_culture_pccontrib_plot <- matched1_culture_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_culture) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("matched1_culture_pccontrib_plot.png")

## Matched group 2 ## 

matched2_culture_pccontrib <- data.table(Contribution = matched2_culture_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched2_culture_pca$var$contrib),
         matched1_order = matched1_culture_pccontrib$Contribution)

matched2_culture_pcdirect <- data.table(Coord =  matched2_culture_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched2_culture_pca$var$coord))

matched2_negative_dir_culture <- which(matched2_culture_pcdirect$Coord < 0)
matched2_culture_pccontrib$Contribution[matched2_negative_dir_culture] <- -1 * matched2_culture_pccontrib$Contribution[matched2_negative_dir_culture]

(matched2_culture_pccontrib_plot <- matched2_culture_pccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_culture) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("matched2_culture_pccontrib_plot.png")

## Home ## 

## Matched group 1 ## 

matched1_home_pccontrib <- data.table(Contribution = matched1_home_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched1_home_pca$var$contrib))

matched1_home_pcdirect <- data.table(Coord = matched1_home_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched1_home_pca$var$coord))

matched1_negative_dir_home <- which(matched1_home_pcdirect$Coord < 0)
matched1_home_pccontrib$Contribution[matched1_negative_dir_home] <- -1 * matched1_home_pccontrib$Contribution[matched1_negative_dir_home]

(matched1_home_pccontrib_plot <- matched1_home_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_home) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("matched1_home_pccontrib_plot.png")

## matched group 2 ## 

matched2_home_pccontrib <- data.table(Contribution = matched2_home_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched2_home_pca$var$contrib),
         matched1_order = matched1_home_pccontrib$Contribution)

matched2_home_pcdirect <- data.table(Coord = matched2_home_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched2_home_pca$var$coord))

matched2_negative_dir_home <- which(matched2_home_pcdirect$Coord < 0)
matched2_home_pccontrib$Contribution[matched2_negative_dir_home] <- -1 * matched2_home_pccontrib$Contribution[matched2_negative_dir_home]

(matched2_home_pccontrib_plot <- matched2_home_pccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_home) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("matched2_home_pccontrib_plot.png")

## school ## 

## matched group 1 ## 

matched1_school_pccontrib <- data.table(Contribution = matched1_school_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched1_school_pca$var$contrib))

matched1_school_pcdirect <- data.table(Coord = matched1_school_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched1_school_pca$var$coord))

matched1_negative_dir_school <- which(matched1_school_pcdirect$Coord < 0)
matched1_school_pccontrib$Contribution[matched1_negative_dir_school] <- -1 * matched1_school_pccontrib$Contribution[matched1_negative_dir_school]

(matched1_school_pccontrib_plot <- matched1_school_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_school) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("matched1_school_pccontrib_plot.png")

## matched group 2 ## 

matched2_school_pccontrib <- data.table(Contribution = matched2_school_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched2_school_pca$var$contrib),
         matched1_order = matched1_school_pccontrib$Contribution)

matched2_school_pcdirect <- data.table(Coord = matched2_school_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched2_school_pca$var$coord))

matched2_negative_dir_school <- which(matched2_school_pcdirect$Coord < 0)
matched2_school_pccontrib$Contribution[matched2_negative_dir_school] <- -1 * matched2_school_pccontrib$Contribution[matched2_negative_dir_school]

(matched2_school_pccontrib_plot <- matched2_school_pccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_school) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("matched2_school_pccontrib_plot.png")

## neighborhood ## 

## matched group 1 ## 

matched1_neighborhood_pccontrib <- data.table(Contribution = matched1_neighborhood_pca$var$contrib[,1])%>%
  mutate(Variable = rownames(matched1_neighborhood_pca$var$contrib))

matched1_neighborhood_pcdirect <- data.table(Coord = matched1_neighborhood_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched1_neighborhood_pca$var$coord))

matched1_negative_dir_neighborhood <- which(matched1_neighborhood_pcdirect$Coord < 0)
matched1_neighborhood_pccontrib$Contribution[matched1_negative_dir_neighborhood] <- -1 * matched1_neighborhood_pccontrib$Contribution[matched1_negative_dir_neighborhood]

(matched1_neighborhood_pccontrib_plot <- matched1_neighborhood_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_neighborhood) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("matched1_neighborhood_pccontrib_plot.png")

## matched group 2 ## 

matched2_neighborhood_pccontrib <- data.table(Contribution = matched2_neighborhood_pca$var$contrib[,1])%>%
  mutate(Variable = rownames(matched2_neighborhood_pca$var$contrib),
         matched1_order = matched1_neighborhood_pccontrib$Contribution)

matched2_neighborhood_pcdirect <- data.table(Coord =  matched2_neighborhood_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched2_neighborhood_pca$var$coord))

matched2_negative_dir_neighborhood <- which(matched2_neighborhood_pcdirect$Coord < 0)
matched2_neighborhood_pccontrib$Contribution[matched2_negative_dir_neighborhood] <- -1 * matched2_neighborhood_pccontrib$Contribution[matched2_negative_dir_neighborhood]

(matched2_neighborhood_pccontrib_plot <- matched2_neighborhood_pccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_neighborhood) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("matched2_neighborhood_pccontrib_plot.png")

## sociodemographics ## 

## matched group 1 ## 

matched1_sociodemo_pccontrib <- data.table(Contribution = matched1_sociodemo_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched1_sociodemo_pca$var$contrib))

matched1_sociodemo_pcdirect <- data.table(Coord = matched1_sociodemo_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched1_sociodemo_pca$var$coord))

matched1_negative_dir_sociodemo <- which(matched1_sociodemo_pcdirect$Coord < 0)
matched1_sociodemo_pccontrib$Contribution[matched1_negative_dir_sociodemo] <- -1 * matched1_sociodemo_pccontrib$Contribution[matched1_negative_dir_sociodemo]

(matched1_sociodemo_pccontrib_plot <- matched1_sociodemo_pccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_sociodem) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme)

ggsave("matched1_sociodemo_pccontrib_plot.png")

## matched group 2 ## 

matched2_sociodemo_pccontrib <- data.table(Contribution = matched2_sociodemo_pca$var$contrib[,1]) %>%
  mutate(Variable = rownames(matched2_sociodemo_pca$var$contrib),
         matched1_order = matched1_sociodemo_pccontrib$Contribution)

matched2_sociodemo_pcdirect <- data.table(Coord = matched2_sociodemo_pca$var$coord[,1]) %>%
  mutate(Variable = rownames(matched2_sociodemo_pca$var$coord))

matched2_negative_dir_sociodemo <- which(matched2_sociodemo_pcdirect$Coord < 0)
matched2_sociodemo_pccontrib$Contribution[matched2_negative_dir_sociodemo] <- -1 * matched2_sociodemo_pccontrib$Contribution[matched2_negative_dir_sociodemo]

(matched2_sociodemo_pccontrib_plot <- matched2_sociodemo_pccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_sociodem) +
    scale_y_continuous(limits = c(-35, 45)) +
    coord_flip() +
    pc_theme )

ggsave("matched2_sociodemo_pccontrib_plot.png")

## Combine plots for matched groups ## 

(matchedplot_pcfig <- grid.arrange(matched1_sociodemo_pccontrib_plot,
                                   matched2_sociodemo_pccontrib_plot,
                                   matched1_home_pccontrib_plot,
                                   matched2_home_pccontrib_plot,
                                   matched1_school_pccontrib_plot,
                                   matched2_school_pccontrib_plot,
                                   matched1_neighborhood_pccontrib_plot,
                                   matched2_neighborhood_pccontrib_plot,
                                   matched1_culture_pccontrib_plot,
                                   matched2_culture_pccontrib_plot,
                                   heights = c(1, 1, 1, 3, 1),
                                   nrow = 5))
ggsave('matchedplot_pcfig.png', matchedplot_pcfig, height = 17, width = 17)

## Compare Variable contributions between release split and matched group PCs ## 

## sociodemographics ## 

relmatch_pccontrib <- data.table(sociodemo_rel1 = rel1_sociodemo_pccontrib$Contribution,
                                 sociodemo_rel2 = rel2_sociodemo_pccontrib$Contribution,
                                 sociodemo_matched1 = matched1_sociodemo_pccontrib$Contribution,
                                 sociodemo_matched2 = matched2_sociodemo_pccontrib$Contribution,
                                 home_rel1 = rel1_home_pccontrib$Contribution,
                                 home_rel2 = rel2_home_pccontrib$Contribution,
                                 home_matched1 = matched1_home_pccontrib$Contribution,
                                 home_matched2 = matched2_home_pccontrib$Contribution,
                                 school_rel1 = rel1_school_pccontrib$Contribution,
                                 school_rel2 = rel2_school_pccontrib$Contribution,
                                 school_matched1 = matched1_school_pccontrib$Contribution,
                                 school_matched2 = matched2_school_pccontrib$Contribution,
                                 neighborhood_rel1 = rel1_neighborhood_pccontrib$Contribution,
                                 neighborhood_rel2 = rel2_neighborhood_pccontrib$Contribution,
                                 neighborhood_matched1 = matched1_neighborhood_pccontrib$Contribution,
                                 neighborhood_matched2 = matched2_neighborhood_pccontrib$Contribution,
                                 culture_rel1 = rel1_culture_pccontrib$Contribution,
                                 culture_rel2 = rel2_culture_pccontrib$Contribution,
                                 culture_matched1 = matched1_culture_pccontrib$Contribution,
                                 culture_matched2 = matched2_culture_pccontrib$Contribution)

(relmatch_sociodemo_corr <- cor(relmatch_pccontrib[,1:4], method = "pearson"))
(relmatch_home_corr <- cor(relmatch_pccontrib[,5:8], method = "pearson"))
(relmatch_school_corr <- cor(relmatch_pccontrib[,9:12], method = "pearson"))
(relmatch_neighborhood_corr <- cor(relmatch_pccontrib[,13:16], method = "pearson"))
(relmatch_culture_corr <- cor(relmatch_pccontrib[,17:20], method = "pearson"))

## p-values for correlations ##

(relmatch_sociodemo_corrpval <- cor.mtest(relmatch_pccontrib[,1:4], method = "pearson"))
(relmatch_home_corrpval <- cor.mtest(relmatch_pccontrib[,5:8], method = "pearson"))
(relmatch_school_corrpval <- cor.mtest(relmatch_pccontrib[,9:12], method = "pearson"))
(relmatch_neighborhood_corrpval <- cor.mtest(relmatch_pccontrib[,13:16], method = "pearson"))
(relmatch_culture_corrpval <- cor.mtest(relmatch_pccontrib[,17:20], method = "pearson"))

rownames(relmatch_sociodemo_corr) <- rep("", 4)
colnames(relmatch_sociodemo_corr) <- rep("", 4)
rownames(relmatch_home_corr) <- rep("", 4)
colnames(relmatch_home_corr) <- rep("", 4)
rownames(relmatch_school_corr) <- rep("", 4)
colnames(relmatch_school_corr) <- rep("", 4)
rownames(relmatch_neighborhood_corr) <- rep("", 4)
colnames(relmatch_neighborhood_corr) <- rep("", 4)
rownames(relmatch_culture_corr) <- rep("", 4)
colnames(relmatch_culture_corr) <- rep("", 4)

corr_color <- rev(brewer.pal(4, "RdBu"))

corr_color2 <- colorRampPalette(c(corr_color[1:2], "white", corr_color[3:4]))

png("relmatch_sociodemo_corrplot.png", res = 700, width = 1600, height = 1200)
(relmatch_sociodemo_corrplot <- corrplot(relmatch_sociodemo_corr, method = 'color', col = corr_color2(100), type = "lower",
                                         diag = FALSE, cl.pos = "n"))
dev.off()

png("relmatch_home_corrplot.png", res = 700, width = 1600, height = 1200)
(relmatch_sociodemo_corrplot <- corrplot(relmatch_home_corr, method = 'color', col = corr_color2(100), type = "lower",
                                         diag = FALSE, cl.pos = "n"))
dev.off()

png("relmatch_school_corrplot.png", res = 700, width = 1600, height = 1200)
(relmatch_sociodemo_corrplot <- corrplot(relmatch_school_corr, method = 'color', col = corr_color2(100), type = "lower",
                                         diag = FALSE, cl.pos = "n"))
dev.off()

png("relmatch_neighborhood_corrplot.png", res = 700, width = 1600, height = 1200)
(relmatch_sociodemo_corrplot <- corrplot(relmatch_neighborhood_corr, method = 'color', col = corr_color2(100), type = "lower",
                                         diag = FALSE, cl.pos = "n"))
dev.off()

png("relmatch_culture_corrplot.png", res = 700, width = 1600, height = 1200)
(relmatch_sociodemo_corrplot <- corrplot(relmatch_culture_corr, method = 'color', col = corr_color2(100), type = "lower",
                                         diag = FALSE, cl.pos = "n"))
dev.off()

## Compare PC loadings with and without regressed covariates ##

## release 1 ## 

(home1_covar_corr <- cor.test(first_release_data$home_PC1, first_release_data$home_PC1_resid, method = "spearman"))
(school1_covar_corr <- cor.test(first_release_data$school_PC1, first_release_data$school_PC1_resid, method = "spearman"))
(neighborhood1_covar_corr <- cor.test(first_release_data$neighborhood_PC1, first_release_data$neighborhood_PC1_resid, method = "spearman"))
(culture1_covar_corr <- cor.test(first_release_data$culture_PC1, first_release_data$culture_PC1_resid, method = "spearman"))
(sociodemo1_covar_corr <- cor.test(first_release_data$sociodemo_PC1, first_release_data$sociodemo_PC1_resid, method = "spearman"))
(genabil1_covar_corr <- cor.test(first_release_data$Neurocog_BPPC1_GenAbility, first_release_data$Neurocog_BPPC1_GenAbility_resid, method = "spearman"))
(execfunc1_covar_corr <- cor.test(first_release_data$Neurocog_BPPC2_ExecFunc, first_release_data$Neurocog_BPPC2_ExecFunc_resid, method = "spearman"))
(learnmem1_covar_corr <- cor.test(first_release_data$Neurocog_BPPC3_LearnMem, first_release_data$Neurocog_BPPC3_LearnMem_resid, method = "spearman"))

## release 2 ## 

(home2_covar_corr <- cor.test(second_release_data$home_PC1, second_release_data$home_PC1_resid, method = "spearman"))
(school2_covar_corr <- cor.test(second_release_data$school_PC1, second_release_data$school_PC1_resid, method = "spearman"))
(neighborhood2_covar_corr <- cor.test(second_release_data$neighborhood_PC1, second_release_data$neighborhood_PC1_resid, method = "spearman"))
(culture2_covar_corr <- cor.test(second_release_data$culture_PC1, second_release_data$culture_PC1_resid, method = "spearman"))
(sociodemo2_covar_corr <- cor.test(second_release_data$sociodemo_PC1, second_release_data$sociodemo_PC1_resid, method = "spearman"))
(genabil2_covar_corr <- cor.test(second_release_data$Neurocog_BPPC1_GenAbility, second_release_data$Neurocog_BPPC1_GenAbility_resid, method = "spearman"))
(execfunc2_covar_corr <- cor.test(second_release_data$Neurocog_BPPC2_ExecFunc, second_release_data$Neurocog_BPPC2_ExecFunc_resid, method = "spearman"))
(learnmem2_covar_corr <- cor.test(second_release_data$Neurocog_BPPC3_LearnMem, second_release_data$Neurocog_BPPC3_LearnMem_resid, method = "spearman"))

## Compare BPPCA and PCA results for neurocognitive scores in complete case subsamples ## 

## Cor test for Neurocog PC1 and Neurocog BPPC1 ## 

## Release 1.1 ## 

(neurocog1_rel1_corr <- cor.test(first_release_data$neurocog_PC1, first_release_data$Neurocog_BPPC1_GenAbility, method = "spearman"))

## Release 2.0.1 ## 

(neurocog1_rel2_corr <- cor.test(second_release_data$neurocog_PC1, second_release_data$Neurocog_BPPC1_GenAbility, method = "spearman"))

## Cor test for Neurocog PC2 and Neurocog BPPC2 ## 

## Release 1.1 ## 

(neurocog2_rel1_corr <- cor.test(first_release_data$neurocog_PC2, first_release_data$Neurocog_BPPC2_ExecFunc , method = "spearman"))

## Release 2.0.1 ## 

(neurocog2_rel2_corr <- cor.test(second_release_data$neurocog_PC2, second_release_data$Neurocog_BPPC2_ExecFunc , method = "spearman"))

## Cor test for Neurocog PC3 and Neurocog BPPC3 ## 

## Release 1.1 ## 

(neurocog3_rel1_corr <- cor.test(first_release_data$neurocog_PC3, first_release_data$Neurocog_BPPC3_LearnMem , method = "spearman"))

## Release 2.0.1 ## 

(neurocog3_rel2_corr <- cor.test(second_release_data$neurocog_PC3, second_release_data$Neurocog_BPPC3_LearnMem , method = "spearman"))

## Compare variable loadings across site splits ## 

## culture ## 

site_culture_corr <- diag(cor(site1_culture_varload, site2_culture_varload, method = "pearson"))
mean(site_culture_corr)
sd(site_culture_corr)
max(site_culture_corr)
min(site_culture_corr)

## home ## 

site_home_corr <- diag(cor(site1_home_varload, site2_home_varload, method = "pearson"))
mean(site_home_corr)
sd(site_home_corr)
max(site_home_corr)
min(site_home_corr)

## school ## 

site_school_corr <- diag(cor(site1_school_varload, site2_school_varload, method = "pearson"))
mean(site_school_corr)
sd(site_school_corr)
max(site_school_corr)
min(site_school_corr)

## neighborhood ## 

site_neighborhood_corr <- diag(cor(site1_neighborhood_varload, site2_neighborhood_varload, method = "pearson"))
mean(site_neighborhood_corr)
sd(site_neighborhood_corr)
max(site_neighborhood_corr)
min(site_neighborhood_corr)

## sociodemo ## 

site_sociodemo_corr <- diag(cor(site1_sociodemo_varload, site2_sociodemo_varload, method = "pearson"))
mean(site_sociodemo_corr)
sd(site_sociodemo_corr)
max(site_sociodemo_corr)
min(site_sociodemo_corr)

## distributions of corr coef ## 

(histogram(site_neighborhood_corr))
(histogram(site_school_corr))
(histogram(site_culture_corr))
(histogram(site_home_corr))
(histogram(site_sociodemo_corr))

site_corr_dist <- data.frame(site_neighborhood_corr,site_school_corr,site_culture_corr,site_home_corr,site_sociodemo_corr )

site_corr_distplots <- list()

(site_corr_distplots[[1]] <- site_corr_dist %>%
    ggplot(aes(x = site_sociodemo_corr)) + 
    geom_histogram(bins = 40, fill = color_sociodem, color = "black") +
    theme_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 15, color = "black")))

(site_corr_distplots[[2]] <- site_corr_dist %>%
    ggplot(aes(x = site_home_corr)) + 
    geom_histogram(bins = 40, fill = color_home, color = "black") +
    theme_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 15, color = "black")))

(site_corr_distplots[[3]] <- site_corr_dist %>%
    ggplot(aes(x = site_school_corr)) + 
    geom_histogram(bins = 40, fill = color_school, color = "black") +
    theme_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 15, color = "black")))

(site_corr_distplots[[4]] <- site_corr_dist %>%
    ggplot(aes(x = site_neighborhood_corr)) + 
    geom_histogram(bins = 40, fill = color_neighborhood, color = "black") +
    theme_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 15, color = "black")))

(site_corr_distplots[[5]] <- site_corr_dist %>%
    ggplot(aes(x = site_culture_corr)) + 
    geom_histogram(bins = 40, fill = color_culture, color = "black") +
    theme_blank() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 15, color = "black")))

plot_margin <- theme(plot.margin = unit(c(1,2,2,2), "cm"))
(site_corrdist_fig <- grid.arrange(grobs = lapply(site_corr_distplots, "+", plot_margin),
                                   heights = c(1, 1, 1),
                                   nrow = 3))

ggsave('site_corrdist_fig.png', site_corrdist_fig, height = 17, width = 17)

## Plot correlation matrix with enviro/sociodemo vars and cog BPPCAs ## 

## Release 1.1 and Release 2.0.1  corr plot ## 

rel1var_corr <- cor(first_release_data[, c(6:12, 20:56, 13:19, 66:68)])
rel2var_corr <- cor(second_release_data[, c(6:12, 20:56, 13:19, 66:68)])

relcomb_corr <- matrix(NA, nrow = 54, ncol = 54)
relcomb_corr[lower.tri(relcomb_corr)] <- rel1var_corr[lower.tri(rel1var_corr)]
relcomb_corr[upper.tri(relcomb_corr)] <- rel2var_corr[upper.tri(rel2var_corr)]
diag(relcomb_corr) <-  0

rel1var_corr_pvals <- cor.mtest(first_release_data[, c(6:12, 20:56, 13:19, 66:68)], method = "spearman")
rel2var_corr_pvals <- cor.mtest(second_release_data[,  c(6:12, 20:56, 13:19, 66:68)], method = "spearman")

relcom_corr_pvals <- matrix(NA, nrow = 54, ncol = 54)
relcom_corr_pvals[lower.tri(relcom_corr_pvals)] <- rel1var_corr_pvals$p[lower.tri(rel1var_corr_pvals$p)]
relcom_corr_pvals[upper.tri(relcom_corr_pvals)] <- rel2var_corr_pvals$p[upper.tri(rel2var_corr_pvals$p)]
diag(relcom_corr_pvals) <- 0

rownames(relcomb_corr) <- c(rep("",54))
colnames(relcomb_corr) <- c(rep(" ", 54))

png("relcom_corrplot.png", res = 700, width = 16000, height = 12000)
(variable_corrplot <- corrplot(relcomb_corr, method = 'color', col = corr_color2(100), type = "full", 
                               diag = FALSE, cl.pos = "r", cl.cex = 1.8, tl.cex = 1.8, tl.col = "black", tl.srt = 0))
dev.off()

png("sociodemo_corrs.png", res = 700, width = 16000, height = 12000)
(foo <- corrplot(relcomb_corr[1:7,1:7], method = 'color', col = corr_color2(100), type = "full", addCoef.col = "black",
                 diag = FALSE, cl.pos = "r", cl.cex = 1.8, tl.cex = 1.8, tl.col = "black", tl.srt = 0))
dev.off()

## Project release 1 and release 2 participants through corresponding comparison PC (e.g. rel 1 culture project through rel 2 culture) ## 

## project release 1.1 participants through release 2.0.1 PCs ##

## home ##

rel1_projrel2_home <- predict.PCA(rel2_home_pca, first_release_data[HomeEnvVar_list$V1])

## school ##

rel1_projrel2_school <- predict.PCA(rel2_school_pca, first_release_data[SchoolEnvVar_list$V1])

## neigborhood ##

rel1_projrel2_neighborhood <- predict.PCA(rel2_neighborhood_pca, first_release_data[NeighborhoodEnvVar_list$V1])

## culture ##

rel1_projrel2_culture <- predict.PCA(rel2_culture_pca, first_release_data[CultureEnvVar_list$V1])

## sociodemo ##

rel1_projrel2_sociodemo <- predict.PCA(rel2_sociodemo_pca, first_release_data[SociodemoVar_list$V1])


## project release 2.0.1 participants throgh release 1.1 PCs ##

## home ##

rel2_projrel1_home <- predict.PCA(rel1_home_pca, second_release_data[HomeEnvVar_list$V1])

## school ##

rel2_projrel1_school <- predict.PCA(rel1_school_pca, second_release_data[SchoolEnvVar_list$V1])

## neigborhood ##

rel2_projrel1_neighborhood <- predict.PCA(rel1_neighborhood_pca, second_release_data[NeighborhoodEnvVar_list$V1])

## culture ##

rel2_projrel1_culture <- predict.PCA(rel1_culture_pca, second_release_data[CultureEnvVar_list$V1])

## sociodemo ##

rel2_projrel1_sociodemo <- predict.PCA(rel1_sociodemo_pca, second_release_data[SociodemoVar_list$V1])

## combine projected and original participant weights ##

rel1_projrel2_PC <- data.table(home_PC1 = rel1_projrel2_home$coord[,1],
                               school_PC1 = rel1_projrel2_school$coord[,1],
                               neighborhood_PC1 = -1 * rel1_projrel2_neighborhood$coord[,1],
                               culture_PC1 = -1 * rel1_projrel2_culture$coord[,1],
                               sociodemo_PC1 = rel1_projrel2_sociodemo$coord[,1],
                               Neurocog_BPPC1_GenAbility = first_release_bppca$Neurocog_BPPC1_GenAbility,
                               Neurocog_BPPC2_ExecFunc = first_release_bppca$Neurocog_BPPC2_ExecFunc,
                               Neurocog_BPPC3_LearnMem = first_release_bppca$Neurocog_BPPC3_LearnMem,
                               Age = first_release_data$Age,
                               Male = first_release_data$Male)

rel2_projrel1_PC <- data.table(home_PC1 = rel2_projrel1_home$coord[,1],
                               school_PC1 = rel2_projrel1_school$coord[,1],
                               neighborhood_PC1 = -1 * rel2_projrel1_neighborhood$coord[,1],
                               culture_PC1 = -1 * rel2_projrel1_culture$coord[,1],
                               sociodemo_PC1 = rel2_projrel1_sociodemo$coord[,1],
                               Neurocog_BPPC1_GenAbility = second_release_bppca$Neurocog_BPPC1_GenAbility,
                               Neurocog_BPPC2_ExecFunc = second_release_bppca$Neurocog_BPPC2_ExecFunc, 
                               Neurocog_BPPC3_LearnMem = second_release_bppca$Neurocog_BPPC3_LearnMem,
                               Age = second_release_data$Age,
                               Male = second_release_data$Male)

rel1_PC <- first_release_data %>%
  dplyr::select(sociodemo_PC1_resid,
                home_PC1_resid,
                school_PC1_resid,
                neighborhood_PC1_resid,
                culture_PC1_resid,
                Neurocog_BPPC1_GenAbility_resid,
                Neurocog_BPPC2_ExecFunc_resid,
                Neurocog_BPPC3_LearnMem_resid)

rel2_PC <- second_release_data %>%
  dplyr::select(sociodemo_PC1_resid,
                home_PC1_resid,
                school_PC1_resid,
                neighborhood_PC1_resid,
                culture_PC1_resid,
                Neurocog_BPPC1_GenAbility_resid,
                Neurocog_BPPC2_ExecFunc_resid,
                Neurocog_BPPC3_LearnMem_resid)

## Regress out age and sex ##

## Combined Release 1 PCs ## 

X <- as.matrix(rel1_projrel2_PC[, c("Age", "Male")])

Y <- as.matrix(rel1_projrel2_PC[, c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1",
                                    "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = rel1_projrel2_PC)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

rel1_projrel2_PC <- rel1_projrel2_PC %>%
  mutate(sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"]) 

## Combined Release 2 PCs ##

X <- as.matrix(rel2_projrel1_PC[, c("Age", "Male")])

Y <- as.matrix(rel2_projrel1_PC[, c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1",
                                    "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = rel2_projrel1_PC)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

rel2_projrel1_PC <- rel2_projrel1_PC %>%
  mutate(sociodemo_PC1_resid = CombMeasures_allcog_resid[,"sociodemo_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid[,"neighborhood_PC1"],
         culture_PC1_resid = CombMeasures_allcog_resid[,"culture_PC1"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid[,"Neurocog_BPPC3_LearnMem"]) 

## Correlation matrix for Release 1 and 2 PCs ## 

## Correlate combined PCs across Rel 1 and Rel 2 ## 

rel1_origproj_corr <- cor(rel1_PC, rel1_projrel2_PC[, c(11:18)], method = "spearman") %>%
  diag()

rel2_origproj_corr <- cor(rel2_PC, rel2_projrel1_PC[, c(11:18)], method = "spearman") %>%
  diag()

## Release 1.1 and Release 2.0.1  corr plot ## 

rel1pc_corr <- cor(first_release_data[, c(81, 78:80, 77, 85:87)])
rel2pc_corr <- cor(second_release_data[, c(81, 78:80, 77, 85:87)])

relcombpc_corr <- matrix(NA, nrow = 8, ncol = 8)
relcombpc_corr[lower.tri(relcombpc_corr)] <- rel1pc_corr[lower.tri(rel1pc_corr)]
relcombpc_corr[upper.tri(relcombpc_corr)] <- rel2pc_corr[upper.tri(rel2pc_corr)]
diag(relcombpc_corr) <- rel_projcomb_corr  

rel1pc_corr_pvals <- cor.mtest(first_release_data[, c(81, 78:80, 77, 85:87)], method = "spearman")
rel2pc_corr_pvals <- cor.mtest(second_release_data[,  c(81, 78:80, 77, 85:87)], method = "spearman")

relcompc_corr_pvals <- matrix(NA, nrow = 8, ncol = 8)
relcompc_corr_pvals[lower.tri(relcompc_corr_pvals)] <- rel1pc_corr_pvals$p[lower.tri(rel1pc_corr_pvals$p)]
relcompc_corr_pvals[upper.tri(relcompc_corr_pvals)] <- rel2pc_corr_pvals$p[upper.tri(rel2pc_corr_pvals$p)]
diag(relcompc_corr_pvals) <- 0

rownames(relcombpc_corr) <- c(rep("",8))
colnames(relcombpc_corr) <- c(rep("",8))

png("relcompc_corrplot.png", res = 700, width = 16000, height = 12000)
(pc_corrplot <- corrplot(relcombpc_corr, method = 'color', col = corr_color2(100), type = "full", addCoef.col = "black", 
                         number.cex = 3, cl.pos = "r", cl.cex = 1.8, tl.cex = 1.8, tl.col = "black", tl.srt = 0))
corrRect(c(1, 4, 3), lwd = 5)
dev.off()

## 4a. Fit mixed-effects models with study site and family as random effect ## 

## General ability ## 

genabil_rel1_mixmodel_fam <- lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                    culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                    sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                    (1 | family_id), data = first_release_data,
                                  control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

genabil_rel2_mixmodel_fam <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## Executive functioning ## 

execfunc_rel1_mixmodel_fam <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

execfunc_rel2_mixmodel_fam <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## Learning/Memory ## 

learnmem_rel1_mixmodel_fam <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

learnmem_rel2_mixmodel_fam <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## summary of MEMs with random effect for study site & family ## 

summary(genabil_rel1_mixmodel_fam)
summary(genabil_rel2_mixmodel_fam)
summary(execfunc_rel1_mixmodel_fam)
summary(execfunc_rel2_mixmodel_fam)
summary(learnmem_rel1_mixmodel_fam)
summary(learnmem_rel2_mixmodel_fam)

## Get R-square terms for models ## 

(rsq_genabil_rel1_mixmodel_fam <-  r.squaredGLMM(genabil_rel1_mixmodel_fam))
(rsq_genabil_rel2_mixmodel_fam <-  r.squaredGLMM(genabil_rel2_mixmodel_fam))
(rsq_execfunc_rel1_mixmodel_fam <- r.squaredGLMM(execfunc_rel1_mixmodel_fam))
(rsq_execfunc_rel2_mixmodel_fam <- r.squaredGLMM(execfunc_rel2_mixmodel_fam))
(rsq_learnmem_rel1_mixmodel_fam <- r.squaredGLMM(learnmem_rel1_mixmodel_fam))
(rsq_learnmem_rel2_mixmodel_fam <- r.squaredGLMM(learnmem_rel2_mixmodel_fam))

## Plot neighborhood x sociodemo interaction ## 

genabil_rel1_inter <- allEffects(genabil_rel1_mixmodel_fam)
genabil_rel2_inter <- allEffects(genabil_rel2_mixmodel_fam)

plot(genabil_rel1_inter, multiline = TRUE, confit = TRUE, ci.style = "bars",
     main = "General Ability - Release 1 Interactions")


plot(genabil_rel2_inter, multiline = TRUE, confit = TRUE, ci.style = "bars",
     main = "General Ability - Release 2 Interactions")

(genabil_rel1_interplot <- plot_model(genabil_rel1_mixmodel_fam, type = "pred", 
                                      terms = c("sociodemo_PC1_resid [-6:6]", "neighborhood_PC1_resid")) +
    theme_blank() +
    xlab("Sociodemographic Component Score") +
    ylab("General Cognitive Ability Component Score") +
    ggtitle("Discovery Sample: Predicted Values of General Cognitive Ability") +
    scale_colour_discrete(name = "Neighborhood Env. Component Score", 
                          labels = c("-1 SD (Less Enriched)", "Mean", "+1 SD (More Enriched)")) + 
    theme(legend.position = "none"))

(genabil_rel2_interplot_legend <- plot_model(genabil_rel2_mixmodel_fam, type = "eff", 
                                             terms = c("sociodemo_PC1_resid [-6:6]", "neighborhood_PC1_resid")) +
    theme_blank() +
    xlab("Sociodemographic Component Score") +
    ylab("General Cognitive Ability Component Score") +
    ggtitle("Replication Sample: Predicted Values of General Cognitive Ability") +
    scale_colour_discrete(name = "Neighborhood Env. Component Score",
                          labels = c("-1 SD (Less Enriched)", "Mean", "+1 SD (More Enriched)")) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)))

(genabil_rel2_interplot <- plot_model(genabil_rel2_mixmodel_fam, type = "pred", terms = c("sociodemo_PC1_resid [-6:6]", "neighborhood_PC1_resid")) +
    theme_blank() +
    xlab("Sociodemographic Component Score") +
    ylab("General Cognitive Ability Component Score") +
    ggtitle("Replication Sample: Predicted Values of General Cognitive Ability") +
    scale_colour_discrete(name = "Neighborhood Env. Component Score",
                          labels = c("-1 SD (Less Enriched)", "Mean", "+1 SD (More Enriched)")) + 
    theme(legend.position = "none",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)))

plot_model(genabil_rel1_mixmodel_fam, type = "int", terms = c("sociodemo_PC1_resid:neighborhood_PC1_resid")) 

plot_model(genabil_rel2_mixmodel_fam, type = "int", terms = c("sociodemo_PC1_resid", "neighborhood_PC1_resid")) 

## 4b. Fit mixed-effects models for matched groups with study site and family as random effect ## 

## General ability ## 

genabil_matched1_mixmodel_fam <- lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                        culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                        sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                        (1 | family_id), data = matched1_data)

genabil_matched2_mixmodel_fam <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                         culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                         sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                         (1 | family_id), data = matched2_data)

## Executive functioning ## 

execfunc_matched1_mixmodel_fam <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                         culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                         sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                         (1 | family_id), data = matched1_data)

execfunc_matched2_mixmodel_fam <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                         culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                         sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                         (1 | family_id), data = matched2_data)

## Learning/Memory ## 

learnmem_matched1_mixmodel_fam <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                         culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                         sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                         (1 | family_id), data = matched1_data)

learnmem_matched2_mixmodel_fam <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                         culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                         sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                         (1 | family_id), data = matched2_data)

##  summary of matched group MEMs with random effects for study site and family ## 

summary(genabil_matched1_mixmodel_fam)
summary(genabil_matched2_mixmodel_fam)
summary(execfunc_matched1_mixmodel_fam)
summary(execfunc_matched2_mixmodel_fam)
summary(learnmem_matched1_mixmodel_fam)
summary(learnmem_matched2_mixmodel_fam)

## Get R-square terms for matched group models ## 

(rsq_genabil_matched1_mixmodel_fam <-  r.squaredGLMM(genabil_matched1_mixmodel_fam))
(rsq_genabil_matched2_mixmodel_fam <-  r.squaredGLMM(genabil_matched2_mixmodel_fam))
(rsq_execfunc_matched1_mixmodel_fam <- r.squaredGLMM(execfunc_matched1_mixmodel_fam))
(rsq_execfunc_matched2_mixmodel_fam <- r.squaredGLMM(execfunc_matched2_mixmodel_fam))
(rsq_learnmem_matched1_mixmodel_fam <- r.squaredGLMM(learnmem_matched1_mixmodel_fam))
(rsq_learnmem_matched2_mixmodel_fam <- r.squaredGLMM(learnmem_matched2_mixmodel_fam))

## 5. Supplementary figures ## 

## Plot hisotgrams for environment, sociodemographic, and cognitive variables ## 

race_eth <- abcd_data %>%
  dplyr::select(subjectid, race_ethnicity) %>% 
  mutate(race_ethnicity = factor(race_ethnicity, 
                                 levels = c("White", "Hispanic", "Black", "Other", "Asian"),
                                 labels = c("Non-Hispanic White", "Hispanic", "Non-Hispanic Black", "Other", "Non-Hispanic Asian")))

figure_data <- abcd_data_final %>%
  mutate(release = ifelse(subject_id %in% first_release_subjects[[1]], 1, 2)) %>%
  left_join(race_eth, by = c("subject_id" = "subjectid"))

rel1_figure_data <- figure_data %>%
  filter(release == 1)

rel2_figure_data <- figure_data %>%
  filter(release == 2)

cog_variance <- data.frame(tbpicvocab = var(rel1_figure_data$NIHtb_PicVocab),
                           tbflanker = var(rel1_figure_data$NIHtb_Flanker),
                           tblist = var(rel1_figure_data$NIHtb_List),
                           tbcardsort = var(rel1_figure_data$NIHtb_CardSort), 
                           tbpattern = var(rel1_figure_data$NIHtb_Pattern), 
                           tbpicture = var(rel1_figure_data$NIHtb_Picture),
                           tbreading = var(rel1_figure_data$NIHtb_Reading),
                           ravlt = var(rel1_figure_data$RAVLT),
                           littleman = var(rel1_figure_data$LMT))

cog_variance[2,1] <- var(rel2_figure_data$NIHtb_PicVocab)
cog_variance[2,2] <- var(rel2_figure_data$NIHtb_Flanker)
cog_variance[2,3] <-var(rel2_figure_data$NIHtb_List)
cog_variance[2,4] <- var(rel2_figure_data$NIHtb_CardSort)
cog_variance[2,5] <-var(rel2_figure_data$NIHtb_Pattern)
cog_variance[2,6] <-var(rel2_figure_data$NIHtb_Picture)
cog_variance[2,7] <- var(rel2_figure_data$NIHtb_Reading)
cog_variance[2,8] <-var(rel2_figure_data$RAVLT)
cog_variance[2,9] <-var(rel2_figure_data$LMT)

## Sociodemo histograms ## 

sociodemo_figs <- list()
plot_margin <- theme(plot.margin = unit(c(1,1,1,1), "cm"))

## Race/Ethnicity ## 

(sociodemo_figs[[1]] <- figure_data %>%
    ggplot(., aes(race_ethnicity, group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = .5) +
    theme_classic() +
    scale_fill_manual(values = c("#889004", "#EBF094")) +
    ylim(0, 3000) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

## household income ## 

(sociodemo_figs[[2]] <- figure_data %>%
    ggplot(., aes(factor(Income), group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = .5) +
    theme_classic() +
    scale_fill_manual(values = c("#889004", "#EBF094")) +
    ylim(0, 2000) +
    ggtitle("") + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

## caregiver education ## 

(sociodemo_figs[[3]] <- figure_data %>%
    ggplot(., aes(factor(Education), group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = .5) +
    theme_classic() +
    scale_fill_manual(values = c("#889004", "#EBF094")) +
    ylim(0, 2000) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(sociodemo_histfig <- grid.arrange(grobs = lapply(sociodemo_figs, "+", plot_margin),
                                   heights = c(1, 1, 1),
                                   nrow = 3))
ggsave('sociodemo_hist.png', sociodemo_histfig, height = 17, width = 10)

## Environment histograms ## 

## home ## 

home_figs <- list()
plot_margin <- theme(plot.margin = unit(c(1,2,2,2), "cm"))

(home_figs[[1]] <- figure_data %>%
    ggplot(., aes(factor(FamConflict_P), group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = .5) +
    theme_classic() +
    scale_fill_manual(values = c("#5EAC6B", "#B1D8B7")) +
    ggtitle("")+
    ylim(0, 1150) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(home_figs[[2]] <- figure_data %>%
    ggplot(., aes(FamConflict_Y, group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#5EAC6B", "#B1D8B7")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(home_figs[[3]] <- figure_data %>%
    ggplot(., aes(ParentAcceptance, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = .5) +
    theme_classic() +
    scale_fill_manual(values = c("#5EAC6B", "#B1D8B7")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(home_figs[[4]] <- figure_data %>%
    ggplot(., aes(ParentMonitoring, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = .5) +
    theme_classic() +
    scale_fill_manual(values = c("#5EAC6B", "#B1D8B7")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(home_histfig <- grid.arrange(grobs = lapply(home_figs, "+", plot_margin),
                              heights = c(1, 1),
                              nrow = 2))
ggsave('home_hist.png', home_histfig, height = 17, width = 20)

## school ##
school_figs <- list()
plot_margin <- theme(plot.margin = unit(c(1,1,1,1), "cm"))

(school_figs[[1]] <- figure_data %>%
    ggplot(., aes(SchoolDisengagement, group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#58AD18", "#7CC644")) +
    ggtitle("") +
    ylim(0, 1500) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(school_figs[[2]] <- figure_data %>%
    ggplot(., aes(SchoolInvolvement, group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#58AD18", "#7CC644")) +
    ggtitle("" ) +
    ylim(0, 1050) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(school_figs[[3]] <- figure_data %>%
    ggplot(., aes(SchoolEnvironment, group = factor(release), fill = factor(release))) +
    geom_bar(stat = "count", position = "identity", color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#58AD18", "#7CC644")) +
    ggtitle("") +
    ylim(0, 850) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(school_histfig <- grid.arrange(grobs = lapply(school_figs, "+", plot_margin),
                                heights = c(1, 1, 1),
                                nrow = 3))
ggsave('school_hist.png', school_histfig, height = 17, width = 10)

## neighborhood ##

neighborhood_figs <- list()
plot_margin <- theme(plot.margin = unit(c(1,2,2,2), "cm"))

(neighborhood_figs[[1]] <- figure_data %>%
    ggplot(., aes(ADI_PercFamPoverty, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[2]] <- figure_data %>%
    ggplot(., aes(ADI_PercPovBelow138, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[3]] <- figure_data %>%
    ggplot(., aes(ADI_PercSingle, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[4]] <- figure_data %>%
    ggplot(., aes(ADI_IncDisparityIdx, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5 ) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[5]] <- figure_data %>%
    ggplot(., aes(ADI_MedianFamInc, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[6]] <- figure_data %>%
    ggplot(., aes(ADI_PercHOwner, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[7]] <- figure_data %>%
    ggplot(., aes(ADI_Edu_HSDip, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[8]] <- figure_data %>%
    ggplot(., aes(ADI_Edu_LTHS, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_figs[[9]] <- figure_data %>%
    ggplot(., aes(ADI_PercUnemp, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#12812C", "#16A637")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(neighborhood_histfig <- grid.arrange(grobs = lapply(neighborhood_figs, "+", plot_margin),
                                      heights = c(1, 1, 1),
                                      nrow = 3))
ggsave('neighborhood_hist.png', neighborhood_histfig, height = 17, width = 17)

## culture ##

culture_figs <- list()
plot_margin <- theme(plot.margin = unit(c(1,2,2,2), "cm"))

(culture_figs[[1]] <- figure_data %>%
    ggplot(., aes(EthIdentExploration, group = factor(release), fill = factor(release))) +
    geom_density( color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0D5126", "#116530")) +
    ggtitle("")  +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs[[2]] <- figure_data %>%
    ggplot(., aes(EthIdentCommitment, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0D5126", "#116530")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs[[3]] <- figure_data %>%
    ggplot(., aes(FamObligation, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0D5126", "#116530")) +
    ggtitle("")  +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs[[4]] <- figure_data %>%
    ggplot(., aes(FamReferent, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0D5126", "#116530")) +
    ggtitle("")  +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs[[5]] <- figure_data %>%
    ggplot(., aes(FamSupport, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0D5126", "#116530")) +
    ggtitle("")  +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs[[6]] <- figure_data %>%
    ggplot(., aes(IndepSelfReliance, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values =c("#0D5126", "#116530")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs[[7]] <- figure_data %>%
    ggplot(., aes(Religion, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0D5126", "#116530")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(culture_figs_arrange <- grid.arrange(grobs = lapply(culture_figs, "+", plot_margin), 
                                      heights = c(1, 1, 1, 1),
                                      nrow = 4) )
ggsave('culture_hist.png', culture_figs_arrange, height = 17, width = 10)

## cog histograms ##

cog_figs <- list()
plot_margin <- theme(plot.margin = unit(c(1,2,2,2), "cm"))

(cog_figs[[1]] <- figure_data %>%
    ggplot(., aes(NIHtb_PicVocab, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[2]] <- figure_data %>%
    ggplot(., aes(NIHtb_Flanker, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[7]] <- figure_data %>%
    ggplot(., aes(NIHtb_List, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[4]] <- figure_data %>%
    ggplot(., aes(NIHtb_CardSort, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[5]] <- figure_data %>%
    ggplot(., aes(NIHtb_Pattern, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[6]] <- figure_data %>%
    ggplot(., aes(NIHtb_Picture, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[3]] <- figure_data %>%
    ggplot(., aes(NIHtb_Reading, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[8]] <- figure_data %>%
    ggplot(., aes(RAVLT, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs[[9]] <- figure_data %>%
    ggplot(., aes(LMT, group = factor(release), fill = factor(release))) +
    geom_density(color = "black", alpha = 0.5) +
    theme_classic() +
    scale_fill_manual(values = c("#0884AF", "#0FADE4")) +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(6, "pt"),
          legend.position = "none"))

(cog_figs_arrange <- grid.arrange(grobs = lapply(cog_figs, "+", plot_margin), 
                                  heights = c(1, 1, 1),
                                  nrow = 3) )
ggsave('cog_hist.png', cog_figs_arrange, height = 17, width = 17)

## Sociodemographic and Environment component variable scores ## 

home_pc_varload <- data.table(variable = rownames(rel2_home_pca$var$coord), 
                              rel1 = rel1_home_pca$var$coord[,1],
                              rel2 = rel2_home_pca$var$coord[,1])

school_pc_varload <- data.table(variable = rownames(rel2_school_pca$var$coord), 
                                rel1 = rel1_school_pca$var$coord[,1],
                                rel2 = rel2_school_pca$var$coord[,1])

neighborhood_pc_varload <- data.table(variable = rownames(rel2_neighborhood_pca$var$coord), 
                                      rel1 = rel1_neighborhood_pca$var$coord[,1],
                                      rel2 = rel2_neighborhood_pca$var$coord[,1])

culture_pc_varload <- data.table(variable = rownames(rel2_culture_pca$var$coord), 
                                 rel1 = rel1_culture_pca$var$coord[,1],
                                 rel2 = rel2_culture_pca$var$coord[,1])

sociodemo_pc_varload <- data.table(variable = rownames(rel2_sociodemo_pca$var$coord), 
                                   rel1 = rel1_sociodemo_pca$var$coord[,1],
                                   rel2 = rel2_sociodemo_pca$var$coord[,1])

write.csv(home_pc_varload, "varload_homepc.csv")
write.csv(school_pc_varload, "varload_schoolpc.csv")
write.csv(neighborhood_pc_varload, "varload_neighborhoodpc.csv")
write.csv(culture_pc_varload, "varload_culturepc.csv")
write.csv(sociodemo_pc_varload, "varload_sociodemopc.csv")

## Neighborhood x sociodemographic interaction plot with observed values ##

## create factor for levels of neighborhood PC ## 
## release 1 ##

first_release_data <- first_release_data %>%
  mutate(neighborhood_level = NA)

mean_neighborhood_rel1 <- mean(first_release_data$neighborhood_PC1_resid)
sd_neighborhood_rel1 <- sd(first_release_data$neighborhood_PC1_resid)

first_release_data[first_release_data[,"neighborhood_PC1_resid"] >= mean_neighborhood_rel1 + sd_neighborhood_rel1, "neighborhood_level"] <- "high"
first_release_data[first_release_data[,"neighborhood_PC1_resid"] <= mean_neighborhood_rel1 - sd_neighborhood_rel1, "neighborhood_level"] <- "low"
first_release_data[is.na(first_release_data$neighborhood_level), "neighborhood_level"] <- "mean"

## release 2 ## 

second_release_data <- second_release_data %>%
  mutate(neighborhood_level = NA)

mean_neighborhood_rel2 <- mean(second_release_data$neighborhood_PC1_resid)
sd_neighborhood_rel2 <- sd(second_release_data$neighborhood_PC1_resid)

second_release_data[second_release_data[,"neighborhood_PC1_resid"] >= mean_neighborhood_rel2 + sd_neighborhood_rel2, "neighborhood_level"] <- "high"
second_release_data[second_release_data[,"neighborhood_PC1_resid"] <= mean_neighborhood_rel2 - sd_neighborhood_rel2, "neighborhood_level"] <- "low"
second_release_data[is.na(second_release_data$neighborhood_level), "neighborhood_level"] <- "mean"
factor(second_release_data$neighborhood_level, levels = c("low", "mean", "high"))

(genabil_rel1_inter_obsplot <- ggplot() +
    theme_blank() +
    geom_point(data = first_release_data, mapping = aes(x = sociodemo_PC1_resid, y = Neurocog_BPPC1_GenAbility_resid, color = factor(neighborhood_level)), 
               alpha = 0.5, size = 0.5) +
    xlab("Sociodemographic Component Score") +
    ylab("General Cognitive Ability Component Score") +
    ggtitle("Discovery Sample: Observed Values of General Cognitive Ability") +
    scale_colour_discrete(name = "Neighborhood Env. Component Score", 
                          labels = c("-1 SD (Less Enriched)", "Mean", "+1 SD (More Enriched)")) + 
    scale_color_manual(values = c("blue", "red", "green3")) +
    xlim(-5,3) +
    theme(legend.position = "none"))

(genabil_rel2_inter_obsplot <- ggplot() + 
    theme_blank() +
    geom_point(data = second_release_data, mapping = aes(x = sociodemo_PC1_resid, y = Neurocog_BPPC1_GenAbility_resid, color = factor(neighborhood_level)), 
               alpha = 0.5, size = 0.5) +
    xlab("Sociodemographic Component Score") +
    ylab("General Cognitive Ability Component Score") +
    ggtitle("Replication Sample: Observed Values of General Cognitive Ability") +
    scale_colour_discrete(name = "Neighborhood Env. Component Score", 
                          labels = c("-1 SD (Less Enriched)", "Mean", "+1 SD (More Enriched)")) + 
    scale_color_manual(values = c("blue", "red", "green3")) +
    xlim(-5,3) +
    theme(legend.position = "none"))

(genabil_rel2_inter_obsplot_legend <- ggplot() + 
    theme_blank() +
    geom_point(data = second_release_data, mapping = aes(x = sociodemo_PC1_resid, y = Neurocog_BPPC1_GenAbility_resid, color = factor(neighborhood_level)), 
               alpha = 0.5, size = 0.5) +
    xlab("Sociodemographic Component Score") +
    ylab("General Cognitive Ability Component Score") +
    ggtitle("Replication Sample: Observed Values of General Cognitive Ability") +
    
    scale_color_manual(values = c("blue", "red", "green3")) +
    xlim(-5,3) +
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(override.aes = list(size=4))))

## Redo analyses with log transformed income variable ## 

first_release_data_loginc <- abcd_data_final_loginc %>%
  filter(subject_id %in% first_release_subjects[[1]]) %>%
  standardize()

first_release_bppca_loginc <- abcd_data_final_loginc %>%
  filter(subject_id %in% first_release_subjects[[1]]) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

second_release_data_loginc <- abcd_data_final_loginc %>%
  filter(subject_id %nin% first_release_subjects[[1]]) %>%
  standardize()

second_release_bppca_loginc <- abcd_data_final_loginc %>%
  filter(subject_id %nin% first_release_subjects[[1]]) %>%
  dplyr::select(Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc, Neurocog_BPPC3_LearnMem)

## Release 1 and Release 2 group PCAs ## 

## Release 1 ## 
rel1_culture_pca_loginc     <- PCA(first_release_data_loginc[CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_home_pca_loginc        <- PCA(first_release_data_loginc[HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_school_pca_loginc      <- PCA(first_release_data_loginc[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_neighborhood_pca_loginc <- PCA(first_release_data_loginc[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_sociodemo_pca_loginc    <- PCA(first_release_data_loginc[SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel1_neurocog_pca_loginc    <- PCA(first_release_data_loginc[NeurocogVar_list$V1], scale.unit = TRUE, graph = TRUE)

first_release_data_loginc <- first_release_data_loginc %>%
  mutate(culture_PC1 = -1 * rel1_culture_pca_loginc$ind$coord[,1],
         home_PC1 = rel1_home_pca_loginc$ind$coord[,1],
         school_PC1 = rel1_school_pca_loginc$ind$coord[,1],
         neighborhood_PC1 = -1 * rel1_neighborhood_pca_loginc$ind$coord[,1],
         sociodemo_PC1 = rel1_sociodemo_pca_loginc$ind$coord[,1],
         neurocog_PC1 = rel1_neurocog_pca_loginc$ind$coord[,1],
         neurocog_PC2 = rel1_neurocog_pca_loginc$ind$coord[,2],
         neurocog_PC3 = rel1_neurocog_pca_loginc$ind$coord[,3], 
         Neurocog_BPPC1_GenAbility = first_release_bppca_loginc$Neurocog_BPPC1_GenAbility,
         Neurocog_BPPC2_ExecFunc = first_release_bppca_loginc$Neurocog_BPPC2_ExecFunc,
         Neurocog_BPPC3_LearnMem = first_release_bppca_loginc$Neurocog_BPPC3_LearnMem)

## Release 2 ##

rel2_culture_pca_loginc     <- PCA(second_release_data_loginc[CultureEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_home_pca_loginc        <- PCA(second_release_data_loginc[HomeEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_school_pca_loginc      <- PCA(second_release_data_loginc[SchoolEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_neighborhood_pca_loginc <- PCA(second_release_data_loginc[NeighborhoodEnvVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_sociodemo_pca_loginc    <- PCA(second_release_data_loginc[SociodemoVar_list$V1], scale.unit = TRUE, graph = TRUE)
rel2_neurocog_pca_loginc         <- PCA(second_release_data_loginc[NeurocogVar_list$V1], scale.unit = TRUE, graph = TRUE)

second_release_data_loginc <- second_release_data_loginc %>%
  mutate(culture_PC1 = -1 * rel2_culture_pca_loginc$ind$coord[,1],
         home_PC1 = rel2_home_pca_loginc$ind$coord[,1],
         school_PC1 = rel2_school_pca_loginc$ind$coord[,1],
         neighborhood_PC1 = -1 * rel2_neighborhood_pca_loginc$ind$coord[,1],
         sociodemo_PC1 = rel2_sociodemo_pca_loginc$ind$coord[,1],
         neurocog_PC1 = rel2_neurocog_pca_loginc$ind$coord[,1],
         neurocog_PC2 = rel2_neurocog_pca_loginc$ind$coord[,2],
         neurocog_PC3 = rel2_neurocog_pca_loginc$ind$coord[,3],
         Neurocog_BPPC1_GenAbility = second_release_bppca_loginc$Neurocog_BPPC1_GenAbility,
         Neurocog_BPPC2_ExecFunc = second_release_bppca_loginc$Neurocog_BPPC2_ExecFunc,
         Neurocog_BPPC3_LearnMem = second_release_bppca_loginc$Neurocog_BPPC3_LearnMem)

## first release resids ## 

X <- as.matrix(first_release_data_loginc[c("Age",
                                           "Male")])

Y <- as.matrix(first_release_data_loginc[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                                           "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = first_release_data_loginc)
CombMeasures_allcog_resid_loginc <- lm.ControlDem$residuals

first_release_data_loginc <- first_release_data_loginc %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid_loginc[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid_loginc[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid_loginc[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid_loginc[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid_loginc[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid_loginc[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid_loginc[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid_loginc[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid_loginc[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid_loginc[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid_loginc[,"Neurocog_BPPC3_LearnMem"])

first_release_resid_loginc <- first_release_data_loginc %>%
  dplyr::select(subject_id, abcd_site, culture_PC1_resid, home_PC1_resid, school_PC1_resid, neighborhood_PC1_resid, 
                sociodemo_PC1_resid, Neurocog_BPPC1_GenAbility_resid, Neurocog_BPPC2_ExecFunc_resid, Neurocog_BPPC3_LearnMem_resid)

## second release resids ## 

X <- as.matrix(second_release_data_loginc[c("Age",
                                            "Male")])

Y <- as.matrix(second_release_data_loginc[c("culture_PC1", "home_PC1", "school_PC1", "neighborhood_PC1", "sociodemo_PC1","neurocog_PC1", "neurocog_PC2",
                                            "neurocog_PC3", "Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc", "Neurocog_BPPC3_LearnMem")])

lm.ControlDem <- lm(Y ~ X, data = second_release_data_loginc)
CombMeasures_allcog_resid_loginc <- lm.ControlDem$residuals

second_release_data_loginc <- second_release_data_loginc %>%
  mutate(culture_PC1_resid = CombMeasures_allcog_resid_loginc[,"culture_PC1"],
         home_PC1_resid = CombMeasures_allcog_resid_loginc[,"home_PC1"],
         school_PC1_resid = CombMeasures_allcog_resid_loginc[,"school_PC1"],
         neighborhood_PC1_resid = CombMeasures_allcog_resid_loginc[,"neighborhood_PC1"],
         sociodemo_PC1_resid = CombMeasures_allcog_resid_loginc[,"sociodemo_PC1"],
         neurocog_PC1_resid = CombMeasures_allcog_resid_loginc[,"neurocog_PC1"],
         neurocog_PC2_resid = CombMeasures_allcog_resid_loginc[,"neurocog_PC2"],
         neurocog_PC3_resid = CombMeasures_allcog_resid_loginc[,"neurocog_PC3"],
         Neurocog_BPPC1_GenAbility_resid = CombMeasures_allcog_resid_loginc[,"Neurocog_BPPC1_GenAbility"],
         Neurocog_BPPC2_ExecFunc_resid = CombMeasures_allcog_resid_loginc[,"Neurocog_BPPC2_ExecFunc"],
         Neurocog_BPPC3_LearnMem_resid = CombMeasures_allcog_resid_loginc[,"Neurocog_BPPC3_LearnMem"])

second_release_resid_loginc <- second_release_data_loginc %>%
  dplyr::select(subject_id, abcd_site, culture_PC1_resid, home_PC1_resid, school_PC1_resid, neighborhood_PC1_resid, 
                sociodemo_PC1_resid, Neurocog_BPPC1_GenAbility_resid, Neurocog_BPPC2_ExecFunc_resid, Neurocog_BPPC3_LearnMem_resid)

## release 1 participant loading comparison for standard income var and log transformed ## 
(rel1_sociodemo_logcor <- cor.test(first_release_data$sociodemo_PC1_resid, first_release_data_loginc$sociodemo_PC1_resid, method = "spearman"))

## release 2 participant loading comparison for standard income var and log transformed ##

(rel2_sociodemo_logcor <- cor.test(second_release_data$sociodemo_PC1_resid, second_release_data_loginc$sociodemo_PC1_resid, method = "spearman"))

## Mixed effects models with log transformed income ## 

## General ability ## 

genabil_rel1_mixmodel_loginc <- lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                       (1 | family_id), data = first_release_data_loginc,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

genabil_rel2_mixmodel_loginc <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                        culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                        sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                        (1 | family_id), data = second_release_data_loginc,
                                      control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## Executive functioning ## 

execfunc_rel1_mixmodel_loginc <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                        culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                        sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                        (1 | family_id), data = first_release_data_loginc,
                                      control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

execfunc_rel2_mixmodel_loginc <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                        culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                        sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                        (1 | family_id), data = second_release_data_loginc,
                                      control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## Learning/Memory ## 

learnmem_rel1_mixmodel_loginc <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                        culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                        sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                        (1 | family_id), data = first_release_data_loginc,
                                      control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

learnmem_rel2_mixmodel_loginc <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                        culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                        sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                        (1 | family_id), data = second_release_data_loginc,
                                      control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## log income mixed model summaries ## 

summary(genabil_rel1_mixmodel_loginc)
summary(genabil_rel2_mixmodel_loginc)
summary(execfunc_rel1_mixmodel_loginc)
summary(execfunc_rel2_mixmodel_loginc)
summary(learnmem_rel1_mixmodel_loginc)
summary(learnmem_rel2_mixmodel_loginc)

## r-squares for ln(income) models ## 

(rsq_genabil_rel1_mixmodel_loginc <-  r.squaredGLMM(genabil_rel1_mixmodel_loginc))
(rsq_genabil_rel2_mixmodel_loginc <-  r.squaredGLMM(genabil_rel2_mixmodel_loginc))

(rsq_execfunc_rel1_mixmodel_loginc <- r.squaredGLMM(execfunc_rel1_mixmodel_loginc))
(rsq_execfunc_rel2_mixmodel_loginc <- r.squaredGLMM(execfunc_rel2_mixmodel_loginc))

(rsq_learnmem_rel1_mixmodel_loginc <- r.squaredGLMM(learnmem_rel1_mixmodel_loginc))
(rsq_learnmem_rel2_mixmodel_loginc <- r.squaredGLMM(learnmem_rel2_mixmodel_loginc))

## Comparing participant weights from PCs generated on full variable set & reliability-selective variable sets ## 

## create variable lists with reliability ## 

reliable_homevars <- data.table(c("FamConflict_P", "FamConflict_Y"), "Home")
reliable_schoolvars <- data.table(c("SchoolEnvironment", "SchoolInvolvement"), "School")
reliable_culturevars <- data.table(c("EthIdentExploration", "EthIdentCommitment", "FamObligation", "FamReferent", 
                                     "FamSupport", "Religion"), "Culture")

## release 1.1 reliability-selective PCs ## 

rel1_home_relpca     <- PCA(first_release_data[reliable_homevars$V1], scale.unit = TRUE)
rel1_school_relpca     <- PCA(first_release_data[reliable_schoolvars$V1], scale.unit = TRUE)
rel1_culture_relpca     <- PCA(first_release_data[reliable_culturevars$V1], scale.unit = TRUE)

reliable_data_rel1 <- first_release_data %>%
  mutate(reliable_culture = -1 * rel1_culture_relpca$ind$coord[,1],
         reliable_home = rel1_home_relpca$ind$coord[,1],
         reliable_school = rel1_school_relpca$ind$coord[,1])

## release 2.0.1 reliability-selective PCs ## 

rel2_home_relpca     <- PCA(second_release_data[reliable_homevars$V1], scale.unit = TRUE)
rel2_school_relpca     <- PCA(second_release_data[reliable_schoolvars$V1], scale.unit = TRUE)
rel2_culture_relpca     <- PCA(second_release_data[reliable_culturevars$V1], scale.unit = TRUE)

reliable_data_rel2 <- second_release_data %>%
  mutate(reliable_culture = -1 * rel2_culture_relpca$ind$coord[,1],
         reliable_home = rel2_home_relpca$ind$coord[,1],
         reliable_school = rel2_school_relpca$ind$coord[,1])

## matched group 1 reliability-selective PCs ## 

matched1_home_relpca     <- PCA(matched1_data[reliable_homevars$V1], scale.unit = TRUE)
matched1_school_relpca     <- PCA(matched1_data[reliable_schoolvars$V1], scale.unit = TRUE)
matched1_culture_relpca     <- PCA(matched1_data[reliable_culturevars$V1], scale.unit = TRUE)

reliable_data_matched1 <- matched1_data %>%
  mutate(reliable_culture = -1 * matched1_culture_relpca$ind$coord[,1],
         reliable_home = matched1_home_relpca$ind$coord[,1],
         reliable_school = matched1_school_relpca$ind$coord[,1])

## matched group 2 reliability-selective PCs ## 

matched2_home_relpca     <- PCA(matched2_data[reliable_homevars$V1], scale.unit = TRUE)
matched2_school_relpca     <- PCA(matched2_data[reliable_schoolvars$V1], scale.unit = TRUE)
matched2_culture_relpca     <- PCA(matched2_data[reliable_culturevars$V1], scale.unit = TRUE)

reliable_data_matched2 <- matched2_data %>%
  mutate(reliable_culture = -1 * matched2_culture_relpca$ind$coord[,1],
         reliable_home = matched2_home_relpca$ind$coord[,1],
         reliable_school = matched2_school_relpca$ind$coord[,1])

## Plot Variable contributions for reliability-selective PCs ## 

## Home ## 

## Rel 1 ## 

rel1_home_relpccontrib <- data.table(Contribution = rel1_home_relpca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_home_relpca$var$contrib))

rel1_home_relpcdirect <- data.table(Coord = rel1_home_relpca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_home_relpca$var$coord))

rel1_negative_dir_relhome <- which(rel1_home_relpcdirect$Coord < 0)
rel1_home_relpccontrib$Contribution[rel1_negative_dir_relhome] <- -1 * rel1_home_relpccontrib$Contribution[rel1_negative_dir_relhome]

(rel1_home_relpccontrib_plot <- rel1_home_relpccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_home) +
    scale_y_continuous(limits = c(-10, 51)) +
    coord_flip() +
    pc_theme)

## Rel 2 ## 

rel2_home_relpccontrib <- data.table(Contribution = rel2_home_relpca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_home_relpca$var$contrib),
         matched1_order = rel1_home_relpccontrib$Contribution)

rel2_home_relpcdirect <- data.table(Coord = rel2_home_relpca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_home_relpca$var$coord))

rel2_negative_dir_relhome <- which(rel2_home_relpcdirect$Coord < 0)
rel2_home_relpccontrib$Contribution[rel2_negative_dir_relhome] <- -1 * rel2_home_relpccontrib$Contribution[rel2_negative_dir_relhome]

(rel2_home_relpccontrib_plot <- rel2_home_relpccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_home) +
    scale_y_continuous(limits = c(-10, 51)) +
    coord_flip() +
    pc_theme )

## School ## 

## Rel 1 ## 

rel1_school_relpccontrib <- data.table(Contribution = rel1_school_relpca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_school_relpca$var$contrib))

rel1_school_relpcdirect <- data.table(Coord = rel1_school_relpca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_school_relpca$var$coord))

rel1_negative_dir_relschool <- which(rel1_school_relpcdirect$Coord < 0)
rel1_school_relpccontrib$Contribution[rel1_negative_dir_relschool] <- -1 * rel1_school_relpccontrib$Contribution[rel1_negative_dir_relschool]

(rel1_school_relpccontrib_plot <- rel1_school_relpccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_school) +
    scale_y_continuous(limits = c(-10, 51)) +
    coord_flip() +
    pc_theme)

## Rel 2 ## 

rel2_school_relpccontrib <- data.table(Contribution = rel2_school_relpca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_school_relpca$var$contrib),
         matched1_order = rel1_school_relpccontrib$Contribution)

rel2_school_relpcdirect <- data.table(Coord = rel2_school_relpca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_school_relpca$var$coord))

rel2_negative_dir_relschool <- which(rel2_school_relpcdirect$Coord < 0)
rel2_school_relpccontrib$Contribution[rel2_negative_dir_relschool] <- -1 * rel2_school_relpccontrib$Contribution[rel2_negative_dir_relschool]

(rel2_school_relpccontrib_plot <- rel2_school_relpccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_school) +
    scale_y_continuous(limits = c(-10, 51)) +
    coord_flip() +
    pc_theme)

## Culture ## 

## Rel 1 ## 

rel1_culture_relpccontrib <- data.table(Contribution = rel1_culture_relpca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel1_culture_relpca$var$contrib))

rel1_culture_relpcdirect <- data.table(Coord = rel1_culture_relpca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel1_culture_relpca$var$coord))

rel1_negative_dir_relculture <- which(rel1_culture_relpcdirect$Coord < 0)
rel1_culture_relpccontrib$Contribution[rel1_negative_dir_relculture] <- -1 * rel1_culture_relpccontrib$Contribution[rel1_negative_dir_relculture]

(rel1_culture_relpccontrib_plot <- rel1_culture_relpccontrib %>%
    ggplot(aes(x = reorder(Variable, Contribution), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_culture) +
    scale_y_continuous(limits = c(-10, 51)) +
    coord_flip() +
    pc_theme)

## Rel 2 ## 

rel2_culture_relpccontrib <- data.table(Contribution = rel2_culture_relpca$var$contrib[,1]) %>%
  mutate(Variable = rownames(rel2_culture_relpca$var$contrib),
         matched1_order = rel1_culture_relpccontrib$Contribution)

rel2_culture_relpcdirect <- data.table(Coord = rel2_culture_relpca$var$coord[,1]) %>%
  mutate(Variable = rownames(rel2_culture_relpca$var$coord))

rel2_negative_dir_relculture <- which(rel2_culture_relpcdirect$Coord < 0)
rel2_culture_relpccontrib$Contribution[rel2_negative_dir_relculture] <- -1 * rel2_culture_relpccontrib$Contribution[rel2_negative_dir_relculture]

(rel2_culture_relpccontrib_plot <- rel2_culture_relpccontrib %>%
    ggplot(aes(x = reorder(Variable, matched1_order), y = Contribution)) +
    geom_bar(stat = "identity", position = "dodge", fill = color_culture) +
    scale_y_continuous(limits = c(-10, 51)) +
    coord_flip() +
    pc_theme )

## Combine reliability selective PC plots ## 

(reliable_pcfig <- grid.arrange(rel1_home_relpccontrib_plot,
                                rel2_home_relpccontrib_plot,
                                rel1_school_relpccontrib_plot,
                                rel2_school_relpccontrib_plot,
                                rel1_culture_relpccontrib_plot,
                                rel2_culture_relpccontrib_plot,
                                heights = c(1, 1, 1),
                                nrow = 3))
ggsave('reliable_pcfig.png', reliable_pcfig, height = 10.5, width = 14)

## Residualize age and sex from participant weights ## 

## release 1 ## 

X <- as.matrix(reliable_data_rel1[c("Age",
                                    "Male")])

Y <- as.matrix(reliable_data_rel1[c("reliable_culture", "reliable_home", "reliable_school")])

lm.ControlDem <- lm(Y ~ X, data = reliable_data_rel1)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

reliable_data_rel1 <- reliable_data_rel1 %>%
  mutate(reliable_home_resid = -1 *CombMeasures_allcog_resid[,"reliable_home"],
         reliable_school_resid = CombMeasures_allcog_resid[,"reliable_school"],
         reliable_culture_resid = CombMeasures_allcog_resid[,"reliable_culture"])

## release 2 ## 

X <- as.matrix(reliable_data_rel2[c("Age",
                                    "Male")])

Y <- as.matrix(reliable_data_rel2[c("reliable_culture", "reliable_home", "reliable_school")])

lm.ControlDem <- lm(Y ~ X, data = reliable_data_rel2)
CombMeasures_allcog_resid <- lm.ControlDem$residuals

reliable_data_rel2 <- reliable_data_rel2 %>%
  mutate(reliable_home_resid = -1 *CombMeasures_allcog_resid[,"reliable_home"],
         reliable_school_resid = CombMeasures_allcog_resid[,"reliable_school"],
         reliable_culture_resid = CombMeasures_allcog_resid[,"reliable_culture"])

## Correlate participant weights ## 

(reliable_culturecorr_rel1 <- cor.test(reliable_data_rel1$reliable_culture_resid, first_release_data$culture_PC1_resid, method = "spearman"))
(reliable_culturecorr_rel2 <- cor.test(reliable_data_rel2$reliable_culture_resid, second_release_data$culture_PC1_resid, method = "spearman"))

(reliable_homecorr_rel1 <- cor.test(reliable_data_rel1$reliable_home_resid, first_release_data$home_PC1_resid, method = "spearman"))
(reliable_homecorr_rel2 <- cor.test(reliable_data_rel2$reliable_home_resid, second_release_data$home_PC1_resid, method = "spearman"))

(reliable_schoolcorr_rel1 <- cor.test(reliable_data_rel1$reliable_school_resid, first_release_data$school_PC1_resid, method = "spearman"))
(reliable_schoolcorr_rel2 <- cor.test(reliable_data_rel2$reliable_school_resid, second_release_data$school_PC1_resid, method = "spearman"))

## Test the value-added for each environmental factor in the mixed-effects models ## 

## General ability ## 

genabil_rel1_null <- lmer(Neurocog_BPPC1_GenAbility_resid ~ (1 | abcd_site) + 
                            (1 | family_id), data = first_release_data,
                          control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

genabil_rel2_null <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ (1 | abcd_site) + 
                             (1 | family_id), data = second_release_data,
                           control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## home removed ## 

genabil_rel1_homeremmodel <- lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                    culture_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                    sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                    (1 | family_id), data = first_release_data,
                                  control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

genabil_rel2_homeremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## school removed ## 

genabil_rel1_schoolremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                       (1 | family_id), data = first_release_data,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


genabil_rel2_schoolremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                       (1 | family_id), data = second_release_data,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## neighborhood removed ## 

genabil_rel1_neighremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid +
                                      culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                      sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                      (1 | family_id), data = first_release_data,
                                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

genabil_rel2_neighremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid +
                                      culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                      sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                      (1 | family_id), data = second_release_data,
                                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


## culture removed ## 

genabil_rel1_cultremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

genabil_rel2_cultremmodel <-  lmer(Neurocog_BPPC1_GenAbility_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## Executive functioning ## 

execfunc_rel1_null <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ (1 | abcd_site) +
                             (1 | family_id), data = first_release_data,
                           control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


execfunc_rel2_null <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ (1 | abcd_site) +
                             (1 | family_id), data = second_release_data,
                           control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


## home removed ## 

execfunc_rel1_homeremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


execfunc_rel2_homeremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## school removed ## 

execfunc_rel1_schoolremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                       (1 | family_id), data = first_release_data,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


execfunc_rel2_schoolremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                       (1 | family_id), data = second_release_data,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
## neighborhood removed ## 

execfunc_rel1_neighremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid +
                                      culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                      sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                      (1 | family_id), data = first_release_data,
                                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


execfunc_rel2_neighremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid +
                                      culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                      sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                      (1 | family_id), data = second_release_data,
                                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## culture removed ## 

execfunc_rel1_cultremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


execfunc_rel2_cultremmodel <- lmer(Neurocog_BPPC2_ExecFunc_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## Learning/Memory ## 

learnmem_rel1_null <- lmer(Neurocog_BPPC3_LearnMem_resid ~ (1 | abcd_site) + 
                             (1 | family_id), data = first_release_data,
                           control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


learnmem_rel2_null <- lmer(Neurocog_BPPC3_LearnMem_resid ~ (1 | abcd_site) +
                             (1 | family_id), data = second_release_data,
                           control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## home removed ## 

learnmem_rel1_homeremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


learnmem_rel2_homeremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     culture_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## school removed ## 

learnmem_rel1_schoolremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                       (1 | family_id), data = first_release_data,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


learnmem_rel2_schoolremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + neighborhood_PC1_resid + 
                                       culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid +
                                       sociodemo_PC1_resid:neighborhood_PC1_resid + sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                       (1 | family_id), data = second_release_data,
                                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## neighborhood removed ## 

learnmem_rel1_neighremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + 
                                      culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                      sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) + 
                                      (1 | family_id), data = first_release_data,
                                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


learnmem_rel2_neighremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid +
                                      culture_PC1_resid + sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                      sociodemo_PC1_resid:culture_PC1_resid + (1 | abcd_site) +
                                      (1 | family_id), data = second_release_data,
                                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## culture removed ## 

learnmem_rel1_cultremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + (1 | abcd_site) + 
                                     (1 | family_id), data = first_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))


learnmem_rel2_cultremmodel <- lmer(Neurocog_BPPC3_LearnMem_resid ~ sociodemo_PC1_resid + home_PC1_resid + school_PC1_resid + neighborhood_PC1_resid + 
                                     sociodemo_PC1_resid:home_PC1_resid + sociodemo_PC1_resid:school_PC1_resid +
                                     sociodemo_PC1_resid:neighborhood_PC1_resid + (1 | abcd_site) +
                                     (1 | family_id), data = second_release_data,
                                   control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

## R-squares for null and reduced models ## 

## General ability ## 

## Release 1 ## 

rsq_genabil_rel1_homerem <- r.squaredGLMM(genabil_rel1_homeremmodel)
rsq_genabil_rel1_schoolrem <- r.squaredGLMM(genabil_rel1_schoolremmodel)
rsq_genabil_rel1_neighrem <- r.squaredGLMM(genabil_rel1_neighremmodel)
rsq_genabil_rel1_cultrem <- r.squaredGLMM(genabil_rel1_cultremmodel)

## Release 2 ## 

rsq_genabil_rel2_homerem <- r.squaredGLMM(genabil_rel2_homeremmodel)
rsq_genabil_rel2_schoolrem <- r.squaredGLMM(genabil_rel2_schoolremmodel)
rsq_genabil_rel2_neighrem <- r.squaredGLMM(genabil_rel2_neighremmodel)
rsq_genabil_rel2_cultrem <- r.squaredGLMM(genabil_rel2_cultremmodel)

## Executive functioning ## 

## Release 1 ## 

rsq_execfunc_rel1_homerem <- r.squaredGLMM(execfunc_rel1_homeremmodel)
rsq_execfunc_rel1_schoolrem <- r.squaredGLMM(execfunc_rel1_schoolremmodel)
rsq_execfunc_rel1_neighrem <- r.squaredGLMM(execfunc_rel1_neighremmodel)
rsq_execfunc_rel1_cultrem <- r.squaredGLMM(execfunc_rel1_cultremmodel)

## Release 2 ## 

rsq_execfunc_rel2_homerem <- r.squaredGLMM(execfunc_rel2_homeremmodel)
rsq_execfunc_rel2_schoolrem <- r.squaredGLMM(execfunc_rel2_schoolremmodel)
rsq_execfunc_rel2_neighrem <- r.squaredGLMM(execfunc_rel2_neighremmodel)
rsq_execfunc_rel2_cultrem <- r.squaredGLMM(execfunc_rel2_cultremmodel)

## Learning/Memory ## 

## Release 1 ## 

rsq_learnmem_rel1_homerem <- r.squaredGLMM(learnmem_rel1_homeremmodel)
rsq_learnmem_rel1_schoolrem <- r.squaredGLMM(learnmem_rel1_schoolremmodel)
rsq_learnmem_rel1_neighrem <- r.squaredGLMM(learnmem_rel1_neighremmodel)
rsq_learnmem_rel1_cultrem <- r.squaredGLMM(learnmem_rel1_cultremmodel)

## Release 2 ## 

rsq_learnmem_rel2_homerem <- r.squaredGLMM(learnmem_rel2_homeremmodel)
rsq_learnmem_rel2_schoolrem <- r.squaredGLMM(learnmem_rel2_schoolremmodel)
rsq_learnmem_rel2_neighrem <- r.squaredGLMM(learnmem_rel2_neighremmodel)
rsq_learnmem_rel2_cultrem <- r.squaredGLMM(learnmem_rel2_cultremmodel)

## Calculate the difference in R-square for full and reduced models ## 

## General ability ## 

## Release 1 ## 

rsq_genabil_rel1_mixmodel_fam[1,1] - rsq_genabil_rel1_homerem[1,1]
rsq_genabil_rel1_mixmodel_fam[1,1] - rsq_genabil_rel1_schoolrem[1,1]
rsq_genabil_rel1_mixmodel_fam[1,1] - rsq_genabil_rel1_neighrem[1,1]
rsq_genabil_rel1_mixmodel_fam[1,1] - rsq_genabil_rel1_cultrem[1,1]

## Release 2 ## 

rsq_genabil_rel2_mixmodel_fam[1,1] - rsq_genabil_rel2_homerem[1,1]
rsq_genabil_rel2_mixmodel_fam[1,1] - rsq_genabil_rel2_schoolrem[1,2]
rsq_genabil_rel2_mixmodel_fam[1,1] - rsq_genabil_rel2_neighrem[1,1]
rsq_genabil_rel2_mixmodel_fam[1,1] - rsq_genabil_rel2_cultrem[1,1]

## Executive functioning ## 

## Release 1 ## 

rsq_execfunc_rel1_mixmodel_fam[1,1] - rsq_execfunc_rel1_homerem[1,1]
rsq_execfunc_rel1_mixmodel_fam[1,1] - rsq_execfunc_rel1_schoolrem[1,1]
rsq_execfunc_rel1_mixmodel_fam[1,1] - rsq_execfunc_rel1_neighrem[1,1]
rsq_execfunc_rel1_mixmodel_fam[1,1] - rsq_execfunc_rel1_cultrem[1,1]

## Release 2 ##

rsq_execfunc_rel2_mixmodel_fam[1,1] - rsq_execfunc_rel2_homerem[1,1]
rsq_execfunc_rel2_mixmodel_fam[1,1] - rsq_execfunc_rel2_schoolrem[1,1]
rsq_execfunc_rel2_mixmodel_fam[1,1] - rsq_execfunc_rel2_neighrem[1,1]
rsq_execfunc_rel2_mixmodel_fam[1,1] - rsq_execfunc_rel2_cultrem[1,1]

## Learning/Memory ## 

## Release 1 ## 

rsq_learnmem_rel1_mixmodel_fam[1,1] - rsq_learnmem_rel1_homerem[1,1]
rsq_learnmem_rel1_mixmodel_fam[1,1] - rsq_learnmem_rel1_schoolrem[1,1]
rsq_learnmem_rel1_mixmodel_fam[1,1] - rsq_learnmem_rel1_neighrem[1,1]
rsq_learnmem_rel1_mixmodel_fam[1,1] - rsq_learnmem_rel1_cultrem[1,1]

## Release 2 ## 

rsq_learnmem_rel2_mixmodel_fam[1,1] - rsq_learnmem_rel2_homerem[1,1]
rsq_learnmem_rel2_mixmodel_fam[1,1] - rsq_learnmem_rel2_schoolrem[1,1]
rsq_learnmem_rel2_mixmodel_fam[1,1] - rsq_learnmem_rel2_neighrem[1,1]
rsq_learnmem_rel2_mixmodel_fam[1,1] - rsq_learnmem_rel2_cultrem[1,1]

## Summary reports for reduced models ## 

## General Ability ## 

## release 1 ## 

summary(genabil_rel1_homeremmodel)
summary(genabil_rel1_schoolremmodel)
summary(genabil_rel1_neighremmodel)
summary(genabil_rel1_cultremmodel)
summary(genabil_rel1_sociodemmodel)

## release 2 ## 

summary(genabil_rel2_homeremmodel)
summary(genabil_rel2_schoolremmodel)
summary(genabil_rel2_neighremmodel)
summary(genabil_rel2_cultremmodel)
summary(genabil_rel2_sociodemmodel)

## Executive Function ## 

## Release 1 ## 

summary(execfunc_rel1_homeremmodel)
summary(execfunc_rel1_schoolremmodel)
summary(execfunc_rel1_neighremmodel)
summary(execfunc_rel1_cultremmodel)
summary(execfunc_rel1_sociodemmodel)

## Release 2 ## 

summary(execfunc_rel2_homeremmodel)
summary(execfunc_rel2_schoolremmodel)
summary(execfunc_rel2_neighremmodel)
summary(execfunc_rel2_cultremmodel)
summary(learnmem_rel1_sociodemmodel)

## Learning/Memory ## 

## Release1 ## 

summary(learnmem_rel1_homeremmodel)
summary(learnmem_rel1_schoolremmodel)
summary(learnmem_rel1_neighremmodel)
summary(learnmem_rel1_cultremmodel)
summary(learnmem_rel1_sociodemmodel)

## Release 2 ## 

summary(learnmem_rel2_homeremmodel)
summary(learnmem_rel2_schoolremmodel)
summary(learnmem_rel2_neighremmodel)
summary(learnmem_rel2_cultremmodel)
summary(learnmem_rel2_sociodemmodel)

## Estimate effect  size for release and match group mixed effects models ## 

## Release 1.1/2.0.1 Discovery and Replication Mixed Effects Models ## 

genabil_rel1_effect <- modelTest(genabil_rel1_mixmodel_fam) 
genabil_rel2_effect <- modelTest(genabil_rel2_mixmodel_fam) 
execfunc_rel1_effect <- modelTest(execfunc_rel1_mixmodel_fam) 
execfunc_rel2_effect <- modelTest(execfunc_rel2_mixmodel_fam) 
learnmem_rel1_effect <- modelTest(learnmem_rel1_mixmodel_fam) 
learnmem_rel2_effect <- modelTest(learnmem_rel2_mixmodel_fam) 

## ln(income) moodels for discovery and replication samples ## 

genabil_loginc_rel1_effect <- modelTest(genabil_rel1_mixmodel_loginc) 
genabil_loginc_rel2_effect <- modelTest(genabil_rel2_mixmodel_loginc) 
execfunc_loginc_rel1_effect <- modelTest(execfunc_rel1_mixmodel_loginc) 
execfunc_loginc_rel2_effect <- modelTest(execfunc_rel2_mixmodel_loginc) 
learnmem_loginc_rel1_effect <- modelTest(learnmem_rel1_mixmodel_loginc) 
learnmem_loginc_rel2_effect <- modelTest(learnmem_rel2_mixmodel_loginc) 

## Matched group mixed effects models ## 

genabil_matched1_effect <- modelTest(genabil_matched1_mixmodel_fam) 
genabil_matched2_effect <- modelTest(genabil_matched2_mixmodel_fam) 
execfunc_matched1_effect <- modelTest(execfunc_matched1_mixmodel_fam) 
execfunc_matched2_effect <- modelTest(execfunc_matched2_mixmodel_fam) 
learnmem_matched1_effect <- modelTest(learnmem_matched1_mixmodel_fam) 
learnmem_matched2_effect <- modelTest(learnmem_matched2_mixmodel_fam) 