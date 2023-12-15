library(dplyr)
library(stringr)
library(MatchIt)
library(plyr)

################## User choices ###################################

#Set the baseline group
baseline = 'HC'
baseline = 'non-cancer'

# Disease list --> Choose diseases and their corresponding ICD10 code
disease_list <- list(
  'AllCancers' = paste('C', c(1:97), sep = ''),
  'Breast' = c('C50'),
  'Prostate' = c('C61'),
  'Colorectal' = c('C18', 'C20'),
  'Malignant_melanoma' = c('C43'),
  'Lymphoma' = c('C82', 'C83', 'C84', 'C85', 'C86'),
  'Cervical' = c('C53'),
  'Uterine' = c('C54'),
  'Testicular' = c('C62'),
  'Ovarian' = c('C56'),
  'Bladder' = c('C67'),
  'Leukaemia' = c('C91', 'C92', 'C93', 'C94', 'C95'),
  'Kidney' = c('C64'),
  'Thyroid' = c('C73'),
  'Lung' = c('C34')
  
)

# omics could be either genomics, proteomics of proteomics
omics = 'genomics'

# output folder
output_folder = '...'

# input folder
input_folder = '...'



################## Helping functions ###################################
# extract patients with certain ICD code ####
ukb_icd_subset_by_ICD = function (data, icd.code, icd.version = 10) 
{
  ukb_case <- data %>% dplyr::select(matches(paste("^diagnoses.*icd", 
                                                   icd.version, sep = ""))) %>% 
    purrr::map_df(~grepl(icd.code, ., perl = TRUE)) %>% rowSums() > 0
  data_subset <- data[ukb_case,]
  return(data_subset)
}



################## Read data ###################################
## Choose either genomics, proteomics and metabolomics

if(omics == 'genomics'){
  omics_path = paste(input_folder, '/EUROS_PRS_9_disease_send.csv', sep = '')
  output_file = 'all'
}

if(omics == 'proteomics'){
  omics_path = paste(input_folder, '/Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt', sep = '')
  output_file = 'prot' #is used in saving results
}

if(omics == 'metabolomics'){
  omics_path = paste(input_folder, '/nmr_biomarker_data_RemovedTechVariation_secondPHASE_FilteredDuplicateV2.csv', sep = '')
  output_file = 'met'  #is used in saving results
}


if(output_file == 'clinical'){
  PID_omics = my_ukb_data$eid
}


#load all UKBB patients
load(paste(input_folder, '/ukb672643.rda', sep = ''))

if(output_file == 'all'){
  PID_omics = read.table(file = omics_path, sep=';', header = TRUE, row.names = 'eid')
  PID_omics = rownames(PID_omics)
} else if(output_file == 'met') {
  PID_omics = read.table(file = omics_path, sep=',', header = TRUE, fill = TRUE, row.names = 'eid')
  PID_omics = rownames(PID_omics)
} else if(output_file == 'prot'){
  PID_omics = read.table(file = omics_path, sep='\t', header = TRUE, fill = TRUE, row.names = 'PID')
  PID_omics = rownames(PID_omics)
}



#clinical with columns only relevant for finding healthy controls
input_file <- paste(input_folder, '/Decoded_ukb_Matrix_Firoj_ControlHealthyICD9and10codes_only.csv', sep = '')
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# A list of patients that withdrawn from the UKBB
withdrawn_patients = read.table(paste(input_folder, '/withdrawn_patients.csv', sep = ''))
my_ukb_data = my_ukb_data[!my_ukb_data$eid %in% withdrawn_patients$V1,]
data = data[!data$eid %in% withdrawn_patients$V1,]
PID_omics = PID_omics[!PID_omics %in% withdrawn_patients$V1]



################## Preprocess data ###################################
#subset to only those with omics
my_ukb_data<- my_ukb_data[my_ukb_data$eid %in% PID_omics, ]


#find healthy controls 
if(baseline == 'HC'){
  rows_with_all_na <- data[rowSums(is.na(data[, -1])) == (ncol(data) - 1), ]
  Prot_healthy_controls = rows_with_all_na[rows_with_all_na$eid %in% PID_omics,]
  ukb_healthy_controls = my_ukb_data[my_ukb_data$eid %in% Prot_healthy_controls$eid,]
}



#find nonCancer
if(baseline == 'non-cancer'){
  patient_group = read.csv('/home/marsm05/Downloads/UKBB_eid_group.csv', sep = ';', row.names = 'X')
  my_ukb_data_nonCancer = my_ukb_data[my_ukb_data$eid %in% patient_group[patient_group[,'group'] %in% c('OD', 'DF'),]$eid,]
  my_ukb_data_nonCancer = my_ukb_data_nonCancer[!my_ukb_data_nonCancer$eid %in% withdrawn_patients$V1,]
  Prot_healthy_controls = my_ukb_data_nonCancer[my_ukb_data_nonCancer$eid %in% PID_omics,]
  ukb_healthy_controls = my_ukb_data[my_ukb_data$eid %in% Prot_healthy_controls$eid,]
}


################## Choose all patients with the diagnosis ###################################



icd.version = '10'

# Repeat for all diseases
for(j in 1:length(disease_list)){
  disease = names(disease_list[j])
  icd.code.list = disease_list[[j]]
  
  ukb_disease_subset_list = list()
  ukb_case_list = list()
  for(icd.code in icd.code.list){
    
    # all patients with a given icd.code
    patients_with_disease = ukb_icd_subset_by_ICD(my_ukb_data, icd.version = icd.version, icd.code = icd.code)
    
    # the below is to remember which columns corresponded to the disease and is later used for creating ToBeSick group
    ukb_case_subset <- patients_with_disease %>% dplyr::select(matches("40006")) %>%
      purrr::map_df(~grepl(icd.code, ., perl = TRUE))
    
    ukb_case_subset = cbind(eid = patients_with_disease$eid,ukb_case_subset)
    
    ukb_case_list = append(ukb_case_list, list(ukb_case_subset))
    
  }
  
  ukb_case = bind_rows(ukb_case_list)
  
  # There might be duplicates in ukb_case (if for example patient has both M05 and M06 in RA)
  # The below code takes only unique values and for each diagnosis returns TRUE if it belongs to either M05 or M06
  coln = colnames(ukb_case)
  ukb_case = ddply(ukb_case, .(eid), function(x) return(as.logical(colSums(x[,2:dim(ukb_case)[2]]))))
  colnames(ukb_case) = coln
  rownames(ukb_case) = ukb_case$eid
  ukb_case = ukb_case[,-1]
  
  ukb_disease_subset = my_ukb_data[my_ukb_data$eid %in% rownames(ukb_case),]
  
  ################## Divide patients between sick and ToBeSick  ###################################
  
  
  # for each patient, we find all disease related columns, then we map them into corresponding column names for diagnosis time,
  # then remove all columns with secondary disease and select the earliest diagnosis time
  earliest_diagnosis = c()
  for(i in 1:dim(ukb_disease_subset)[1]){
    diagnoses = colnames(ukb_case)[unlist(ukb_case[i,])]
    diagnoses = gsub("type_of_cancer_icd10_f40006", "date_of_cancer_diagnosis_f40005", diagnoses) # This is because of the codes for diagnosis and time of diagnosis is different
    diagnoses_times = ukb_disease_subset[i,diagnoses]
    if(length(diagnoses_times) == 1){
      earliest_diagnosis = c(earliest_diagnosis, as.character(diagnoses_times))
    }
    else{
      earliest_diagnosis = c(earliest_diagnosis, min(as.vector(t(diagnoses_times))))
    }
  }
  
  # Compare sampling time and diagnosis time and divide the patients
  days_diagnosis_after_sample = as.Date(earliest_diagnosis) - ukb_disease_subset[,'date_of_attending_assessment_centre_f53_0_0'] 
  ukb_disease_sick = ukb_disease_subset[days_diagnosis_after_sample < 0,]
  ukb_disease_ToBe_sick = ukb_disease_subset[days_diagnosis_after_sample > 0,]
  
  
  ################## Prepare files for Matching based on age and sex  ###################################
  clinical_variables = c('eid', 'age_at_recruitment_f21022_0_0', 'sex_f31_0_0', 'hand_grip_strength_left_f46_0_0', 'hand_grip_strength_right_f47_0_0', 
                         'smoking_status_f20116_0_0', 'alcohol_drinker_status_f20117_0_0',
                         'diabetes_diagnosed_by_doctor_f2443_0_0', 'diastolic_blood_pressure_automated_reading_f4079_0_0',
                         'systolic_blood_pressure_automated_reading_f4080_0_0', 'body_mass_index_bmi_f21001_0_0', 
                         'haemoglobin_concentration_f30020_0_0',
                         'haematocrit_percentage_f30030_0_0', 'mean_corpuscular_volume_f30040_0_0', 'mean_corpuscular_haemoglobin_f30050_0_0',
                         'red_blood_cell_erythrocyte_distribution_width_f30070_0_0', 'platelet_crit_f30090_0_0', 'mean_platelet_thrombocyte_volume_f30100_0_0',
                         'platelet_distribution_width_f30110_0_0', 'platelet_distribution_width_f30110_0_0', 'lymphocyte_percentage_f30180_0_0', 
                         'monocyte_percentage_f30190_0_0', 'neutrophill_percentage_f30200_0_0', 'basophill_percentage_f30220_0_0', 'reticulocyte_percentage_f30240_0_0',
                         'mean_reticulocyte_volume_f30260_0_0', 'immature_reticulocyte_fraction_f30280_0_0', 'high_light_scatter_reticulocyte_percentage_f30290_0_0',
                         'creatinine_enzymatic_in_urine_f30510_0_0', 'potassium_in_urine_f30520_0_0', 'sodium_in_urine_f30530_0_0', 'alkaline_phosphatase_f30610_0_0', 
                         'alanine_aminotransferase_f30620_0_0', 'apolipoprotein_b_f30640_0_0', 'aspartate_aminotransferase_f30650_0_0', 'urea_f30670_0_0',
                         'cholesterol_f30690_0_0', 'creatinine_f30700_0_0', 'creactive_protein_f30710_0_0', 
                         'gamma_glutamyltransferase_f30730_0_0', 'glycated_haemoglobin_hba1c_f30750_0_0', 'igf1_f30770_0_0', 'ldl_direct_f30780_0_0',
                         'total_bilirubin_f30840_0_0',  'triglycerides_f30870_0_0', 'urate_f30880_0_0')
  
  
  ukb_disease_sick = ukb_disease_sick[clinical_variables]
  ukb_healthy_controls = ukb_healthy_controls[clinical_variables]
  ukb_disease_sick['group'] = rep('Prevalent', dim(ukb_disease_sick)[1])
  ukb_healthy_controls['group'] = rep('HC', dim(ukb_healthy_controls)[1])
  ukb_healthy_controls = na.omit(ukb_healthy_controls) #necessary for Genomics
  
  
  ################## Match healthy to sick  ###################################
  
  #using MatchIt function that pairs sick with healthy based on age and sex
  
  summary = rbind(ukb_disease_sick, ukb_healthy_controls)
  
  a = matchit(group ~ age_at_recruitment_f21022_0_0, data = summary, method = "nearest", exact = ~as.factor(sex_f31_0_0), 
              distance = "euclidean",ratio = 1)
  
  # paired controls
  ukb_healthy_controls_PairedToPrevalent = ukb_healthy_controls[rownames(ukb_healthy_controls) %in% as.vector(a$match.matrix),]
  

  ################## Output files  ###################################
  
  file_path= paste(output_folder, "/ukb_", disease, "_", output_file, "_prevalent.csv", sep = '')
  write.csv(ukb_disease_sick, file = file_path, row.names = FALSE)

  
  file_path= paste(utput_folder, "/ukb_", disease, "_", output_file, "_", baseline ,"_PairedToPrevalent.csv", sep = '')
  write.csv(ukb_healthy_controls_PairedToPrevalent, file = file_path, row.names = FALSE)

}


