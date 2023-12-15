library(ukbtools)
library(dplyr)
library(ukbnmr)
library(data.table)

decoded <- fread("/data/sharedData/UK_BIOBANK_DATA/Download_Data/Metabolomics_Data/SecondPhaseData/ukb674198.csv")


nmr <- extract_biomarkers(decoded)
biomarker_qc_flags <- extract_biomarker_qc_flags(decoded)
sample_qc_flags <- extract_sample_qc_flags(decoded)

fwrite(biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv") ##saving QC flags
#####Remove technical variation from the DATA BY USING UKBNMR package####
fwrite(nmr, file="/data/sharedData/UK_BIOBANK_DATA/Download_Data/Metabolomics_Data/SecondPhaseData/nmr_biomarker_data_secondPhase_RemovedTechVar.csv") ##data for further processing
fwrite(biomarker_qc_flags, file="/data/sharedData/UK_BIOBANK_DATA/Download_Data/Metabolomics_Data/SecondPhaseData/nmr_biomarker_qc_flags_secondPhase.csv")
fwrite(sample_qc_flags, file="/data/sharedData/UK_BIOBANK_DATA/Download_Data/Metabolomics_Data/SecondPhaseData/nmr_sample_qc_flags.csv")

MID_metabo = read.csv(file = '/data/sharedData/UK_BIOBANK_DATA/Download_Data/Metabolomics_Data/SecondPhaseData/nmr_biomarker_data_secondPhase_RemovedTechVar.csv')

# use visit index colum in the data file to filter out repeated assesment occurance and keep only first assesment data
MID_metabo <- MID_metabo %>%
  group_by(eid) %>%
  filter(!(visit_index == 1 & any(visit_index == 0)))

# Report duplicated 'eid' values
duplicated_eids <- MID_metabo %>%
  filter(duplicated(eid) | duplicated(eid, fromLast = TRUE)) %>%
  pull(eid)
     if (length(duplicated_eids) > 0) {
  		cat("Duplicated eids are:", duplicated_eids, "\n")
		} else {
  			cat("No duplicates found.\n")
			}
MID_metabo <- MID_metabo %>%
  distinct(eid, .keep_all = TRUE) %>%
  select(-visit_index)  # Remove the visit_index column, as this colum is not necessary anymore and not part of the metabolomics data

save_path <- "/data/sharedData/UK_BIOBANK_DATA/Download_Data/Metabolomics_Data/SecondPhaseData/nmr_biomarker_data_RemovedTechVariation_secondPHASE_FilteredDuplicateV2.csv" ###Metabolomics data for downstream processing
write.csv(MID_metabo, file = save_path, row.names = FALSE)

















