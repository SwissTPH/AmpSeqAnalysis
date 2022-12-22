##########################
# Process AmpSeq results - calcuate MOI and haplotype frequencies across replicates
#
# created 12.12.2022
# monica.golumbeanu@unibas.ch
##########################

library(dplyr)

process_raw_data = function(hap_tab) {
  
  # Ensure the table is a data frame and remove the empty column
  hap_tab = as.data.frame(hap_tab)
  hap_tab$NA. = NULL
  
  # Remove the rows that do not correpsond to a haplotype
  hap_tab = hap_tab %>% dplyr::filter(!(Haplotype %in% c("Noise", "Singelton", "Chimera", "Indels")))
  # Order the rows by sample name
  hap_tab = hap_tab %>% dplyr::arrange(SampleName)
  # Identify the replicate number and calculate MOI for each sample
  hap_tab = hap_tab %>% 
                dplyr::group_by(SampleName) %>% 
                dplyr::mutate(Replicate = match(SampleID, unique(SampleID)),
                              MOI = length(unique(Haplotype)))
  # Calculate the "frequency" of each haplotype for each replicate
  # To do so, for each replicate, we divide the number of reads of each haplotype 
  # to the total number of reads of the replicate
  hap_tab = hap_tab %>% 
                dplyr::group_by(Replicate) %>%
                dplyr::mutate(HapFreqRep = Reads/sum(Reads))
  # We calculate the average frequency for each haplotype across the replicates
  hap_tab = hap_tab %>% 
    dplyr::group_by(SampleName, Haplotype) %>%
    dplyr::mutate(HapFreq = mean(HapFreqRep))
  
  return(hap_tab)
}

process_marker = function(marker_name) {
  
  marker_hap_tab = read.table(paste0("~/GitRepos/ampseq_pipeline/results_patient_samples/finalTabHaplotype_", 
                                     marker_name, ".txt"), 
                              stringsAsFactors = FALSE, header = TRUE)
  
  csp_processed = process_raw_data(marker_hap_tab)
  
  write.table(csp_processed, paste0("~/GitRepos/ampseq_pipeline/results_patient_samples/processedTabHaplotype_", 
                                    marker_name, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# Process data for all markers
# process_marker("csp")
# process_marker("cpp")
# process_marker("cpmp")
process_marker("ama1-D3")
# process_marker("msp7")
