###########################
# Script for identifying microsatellite peaks
#
# INPUT: - microsatellite_peaks.xlsx = excel table with the following columns:
#           Sample ID, Size, Height, Marker, Final peaks (optional)
#        - microsatellite_rules.xlsx = excel table with the marker rules and following columns:
#            Marker, Abs_Min_Height, Perc_Highest, Perc_Minus3, Perc_Plus3, Perc_MinusPlus3, Perc_Plus6
#        - the contents of the Marker columns should match between the two files 
#               (i.e., there should be a rule in microsatellite_rules.xlsx 
#                       for all the markers in the microsatellite_peaks.xlsx)
#
# OUTPUT: - excel table with the following columns:
#           Sample ID, Size, Height, Marker, Final peaks (optional), Detected_peaks, Signal
#
# created 09.10.2022
# monica.golumbeanu@unibas.ch
##########################


# Command to install packages (to be run only once):
# install.packages("readxl") 
# install.packages("dplyr")

# For loading the packages:
library(readxl)
library(dplyr)

# Function which contains the criteria for detecting peaks
# INPUT: sample_marker_data = subset of the raw data which corresponds to the marker of interest
#        calling_rule = rule to be applied for detecting the peaks
# OUTPUT: table with the detection results
call_peaks = function(sample_marker_data, calling_rule) {
  print("Applying microsatellite algorithm ...")
  
  sample_marker_data = as.data.frame(sample_marker_data)
  # Create results column and initialize
  sample_marker_data$Detected_peaks = "NO"
  sample_marker_data$ID = c(1:nrow(sample_marker_data))
  
  # Identify signal from background
  if ((calling_rule$Abs_Min_Height > 0) == TRUE) {
    background_noise = calling_rule$Abs_Min_Height
  } else {
    background_noise = max(sample_marker_data$Height)*calling_rule$Perc_Highest/100
  }
  print(paste("The level of background noise for marker ", calling_rule$Marker, "is: ", background_noise))
  sample_marker_data$Signal = as.integer(sample_marker_data$Height > background_noise)
  
  # Initially all peaks where there is signal are marked as detected
  # They will be dropped as soon as any of the rules does not apply
  sample_marker_data[which(sample_marker_data$Signal == 1), "Detected_peaks"] = "YES"

  # Initialize the entries with results of the detection rules
  sample_marker_data$Exclusion_reason = ""
  sample_marker_data[which(sample_marker_data$Signal == 0), "Exclusion_reason"] = "Background;"
  
  # Round up the sizes of the fragments
  sample_marker_data$Size = round(sample_marker_data$Size)
  
  # Apply each rule to each signal position
  for (i in which(sample_marker_data$Signal == 1)) {
    
    # Extract the fragment size
    fragment_size = sample_marker_data[i, "Size"]
    # Identify the -3, +3 and +6 peaks
    minus_3_peak = sample_marker_data[which(sample_marker_data$Size == (fragment_size - 3)), ]
    plus_3_peak = sample_marker_data[which(sample_marker_data$Size == (fragment_size + 3)), ]
    plus_6_peak = sample_marker_data[which(sample_marker_data$Size == (fragment_size + 6)), ]
    
    # # Just for testing the cutoff rules on individual peaks
    # if (sample_marker_data[i, "Size"] == 171 &&	sample_marker_data[i, "Height"] == 1174) {
    #   print("here")
    # }
      
    # Check for the -3 peak rule
    if(!is.na(calling_rule$Perc_Minus3) && 
       nrow(minus_3_peak) > 0 &&
       any(sample_marker_data[i, "Height"] <= minus_3_peak[, "Height"]*calling_rule$Perc_Minus3/100)) 
    {
      sample_marker_data[i, "Detected_peaks"] = "NO"
      sample_marker_data[i, "Exclusion_reason"] = paste(sample_marker_data[i, "Exclusion_reason"], "Minus3;")
    }
    
    # Check for the +3 peak rule
    if(!is.na(calling_rule$Perc_Plus3) &&
       nrow(plus_3_peak) > 0 &&
       any(sample_marker_data[i, "Height"] <= plus_3_peak[, "Height"]*calling_rule$Perc_Plus3/100))
    {
          sample_marker_data[i, "Detected_peaks"] = "NO"
          sample_marker_data[i, "Exclusion_reason"] = paste(sample_marker_data[i, "Exclusion_reason"], "Plus3;")
    }
    
    # Check for the +6 peak rule
    if(!is.na(calling_rule$Perc_Plus6) && 
       nrow(plus_6_peak) > 0 &&
       any(sample_marker_data[i, "Height"] <= plus_6_peak[, "Height"]*calling_rule$Perc_Plus6/100)) 
    {
          sample_marker_data[i, "Detected_peaks"] = "NO"
          sample_marker_data[i, "Exclusion_reason"] = paste(sample_marker_data[i, "Exclusion_reason"], "Plus6;")
    }
    
    # Check for the +/-3 peak rule
    if(!is.na(calling_rule$Perc_MinusPlus3) && 
       nrow(minus_3_peak) > 0 &&
       nrow(plus_3_peak) > 0 &&
       any(sample_marker_data[i, "Height"] <= (plus_3_peak[, "Height"] + minus_3_peak[, "Height"])*calling_rule$Perc_MinusPlus3/100)) 
    {
          sample_marker_data[i, "Detected_peaks"] = "NO"
          sample_marker_data[i, "Exclusion_reason"] = paste(sample_marker_data[i, "Exclusion_reason"], "PlusMinus3;")
    }
  } 
  
  # Filter the stutter peaks (consecutively detected peaks with fragment size difference <= 2)
  for (i in which(sample_marker_data$Detected_peaks == "YES")) {
    id_same_fragment = which(sample_marker_data$Size >= sample_marker_data[i, "Size"] & 
                               sample_marker_data$Size < sample_marker_data[i, "Size"] + 3 &
                               sample_marker_data$ID > sample_marker_data[i, "ID"])
    same_fragment = sample_marker_data[id_same_fragment, ]
    
    # # Just for testing the stutter peak filtering on individual peaks
    # if (sample_marker_data[i, "Size"] == 168 &&	sample_marker_data[i, "Height"] == 32414) {
    #   print("here")
    # }

    if (nrow(same_fragment) > 0) {
      if(any(same_fragment$Height >= sample_marker_data[i, "Height"] &
            same_fragment$Detected_peaks == "YES")) {
              sample_marker_data[i, "Detected_peaks"] = "NO"
              sample_marker_data[i, "Exclusion_reason"] = paste(sample_marker_data[i, "Exclusion_reason"], "Size") 
      } else {
              sample_marker_data[id_same_fragment, "Detected_peaks"] = "NO"
              sample_marker_data[id_same_fragment, "Exclusion_reason"] = paste(sample_marker_data[id_same_fragment, "Exclusion_reason"], "Size")
      }
    }  
  }
  
  return(sample_marker_data) 
}

# Main function which calls the detection algorithm
# INPUT: raw_peaks = raw data containing the peak details
#        calling_rules = table with the thresholds for calling criteria
#        result_file = path to the file where the result should be saved
# OUTPUT: the results are stored in result_file
run_detection_algorithm = function(raw_peaks, calling_rules, result_file) {
  # Identify the pairs of markers and rules
  sample_marker_pairs = unique(raw_peaks[, c("Sample_ID", "Marker")])
  colnames(sample_marker_pairs) = c("Sample_ID", "Marker")
  
  # Initialize final list
  final_result = NULL
  
  # Select each sample and marker and call peaks 
  for (i in 1:nrow(sample_marker_pairs)) {
    marker_calling_rule = calling_rules %>% filter(Marker == sample_marker_pairs$Marker[i])
    sample_marker_subtable = raw_peaks %>% 
      filter(Sample_ID == sample_marker_pairs$Sample_ID[i] &
               Marker == sample_marker_pairs$Marker[i])
    if (nrow(sample_marker_subtable) > 0 & nrow(marker_calling_rule) > 0) {
      list_peaks = call_peaks(sample_marker_subtable, marker_calling_rule)
      final_result = rbind.data.frame(final_result, list_peaks) 
    }
  }
  
  # Remove ID column
  final_result$ID = NULL
  
  # Write the results to a table
  write.table(final_result, result_file, sep ="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}

# TO MODIFY THE FILE PATHS ACCORDING TO YOUR SYSTEM:
# Read the raw peak data and the microsatellite calling rules
raw_peaks = read_excel("~/genotyping/test_data/microsatellite_peaks_ratios_PfPK2.xlsx")
calling_rules = read_excel("~/genotyping/test_data/microsatellite_rules_abs_11Nov.xlsx")

# Define the file where results will be saved:
result_file = "~/genotyping/microsatellite_analysis/scripts/results_PfPK2.txt"

# Main command that runs the peak detection algorithm:
run_detection_algorithm(raw_peaks, calling_rules, result_file)


