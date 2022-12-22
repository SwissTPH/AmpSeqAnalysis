#################################
# Script preparing the sample and barcode input files 
# for the analysis of the patient samples
#
# created 09.12.2022
# monica.golumbeanu@unibas.ch
#################################

library(readxl)

sample_desc = as.data.frame(read_excel("~/genotyping/ampSeq_patients/PatSamples_MiSeq_0001_sample_info_new.xlsx"))
sample_desc = sample_desc[which(!is.na(sample_desc$sample_id)), ]

# Prepare the sample file
final_sample_tab = NULL
final_sample_tab$SampleID = sample_desc$amplicon_id
final_sample_tab$BarcodeID_F = sample_desc$fw_name
final_sample_tab$BarcodeID_R = sample_desc$rv_name
final_sample_tab$SampleName = substr(sample_desc$sample_id, 1, str_locate(sample_desc$sample_id, "_Rep") - 1)
final_sample_tab$Replicate = sample_desc$sample_id
final_sample_tab = as.data.frame(final_sample_tab)

write.table(x = final_sample_tab, file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/input_files/sample_file.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Prepare the barcode files
# Forward
barcodes_f = unique(sample_desc[, c("fw_name", "fw_barcode")])
barcodes = ""
for (i in 1:nrow(barcodes_f)) {
  barcodes = paste(barcodes, paste0(">", barcodes_f[i, 1], "\n", barcodes_f[i, 2]), sep = "\n")
}
barcodes = substr(barcodes, 2, nchar(barcodes))
cat(barcodes, file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/input_files/barcode_F.fasta")

# Reverse
barcodes_r = unique(sample_desc[, c("rv_name", "rv_barcode")])
barcodes = ""
for (i in 1:nrow(barcodes_r)) {
  barcodes = paste(barcodes, paste0(">", barcodes_r[i, 1], "\n", barcodes_r[i, 2]), sep = "\n")
}
barcodes = substr(barcodes, 2, nchar(barcodes))
cat(barcodes, file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/input_files/barcode_R.fasta")


