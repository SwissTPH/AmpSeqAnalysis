###########################
# Main script using the HaplotypR package to detect 
# haplotypes in AmpSeq data
# 
# created 28.11.2022
# monica.golumbeanu@unibas.ch
##########################

# Load the necessary packages
library("HaplotypR")
library("ShortRead")

print(Sys.time())
start_time = Sys.time()

# Retrieve inputs
# args = commandArgs(TRUE)
# input_folder = args[1]
# output_folder = args[2]

# For testing:
input_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/input_files/"
output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/output_files2/"

# Create the output folder if it does not already exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

##################
# INPUT FILES
##################
## Construct the file paths of the input files necessary for the workflow and 
# check if the files exist

print("Locating input files...")

# File with the marker sequences (Forward, Reverse, Reference sequence)
primerFile = file.path(paste0(input_folder, "marker_file.txt"))
if(!file.exists(primerFile)) {
  stop("Primer file not found. Make sure that a file called markerFile.txt exists
       in the input folder")
}

# Sample file containing the IDs, forward and reverse barcode IDs, sample name and replicate 
sampleFile = file.path(paste0(input_folder, "sample_file.txt"))
if(!file.exists(primerFile)) {
  stop("Sample file not found. Make sure that a file called sampleFile.txt exists
       in the input folder")
}

# Fasta files with the barcodes ID and genomic sequence for forward and reverse
fnBarcodeF = file.path(paste0(input_folder, "barcode_F.fasta"))
if(!file.exists(fnBarcodeF)) {
  stop("Forward barcode file not found. Make sure that a file called barcode_Fwd.fasta exists
       in the input folder")
}

fnBarcodeR = file.path(paste0(input_folder, "barcode_R.fasta"))
if(!file.exists(fnBarcodeR)) {
  stop("Reverse barcode file not found. Make sure that a file called barcode_Rev.fasta exists
       in the input folder")
}

# All the fastq.gz files containing the raw read sequences
reads = list.files(input_folder, pattern="fastq.gz", full.names = T)
if(length(reads) == 0) {
  stop("Files containing the raw read sequences not found. Make sure that the corresponding fastq.gz files exist
       in the input folder")
}

################
# Demultiplexing
################

print("Demultiplexing by sample...")

# 1) Demultiplex by sample
# OUTPUT: 
# dePlexSample/ folder contains for each sample the corresponding reads (.fastq.gz)
# SummaryDemultiplexSample.txt contains for each sample and barcode pair the number of reads and the paths to the files with the reads

# Create output subdirectory to store the results of demultiplexing by sample
outDeplexSampleDir = file.path(output_folder, "dePlexSample")
if (dir.exists(outDeplexSampleDir)) {
  print("Removing already existing folder for sample demultiplexing.")
  unlink(outDeplexSampleDir, recursive = TRUE)
}
dir.create(outDeplexSampleDir)

# Demultiplex reads by samples
dePlexSample = demultiplexReads(reads[1], reads[2], fnBarcodeF, 
                                fnBarcodeR, outDeplexSampleDir)

# Rename output files to sample files
sampleTab = read.delim(sampleFile, stringsAsFactors = F)
dePlexSample = renameDemultiplexedFiles(sampleTab, dePlexSample)

#Remove samples with NA in SampleID --> phiX and any other non existing barcode combinations
dePlexSample = dePlexSample[!is.na(dePlexSample$SampleID),]

# Save summary table for the demultiplexing by sample operation
write.table(dePlexSample, file.path(output_folder, "SummaryDemultiplexSample.txt"), sep="\t", row.names=F)


print("Demultiplexing by marker ...")

# 2) Demultiplexing by marker 
# OUTPUT: 
# dePlexMarker/ folder contains for each sample and marker the corresponding reads (.fastq.gz)
# SummaryDemultiplexMarker.txt contains for each sample, marker and barcode pair the number of reads and the paths to the files containing the reads

# Create output subdirectory to store the results of demultiplexing by marker
outDeplexMarkerDir = file.path(output_folder, "dePlexMarker")
if (dir.exists(outDeplexMarkerDir)) {
  print("Removing already existing folder for marker demultiplexing.")
  unlink(outDeplexMarkerDir, recursive = TRUE)
} 
dir.create(outDeplexMarkerDir)

# Process each marker
markerTab = read.delim(primerFile, stringsAsFactors = F)
dePlexMarker = demultiplexByMarker(dePlexSample, markerTab, outDeplexMarkerDir)
dePlexMarker = dePlexMarker[!is.na(dePlexMarker$FileR1),]

# Save summary table
write.table(dePlexMarker, file.path(output_folder, "SummaryDemultiplexMarker.txt"), sep="\t", row.names=F)

################
# Merging reads
################

print("Merging reads ...")

## Create output subdirectory for storing the merged reads
outProcFilesDir = file.path(output_folder, "processedReads/")
if (dir.exists(outProcFilesDir)) {
  print("Removing already existing folder for merging reads.")
  unlink(outProcFilesDir, recursive = TRUE)
}
dir.create(outProcFilesDir)

postfix = "_merge"
refSeq = DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) = markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(output_folder, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

procReads = mergeAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), 
                               outProcFilesDir, method = "vsearch")
procReads = cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
procReads = procReads[-which(procReads$SampleName == "NTC"),]  
procReads = procReads[procReads$numRead >= 190,]
write.table(procReads, file.path(output_folder, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F, quote=F)

###############
# Calling SNPs
###############

print("Calling SNPs ...")

# Calculate mismatch rate and call SNPs

# Options
minMMrate <- 0.05
minOccGen <- 2

# process each marker
snpLst <- lapply(markerTab$MarkerID, function(marker){
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]), 
                                            refSeq[marker], 
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                            minCoverage=100L)
  names(seqErrLst) <- procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(output_folder, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate = minMMrate, minReplicate = minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  
  snps <- data.frame(Chr = marker, Pos = potSNP, Ref = snpRef, Alt = "N", stringsAsFactors = F)
  write.table(snps, file = file.path(output_folder, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                      minMMrate*100, minOccGen, marker, postfix)), 
              row.names = F, col.names = T, sep = "\t", quote = F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(output_folder, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                                   minMMrate*100, minOccGen, marker, postfix)), 
      width=1500 , height=600)
  matplot(seqErr, type = "p", pch = 16, cex = 0.4, col = "#00000088", ylim = c(0, 1),
          ylab = "Mismatch Rate", xlab = "Base Position", main = marker, cex.axis = 2, cex.lab = 2)
  abline(v = snps[,"Pos"], lty = 2, col = "grey")
  abline(h = minMMrate, lty = 1, col = "red")
  dev.off()
  
  return(snps)
})
names(snpLst) = markerTab$MarkerID

###############
# Calling Haplotypes
###############

print("Calling haplotypes ...")

# call haplotype options
minCov = 3
detectionLimit = 1/100
minOccHap = 2 
minCovSample = 25

# call final haplotypes
finalTab <- createFinalHaplotypTable(
  outputDir = output_folder, sampleTable = procReads, markerTable = markerTab, referenceSeq = refSeq,
  snpList = snpLst, postfix = postfix, minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample)

write.table(finalTab$csp, file.path(output_folder, "finalTabHaplotype_csp.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$cpmp, file.path(output_folder, "finalTabHaplotype_cpmp.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$`ama1-D3`, file.path(output_folder, "finalTabHaplotype_ama1-D3.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$msp7, file.path(output_folder, "finalTabHaplotype_msp7.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$cpp, file.path(output_folder, "finalTabHaplotype_cpp.txt"), sep="\t", row.names=F, quote=F)

print(Sys.time())
end_time = Sys.time()

print(as.difftime(end_time - start_time, format = "%X", units = "auto", tz = "CET"))

