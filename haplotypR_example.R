###########################
# Example of running HaplotypR on an example dataset
#
# created 21.10.2022
# monica.golumbeanu@unibas.ch
###########################

# Load the necessary packages
library("HaplotypR")
library("ShortRead")

# Define output directory and create it
outputDir = "~/genotyping/exampleHaplotypR/"  
if(!dir.exists(outputDir)) {
  dir.create(outputDir, recursive=T)
}
  
# Copy example files to output directory
file.copy(from = system.file(package="HaplotypR", "extdata"), to = outputDir, recursive = T)
raw_data_file = paste0(outputDir, "extdata/")

# List files example files in output directory
dir(raw_data_file)

# Description of the example data: 414859 reads, 2 markers (csp, cpmp)

##################
# Demultiplexing
#################

### Demultiplexing by sample ###
# set input file path
# File with the marker sequences (Forward, Reverse, Reference sequence)
primerFile = paste0(raw_data_file, "markerFile.txt")
# Sample file containing the IDs, forward and reverse barcode IDs, sample name and replicate 
sampleFile = paste0(raw_data_file, "sampleFile.txt")
# Fasta files with the barcodes ID and genomic sequence for forward and reverse
fnBarcodeF = paste0(raw_data_file, "barcode_Fwd.fasta")
fnBarcodeR = paste0(raw_data_file, "barcode_Rev.fasta")
# All the fastq.gz files containing the read raw sequences
reads = list.files(raw_data_file, pattern="reads", full.names = T)

# create output subdirectory 
outDeplexSample = file.path(outputDir, "dePlexSample")
dir.create(outDeplexSample)

# Demultiplex reads by samples
dePlexSample = demultiplexReads(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample)

# Rename output files to sample files
sampleTab = read.delim(sampleFile, stringsAsFactors=F)
dePlexSample = renameDemultiplexedFiles(sampleTab, dePlexSample)
# 414859 reads -> 1659436 lines in the fastq file

# Save summary table for the demultiplexing by sample operation
write.table(dePlexSample, file.path(outputDir, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)
# dePlexSample contains for each barcode the number of corresponding reads and the generated file

### Demultiplexing by marker and truncate primer sequence ###
# create output subdirectory 
outDeplexMarker = file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
markerTab = read.delim(primerFile, stringsAsFactors=F)
dePlexMarker = demultiplexByMarker(dePlexSample, markerTab, outDeplexMarker)

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)

### Trim to a fixed length ###
# create output subdirectory 
outProcFiles = file.path(outputDir, "processedReads")
dir.create(outProcFiles)

# Trim options for trimming to a fixed length
numNtF = 190
numNtR = 120
postfix = sprintf("_bind%.0f_%.0f", numNtF, numNtR)

# Adjust reference to trim options and save as fasta file
refSeq = as.character(markerTab$ReferenceSequence)
refSeq = DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
names(refSeq) = markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

# Fuse paired reads
procReads = bindAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles, 
                               read1Length = numNtF, read2Length = numNtR)
procReads = cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F)

############
# Call SNPs
############

# Options
minMMrate = 0.5
minOccGen = 2

# process each marker
snpLst = lapply(markerTab$MarkerID, function(marker){
  # Calculate mismatch rate
  seqErrLst = calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]), 
                                            refSeq[marker], 
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                            minCoverage=100L)
  names(seqErrLst) = procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr = do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP = callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef = unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  snps = data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt = "N", stringsAsFactors = F)
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                      minMMrate*100, minOccGen, marker, postfix)), 
              row.names = F, col.names = T, sep="\t", quote = F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                                   minMMrate*100, minOccGen, marker, postfix)), 
      width = 1500 , height = 600)
  matplot(seqErr, type = "p", pch = 16, cex = 0.4, col = "#00000088", ylim = c(0, 1),
          ylab = "Mismatch Rate", xlab = "Base Position", main = marker, cex.axis = 2, cex.lab = 2)
  abline(v = snps[,"Pos"], lty = 2, col = "grey")
  abline(h = minMMrate, lty = 1, col = "red")
  dev.off()
  
  return(snps)
})
names(snpLst) = markerTab$MarkerID

#################
# Call haplotypes
#################

# call haplotype options
minCov = 3
detectionLimit = 1/100
minOccHap = 2
minCovSample = 25

# remove samples without reads
procReads = procReads[procReads$numRead>0,]

# call final haplotypes
finalTab = createFinalHaplotypTable(
  outputDir = outputDir, sampleTable = procReads, markerTable = markerTab, referenceSeq = refSeq,
  snpList = snpLst, postfix = postfix, minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample)
