#!/usr/bin/env Rscript

#library path and load libraries

library.path <- .libPaths("/scicore/home/gagneux/rutlil00/Rpackages")

library("devtools", lib.loc = library.path)
library("git2r", lib.loc = library.path)
library("HaplotypR", lib.loc = library.path)
library("ShortRead", lib.loc = library.path)

outputDir <- "scicore/home/mynaco16/GROUP/analysis/annina_thesis/run2_modified/001_haplotype_calling"  
# Create output directory
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)

primerFile <- "/scicore/home/mynaco16/GROUP/analysis/annina_thesis/run2/001_haplotype_calling/input_files/marker_file.txt"
sampleFile <- "/scicore/home/mynaco16/GROUP/analysis/annina_thesis/run2/001_haplotype_calling/input_files/sample_file_annina_thesis_run2.txt"
fnBarcodeF <- "/scicore/home/mynaco16/GROUP/analysis/annina_thesis/run2/001_haplotype_calling/input_files/barcode_F.fasta"
fnBarcodeR <- "/scicore/home/mynaco16/GROUP/analysis/annina_thesis/run2/001_haplotype_calling/input_files/barcode_R.fasta"
reads <- list.files("/scicore/home/mynaco16/GROUP/analysis/annina_thesis/run2/001_haplotype_calling/input_files", pattern="Pool4", full.names = T)

# create output subdirectory 
outDeplexSample <- file.path(outputDir, "dePlexSample")
dir.create(outDeplexSample)

# demultiplex by samples
dePlexSample <- demultiplexReads(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample)
##total 540000 demultiplexed reads

# rename output files to sample files
sampleTab <- read.delim(sampleFile, stringsAsFactors=F)
dePlexSample <- renameDemultiplexedFiles(sampleTab, dePlexSample)

#subset: remove samples with NA in SampleID --> phiX and any other non existing barcode combination
dePlexSample <- dePlexSample[!is.na(dePlexSample$SampleID),]

# save summary table
write.table(dePlexSample, file.path(outputDir, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)


# Run demultiplex by marker and truncate primer sequence

# create output subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
markerTab <- read.delim(primerFile, stringsAsFactors=F)
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, outDeplexMarker)

# remove samples without sequence reads
dePlexMarker <- dePlexMarker[!is.na(dePlexMarker$FileR1),]

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)


#Method for method work only for overlapping sequence read pair by merging the overlap of the forward and reverse read (using vsearch wrapper).


## create output subdirectory 
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

postfix <- "_merge"
refSeq <- DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

procReads <- mergeAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles)
procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F, quote=F)

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
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
 #browser()
  snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                      minMMrate*100, minOccGen, marker, postfix)), 
              row.names=F, col.names=T, sep="\t", quote=F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                                   minMMrate*100, minOccGen, marker, postfix)), 
      width=1500 , height=600)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
  abline(v=snps[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
  dev.off()
  
  return(snps)
})
names(snpLst) <- markerTab$MarkerID


#Call Haplotypes

# call haplotype options
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 2 
minCovSample <- 25

# call final haplotypes
finalTab <- createFinalHaplotypTable(
  outputDir = outputDir, sampleTable = procReads, markerTable = markerTab, referenceSeq = refSeq,
  snpList = snpLst, postfix = postfix, minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample)

write.table(finalTab$csp, file.path(outputDir, "finalTabHaplotype_csp.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$cpmp, file.path(outputDir, "finalTabHaplotype_cpmp.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$`ama1-D3`, file.path(outputDir, "finalTabHaplotype_ama1-D3.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$msp7, file.path(outputDir, "finalTabHaplotype_msp7.txt"), sep="\t", row.names=F, quote=F)
write.table(finalTab$cpp, file.path(outputDir, "finalTabHaplotype_cpp.txt"), sep="\t", row.names=F, quote=F)
