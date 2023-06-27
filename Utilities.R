#################################################################
# Utilities.R
# Written by Sami Ul Haq
# Date updated: March 16, 2020
# THis file contains numerous functions that help with 
# analysis of cfMeDIP data
#################################################################


library(MEDIPS)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

library(Repitools)
library(DESeq2)

library(GenomicRanges)
library(MEDIPS)
library(gtools)

library(annotatr)
library(ggplot2)



medips.set.creator <- function(fileNames, outputFormat, chr.start=1, chr.end=22) {
  #' medips.set.creator
  #' 
  #' @description 
  #' This function makes medips dataset
  #' using a BAM file and exports it as a WIG file 
  #' default window size is 300bp
  #' 
  #' 
  #' @usage 
  #' medips.set.creator(list.of.bam.files, "rpkm")
  #' medips.set.creator(list.of.bam.files, "counts")
  #' 
  
  if(chr.end == "X") {
    chr.select <- paste0("chr", c(1:22, "X"))
  } else {
    chr.select <- paste0("chr", c(1:22))
  }
  
  
  # reference genome
  BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
  
  # to prevent PCR overamplification to mess with results
  # a max number of stacked reads per genomic position
  uniq <- 1e-3
  
  # reads extended to 300 nucleotides
  extend <- 300
  
  # reads are shifted by a specific # of nucleotides
  shift <- 0
  
  # genome divided into windows of size (ws)
  ws <- 300 #MEDIPS default is 100  

  for(i in fileNames) {
    # create MeDIP set for hESC experimental conditions
    exp_MeDIP_set = MEDIPS.createSet(file=i, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select = chr.select)
    
    wig_name <- paste0( substr(i, 1, nchar(i)-4), ".wig")
    #cat(wig_name, "\n")
    # exporting medip data as WIG file
    MEDIPS.exportWIG(exp_MeDIP_set, file=wig_name, format=outputFormat)
  }
  
}


wig.matrix.maker <- function(fileNames, chrom.start, chrom.end) {
  #'
  #' @description 
  #' This function makes a counts matrix from wig files of 300bp window sizes
  #' using the given range of chromosome
  #' 
  #' @usage 
  #' wig.matrix.maker(list.of.wig.files, 1, 22)
  #' 
  
  genome_300bp <- data.frame(genomeBlocks(BSgenome.Hsapiens.UCSC.hg19, chrs=seqnames(BSgenome.Hsapiens.UCSC.hg19)[chrom.start:chrom.end], width=300))
  genome_300bp$name_format <- paste(genome_300bp$seqnames, genome_300bp$start, genome_300bp$end, sep=".")

  # for any rows that are not 300 nt long (ends of chromosomes) had a wierd 295 width element
  genome_300bp <- subset(genome_300bp, genome_300bp$width == 300)

  # creates an empty matrix
  boilerplate <- matrix(nrow=length(genome_300bp$name_format), ncol=0)
  # the names of the matrix rows are the 300 bp windows
  row.names(boilerplate) <- genome_300bp$name_format
  
  for(i in fileNames) {
    cat(i, "\n")
    # wig file from MEDIPS is imported
    wig_file <- import.wig(i)

    normal_wig <- GRanges(wig_file)

    norm_df <- data.frame(normal_wig)

    rm(wig_file)
    rm(normal_wig)

    # the counts are added to the matrix
    boilerplate <- cbind(boilerplate, norm_df$score)
  }
  
  # the columns are given names (the .wig ending is removed)
  colnames(boilerplate) <- substr(fileNames, 1, nchar(fileNames)-4)
  
  return(boilerplate)
}


medips.qc <- function(filename, chromosomes) {
  #'
  #' @description 
  #' This function performs QC metrics using the MEDIPS library
  #' including: saturation analysis & CPG coverage (pie & histogram)
  #' 
  #' This function requires to state what chromosomes the analysis will
  #' be done on i.e. for chromosomes 1 to 22, pass '1:22'
  #' 
  #' @usage 
  #' medips.qc(name.of.bam.file, 1:22)
  #' 
  
  cat("\nWorking with file: ", filename, "\n")
  
  # reference genome
  BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
  
  # to prevent PCR overamplification to mess with results
  # a max number of stacked reads per genomic position
  uniq <- 1e-3
  
  # reads extended to 300 nucleotides
  extend <- 300
  
  # reads are shifted by a specific # of nucleotides
  shift <- 0
  
  # genome divided into windows of size (ws)
  ws <- 300 #MEDIPS default is 100
  
  # Saturation analysis
  # tests whether given mapped reads enough to generate
  # saturated/reproducible coverage profile of ref genome
  # 
  # How: breaks dataset into two groups
  # that is further divided into nit subgroups with ws window_size
  # Iterated where for each subset, more reads are considered = genomic coverages more similar?
  # tests the BAM file
  
  # chromosomes 20:22
  sat_anal <- MEDIPS.saturation(file=filename, BSgenome=BSgenome, uniq=uniq, extend=extend, shift=shift, window_size=ws, chr.select = paste0("chr", chromosomes), nit=10, empty_bins=TRUE, rank=FALSE)
  
  # $maxEstCor = saves the estimated correlation (pearson)
  # $maxTruCor = saves the actual correlation of bam file (pearson)
  
  png(filename = paste0(filename, "_saturation_analysis", ".png"))
  #plots the estimated vs actual correlation of the saturation analysis
  sat_anal_graph <- MEDIPS.plotSaturation(sat_anal)
  dev.off()
  
  # Sequence Pattern Coverage
  # tests # of CpGs covered by the short reads
  # depth of CpG coverage is also tested
  cpg_coverage <- MEDIPS.seqCoverage(file=filename, pattern="CG", BSgenome=BSgenome, chr.select=paste0("chr", chromosomes), extend=extend, shift=shift, uniq=uniq)

  
  # visualizing coverage with pie chart
  png(filename=paste0(filename, "_sequence_coverage", ".png"))
  corr_graph <- MEDIPS.plotSeqCoverage(seqCoverageObj=cpg_coverage, type="pie", cov.level = c(0,1, 2, 3, 4, 5))
  dev.off()
  
  # visualizing coverage with histogram
  png(filename=paste0(filename, "_sequence_coverage_histogram", ".png"))
  corr_hist <- MEDIPS.plotSeqCoverage(seqCoverageObj=cpg_coverage, type="hist", t = 15, main="Sequence pattern coverage, Tumour")
  dev.off()
  
}

medips.cpg.enrichment.reference.genome <- function(BSgenome="BSgenome.Hsapiens.UCSC.hg19") {
  #'
  #' @description 
  #' This function calculates the reference hg19 human genome 
  #' RelH and GoGe values
  #' 
  
  if (is.null(BSgenome)) {
    stop("Must specify a BSgenome library.")
  }
  
  # saves the genomic database
  dataset = get(ls(paste("package:", BSgenome, sep = "")))
  
  cat(paste("Calculating CpG density for the reference genome",BSgenome, "...\n", sep = " "))

  # PDict is a container for storing a preprocessed dictionary of DNA patterns that can
  #later be passed to the matchPDict function for fast matching against a reference sequence
  CG <- DNAStringSet("CG")
  pdict0 <- PDict(CG)

  # makes an object on which bsapply (apply function for BSGenome) can be used
  # dictionary Of Patterns (CG) Against A Reference genome (hg19) using countPDict (type of matchPDict function)
  params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))

  #genome.CG = the numbe of CpGs within the reference genome
  genome.CG = sum(bsapply(params, pdict = pdict0))

  #cat("number of CG in reference genome: ", genome.CG, "\n")

  # counts the frequency of an alphabet in the reference genome
  params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify = TRUE)
  alphabet = bsapply(params)

  # total bases in entire genome
  genome.l = sum(as.numeric(alphabet))

  #cat("total number of bases in entire genome: ", genome.l, "\n")

  # index 1 = A   2 = C   3 = G   4 = T

  #regions.C = the number of Cs within the regions
  genome.C = as.numeric(sum(alphabet[2, ]))
  #regions.G = the number of Gs within the regions
  genome.G = as.numeric(sum(alphabet[3, ]))

  #cat("Total # of C's in ref genome: ", genome.C, "\n")
  #cat("Total # of G's in ref genome: ", genome.G, "\n")


  #genome.relH = the relative frequency of CpGs within the reference genome
  genome.relH = genome.CG/genome.l * 100

  #cat("genome.relH: ", genome.relH, "\n")

  #genome.GoGe = the observed/expected ratio of CpGs within the reference genome
  genome.GoGe = (genome.CG * genome.l)/(genome.C * genome.G)
  
  #cat("genome.GoGe: ", genome.GoGe, "\n")
  
  # garbage collection function = frees up memory for system
  #gc()
  
  vals.to.return <- list(genome.relH, genome.GoGe)
  return(vals.to.return)
}




medips.cpg.enrichment <- function(file, BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select, uniq=1e-3, extend=300, shift=0, paired=TRUE, genome.relH, genome.GoGe) {
  #'
  #' @description 
  #' This function does CPG enrichment analysis QC metrics on MEDIP BAM files
  #' using functions from the MEDIPS package. It specifically calculates:
  #' 1. enrichment.score.relH
  #' 2. enrichment.score.GoGe
  #' 
  #' @usage 
  #' medips.cpg.enrichment(bam.file, chromsomes.range)
  #' 
  
  
  if (is.null(BSgenome)) {
    stop("Must specify a BSgenome library.")
  }
  
  fileName = unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
  path = paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/")))) - 1], collapse = "/")
  if (path == "") {
    path = getwd()
  }
  if (!fileName %in% dir(path)) {
    stop(paste("File", fileName, " not found in", path, 
               sep = " "))
  }
  
  # saves the genomic database
  dataset = get(ls(paste("package:", BSgenome, sep = "")))
  
  # GRange.Reads saves the passed file name into a GRange object
  if (!paired) {
    GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)
  } else {
    GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)
  }
  
  # gets the number of chromosomes from the passed file GRange object
  if (length(unique(seqlevels(GRange.Reads))) > 1) {
    # mixedsort arranges the chromsomes names chr1 chr3 chr2 => chr1 chr2 chr3 (sequentially orders them)
    chromosomes = mixedsort(unique(seqlevels(GRange.Reads)))
  }
  
  if (length(unique(seqlevels(GRange.Reads))) == 1) {
    chromosomes = unique(seqlevels(GRange.Reads)) # no need to do mixedsort because only 1 level exists
  }
  
  # shows which reference dataset is being used
  cat(paste("Loading chromosome lengths for ", BSgenome, "...\n",sep = ""))
  
  # gets the length of the chromosomes that are being used
  chr_lengths = as.numeric(seqlengths(dataset)[chromosomes])
  
  # all ranges are restricted (?)
  ranges(GRange.Reads) <- restrict(ranges(GRange.Reads), +1)
  
  # number of chromosomes being used (used to make the number of rows for matrix below)
  total = length(chromosomes)
  
  cat("Calculating CpG density for given regions...\n")
  
  # splits the selected file's GRange object by CHROMOSOMES
  splitted.regions <- split(GRange.Reads, seqnames(GRange.Reads))
  
  # This function calculates the number of C, G, & CG in the selected regions passed into the function
  regions.c.g.cpg.finder <- function(x){
    # x = chromomsal regions
    
    i= which( mixedsort(chromosomes)%in%names(x) )
    # sets the range of the input chromosomal region to the length of that chromsome
    ranges(x) <- restrict(ranges(x),end=chr_lengths[which(chromosomes %in% seqnames(x))])
    
    # gets the associated DNA bases for that range
    y=DNAStringSet(getSeq(dataset, names=seqnames(x), start=start(x), end=end(x), as.character=TRUE))
    
    # finds number of C, G, CG sites as well as width and length of DNA sequence
    c(sum(as.numeric(vcountPattern("CG",y))),
      sum(as.numeric(vcountPattern("C",y))),
      sum(as.numeric(vcountPattern("G",y))),
      sum(as.numeric(width(y))),
      length(y))
  }
  
  # the data is converted to a matrix 
  # each row is chromosome (from input file)
  seq <- matrix(unlist(lapply(splitted.regions, regions.c.g.cpg.finder), use.names=FALSE), ncol=5, nrow=total, byrow=TRUE)
  
  # sums the columns
  Value = colSums(seq)
  
  # calculates how many bases were not processed
  unused = length(GRange.Reads) - Value[5]
  if (unused != 0)
    cat(unused, "unused sequences, limits out of range\n")
  
  # saves the total number of CG C and G in the regions (from the associated file)
  regions.CG = Value[1]
  regions.C = Value[2]
  regions.G = Value[3]
  all.genomic = Value[4]
  
  # regions.relH = the relative frequency of CpGs within the regions
  regions.relH = as.numeric(regions.CG)/as.numeric(all.genomic) *  100
  
  #cat("regions.relH: ", regions.relH, "\n")
  
  #regions.GoGe	= the observed/expected ratio of CpGs within the regions
  regions.GoGe = (as.numeric(regions.CG) * as.numeric(all.genomic)) / (as.numeric(regions.C) * as.numeric(regions.G))
  
  #cat("regions.GoGe: ", regions.GoGe, "\n")
  
  
  
  
  # 
  # cat(paste("Calculating CpG density for the reference genome",BSgenome, "...\n", sep = " "))
  # 
  # # PDict is a container for storing a preprocessed dictionary of DNA patterns that can 
  # #later be passed to the matchPDict function for fast matching against a reference sequence
  # CG <- DNAStringSet("CG")
  # pdict0 <- PDict(CG)
  # 
  # # makes an object on which bsapply (apply function for BSGenome) can be used
  # # dictionary Of Patterns (CG) Against A Reference genome (hg19) using countPDict (type of matchPDict function)
  # params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
  # 
  # #genome.CG = the numbe of CpGs within the reference genome
  # genome.CG = sum(bsapply(params, pdict = pdict0))
  # 
  # #cat("number of CG in reference genome: ", genome.CG, "\n")
  # 
  # # counts the frequency of an alphabet in the reference genome
  # params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify = TRUE)
  # alphabet = bsapply(params) 
  # 
  # # total bases in entire genome
  # genome.l = sum(as.numeric(alphabet))
  # 
  # #cat("total number of bases in entire genome: ", genome.l, "\n")
  # 
  # # index 1 = A   2 = C   3 = G   4 = T
  # 
  # #regions.C = the number of Cs within the regions
  # genome.C = as.numeric(sum(alphabet[2, ]))
  # #regions.G = the number of Gs within the regions
  # genome.G = as.numeric(sum(alphabet[3, ]))
  # 
  # #cat("Total # of C's in ref genome: ", genome.C, "\n")
  # #cat("Total # of G's in ref genome: ", genome.G, "\n")
  # 
  # 
  # #genome.relH = the relative frequency of CpGs within the reference genome
  # genome.relH = genome.CG/genome.l * 100
  # 
  # #cat("genome.relH: ", genome.relH, "\n")
  # 
  # #genome.GoGe = the observed/expected ratio of CpGs within the reference genome
  # genome.GoGe = (genome.CG * genome.l)/(genome.C * genome.G)
  
  #cat("genome.GoGe: ", genome.GoGe, "\n")
  
  
  
  
  
  #enrichment.score.relH = regions.relH/genome.relH
  enrichment.score.relH = regions.relH/genome.relH
  
  # enrichment.score.GoGe	= regions.GoGe/genome.GoGe
  enrichment.score.GoGe = regions.GoGe/genome.GoGe
  
  # garbage collection function = frees up memory for system
  gc()
  
  # # puts the values to return in a convenient package
  # values.to.return <- list(genome = BSgenome, regions.CG = regions.CG,
  #                          regions.C = regions.C, regions.G = regions.G, regions.relH = regions.relH,
  #                          regions.GoGe = regions.GoGe, genome.C = genome.C, genome.G = genome.G,
  #                          genome.CG = genome.CG, genome.relH = genome.relH, genome.GoGe = genome.GoGe,
  #                          enrichment.score.relH = enrichment.score.relH, enrichment.score.GoGe = enrichment.score.GoGe)
  # 
  cat("\n\nWorking with file: ", file, "\n")
  cat("\nEnrichment.score.relH: ", enrichment.score.relH, "\n")
  cat("\nEnrichment.score.GoGe: ", enrichment.score.GoGe, "\n\n")
  
  # for successful MeDIPS: enrichment.score.relH value >3 and enrichment.score.GoGe value >1.7.
}


medips.methylation.profile.maker <- function(bams, chromosomes.used) {
  #' 
  #' @description 
  #' This function takes in a bunch of bam files and 
  #' makes a the methylation profile (makes a MEDIPS set)
  #' 
  #' @usage 
  #' medips.methylation.profile.maker(list.of.bam.files, 1:22)
  #' 
  
  # reference genome
  BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
  
  # to prevent PCR overamplification to mess with results
  # a max number of stacked reads per genomic position
  uniq <- 1e-3
  
  # reads extended to 300 nucleotides
  extend <- 300
  
  # reads are shifted by a specific # of nucleotides
  shift <- 0
  
  # genome divided into windows of size (ws)
  ws <- 300 #MEDIPS default is 100
  
  # holds the bam file MEDIP data
  meth.profile.set <- c()
  
  for(i in bams) {
    # each file is turned into a medips set and added to the holder vector
    currently.processed <- MEDIPS.createSet(file=i, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select = paste0("chr", chromosomes.used))
    meth.profile.set <- c(meth.profile.set, currently.processed)
  }
  
  # methylation profile made
  methylation.profile = MEDIPS.meth(MSet1=meth.profile.set)
  
  return(methylation.profile)
}


dmr.analysis.medips <- function(set.1.bams, set.2.bams, chromosomes.used, significance.threshold=0.05) {
  #' 
  #' @description 
  #' This function analyzes for DMRs in two different sets of bam files
  #' using the MEDIPS method
  #' 
  #' @usage 
  #' dmr.analysis.medips(set.of.1st.bam.files, set.of.2nd.bam.files, chromosomes.used, sig.cutoff)
  #' 
  
  # reference genome
  BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
  # to prevent PCR overamplification to mess with results
  # a max number of stacked reads per genomic position
  uniq <- 1e-3
  # reads extended to 300 nucleotides
  extend <- 300
  # reads are shifted by a specific # of nucleotides
  shift <- 0
  # genome divided into windows of size (ws)
  ws <- 300 #MEDIPS default is 100
  
  # holds the MEDIPS set for set 1
  set.1 <- c()
  
  # holds the MEDIPS set for set 2
  set.2 <- c()
  
  for(i in set.1.bams) {
    # each file is turned into a medips set and added to the holder vector
    currently.processed.set1 <- MEDIPS.createSet(file=i, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select = paste0("chr", chromosomes.used))
    set.1 <- c(set.1, currently.processed.set1)
  }
  
  for(j in set.2.bams) {
    # each file is turned into a medips set and added to the holder vector
    currently.processed.set2 <- MEDIPS.createSet(file=j, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select = paste0("chr", chromosomes.used))
    set.2 <- c(set.2, currently.processed.set2)
  }
  
  
  ######################
  #Revise for this function
  
  # create coupling set for normalizing CpG density
  CS = MEDIPS.couplingVector(pattern = "CG", refObj = set.1)
  
  # calculating methylation profiles and DMRs
  DMR_analysis = MEDIPS.meth(MSet1=set.1, MSet2=set.2, CSet = CS, p.adj = "bonferroni", diff.method = "edgeR", MeDIP = T, CNV = F, minRowSum = 10)
  
  #selects significant DMR windows
  DMR_significants = MEDIPS.selectSig(results=DMR_analysis, p.value=significance.threshold, adj=T, ratio=NULL, bg.counts=NULL, CNV=F)
  
  # the significant windows is returned as a list
  return(DMR_significants)
  
}



