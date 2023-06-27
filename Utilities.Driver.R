source("Utilities.R")

# this gets all the bam files of interest
bams <- Sys.glob("*.bam")

chromosomes.used <- 1:22
chr.select <- paste0("chr", chromosomes.used)


#################################################
# this does qc analysis on the BAM files

# genome RelH & GoGe values are calculated
genome.ref.values <- medips.cpg.enrichment.reference.genome()
genome.relH <- genome.ref.values[[1]]
genome.GoGe <- genome.ref.values[[2]]


# the MEDIPS specific QCs are done
for (bam.file in bams) {
  medips.qc(bam.file, chromosomes.used)
  medips.cpg.enrichment(bam.file, BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select = chr.select, genome.relH = genome.relH, genome.GoGe = genome.GoGe)

}


#################################################
# this part makes the count matrix from bam files and
# outputs as WIG files

medips.set.creator(bams, "count")

# uncomment this for rpkm format
# medips.set.creator(bams, "rpkm")

# all the WIG files that were created
wigs <- Sys.glob("*.wig")
# the matrix is made for chromosomes 1 to 22
matri <- wig.matrix.maker(wigs, 1, 22)



