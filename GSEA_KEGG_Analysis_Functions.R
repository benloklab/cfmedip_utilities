##################################################################################
# This is a Utility script with a bunch of KEGG, GSEA, GO, analysis functions



## Helpful link: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#go-enrichment-analysis

#library(clusterProfiler)

library(GenomicRanges)

library(annotatr)
library(org.Hs.eg.db)

#library(fgsea)


medip.regions.to.GRange.maker <- function(list.of.regions) {
  #' 
  #' medip.regions.to.GRange.maker
  #' 
  #' @description 
  #' This function takes as input a character vector of windows
  #' example. c("chr.1.100.400", "chr10.200.500)
  #' and converts it into a GRange
  #' 
  #' 
  
  
  chromosome.names <- sapply(list.of.regions, FUN=function(x) {
    strsplit(x, ".", fixed=TRUE)[[1]][1]
  })
  start.region <- sapply(list.of.regions, FUN=function(x) {
    strsplit(x, ".", fixed=TRUE)[[1]][2]
  })
  end.region <- sapply(list.of.regions, FUN=function(x) {
    strsplit(x, ".", fixed=TRUE)[[1]][3]
  })
  
  regions.to.data.frame <- data.frame(chr=chromosome.names, start=start.region, end=end.region)
  regions.to.GRange <- makeGRangesFromDataFrame(regions.to.data.frame)
  
  return(regions.to.GRange)
  
}


medip.Grange.to.gene.symbols <- function(medip.ranges.Grange) {
  #' 
  #' medip.Grange.to.gene.symbols
  #' 
  #' @description 
  #' This function takes as input a GRange of medip regions
  #' i.e. Chr1 100 - 400bp
  #' 
  #' and then finds the genes that are in those regions
  #' NOTE: this is only mapping to gene promoter regions, cpg islands
  #' 
  #' @usage 
  #' 
  #' To use this, first do:
  #' GRANGE.of.MeDIP.Regions <- medip.regions.to.GRange.maker(gene.windows)
  #' medip.Grange.to.gene.symbols(GRANGE.of.MeDIP.Regions)
  #' 
  #' 
  
  # the types of annotations that will be downloaded from reference
  #annots = c("hg19_genes_promoters", "hg19_cpg_islands", "hg19_basicgenes")
  annots = c("hg19_basicgenes")
  annotations = build_annotations(genome = 'hg19', annotations = annots)
  
  overlapping.regions <- subsetByOverlaps(annotations, medip.ranges.Grange)
  gene.symbols.from.overlaps <- overlapping.regions$symbol
  
  #   removes NAs
  gene.symbols.from.overlaps <- gene.symbols.from.overlaps[!is.na(gene.symbols.from.overlaps)]
  # removes duplicated regions
  gene.symbols.from.overlaps <- unique(gene.symbols.from.overlaps)
  return(gene.symbols.from.overlaps)
  
}


kegg.analysis <- function(list.of.gene.symbols) {
  #'
  #' kegg.analysis
  #' 
  #' @description 
  #' This function takes as input a vector of gene symbols,
  #' It then converts them to UNIPROT IDs
  #' and then runs the KEGG analysis on them
  #' 
  #' This function returns an enrichKEGG object
  #' 
  #' @usage 
  #' 
  #' kegg.object <- kegg.analysis(vector.of.gene.symbols)
  #' 
  #' look at the top pathways using:
  #' head(kegg.object, n=20)
  #' 
  #' look at a pathway in greater detail:
  #' browseKEGG(kegg.object, "hsa_ID_from_KEGG_object")
  #'
  
  ### KEGG Pathway Analysis
  
  # maps gene symbols to ENTREZID
  converted.to.entrez.ids <- bitr(list.of.gene.symbols, fromType = "SYMBOL", toType = c("ENTREZID", "UNIPROT"), OrgDb = org.Hs.eg.db)
  
  # maps the ENTREZ IDs
  KEGG.enrichment.analysis <- enrichKEGG(gene=converted.to.entrez.ids$UNIPROT , organism='hsa', keyType="uniprot")
  
  head(KEGG.enrichment.analysis, n=20)
  
  return(KEGG.enrichment.analysis)
  
}






GSEA.on.medip <- function(deseq2.results.object, pathways.reference.file="msigdb.v7.2.symbols.gmt") {
  #' 
  #' GSEA.on.medip
  #' 
  #' @description 
  #' This function takes as input, a deseq2 results object. The following happens:
  #' 1. NAs are omitted
  #' 2. dataset is sorted from highest LogFC to lowest
  #' 3. the regions are converted to gene symbols
  #' 4. analysis is done on a passed .gmt file (default = Hallmark gene set from MSigDB)
  #' 
  #' the function returns an fgsea object AND a ranked.genes.stats (ranked genes with gene symbols) AND the pathways reference file
  #' 
  #' @usage 
  #' 
  #' gsea.object <- GSEA.on.medip(deseq2.results.object)
  #' 
  #' # results are ordered by NES
  #' # NES = normalized enrichment score (enrichment score normalized to mean enrichment of random samples of the same size)
  #' # p-adj = Benjamini-Hochberg adjusted p-value
  #' # leadingEdge = vector with indexes of leading edge genes that drive the enrichment
  #' head(gsea.object[[1]][order(padj, -abs(NES)), ], n=10)
  #' 
  #' # one of the results is plotted using an enrichment plot
  #' plotEnrichment(gsea.object[[3]][["MIKKELSEN_MCV6_HCP_WITH_H3K27ME3"]], gsea.object[[2]])
  #' 
  #' 
  #' 
  #' 
  
  
  # converts Deseq2 results to data frame
  dsq.results <- data.frame(deseq2.results.object)
  # removes NAs
  dsq.results <- na.omit(dsq.results)
  
  barplot(dsq.results$log2FoldChange)
  
  
  # dataset is sorted from HIGHEST LogFC to LOWEST LogFC
  # THis step is CRITICAL for GSEA
  dsq.results <- dsq.results[order(-dsq.results$log2FoldChange),]
  
  # the names of the genes (Gene Symbols are extracted)
  ranked.genes.list <- medip.regions.to.GRange.maker(rownames(dsq.results))
  ranked.genes.list <- medip.Grange.to.gene.symbols(ranked.genes.list)
  
  # The DESEQ2 gene stats are extracted
  ranked.genes.stats <- dsq.results$stat
  # these gene stats are named with the Gene Symbols
  # THis is important because GSEA is done on the gene stats
  names(ranked.genes.stats) <- ranked.genes.list
  
  # plotting LogFCs
  barplot(dsq.results$log2FoldChange)
  
  # By default, using the Hallmark gene set from MSigDB. (msigdb.v7.2.symbols.gmt)
  # This file has been downloaded from Broad Institute (https://www.gsea-msigdb.org/gsea/downloads.jsp)
  pathways.hallmark <- gmtPathways( pathways.reference.file )
  
  fgseaRes <- fgsea(pathways.hallmark, ranked.genes.stats, minSize=15, maxSize = 500, nperm=1000)
  
  head(fgseaRes[order(padj, -abs(NES)), ], n=10)
  
  ggplot(fgseaRes[1:20], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") +
    theme_minimal()
  
  fgseaRes <- fgseaRes[order(fgseaRes$padj, -fgseaRes$NES), ]
  
  plotGseaTable(pathways.hallmark[fgseaRes$pathway[1:10]],
                ranked.genes.stats,
                fgseaRes,
                gseaParam = 0.5)
  
  
  return(list(fgseaRes, ranked.genes.stats, pathways.hallmark))
  
}




windows.overlapping.with.features <- function(grange.of.interest, grange.with.features.of.interest) {
  #' This function finds windows in one GRange that overlap with another GRange
  #' @param grange.of.interest
  #' @param grange.with.features.of.interest
  #' @usage windows.overlapping.with.features(GRange.1, GRange.2)
  #' 
  
  # finds the windows that overlap with grange.with.features.of.interest
  feature.windows <- subsetByOverlaps(grange.of.interest, grange.with.features.of.interest)
  feature.windows <- data.frame(feature.windows)
  
  # converts the chrm name, start, and end of each window to a string for comparison to the original matrix
  feature.windows.only.windows <- paste0(feature.windows$seqnames, ".", feature.windows$start, ".", feature.windows$end)
  return(feature.windows.only.windows)
}



genomic.window.annotator <- function(queryGrange, subjectGrange, removeDuplicateRegions=TRUE) {
  #' This function finds finds the overlapping windows between
  #' two GRange objects. This function also combines the metadata
  #' from both GRanges.
  #' Note: this function removes duplicated regions. Only unique
  #' regions are kept!
  #' @usage genomic.window.annotator(grange.Of.MeDIP.Hits, grange.of.ENSEMBL.annotations)
  
  m <- findOverlaps(query = queryGrange, subject = subjectGrange)
  rearranged.Grange <- queryGrange[queryHits(m)]
  mcols(rearranged.Grange) <- cbind.data.frame(mcols(rearranged.Grange),
                                               mcols(subjectGrange[subjectHits(m)]))
  
  if(removeDuplicateRegions == FALSE) {
    return(rearranged.Grange)
  } else {
    # only the unique elements are kept 
    rearranged.Grange <- rearranged.Grange[isUnique(rearranged.Grange)]
        return(rearranged.Grange)
    
  }
}



