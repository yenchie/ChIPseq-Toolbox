# peakAnnt.Gm189.v.0.1.R
# author: YCW
# date: 20231115
# Description: run DiffBind analysis
# Version: GenomicFeatures_1.46.5; ChIPseeker_1.30.3
# Input:
#    1. samplesheet: 
#         column: "SampleID"   "Factor"     "bamReads"   "Condition"  "ControlID" 
#                 "bamControl" "Peaks"      "Tissue"     "PeakCaller" "Replicate" 
#                 "Treatment"
#
# Output:
#    1. DiffBind_normalized.RData: data that done dba.counting and normalized
#    2. dba.report.csv: differentail bind table
#    3. plot:
#           - corraletion matrix
#           - heatmap
#           - PCA plot
#           - MA plot
#           - volcano plot
#           - affinity box plot
# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf

rm(list=ls())
library(dplyr)
library(GenomicFeatures)
library(ChIPseeker)
#sessionInfo()


args = commandArgs(trailingOnly=TRUE)
cat("
          #USAGE:Rscript --vanilla /bcst/JYL/YCWang/Rscript/peakAnnt.Gm189.v.0.1.R [arg1] [arg2]\n
          #Arg1: input folder path, ex. /bcst/JYL/YCWang/testing/ATAC_test/peakannt/bed \n
          #Arg2: output file path ex. /bcst/JYL/YCWang/testing/ATAC_test/peakannt \n
          ")

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("2 arguments must be supplied (input).n", call.=FALSE)
}

#----------------------
setwd(args[2])

### usage: ###
# makeTxDbFromGFF(file,
#                 format=c("auto", "gff3", "gtf"),
#                 dataSource=NA,
#                 organism=NA,
#                 taxonomyId=NA,
#                 circ_seqs=NULL,
#                 chrominfo=NULL,
#                 miRBaseBuild=NA,
#                 metadata=NULL,
#                 dbxrefTag)


## run-----
# txdb <- makeTxDbFromGFF(file="/bcst/JYL/JYL_qnap/db/Gm/Gm_v_1.1/Gmax189_gene_exons.gff3",
#                         format = "gff3",
#                         dataSource="http://www.phytozome.org/",
#                         organism="Glycine max")
# 
# saveDb(txdb, "Gm189")
txdb <- loadDb("/bcst/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/TxDb/Gm189")
#gene<-genes(txdb)%>%as.data.frame()

datapath<- args[1]
bedfile = datapath%>%list.files(full.names = T)
peak <- readPeakFile(bedfile)

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb = txdb )


# Output annotation
dfGRanges = data.frame(as.GRanges(peakAnno))
write.table(dfGRanges,file = "diffPeakAnnotation.txt",sep="\t",row.names=F,quote = F)

# Visualize Genomic Annotation
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)

upsetplot(peakAnno, vennpie=TRUE)

# Visualize distribution of TF-binding loci relative to TSS

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")



# ChIP peak annotation comparision
# Profile of several ChIP peak data binding to TSS region
# Average profiles
files = list.files(datapath, full.names = T)[13:14]

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

# ChIP peak annotation comparision
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
