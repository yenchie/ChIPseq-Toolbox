# annoPeak
##' Class "annoPeak"
##' This class represents the output of ChIPseeker Annotation
##'
##'
##' @name annoPeak-class
##' @aliases annoPeak-class
##' show,annoPeak-method vennpie,annoPeak-method
##' plotDistToTSS,annoPeak-method plotAnnoBar,annoPeak-method
##' plotAnnoPie,annoPeak-method upsetplot,annoPeak-method
##' subset,annoPeak-method
##'
##' @docType class
##' @slot anno annotation
##' @slot tssRegion TSS region
##' @slot level transcript or gene
##' @slot hasGenomicAnnotation logical
##' @slot detailGenomicAnnotation Genomic Annotation in detail
##' @slot annoStat annotation statistics
##' @slot peakNum number of peaks
##' @exportClass annoPeak
##' @author YCW
##' @keywords classes

setClass("annoPeak",
    representation = representation(
        annoGene = "character",
        gr.ref = "GRanges",
        gr.peak = "GRanges",
        promoterRegion = "numeric",
        Annotation = "data.frame",
        annoStat = "data.frame",
        peakNum = "numeric"
    )
)

setMethod("show", "annoPeak", function(object) {
    cat("annoPeak Object\n")
    cat("all peak number:", object@peakNum, "\n")
    cat("annotated gene number:", length(object@annoGene), "\n")
    print(object@annoStat)
})
