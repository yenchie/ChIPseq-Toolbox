## peakannt.chipseeker.fun.v.1.R
# author: YCW
# date:2024.03.19

library(dplyr)
library(stringr)
library(GenomicFeatures)
library(ChIPseeker)
library(DOSE)
library(ggplot2)
library(ggpubr)
library(forcats)
library(tibble)

rm(list = ls())

## creat your references ----
## example:
# seqlength <- read.delim("./ref.genome.faCount", stringsAsFactors = F, header = T) %>%
#   dplyr::select(X.seq, len) %>%
#   filter(!X.seq %>% str_detect("total")) # input the faCount of your reference

# txdb <- makeTxDbFromGFF(
#   file = "./ref.genome_gene_exons.gff3", # revise for your ref.
#   format = "gff3",
#   dataSource = "http://www.phytozome.org/",
#   organism = "Glycine max", # revise for your species
#   chrominfo = Seqinfo(
#     seqnames = seqlength$X.seq,
#     seqlengths = seqlength$len,
#     genome = "Gm189"
#   )
# ) # revise for your ref.
# saveDb(txdb, "./ref.genome")

## function ----
getCurrentFileLocation <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

source.path <- getCurrentFileLocation()
source(file.path(source.path %>% list.files("plotAnntBar.fun.v.1.R", full.names = T)))


annt.gene.fun <- function(x) {
  d1 <- as.data.frame(x)
  print(nrow(d1))
  d1 <- d1 %>% filter(!annotation %>% str_detect("Distal Intergenic"))
  print(nrow(d1))
  print(d1$annotation %>% str_extract("^\\w+") %>% unique())
  print(length(d1$geneId %>% unique()))
  return(d1$geneId %>% unique())
}


## core fun -----

txdb <- loadDb("/bcst/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/TxDb/Gm189")


print("files formate: col1: path; col2: ID # names of peaksets")

peak.Anno.fun() <- function(files, outpath, date) {
  print(outpath)

  # beds.path <- "./bed"
  # date <- " "

  ## input -----

  # files <- data.frame(
  #   path = beds.path %>% list.files("bed$", full.names = T)
  # ) %>%
  #   mutate(ID = path %>% as.character() %>% basename()) # names for peakset

  grl1 <- NULL
  grl <- GRangesList()
  for (id in 1:unique(files$ID)) {
    gr <- read.table(files$path[bedfile$ID == id], header = FALSE)
    colnames(gr)[1:3] <- c("seqnames", "start", "end")
    gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE, starts.in.df.are.0based = T)
    ID <- paste(files$id[bedfile])
    grl[[ID]] <- gr
  }

  head(grl)
  length(grl)
  names(grl)
  grl1 <- grl

  # plotting -----
  print(file.path(outpath, paste0("ChIPseeker.plot.", ".", date, ".pdf")))
  pdf(file.path(outpath, paste0("ChIPseeker.plot.", ".", date, ".pdf")))
  promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 0) # signed your definition of promoter
  options(ChIPseeker.downstreamDistance = 0)
  peakAnnoList <- lapply(grl1, annotatePeak,
    TxDb = txdb,
    tssRegion = c(-1000, -1), # signed your definition for TSS
    verbose = TRUE,
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c(
      "Promoter", "5UTR", "3UTR", "Exon", "Intron",
      "Downstream", "Intergenic"
    ),
    level = "gene", sameStrand = F
  )
  return(peakAnnoList)

  makeVennDiagram(grl1,
    NameOfPeaks = c(names(grl1)),
    scaled = FALSE, euler.d = FALSE, totalTest = lengths(grl1) %>% max(),
    connectedPeaks = "keepAll",
    fill = c("white", "white", "white"), # circle fill color
    col = c("#D55E00", "#0072B2", "#009E73"), # circle border color
    cat.col = c("#D55E00", "#0072B2", "#009E73")
  )


  print(plotAnntBar(peakAnnoList))

  # print(plotAnnoBar(peakAnnoList))
  print(plotDistToTSS(peakAnnoList))
  # gene annotation intersection ------------

  genes <- lapply(peakAnnoList, function(i) {
    annt.gene.fun(i)
  })

  gene.table <- NULL
  for (l in names(genes)) {
    x <- NULL
    x <- genes[[l]]
    head(x)
    d1 <- data.frame(geneID = x %>% unlist()) %>%
      mutate(group = l)
    gene.table <- gene.table %>% bind_rows(d1)
  }
  head(gene.table)
  print(file.path(outpath, paste0("ChIPseeker.genetable.", ".", date, ".txt")))
  write.table(gene.table, file.path(outpath, paste0("ChIPseeker.genetable.", ".", date, ".txt")), quote = F, row.names = F, sep = "\t")


  # Count the number of variables in the list
  num_variables <- sapply(genes, function(x) length(unlist(x)))
  cat(num_variables)
  element_lengths <- lengths(genes)
  result_df <- data.frame(method = names(genes), n = element_lengths) %>%
    mutate(Group = g)
  res <- res %>% bind_rows(result_df)

  tryCatch(
    {
      print(vennplot(genes, by = "Vennerable"))
    },
    error = function(e) {
      cat("An error occurred:", conditionMessage(e), "\n")
    }
  )
  tryCatch(
    {
      print(vennplot(genes, by = "gplots"))
    },
    error = function(e) {
      cat("An error occurred:", conditionMessage(e), "\n")
    }
  )

  Vgene <- Venn(genes)
  print(plot(Vgene, doWeights = F))

  pcentFun <- function(x) {
    100 * (x / sum(x))
  }
  filtered_fin_list <- genes
  V <- Vgene

  tryCatch(
    {
      if (length(filtered_fin_list) <= 3) {
        VennList <- compute.Venn(V, doWeights = F, type = "circles")
      } else if (length(filtered_fin_list) == 4) {
        VennList <- compute.Venn(V, doWeights = F, type = "ellipses")
      } else {
        VennList <- compute.Venn(V, doWeights = F, type = "ChowRuskey")
      }
    },
    error = function(e) {
      VennList <- compute.Venn(V, doWeights = F, type = "AWFE")
      # cat("An error occurred:", conditionMessage(e), "\n")
    }
  )
  tryCatch(
    {
      Weight <- V@IndicatorWeight %>%
        as.data.frame() %>%
        rownames_to_column()
      Weight_sorted <- Weight[order(match(Weight$rowname, VennList@FaceLabels$FaceName)), ]
      areas <- Weight_sorted$.Weight
      names(areas) <- Weight_sorted$rowname
      areasPcent <- round(pcentFun(areas), digits = 2)
      VennList@FaceLabels$Signature <- paste0(areas, "\n", areasPcent %>% round(1), "%")

      # plot(VennList, show = list(FaceText = c("signature"), DarkMatter = F))
      gp <- VennThemes(VennList, increasingLineWidth = F, colourAlgorithm = "signature")
      modify_vector <- function(element, option, value) {
        if (option %in% names(element)) {
          element[[option]] <- value
        }
        return(element)
      }

      gp[["SetText"]] <- lapply(gp[["SetText"]], function(x) modify_vector(x, "fontsize", 12))
      gp[["FaceText"]] <- lapply(gp[["FaceText"]], function(x) modify_vector(x, "fontsize", 12))

      gridExtra::grid.arrange(
        grid::grid.grabExpr(height = 4, width = 6, plot(
          VennList,
          show = list(
            FaceText = c("signature"),
            DarkMatter = F
          ),
          gp = gp
        )),
        top = textGrob(
          paste(g, "Intersection of Gene Sets in Venn"),
          gp = gpar(fontsize = 15)
        )
      )
    },
    error = function(e) {
      cat("An error occurred:", conditionMessage(e), "\n")
    }
  )

  dev.off()
}
## fin---
