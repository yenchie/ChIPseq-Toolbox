## peakannt.chipseeker.fun.v.1.R
# author: YCW
# date:2024.03.19
# revised: 2024.05.13
# reviser: YCW

# dependency:
# beds2grl.fun.R
# plotAnntBar.fun.v.1.R
# txdb

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
# getCurrentFileLocation <- function() {
#   this_file <- commandArgs() %>%
#     tibble::enframe(name = NULL) %>%
#     tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
#     dplyr::filter(key == "--file") %>%
#     dplyr::pull(value)
#   if (length(this_file) == 0) {
#     this_file <- rstudioapi::getSourceEditorContext()$path
#   }
#   return(dirname(this_file))
# }

# source.path <- getCurrentFileLocation()
# source(file.path(source.path %>% list.files("plotAnntBar.fun.v.1.R", full.names = T)))


annt.gene.fun <- function(x) {
  d1 <- as.data.frame(x)
  print(nrow(d1))
  d1 <- d1 %>% filter(!annotation %>% str_detect("Distal Intergenic"))
  print(nrow(d1))
  print(d1$annotation %>% str_extract("^\\w+") %>% unique())
  print(length(d1$geneId %>% unique()))
  head(d1)
  d2 <- data.frame(geneId = c(d1$geneId, d1$flank_geneIds)) %>%
    separate_longer_delim(geneId, delim = ";") %>%
    distinct() %>%
    na.omit()
  print(nrow(d2))
  head(d2)
  return(d2$geneId %>% unique())
}


## core fun -----

print("files formate: col1: path; col2: ID # names of peaksets")
print("0-base coordinate system in use \n please check your input files format")

peak.Anno.fun <- function(files = NULL, grl = NULL, outpath = "./", date = date, plot = F, output.gene.table = T, Txdb = Txdb) {
  txdb <- loadDb(Txdb)
  print(outpath)

  if (is.null(grl)) {
    grl <- beds2grl(files)
  } else {
    print("grl as input")
  }

  head(grl)
  length(grl)
  names(grl)

  promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 0) # signed your definition of promoter
  options(ChIPseeker.downstreamDistance = 0)

  grl1 <- grl[which(lengths(grl) != 0)]
  head(grl1)
  length(grl1)
  names(grl1)

  peakAnnoList <- lapply(grl1, annotatePeak,
    TxDb = txdb,
    tssRegion = c(-1000, -1), # signed your definition for TSS
    verbose = TRUE,
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c(
      "Promoter", "5UTR", "3UTR", "Exon", "Intron",
      "Downstream", "Intergenic"
    ),
    level = "gene", sameStrand = F, addFlankGeneInfo = TRUE, flankDistance = 0
  )
  genes <- lapply(peakAnnoList, function(i) {
    annt.gene.fun(i)
  })
  head(genes)
  length(genes)
  names(genes)

  if (output.gene.table == T) {
    gene.table <- gene.table.0 <- NULL
    for (l in names(grl)) {
      d1 <- NULL
      x <- NULL
      if (l %in% names(genes)) {
        x <- genes[[l]]
        head(x)
        d1 <- data.frame(geneID = x %>% unlist()) %>%
          mutate(group = l)
        print(nrow(d1))
        gene.table <- gene.table %>% bind_rows(d1)

        write.table(d1$geneID, file.path(outpath, paste0(l, "ChIPseeker.geneAnno.", date, ".gene.id")), quote = F, row.names = F, sep = "\t", col.names = F)
      } else {
        d1 <- NULL
        print(paste(id, "none annotated genes, output empty file"))
        write.table(d1, file.path(outpath, paste0(l, ".ChIPseeker.geneAnno.", date, ".gene.id")), quote = F, row.names = F, sep = "\t", col.names = F)
        gene.table.0 <- gene.table.0 %>% bind_rows(data.frame(group = l, n = 0))
      }
    }
    print(file.path(outpath, paste0(l, ".ChIPseeker.geneAnno.", date, ".gene.id")))

    head(gene.table)
    gene.count <- gene.table %>%
      dplyr::count(group) %>%
      bind_rows(gene.table.0) %>%
      dplyr::rename(count = n)

    print(file.path(outpath, paste0("ChIPseeker.gene.number.", date, ".txt")))
    write.table(gene.count, file.path(outpath, paste0("ChIPseeker.gene.number.", date, ".txt")), quote = F, row.names = F, sep = "\t")
  }

  if (plot == T) {
    # plotting -----
    print(file.path(outpath, paste0("ChIPseeker.plot.", date, ".pdf")))
    pdf(file.path(outpath, paste0("ChIPseeker.plot.", date, ".pdf")))

    tryCatch(
      {
        makeVennDiagram(grl1,
          NameOfPeaks = c(names(grl1)),
          scaled = FALSE, euler.d = FALSE, totalTest = lengths(grl1) %>% max(),
          connectedPeaks = "keepAll",
          fill = rainbow(length(grl1)), # circle fill color
          col = rainbow(length(grl1)), # circle border color
          cat.col = rainbow(length(grl1))
        )
      },
      error = function(e) {
        cat("An error occurred:", conditionMessage(e), "\n")
      }
    )

    print(plotAnntBar(peakAnnoList))

    # print(plotAnnoBar(peakAnnoList))
    print(plotDistToTSS(peakAnnoList))

    # gene annotation intersection ------------
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
            paste("Intersection of Gene Sets in Venn"),
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

  return(peakAnnoList)
}
## fin---
