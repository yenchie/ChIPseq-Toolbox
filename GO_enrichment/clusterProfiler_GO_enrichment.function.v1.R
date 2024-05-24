# clusterProfiler_GO_enrichment.function.v1.R

# revised date: 2024.04.16
# add DAG plot
# 2024.02.06
# - revised if didn'tenrichment any result, print "no any GO terms be enrichment"

require(clusterProfiler)
require(topGO)
require(enrichplot)
require(dplyr)
require(stringr)
require(tidyr)
require(GO.db)
require(ggplot2)

GO_enricher <- function(datapath, outpath, cutoff, GO.path, termfile.path) {
  GO_ALL <- read.table(GO.path, sep = "\t", stringsAsFactors = FALSE)
  head(GO_ALL)
  colnames(GO_ALL) <- c("geneID", "ID", "Ontology")

  GO_BP <- GO_ALL %>% filter(Ontology == "BP")
  GO_MF <- GO_ALL %>% filter(Ontology == "MF")
  GO_CC <- GO_ALL %>% filter(Ontology == "CC")

  termfile <- read.table(termfile.path, sep = "\t", stringsAsFactors = FALSE, header = FALSE, fill = TRUE, quote = "")
  colnames(termfile) <- c("ID", "Description", "Ontology")
  termfile %>% head()
  termfile %>% nrow()
  length(termfile$V1 %>% unique()) %>% print()
  which(!GO_ALL$ID %in% termfile$ID) %>%
    length() %>%
    print()

  # query GO term description
  # query.GO <- GO_ALL %>%
  #   left_join(termfile) %>%
  #   filter(Description %>% is.na()) %>%
  #   dplyr::select(ID) %>%
  #   distinct()
  # term <- list()
  # for (go in query.GO$GO) {
  #   if (GOTERM[[go]] %>% is.null()) {
  #     print(paste("non description", go))
  #     cat("--------------------------------------\n")
  #   } else {
  #     term[[go]] <- GOTERM[[go]]
  #     print(GOTERM[[go]])
  #     cat("--------------------------------------\n")
  #   }
  # }

  # termfile.1 <- data.frame(
  #   ID = sapply(term, function(x) x@GOID),
  #   Description = sapply(term, function(x) x@Term)
  # )
  # termfile <- termfile %>% bind_rows(termfile.1)
  # GO_ALL <- GO_ALL %>%
  #   left_join(termfile)



  datapath <- file.path(datapath)
  print(paste("read id files in the folder:", datapath))
  print(paste("please notice that file name should be gene.id|TF.id in the end"))
  filelist <- list.files(file.path(datapath), pattern = "id$", full.names = T)
  print(paste("files in", datapath, "imported"))
  print(filelist)

  # create output path
  dir.create(file.path(datapath, "output", "plot"), recursive = T, mode = "0770")
  dir.create(file.path(datapath, "output", "GOterm"), recursive = T, mode = "0770")
  # file=1
  for (file in 1:length(filelist)) {
    gene.id <- read.delim(filelist[file], stringsAsFactors = F, header = F)$V1 %>% as.vector()
    head(gene.id)
    name <- filelist[file] %>%
      as.character() %>%
      basename() %>%
      str_remove(".id$")
    print(name)
    print(unique(GO_ALL$Ontology))

    # go="BP"
    for (go in unique(GO_ALL$Ontology)) {
      print(paste("start:", name, go))

      TERM2GENE <- GO_ALL %>%
        dplyr::filter(Ontology == go)
      head(TERM2GENE)
      annt.size <- count(TERM2GENE, ID)
      head(annt.size)
      print(paste("total input genes:", gene.id %>% length()))
      print(paste("gene in GOterm:", go, which(gene.id %in% TERM2GENE$geneID) %>% length()))
      # count(TERM2GENE, GO)
      # n_distinct(TERM2GENE$gene)
      ego <- enricher(
        gene.id,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        # universe	background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
        universe = NULL,
        qvalueCutoff = 1,
        minGSSize = 1,
        maxGSSize = max(annt.size$n),
        TERM2GENE = TERM2GENE[c("ID", "geneID")] # ,
        # TERM2NAME = TERM2GENE[c("GO", "descritpion")]
      )
      ego@ontology <- go
      pdf(file.path(datapath, "output", "plot", paste0(name, ".GO.enrichment.DAG.", go, ".pdf")))
      plotGOgraph(ego)
      dev.off()
      ego.t <- ego %>%
        as.data.frame() %>%
        mutate(Ontology = go)
      head(ego.t)
      colnames(ego.t)

      ego.t$qval <- p.adjust(ego.t$pvalue, method = "BH")
      # GeneRatio = genes of interest in the gene set / total genes of interest.
      #             the fraction of differentially expressed genes found in the gene set.
      # BgRatio = size of the geneset / size of all of the unique genes in the collection of genesets

      if (ego.t$qvalue %>% is.na() %>% which() %>% length() == 0 && ego.t$Ontology %>% length() != 0) {
        p1 <-
          dotplot(
            ego,
            color = "qvalue",
            showCategory = length(ego.t$ID %>% unique()),
            font.size = 12,
            title = paste(name, go, sep = "_")
          ) +
          scale_y_discrete(
            labels = function(x) {
              str_wrap(x, width = 40)
            }
          )

        if (length(ego.t$ID %>% unique()) > 0) {
          p2 <- dotplot(
            ego,
            color = "qvalue",
            showCategory = 20,
            font.size = 12,
            title = paste(name, go, sep = "_")
          ) +
            scale_y_discrete(
              labels = function(x) {
                str_wrap(x, width = 40)
              }
            )
        } else {
          print(paste("no GO term enriched in", name, go))

          p2 <-
            dotplot(
              ego,
              color = "qvalue",
              showCategory = length(ego.t$ID %>% unique()),
              font.size = 12,
              title = paste(name, go, sep = "_")
            ) +
            scale_y_discrete(
              labels = function(x) {
                str_wrap(x, width = 40)
              }
            )
        }
      } else {
        print(paste(
          "there is NA in qvalue or No gene set have enough size ,
        that may due to input not enough gene ids",
          name,
          go
        ))
        p1 <- NULL
        p2 <- NULL
      }


      if (length(ego.t$ID %>% unique()) < 20) {
        zoom <- 1
      } else {
        zoom <- length(ego.t$ID %>% unique()) / 20
      }

      print(paste("estimate ploting scale factor:", zoom))


      # create output path
      dir.create(file.path(datapath, "output", "plot"), recursive = T, mode = "0770")
      dir.create(file.path(datapath, "output", "GOterm"), recursive = T, mode = "0770")
      print(paste("output:", file.path(datapath, "output", "plot"), name, go))
      print(paste("output:", file.path(datapath, "output", "GOterm"), name, go))

      ggsave(
        filename = paste(name, go, "dotplot_all_terms.pdf", sep = "_"),
        plot = p1,
        path = file.path(datapath, "output", "plot"),
        width = 10,
        height = 12 * zoom,
        limitsize = FALSE
      )
      ggsave(
        filename = paste(name, go, "dotplot_top20_qval_terms.pdf", sep = "_"),
        plot = p2,
        path = file.path(datapath, "output", "plot"),
        width = 10,
        height = 12
      )

      write.table(ego.t, file.path(datapath, "output", "GOterm", paste(name, go, "clusterprofile.enrichment.txt", sep = "_")), sep = "\t", quote = F, row.names = F, col.names = T, eol = "\n")
    }
  }
}
