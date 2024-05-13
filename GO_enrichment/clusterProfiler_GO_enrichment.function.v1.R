# clusterProfiler_GO_enrichment.function.v1.R

# revised date: 2024.04.16
# add DAG plot
# 2024.02.06
# - revised if didn'tenrichment any result, print "no any GO terms be enrichment"
GO_enricher <- function(data.path, outpath, cutoff) {
  require(clusterProfiler)
  require(topGO)
  require(enrichplot)
  require(dplyr)
  require(stringr)
  require(tidyr)
  require(GO.db)
  require(ggplot2)



  # read all the table.
  featurepath <- "/bcst/JYL/JYL_qnap/Project/Gm/ChIPseq_analysis/LOG/Gmax189/Methylation/Analysis/Other_Analysis/Bivalent_H3K4me3_H3K27me3/GO_enrichment/input/feature"
  print(paste("read annotation files:", featurepath))

  # all.genes<-read.table(file.path(featurepath, "Gmax_189_Gene_Model.lengths.txt"), sep="\t",stringsAsFactors=FALSE, header=FALSE)
  GO_BP <- read.table(file.path(featurepath, "Gmax.soybase.GO_Biological_Process.txt"), sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::rename(GO = V2, gene = V1) %>%
    mutate(Ontology = "GO_BP")
  GO_MF <- read.table(file.path(featurepath, "Gmax.soybase.GO_Molecular_Function.txt"), sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::rename(GO = V2, gene = V1) %>%
    mutate(Ontology = "GO_MF")
  GO_CC <- read.table(file.path(featurepath, "Gmax.soybase.GO_Cellular_Component.txt"), sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::rename(GO = V2, gene = V1) %>%
    mutate(Ontology = "GO_CC")
  termfile <- read.table(file.path(featurepath, "GOterm.txt"), sep = "\t", stringsAsFactors = FALSE, header = FALSE, fill = TRUE, quote = "") %>% dplyr::rename(GO = V1, descritpion = V2, taq = V3)

  GO.all <- bind_rows(GO_BP, GO_MF) %>% bind_rows(GO_CC)
  GO.all$Ontology <- GO.all$Ontology %>% as.character()
  unique(GO.all$GO) %>% length()
  unique(termfile$GO) %>% length()
  # query GO term description


  query.GO <- GO.all %>%
    left_join(termfile) %>%
    dplyr::select(-c(taq, V4)) %>%
    filter(descritpion %>% is.na()) %>%
    dplyr::select(GO) %>%
    distinct()
  term <- list()
  for (go in query.GO$GO) {
    if (GOTERM[[go]] %>% is.null()) {
      print(paste("non description", go))
      cat("--------------------------------------\n")
    } else {
      term[[go]] <- GOTERM[[go]]
      print(GOTERM[[go]])
      cat("--------------------------------------\n")
    }
  }

  termfile.1 <- data.frame(
    GO = sapply(term, function(x) x@GOID),
    descritpion = sapply(term, function(x) x@Term)
  )
  termfile <- termfile %>% bind_rows(termfile.1)
  GO.all <- GO.all %>%
    left_join(termfile) %>%
    dplyr::select(-c(taq, V4))



  data.path <- file.path(data.path)
  print(paste("read id files in the folder:", data.path))
  print(paste("please notice that file name should be gene.id|TF.id in the end"))
  filelist <- list.files(file.path(data.path), pattern = "gene.id|TF.id", full.names = T)
  print(paste("files in", data.path, "imported"))
  # gene.id<- read.delim(filelist[1], stringsAsFactors = F, header = F)$V1%>%as.vector()

  # create output path
  dir.create(file.path(outpath, "output", "plot"), recursive = T, mode = "0770")
  dir.create(file.path(outpath, "output", "GOterm"), recursive = T, mode = "0770")

  for (file in 1:length(filelist)) {
    gene.id <- read.delim(filelist[file], stringsAsFactors = F, header = F)$V1 %>% as.vector()
    name <- filelist[file] %>%
      as.character() %>%
      basename() %>%
      str_remove("_id")

    print(unique(GO.all$Ontology))

    for (go in unique(GO.all$Ontology)) {
      print(paste("start:", name, go))

      TERM2GENE <- GO.all %>%
        dplyr::filter(Ontology == go)
      annt.size <- count(TERM2GENE, GO)
      print(paste("total input genes:", gene.id %>% length()))
      print(paste("gene in GOterm:", go, which(gene.id %in% TERM2GENE$gene) %>% length()))
      # count(TERM2GENE, GO)
      # n_distinct(TERM2GENE$gene)
      ego <- enricher(
        gene.id,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        # universe	background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
        universe = NULL,
        qvalueCutoff = cutoff,
        minGSSize = 1,
        maxGSSize = max(annt.size$n),
        TERM2GENE = TERM2GENE[c("GO", "gene")] # ,
        # TERM2NAME = TERM2GENE[c("GO", "descritpion")]
      )
      ego@ontology <- go %>% str_remove("GO_")
      pdf(file.path(outpath, "output", "plot", paste0(name, ".GO.enrichment.DAG.", go, ".pdf")))
      plotGOgraph(ego)
      dev.off()
      ego.t <- ego %>%
        as.data.frame() %>%
        mutate(Ontology = go)
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
      dir.create(file.path(data.path, "output", "plot"), recursive = T, mode = "0770")
      dir.create(file.path(data.path, "output", "GOterm"), recursive = T, mode = "0770")
      print(paste("output:", file.path(data.path, "output", "plot"), name, go))
      print(paste("output:", file.path(data.path, "output", "GOterm"), name, go))

      ggsave(
        filename = paste(name, go, "dotplot_all_terms.pdf", sep = "_"),
        plot = p1,
        path = file.path(outpath, "output", "plot"),
        width = 10,
        height = 12 * zoom,
        limitsize = FALSE
      )
      ggsave(
        filename = paste(name, go, "dotplot_top20_qval_terms.pdf", sep = "_"),
        plot = p2,
        path = file.path(outpath, "output", "plot"),
        width = 10,
        height = 12
      )

      write.table(ego.t, file.path(outpath, "output", "GOterm", paste(name, go, "clusterprofile.enrichment.txt", sep = "_")), sep = "\t", quote = F, row.names = F, col.names = T, eol = "\n")
    }
  }
}
