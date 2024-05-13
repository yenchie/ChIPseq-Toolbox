#########################################
# GOseq.hyper.FDR.DEG.fun,v.2.R
# revised date: 2024.05.06
# reviser: YCW
###
# revised date: 2024.05.06
# - revised duplicated line to a function: GO.seq.ego(pwf, GO_DB, termfile, q.cut.off= 0.05)
# - added DAG plot
# - revise output table colnames
#   [1] "category" ==> "ID" | "pval" ==> "pvalue" | "numDEInCat" ==> "geneHits" | "numInCat"==> "pathGenes" | "qval"==>"FDR" |
#   [6] "Description"  "Ontology" "entries"
#
## 2024.02.06
# - print "no any GO terms were enriched" if no results were found.
#########################################
# Got this version from Min on 20171108
# I Ran all command lines on 20171108
# I got one text file and three png tree files without problem.
# Three png files are trees of major GO categories (BP, MF and CC).
#
# USAGE:
# When runnig this Rscript:
# (1) five files need to be put on the working space:
# GOterm.txt
# Gene_Model.lengths.txt
# GO_Molecular_Function.txt
# GO_Cellular_Component.txt
# GO_Biological_Process.txt
#
# (2) input folder: files with gene ID
#
########################################
#
# This is for subregions specific expression genes GO terms enrichment.
# three type of GO terms test independent.
# Do hypergeometric and FDR adjust.
########################################
# core function:
# GOseq.ego.fun.v.1.R
########################################

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("goseq", quietly = TRUE)) {
  BiocManager::install("goseq")
}
########################################

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
print(source.path)
print(file.path(source.path %>% list.files("GOseq.ego.fun.v.1.R", full.names = T)))
source(file.path(source.path %>% list.files("GOseq.ego.fun.v.1.R", full.names = T)))

GOseq.hyper.FDR.DEG.fun <- function(datapath) {
  require(goseq)
  require(dplyr)
  require(stringr)
  require(tidyr)
  require(GO.db)
  require(ggplot2)


  # read all the table.
  featurepath <- "/bcst/JYL/JYL_qnap/WHHsieh/GO_enrichment/GO_input"
  all.genes <- read.table(file.path(featurepath, "Ro_v3.genelength.txt"), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  head(all.genes)
  nrow(all.genes)
  length(all.genes$V1 %>% unique())

  GO_BP <- read.table(file.path(featurepath, "Ro_protein_v3_Go47229_BiologicalProcess.txt"), sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(Ontology = "BP")
  head(GO_BP)
  nrow(GO_BP)
  length(GO_BP$V1 %>% unique())
  length(GO_BP$V2 %>% unique())
  GO_MF <- read.table(file.path(featurepath, "Ro_protein_v3_Go47229_MolecularFunction.txt"), sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(Ontology = "MF")
  head(GO_MF)
  nrow(GO_MF)
  length(GO_MF$V1 %>% unique())
  length(GO_MF$V2 %>% unique())
  GO_CC <- read.table(file.path(featurepath, "Ro_protein_v3_Go47229_CellularComponent.txt"), sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(Ontology = "CC")
  head(GO_CC)
  nrow(GO_CC)
  length(GO_CC$V1 %>% unique())
  length(GO_CC$V2 %>% unique())

  GO_ALL <- GO_BP %>%
    bind_rows(GO_MF) %>%
    bind_rows(GO_CC)
  colnames(GO_ALL) <- c("geneID", "ID", "Ontology")
  head(GO_ALL)
  nrow(GO_ALL)
  length(GO_ALL$geneID %>% unique())
  length(GO_ALL$ID %>% unique()) # 5,717

  GO_ALL %>%
    group_by(ID) %>%
    dplyr::count(Ontology)
  length(GO_BP$geneID %>% unique()) + length(GO_MF$geneID %>% unique()) + length(GO_CC$geneID %>% unique())
  length(GO_BP$ID %>% unique()) + length(GO_MF$ID %>% unique()) + length(GO_CC$ID %>% unique())


  type2Ont <- data.frame(
    V3 = c("P", "F", "C"),
    Ontology = c("BP", "MF", "CC")
  )
  termfile <- read.table(file.path(featurepath, "GOterm47229.txt"), sep = "\t", stringsAsFactors = FALSE, header = FALSE, fill = TRUE, quote = "") %>% left_join(type2Ont)
  colnames(termfile) <- c("ID", "Description", "Ont", "note", "Ontology")
  termfile %>% head()
  termfile %>% nrow()
  length(termfile$V1 %>% unique()) # 47,229
  which(!GO_ALL$ID %in% termfile$ID) %>% length()
  which(!(GO_ALL$geneID %>% unique()) %in% (all.genes$V1 %>% unique())) %>% length()

  setdiff(GO_ALL$geneID, all.genes$V1)
  setdiff(all.genes$V1, GO_ALL$geneID)
  print(paste("All annotated genes in GO term list", GO_ALL$geneID %>% unique() %>% length()))
  print(paste("All  GO term number", GO_ALL$ID %>% unique() %>% length()))

  print(paste("All annotated genes in BP GO term list", GO_BP$V1 %>% unique() %>% length()))
  print(paste("All BP GO term number", GO_BP$V2 %>% unique() %>% length()))

  print(paste("All annotated genes in MF GO term list", GO_MF$V1 %>% unique() %>% length()))
  print(paste("All MF GO term number", GO_MF$V2 %>% unique() %>% length()))

  print(paste("All annotated genes in CC GO term list", GO_CC$V1 %>% unique() %>% length()))
  print(paste("All CC GO term number", GO_CC$V2 %>% unique() %>% length()))

  # query.GO <- GO_ALL %>%
  #     left_join(termfile) %>%
  #     dplyr::select(-c(Ont, note)) %>%
  #     distinct() %>%
  #     filter(Description %>% is.na()) %>%
  #     dplyr::select(ID) %>%
  #     distinct()
  #   term <- list()
  #   for (go in query.GO$ID) {
  #     if (GOTERM[[go]] %>% is.null()) {
  #       print(paste("non description", go))
  #       cat("--------------------------------------\n")
  #     } else {
  #       term[[go]] <- GOTERM[[go]]
  #       print(GOTERM[[go]])
  #       cat("--------------------------------------\n")
  #     }
  #   }

  #   termfile.1 <- data.frame(
  #     GO = sapply(term, function(x) x@GOID),
  #     descritpion = sapply(term, function(x) x@Term)
  #   )
  #   termfile <- termfile %>% bind_rows(termfile.1)
  #   GO_ALL <- GO_ALL %>%
  #     left_join(termfile) %>%
  #     dplyr::select(-c(Ont, note))

  ###  genes -----
  assayed.genes <- all.genes$V1
  gene.length <- as.integer(all.genes$V2)
  names(gene.length) <- assayed.genes

  ###################################################################
  # List all the gene file.
  datapath <- datapath
  filelist <- datapath %>%
    list.dirs() %>%
    list.files(pattern = "id$", full.names = T)
  print("read input gene list files:")
  print(filelist)

  # Start analysis. ---------
  dir.create(file.path(datapath, "output", "GOterm"), recursive = T, mode = "0770")
  dir.create(file.path(datapath, "output", "plot"), recursive = T, mode = "0770")

  a <- 0
  for (i in filelist) {
    name <- NA
    a <- a + 1
    print(a)
    de.genes <- read.table(i, sep = "\t", stringsAsFactors = FALSE, header = FALSE) %>% distinct()
    nrow(de.genes)

    gene.vector <- as.integer(assayed.genes %in% de.genes$V1)
    names(gene.vector) <- assayed.genes
    pwf <- nullp(gene.vector, bias.data = gene.length, plot.fit = F)
    head(pwf)
    # plotPWF(pwf,binsize=200)

    ## GO terms enrichment of Biological Process ------

    q.cut.of <- 0.05
    BP <- GO.seq.ego(pwf, GO_BP, termfile = termfile, 0.05)
    outpath <- file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", paste0("_DAG.BP", q.cut.of, ".pdf"))))
    DAG.GOseq.fun(GOseq.result = BP$ego, GO_DB = GO_BP, GO.ontology = "BP", termfile, GOI.list = de.genes$V1, q.cut.off = 0.05, outpath = outpath)

    CC <- GO.seq.ego(pwf, GO_CC, termfile = termfile, 0.05)
    outpath <- file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", paste0("_DAG.CC", q.cut.of, ".pdf"))))
    DAG.GOseq.fun(CC$ego, GO_CC, "CC", termfile, de.genes$V1, 0.05, outpath)

    MF <- GO.seq.ego(pwf, GO_MF, termfile = termfile, 0.05)
    outpath <- file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", paste0("_DAG.MF", q.cut.of, ".pdf"))))
    DAG.GOseq.fun(MF$ego, GO_MF, "MF", termfile, de.genes$V1, 0.05, outpath)

    # mix table into one.
    enriched.go.BP <- BP$fin.result
    head(enriched.go.BP) %>% print()
    enriched.go.CC <- CC$fin.result
    head(enriched.go.CC) %>% print()
    enriched.go.MF <- MF$fin.result
    head(enriched.go.MF) %>% print()


    enriched.go.annot <- NULL
    enriched.go.annot <- rbind(enriched.go.BP, enriched.go.CC, enriched.go.MF)
    head(enriched.go.annot)

    enriched.go.annot <- enriched.go.annot %>%
      dplyr::select(
        "ID",
        "pvalue",
        "geneHits",
        "pathGenes",
        "FDR",
        "term", "Description",
        "Ontology",
        "entries"
      )

    name <- file.path(datapath, "output", "GOterm", paste0(i %>% basename() %>% str_replace("\\.id", "_GOseq.enrichment.txt")))
    print("output enriched result:")
    print(name)
    summary.number <- dplyr::count(enriched.go.annot, Ontology) %>% dplyr::rename(term.number = n)
    print("total  enriched terms:")
    print(summary.number)
    print(enriched.go.annot$Ontology %>% unique())
    write.table(enriched.go.annot, name, sep = "\t", quote = FALSE, row.names = F)
    print(paste("file:", i, "GO term enrichment done"))


    ## plotting -----
    p1 <- NULL
    d1 <- NULL
    if (!is.na(name)) {
      d1 <-
        read.delim(
          name,
          stringsAsFactors = F
        ) %>% arrange(FDR)
      colnames(d1)
      print(d1$Ontology %>% unique())
      p1 <- d1 %>%
        na.omit() %>%
        mutate(log10.qval = -log10(FDR)) %>%
        filter(!log10.qval %>% is.infinite()) %>%
        group_by(Ontology) %>%
        top_n(10, wt = log10.qval) %>%
        ungroup() %>%
        mutate(
          Description = Description %>% str_wrap(width = 40),
          GeneRatio = geneHits * 100 / pathGenes
        ) %>%
        ggplot(aes(
          x = GeneRatio,
          y = reorder(Description, GeneRatio),
          colour = log10.qval %>% round(1),
          size = geneHits
        )) +
        geom_point() +
        scale_size(range = c(6, 16)) +
        scale_color_gradient(low = "blue", high = "red") +
        expand_limits(x = 0) +
        labs(
          x = "GeneRatio (%) ",
          # the fraction of genes-of-interest found in the gene set.
          y = "GO term",
          colour = "-log10.qval",
          size = "Count"
        ) +
        theme(axis.text.y = element_text(face = "bold", size = 12.5, vjust = 0.5)) +
        facet_wrap(~Ontology, scales = "free_y", ncol = 2)

      print(paste(
        "output:",
        file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", "_dotplot_top10_terms.pdf")))
      ))

      ggsave(
        filename = paste0(i %>% basename() %>% str_replace("\\.id", "_dotplot_top10_terms.pdf")),
        plot = p1,
        path = file.path(datapath, "output", "plot"),
        width = 16,
        height = 16
      )

      # complete go term plot
      d2 <- NULL
      for (ont in d1$Ontology %>% unique()) {
        d2 <- d1 %>% # na.omit()%>%
          filter(Ontology == ont) %>%
          mutate(log10.qval = -log10(FDR)) %>%
          filter(!log10.qval %>% is.infinite()) %>%
          mutate(
            Description = Description %>% str_wrap(width = 40),
            GeneRatio = geneHits * 100 / pathGenes
          )

        p2 <- d2 %>%
          ggplot(aes(
            x = GeneRatio,
            y = reorder(Description, GeneRatio),
            colour = log10.qval %>% round(1),
            size = geneHits
          )) +
          geom_point() +
          scale_size(range = c(6, 16)) +
          scale_color_gradient(low = "blue", high = "red") +
          expand_limits(x = 0) +
          labs(
            title = paste("GO:", ont),
            x = "GeneRatio (%) ",
            # the fraction of genes-of-interest found in the gene set.
            y = "GO term",
            colour = "-log10.qval",
            size = "Count"
          ) +
          theme(axis.text.y = element_text(
            face = "bold", size = 12.5, vjust =
              0.5
          ))

        print(paste(
          "output:",
          file.path(datapath, "output", "plot", paste0(ont, "_", i %>% basename() %>% str_replace("\\.id", "_dotplot_all_terms.pdf")))
        ))

        if (length(d2$category %>% unique()) < 20) {
          zoom <- 1
        } else {
          zoom <- length(d2$category %>% unique()) / 20
        }
        ggsave(
          filename = paste0(ont, "_", i %>% basename() %>% str_replace("\\.id", "_dotplot_all_terms.pdf")),
          plot = p2,
          path = file.path(datapath, "output", "plot"),
          width = 10,
          height = 12 * zoom,
          limitsize = FALSE
        )
      }
    } else {
      print(paste("No GO terms were enriched in the file:", i))
    }
  }

  Sys.sleep(5)
}
