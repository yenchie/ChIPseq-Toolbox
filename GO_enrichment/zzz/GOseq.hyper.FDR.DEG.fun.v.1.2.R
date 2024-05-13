#########################################
# GOseq.hyper.FDR.DEG.fun,v.1.2.R
# revised date: 20230914
# reviser: YCW
###
# revised date: 2024.04.17
# v.1.2. add DAG plot
# 2024.02.06
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
# Gmax_189_Gene_Model.lengths.txt
# Gmax.soybase.GO_Molecular_Function.txt
# Gmax.soybase.GO_Cellular_Component.txt
# Gmax.soybase.GO_Biological_Process.txt
#
# (2) input file: a file with gene ID
#
# (3) Need to change the InputFileName in this command line, but keep the "":
# filelist=list.files(pattern="InputFileName")
#
# (4) When running the script, need to put all files in the same directory.
#
########################################
#
# This is for subregions specific expression genes GO terms enrichment.
# three type of GO terms test independent.
# Do hypergeometric and FDR adjust.
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("goseq", quietly = TRUE)) {
  BiocManager::install("goseq")
}
########################################
source("/homes/yenching/YCWang/Rscript/1_Function/CMD_GO_USE_THIS_ONE/jamenrich-enrichdf2er.R")
source("/homes/yenching/YCWang/Rscript/1_Function/CMD_GO_USE_THIS_ONE/jamenrich-import.R")
GOseq.hyper.FDR.DEG.fun <- function(datapath) {
  rm(list = ls())
  require(goseq)
  require(dplyr)
  require(stringr)
  require(tidyr)
  require(GO.db)
  require(ggplot2)

  # read all the table.
  featurepath <- "/bcst/JYL/JYL_qnap/Project/Gm/ChIPseq_analysis/LOG/Gmax189/Methylation/Analysis/Other_Analysis/Bivalent_H3K4me3_H3K27me3/GO_enrichment/input/feature"
  all.genes <- read.table(file.path(featurepath, "Gmax_189_Gene_Model.lengths.txt"), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  GO_BP <- read.table(file.path(featurepath, "Gmax.soybase.GO_Biological_Process.txt"), sep = "\t", stringsAsFactors = FALSE)
  GO_MF <- read.table(file.path(featurepath, "Gmax.soybase.GO_Molecular_Function.txt"), sep = "\t", stringsAsFactors = FALSE)
  GO_CC <- read.table(file.path(featurepath, "Gmax.soybase.GO_Cellular_Component.txt"), sep = "\t", stringsAsFactors = FALSE)
  type2Ont <- data.frame(
    V3 = c("P", "F", "C"),
    Ontology = c("BP", "MF", "CC")
  )
  termfile <- read.table(file.path(featurepath, "GOterm.txt"), sep = "\t", stringsAsFactors = FALSE, header = FALSE, fill = TRUE, quote = "") %>% left_join(type2Ont)

  # x<- GO_BP%>%bind_rows(GO_MF)%>%bind_rows(GO_CC)
  # n_distinct(x$V2)
  # non_ant<- x%>%filter(!V2%in%termfile$V1)
  # length(unique(non_ant$V2))

  #
  assayed.genes <- all.genes$V1
  gene.length <- as.integer(all.genes$V2)
  names(gene.length) <- assayed.genes

  ###################################################################
  ##################################################################
  ### change input file name HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # List all the gene file.
  datapath <- datapath
  filelist <- datapath %>%
    list.dirs() %>%
    list.files(pattern = "id$", full.names = T)
  print("read input gene list files:")
  print(filelist)

  ############################################################
  ##################################################################
  # Just need one input as above line.
  # The following are added by min, and then Min ran it one by one.
  #
  ##### START of lines modified by Min ##########################
  # filelist=list.files(pattern="LT0.05.DMVgene.TF.ID.txt")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.ID.txt")
  # filelist=list.files(pattern="LT0.05.notInDMV.ID.txt")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.GLOB-")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.HRT-")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.COT-")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.EM-")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.GOseq.enrichment.geneNotInBP.txt")
  # filelist=list.files(pattern="LT0.05.DMVgene.5foldUP-from-edgeR-Jungim.GOseq.enrichment.geneNotInBP-MF.txt")
  # filelist=list.files(pattern="Gmax_v1.1_189_geneID.txt")
  # filelist=list.files(pattern="seedonly.LT0.05.DMVgene.ID.txt")

  # filelist=list.files(pattern="random-set1.txt")
  # filelist=list.files(pattern="random-set2.txt")
  # filelist=list.files(pattern="random-set3.txt")
  # filelist=list.files(pattern="random-set4.txt")
  # filelist=list.files(pattern="random-set5.txt")
  ##### END of lines modified by Min ##########################

  # Start analysis.
  # dir.create(file.path(datapath, "output"), recursive = T, mode = "0770")
  dir.create(file.path(datapath, "output", "GOterm"), recursive = T, mode = "0770")
  dir.create(file.path(datapath, "output", "plot"), recursive = T, mode = "0770")

  a <- 0
  for (i in filelist) {
    name <- NA
    a <- a + 1
    print(a)
    de.genes <- read.delim(i, sep = "\t", stringsAsFactors = FALSE, header = FALSE)

    gene.vector <- as.integer(assayed.genes %in% de.genes$V1)
    names(gene.vector) <- assayed.genes

    pwf <- nullp(gene.vector, bias.data = gene.length, plot.fit = F)
    # plotPWF(pwf,binsize=200)
    # GO terms enrichment of Biological Process
    go.wall.BP <- goseq(pwf, gene2cat = GO_BP, method = "Hypergeometric", use_genes_without_cat = T) # Wallenius
    go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue, method = "BH")
    # # use_genes_without_cat test ----
    # go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric",use_genes_without_cat=T)
    # go.wall.BP.pwf <- goseq(pwf,gene2cat=GO_BP, method="Hypergeometric", use_genes_without_cat=F)
    # go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
    # go.wall.BP.pwf$qval <- p.adjust(go.wall.BP.pwf$over_represented_pvalue,method="BH")
    # plot(
    #   log10(go.wall.BP.pwf[, 8]),
    #   log10(go.wall.BP[match(go.wall.BP.pwf[, 1], go.wall.BP[, 1]), 8]),
    #   xlab = "use_genes_without_cat=F, log10(qval)",
    #   ylab = "use_genes_without_cat=T, log10(qval)",
    #   xlim = c(-3, 0),
    #   ylim = c(-3, 0)
    # )
    # abline(0,1,col=3,lty=2)
    # abline(h=-2, col="red",lty=4)
    # abline(v=-2, col="red",lty=4)
    # #x<-log10(0.01)
    #
    #
    # # different method test ----
    # go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric",use_genes_without_cat=T)
    # go.wall.BP.pwf <- goseq(pwf,gene2cat=GO_BP, method="Wallenius", use_genes_without_cat=T)
    # go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
    # go.wall.BP.pwf$qval <- p.adjust(go.wall.BP.pwf$over_represented_pvalue,method="BH")
    # plot(
    #   log10(go.wall.BP.pwf[, 8]),
    #   log10(go.wall.BP[match(go.wall.BP.pwf[, 1], go.wall.BP[, 1]), 8]),
    #   xlab = "Wallenius, log10(qval)",
    #   ylab = "Hypergeometric, log10(qval)",
    #   xlim = c(-3, 0),
    #   ylim = c(-3, 0)
    # )
    # abline(0,1,col=3,lty=2)
    # abline(h=-2, col="red",lty=4)
    # abline(v=-2, col="red",lty=4)
    # #x<-log10(0.01)
    #
    # # different method, different bg test ----
    # go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric",use_genes_without_cat=T)
    # go.wall.BP.pwf <- goseq(pwf,gene2cat=GO_BP, method="Wallenius", use_genes_without_cat=F)
    # go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
    # go.wall.BP.pwf$qval <- p.adjust(go.wall.BP.pwf$over_represented_pvalue,method="BH")
    # plot(
    #   log10(go.wall.BP.pwf[, 8]),
    #   log10(go.wall.BP[match(go.wall.BP.pwf[, 1], go.wall.BP[, 1]), 8]),
    #   xlab = "Wallenius, use_genes_without_cat=F, log10(qval)",
    #   ylab = "Hypergeometric, use_genes_without_cat=T, log10(qval)",
    #   xlim = c(-3, 0),
    #   ylim = c(-3, 0)
    # )
    # abline(0,1,col=3,lty=2)
    # abline(h=-2, col="red",lty=4)
    # abline(v=-2, col="red",lty=4)
    # #x<-log10(0.01)
    ## -----



    # enriched.GO=go.wall.BP$category[p.adjust(go.wall.BP$over_represented_pvalue,method="BH")<.05]
    # print(enriched.GO)
    enriched.go.BP <- subset(go.wall.BP, go.wall.BP$qval < 0.05)

    # if have Go terms qval lower than 0.05, draw a tree.
    if (nrow(enriched.go.BP) > 0) {
      name <- file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", "_GOseq.enrichment.DAG.BP.pdf")))
      pdf(name)
      enriched.go.BP.DAG <- enriched.go.BP %>%
        dplyr::rename(ID = category, pvalue = over_represented_pvalue, Count = numDEInCat) %>%
        mutate(Description = ID, FDR = qval, GeneRatio = paste0(Count, "/", numInCat)) %>%
        na.omit()
      enriched.go.BP_res <- enrichDF2enrichResult(
        enriched.go.BP.DAG,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        keyColname = "ID",
        pathGenes = "numInCat",
        # geneColname = "gene",
        geneHits = "Count",
        # geneRatioColname = "GeneRatio",
        geneDelim = "[,/ ]+",
        geneSep = ",",
        pvalueColname = "FDR",
        descriptionColname = "ID",
        msigdbGmtT = NULL,
        verbose = FALSE
      )
      enriched.go.BP_res@ontology <- "BP"
      enriched.go.BP_res %>% plotGOgraph()
      dev.off()

      goIDs <- enriched.go.BP[, 1]
      pvalue <- enriched.go.BP[, 8]
      # name=paste0(substr(i,1,nchar(i)-7),"GOseq.enrichment.BP.png")
      # pngRes <- getAmigoTree(goIDs=goIDs, pvalues=pvalue, pcolors=c("white","red"), psplit=c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001), filename=name)

      enriched.go.wall.annot.BP <- merge(enriched.go.BP, termfile, by.x = "category", by.y = "V1", all.x = TRUE)

      # for table
      Data1 <- NULL
      for (j in c(1:nrow(enriched.go.wall.annot.BP))) {
        referencesGenes <- GO_BP[GO_BP[, 2] == enriched.go.wall.annot.BP[j, 1], 1]
        NumberOfReferencesGenes <- length(referencesGenes)
        Genes <- intersect(referencesGenes, de.genes[, 1])
        entries <- NULL
        for (k in c(1:length(Genes))) {
          entries <- paste(entries, Genes[k], sep = ",")
        }
        entries <- substr(entries, 2, nchar(entries))
        NumberOfGenes <- length(Genes)
        Data2 <- cbind(NumberOfGenes, NumberOfReferencesGenes, entries)
        Data1 <- rbind(Data1, Data2)
      }
      enriched.go.wall.annot.BP <- cbind(enriched.go.wall.annot.BP, Data1)
      enriched.go.wall.annot.BP <- enriched.go.wall.annot.BP[order(enriched.go.wall.annot.BP$qval), ]
    }
    if (nrow(enriched.go.BP) == 0) {
      enriched.go.wall.annot.BP <- NULL
    }

    # GO terms enrichment of molecular function.
    go.wall.MF <- goseq(pwf, gene2cat = GO_MF, method = "Hypergeometric", use_genes_without_cat = T)
    go.wall.MF$qval <- p.adjust(go.wall.MF$over_represented_pvalue, method = "BH")
    # enriched.GO=go.wall.MF$category[p.adjust(go.wall.MF$over_represented_pvalue,method="BH")<.05]
    # print(enriched.GO)
    enriched.go.MF <- subset(go.wall.MF, go.wall.MF$qval < 0.05)

    # if have Go terms qval lower than 0.05, draw a tree.
    if (nrow(enriched.go.MF) > 0) {
      name <- file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", "_GOseq.enrichment.DAG.MF.pdf")))
      pdf(name)
      enriched.go.MF.DAG <- enriched.go.MF %>%
        dplyr::rename(ID = category, pvalue = over_represented_pvalue, Count = numDEInCat) %>%
        mutate(Description = ID, FDR = qval, GeneRatio = paste0(Count, "/", numInCat)) %>%
        na.omit()
      enriched.go.MF_res <- enrichDF2enrichResult(enriched.go.MF.DAG,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        keyColname = "ID",
        pathGenes = "numInCat",
        # geneColname = "gene",
        geneHits = "Count",
        # geneRatioColname = "GeneRatio",
        geneDelim = "[,/ ]+",
        geneSep = ",",
        pvalueColname = "FDR",
        descriptionColname = "ID",
        msigdbGmtT = NULL,
        verbose = FALSE
      )
      enriched.go.MF_res@ontology <- "MF"
      enriched.go.MF_res %>% plotGOgraph()
      dev.off()

      goIDs <- enriched.go.MF[, 1]
      pvalues <- enriched.go.MF[, 8]
      name <- paste(substr(i, 1, nchar(i) - 7), "GOseq.enrichment.MF", sep = ".")
      # pngRes <- getAmigoTree(goIDs=goIDs,pvalues=pvalues,pcolors=c("white","red"),
      #                        psplit=c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001),filename=name)

      enriched.go.wall.annot.MF <- merge(enriched.go.MF, termfile, by.x = "category", by.y = "V1", all.x = TRUE)


      # for table
      Data1 <- NULL
      for (j in c(1:nrow(enriched.go.wall.annot.MF))) {
        referencesGenes <- GO_MF[GO_MF[, 2] == enriched.go.wall.annot.MF[j, 1], 1]
        NumberOfReferencesGenes <- length(referencesGenes)
        Genes <- intersect(referencesGenes, de.genes[, 1])
        entries <- NULL
        for (k in c(1:length(Genes))) {
          entries <- paste(entries, Genes[k], sep = ",")
        }
        entries <- substr(entries, 2, nchar(entries))
        NumberOfGenes <- length(Genes)
        Data2 <- cbind(NumberOfGenes, NumberOfReferencesGenes, entries)
        Data1 <- rbind(Data1, Data2)
      }
      enriched.go.wall.annot.MF <- cbind(enriched.go.wall.annot.MF, Data1)
      enriched.go.wall.annot.MF <- enriched.go.wall.annot.MF[order(enriched.go.wall.annot.MF$qval), ]
    }
    if (nrow(enriched.go.MF) == 0) {
      enriched.go.wall.annot.MF <- NULL
    }

    # GO terms enrichment of Cellular component.
    go.wall.CC <- goseq(pwf, gene2cat = GO_CC, method = "Hypergeometric", use_genes_without_cat = T)
    go.wall.CC$qval <- p.adjust(go.wall.CC$over_represented_pvalue, method = "BH")
    # enriched.GO=go.wall.CC$category[p.adjust(go.wall.CC$over_represented_pvalue,method="BH")<.05]
    # print(enriched.GO)
    enriched.go.CC <- subset(go.wall.CC, go.wall.CC$qval < 0.05)

    # if have Go terms qval lower than 0.05, draw a tree.
    if (nrow(enriched.go.CC) > 0) {
      name <- file.path(datapath, "output", "plot", paste0(i %>% basename() %>% str_replace("\\.id", "_GOseq.enrichment.DAG.MF.pdf")))
      pdf(name)
      enriched.go.CC.DAG <- enriched.go.CC %>%
        dplyr::rename(ID = category, pvalue = over_represented_pvalue, Count = numDEInCat) %>%
        mutate(Description = ID, FDR = qval, GeneRatio = paste0(Count, "/", numInCat)) %>%
        na.omit()
      enriched.go.CC_res <- enrichDF2enrichResult(enriched.go.CC.DAG,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        keyColname = "ID",
        pathGenes = "numInCat",
        # geneColname = "gene",
        geneHits = "Count",
        # geneRatioColname = "GeneRatio",
        geneDelim = "[,/ ]+",
        geneSep = ",",
        pvalueColname = "FDR",
        descriptionColname = "ID",
        msigdbGmtT = NULL,
        verbose = FALSE
      )
      enriched.go.CC_res@ontology <- "CC"
      enriched.go.CC_res %>% plotGOgraph()
      dev.off()

      goIDs <- enriched.go.CC[, 1]
      pvalues <- enriched.go.CC[, 8]
      name <- paste(substr(i, 1, nchar(i) - 7), "GOseq.enrichment.CC", sep = ".")
      # pngRes <- getAmigoTree(goIDs=goIDs,pvalues=pvalues,pcolors=c("white","red"),
      #                        psplit=c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001),filename=name)

      enriched.go.wall.annot.CC <- merge(enriched.go.CC, termfile, by.x = "category", by.y = "V1", all.x = TRUE)


      # for table
      Data1 <- NULL
      for (j in c(1:nrow(enriched.go.wall.annot.CC))) {
        referencesGenes <- GO_CC[GO_CC[, 2] == enriched.go.wall.annot.CC[j, 1], 1]
        NumberOfReferencesGenes <- length(referencesGenes)
        Genes <- intersect(referencesGenes, de.genes[, 1])
        entries <- NULL
        for (k in c(1:length(Genes))) {
          entries <- paste(entries, Genes[k], sep = ",")
        }
        entries <- substr(entries, 2, nchar(entries))
        NumberOfGenes <- length(Genes)
        Data2 <- cbind(NumberOfGenes, NumberOfReferencesGenes, entries)
        Data1 <- rbind(Data1, Data2)
      }
      enriched.go.wall.annot.CC <- cbind(enriched.go.wall.annot.CC, Data1)
      enriched.go.wall.annot.CC <- enriched.go.wall.annot.CC[order(enriched.go.wall.annot.CC$qval), ]
    }
    if (nrow(enriched.go.CC) == 0) {
      enriched.go.wall.annot.CC <- NULL
    }

    # mix table into one.
    if (nrow(enriched.go.BP) + nrow(enriched.go.CC) + nrow(enriched.go.MF) > 0) {
      enriched.go.wall.annot <- NULL
      enriched.go.wall.annot <- rbind(enriched.go.wall.annot.BP, enriched.go.wall.annot.MF, enriched.go.wall.annot.CC)

      enriched.go.wall.annot <- enriched.go.wall.annot %>%
        dplyr::select(
          "category",
          "over_represented_pvalue",
          "numDEInCat",
          "numInCat",
          "qval",
          "V2",
          "Ontology",
          "entries"
        )

      colnames(enriched.go.wall.annot)[2] <- "pval"
      colnames(enriched.go.wall.annot)[6] <- "Description"
      # colnames(enriched.go.wall.annot)[7]="type"
      name <- file.path(datapath, "output", "GOterm", paste0(i %>% basename() %>% str_replace("\\.id", "_GOseq.enrichment.txt")))
      print("output enriched result:")
      print(name)
      summary.number <- dplyr::count(enriched.go.wall.annot, Ontology) %>% dplyr::rename(term.number = n)
      print("total  enriched terms:")
      print(summary.number)

      write.table(enriched.go.wall.annot, name, sep = "\t", quote = FALSE, row.names = F)
      print(paste("file:", i, "GO term enrichment done"))
    }
    p1 <- NULL
    d1 <- NULL
    if (!is.na(name)) {
      d1 <-
        read.delim(
          name,
          stringsAsFactors = F
        ) %>% arrange(qval)
      colnames(d1)
      p1 <- d1 %>%
        na.omit() %>%
        mutate(log10.qval = -log10(qval)) %>%
        filter(!log10.qval %>% is.infinite()) %>%
        group_by(Ontology) %>%
        top_n(10, wt = log10.qval) %>%
        ungroup() %>%
        mutate(
          Description = Description %>% str_wrap(width = 40),
          GeneRatio = numDEInCat * 100 / numInCat
        ) %>%
        ggplot(aes(
          x = GeneRatio,
          y = reorder(Description, GeneRatio),
          colour = log10.qval %>% round(1),
          size = numDEInCat
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
          mutate(log10.qval = -log10(qval)) %>%
          filter(!log10.qval %>% is.infinite()) %>%
          mutate(
            Description = Description %>% str_wrap(width = 40),
            GeneRatio = numDEInCat * 100 / numInCat
          )

        p2 <- d2 %>%
          ggplot(aes(
            x = GeneRatio,
            y = reorder(Description, GeneRatio),
            colour = log10.qval %>% round(1),
            size = numDEInCat
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
}
