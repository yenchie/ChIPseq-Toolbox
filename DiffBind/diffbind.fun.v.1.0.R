## diffbind.fun.v.1.0.R
# Author: YCW
# Date: 2024/03/12
# Update: 2024.04.09
# Updata log:
# 1. 2024.04.09, version 1.0
#

require(DiffBind)
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)
require(ggpubr)
require(scales)
require(grid)
require(reshape2)
require(purrr)
require(tibble)

dba.plot <- function(data, title) {
  datac <- dba.count(data, peaks = NULL, score = DBA_SCORE_NORMALIZED)
  plot(
    datac,
    sub = title,
    attributes = c(DBA_CONDITION, DBA_REPLICATE),
    ColAttributes = DBA_CONDITION,
    distMethod = "spearman",
    margin = 5
  )

  heatmap.method <- list(
    hclust = c("ward.D2", "complete", "average", "mcquitty"),
    dist = c("euclidean", "manhattan", "correlation", "maximum")
  )
  for (hclust in 1:length(heatmap.method$hclust)) {
    for (dist in 1:length(heatmap.method$dist)) {
      hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
      dba.plotHeatmap(
        datac,
        attributes = c(DBA_CONDITION, DBA_REPLICATE),
        ColAttributes = DBA_CONDITION,
        correlations = FALSE,
        scale = "row",
        margin = 5,
        colScheme = hmap,
        main = paste("heatmap", heatmap.method$dist[dist], heatmap.method$hclust[hclust], sep = "_"),
        distMethod = heatmap.method$dist[dist],
        hclustfun = function(x) {
          hclust(x, method = heatmap.method$hclust[hclust])
        },
        sub = title
      )
    }
  }

  # PCA plot
  dba.plotPCA(datac, sub = title)
}
value.box.plot <- function(data, title, peakset = NULL) {
  if (!is.null(peakset)) {
    print("specific peakset in use")
  } else {
    print("NULL peakset in use")
    peakset <- NULL
  }
  datac <- dba.count(data, peaks = NULL, score = DBA_SCORE_NORMALIZED)
  nor.value <- dba.peakset(datac, bRetrieve = T, DataType = DBA_DATA_FRAME) %>%
    mutate(ID = paste(CHR, START, END, sep = "_"), .before = CHR) %>%
    dplyr::select(-c(CHR, START, END)) %>%
    pivot_longer(
      cols = 2:nrow(datac$samples) + 1,
      names_to = "SampleID",
      values_to = "nor.value"
    ) %>%
    mutate(
      stage = SampleID %>% str_extract("^\\w+"),
      SampleID = SampleID %>% str_remove("\\.H3.*")
    ) %>%
    mutate(stage = factor(stage,
      levels = c("cot", "em", "mm", "lm", "pd1", "dry", "sdlg27h", "sdlg8day")
    )) %>%
    arrange(stage)
  input <- nor.value %>%
    mutate(SampleID = SampleID %>%
      factor(levels = unique(SampleID)))
  print(head(input))

  outlier <- quantile(input$nor.value, c(0.25, 0.75))
  if (outlier[1] - IQR(input$nor.value) * 1.5 >= 0) {
    box.length <- c(
      outlier[1] - IQR(input$nor.value) * 1.5,
      outlier[2] + IQR(input$nor.value) * 1.5
    )
  } else {
    box.length <- c(
      0,
      outlier[2] + IQR(input$nor.value) * 1.5
    )
  }

  p1 <- input %>%
    ggplot() +
    geom_boxplot(aes(
      x = SampleID,
      y = nor.value
    ), outlier.shape = NA) +
    xlab("Sample") +
    ylab("value") +
    scale_y_continuous(limits = box.length) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = title)

  outlier <- quantile((input$nor.value + 1) %>% log2(), c(0.25, 0.75))
  if (outlier[1] - IQR((input$nor.value + 1) %>% log2()) * 1.5 >= 0) {
    box.length.log <- c(
      outlier[1] - IQR((input$nor.value + 1) %>% log2()) * 1.5,
      outlier[2] + IQR((input$nor.value + 1) %>% log2()) * 1.5
    )
  } else {
    box.length.log <- c(
      0,
      outlier[2] + IQR((input$nor.value + 1) %>% log2()) * 1.5
    )
  }
  p2 <- input %>%
    ggplot() +
    geom_boxplot(aes(
      x = SampleID,
      y = (nor.value + 1) %>% log2()
    ), outlier.shape = NA) +
    xlab("Sample") +
    ylab("value") +
    scale_y_continuous(limits = box.length.log) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = title)

  return(
    ggarrange(
      p1,
      p2,
      ncol = 2,
      nrow = 1,
      align = "v",
      font.label = list(
        size = 12,
        color = "black",
        face = "bold",
        family = NULL,
        position = "top"
      )
    )
  )
}
contrast.dba <- function(norm, title) {
  dba.analyz <- dba.analyze(norm, method = DBA_ALL_METHODS)
  #  summary of results
  th <- norm$config$th
  print(th)
  comp.table <- dba.show(dba.analyz,
    bContrasts = T,
    th = th
  )
  print("compare table: ")
  print(comp.table)

  # # output compare table ----
  # if (dir.exists("./comp_table")){
  #   setwd("./comp_table")
  # } else {
  #   dir.create("./comp_table")
  #   setwd("./comp_table")
  # }
  # #  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
  # inALL=NULL
  # pdf(paste("./dba.analyz.venn", title, ".pdf"))
  # for (OL in 1:length(comp.table$Factor)) {
  #   if(comp.table[OL,]$DB.edgeR >0 && comp.table[OL,]$DB.DESeq2 >0 ){
  #     Venn<- dba.plotVenn(dba.analyz, contrast=OL,
  #                         th=th, method=DBA_ALL_METHODS,
  #                         DataType=DBA_DATA_FRAME)
  #     inAll.1 <- data.frame(comp.table[OL,])%>%
  #       mutate(inBoth=Venn$inAll%>%length)
  #   }else{
  #     inAll.1 <- data.frame(comp.table[OL,])%>%
  #       mutate(inBoth=0)
  #   }
  #   inALL<-inALL%>%bind_rows(inAll.1)
  # }
  # dev.off()
  #
  # write.table(inALL,
  #             file = paste0("./comp.table.", title, ".txt"),
  #             sep = "\t", quote = FALSE, row.names = TRUE)
  # setwd(outpath)


  # output report -----
  print("report table outout")

  for (i in 1:nrow(comp.table)) {
    if (dir.exists(file.path(outpath, "DBS_table"))) {
      setwd(file.path(outpath, "DBS_table"))
    } else {
      dir.create(file.path(outpath, "DBS_table"))
      setwd(file.path(outpath, "DBS_table"))
    }

    if (dir.exists(file.path(outpath, "DBS_table", title))) {
      setwd(file.path(outpath, "DBS_table", title))
    } else {
      dir.create(file.path(outpath, "DBS_table", title))
      setwd(file.path(outpath, "DBS_table", title))
    }

    print(paste0(comp.table$Group[i], comp.table$Group2[i]))
    comp1.edgeR <-
      dba.report(dba.analyz,
        method = DBA_EDGER,
        contrast = i,
        bCalled = T,
        bNormalized = T,
        bCounts = T,
        bUsePval = F,
        fold = 0,
        th = 1
      )
    # sum(comp1.edgeR$Fold > 0)
    # sum(comp1.edgeR$Fold < 0)
    # sum(comp1.edgeR$FDR < th)


    comp1.deseq <-
      dba.report(dba.analyz,
        method = DBA_DESEQ2,
        contrast = i,
        bCalled = T,
        bNormalized = T,
        bCounts = T,
        bUsePval = F,
        fold = 0,
        th = 1
      )

    # sum(comp1.edgeR$Fold > 0)
    # sum(comp1.edgeR$Fold < 0)
    # sum(comp1.edgeR$FDR < th)


    # EdgeR
    out <- as.data.frame(comp1.edgeR)
    write.table(
      out,
      file = paste0(title, ".", comp.table$Group[i], "_", comp.table$Group2[i], ".dba.report_edgeR.txt"),
      sep = "\t",
      quote = F,
      row.names = F
    )
    # DESeq2
    out <- as.data.frame(comp1.deseq)
    write.table(
      out,
      file = paste0(title, ".", comp.table$Group[i], "_", comp.table$Group2[i], ".dba.report_DEseq2.txt"),
      sep = "\t",
      quote = F,
      row.names = F
    )
    setwd(outpath)
  }

  # MA plot -----
  # if (dir.exists("./MA_plot")){
  #   setwd("./MA_plot")
  # } else {
  #   dir.create("./MA_plot")
  #   setwd("./MA_plot")
  # }
  # pdf(paste0("./dba.analyz.MA.", title, ".pdf"), width = 18, height = 23)
  # par(mfrow = c(11, 6))
  # for (i in 1:nrow(comp.table)) {
  #   dba.plotMA(
  #     dba.analyz,
  #     method = DBA_DESEQ2,
  #     contrast = i,
  #     bUsePval=F,
  #     th=th,
  #     fold=fold, bSignificant=T,
  #     bNormalized = TRUE, bSmooth= T,
  #     sub = paste("DEseq2", title)
  #   )
  #   dba.plotMA(
  #     dba.analyz,
  #     method = DBA_EDGER,
  #     contrast = i,
  #     bUsePval=F,
  #     th=th,
  #     fold=fold, bSignificant=T,
  #     bNormalized = TRUE, bSmooth= T,
  #     sub = paste("edgeR", title)
  #   )
  # }
  # dev.off()
  setwd(outpath)
}
DiffBind.build.fun <- function(samplesheet_path = samplesheet_path, outpath = outpath, peakset = NULL) {
  if (!is.null(peakset)) {
    print("specific peakset in use")
  } else {
    print("NULL peakset in use, the peakset would use consensus peak set, which Peaks identified in at least two replicates of at least one sample group.")
    peakset <- NULL
  }
  print(peakset)
  # parameter for cut-off
  th <- 0.05
  # MA.fold=1
  level <- c(
    "cot",
    "em",
    "mm",
    "lm",
    "ed",
    "md",
    "ld",
    "dry",
    "sdlg3h",
    "sdlg16h",
    "sdlg27h",
    "s8d"
  )


  setwd(outpath)
  sessionInfo()

  # DiffBind analysis -----
  samples <-
    read.csv(samplesheet_path,
      stringsAsFactors = F
    ) %>%
    # group_by(Condition, Factor) %>%
    # mutate( Replicate=row_number()) %>%
    # ungroup() %>%
    mutate(Condition = Condition %>% factor(levels = level)) %>%
    arrange(Condition)

  head(samples)
  unique(samples$SampleID)
  getwd()
  if (file.exists("./dba_DBA.RData") == T) {
    # Load the DBdata.count.RData -----
    lnames <- load("./dba_DBA.RData")
    DBdata.count <- DBAobject
  } else {
    # Load the sample sheet contains necessary libraries and data. -----
    file_check <- c(
      file.exists(samples$bamReads),
      file.exists(samples$bamControl),
      file.exists(samples$Peaks)
    )
    if (which(file_check == F) %>% length() == 0) {
      cat("all files exist \n")
    }

    # Create a DBA object. -----
    DBdata <-
      dba(
        sampleSheet = samples, # peakCaller="macs",
        attributes = c(DBA_ID, DBA_REPLICATE, DBA_FACTOR, DBA_CONDITION),
        config = list(
          th = th,
          AnalysisMethod = DBA_ALL_METHODS,
          doBlacklist = F,
          doGreylist = F
        )
      )


    # Counting values ------
    # generate values from all score parameters for dba.count.
    # via correlation heatmap, hcluster heatmap, PCA plot and box plot of values between to visulization results of each score.
    if (!is.null(peakset)) {
      print("specific peakset input")
      DBdata.count <- dba.count(DBdata,
        peaks = peakset, bUseSummarizeOverlaps = T,
        summits = TRUE,
        bParallel = T,
        bScaleControl = T,
        bSubControl = T
      )
    } else {
      DBdata$minOverlap %>% print()
      DBdata.1 <- dba.peakset(DBdata, consensus = DBA_CONDITION, minOverlap = 0.66)
      DBdata.1$minOverlap %>% print()
      DBdata.1 <- dba(DBdata.1, mask = DBdata.1$masks$Consensus, minOverlap = 1)
      dba.show(DBdata.1)
      # dba.plotVenn(DBdata.1, DBdata.1$masks$Consensus)
      cons.peaks <- dba.peakset(DBdata.1, bRetrieve = T)
      print("consensus peaks head")
      print(cons.peaks)

      print("summits optional in dba.count is in use")
      # if( samples$Peaks %>% str_detect("narrow")%>% unique() ){
      # print("peak type is narrow peak, summits optional in dba.count is in use")
      DBdata.count <- dba.count(DBdata,
        peaks = cons.peaks, bUseSummarizeOverlaps = T,
        filter = 0,
        summits = T,
        bParallel = T,
        bScaleControl = T,
        bSubControl = T
      )

      # }else if( samples$Peaks %>% str_detect("broad") %>% unique() ){
      #   print("peak type is broad peak, summits optional in dba.count is not in use")
      #   DBdata.count<-dba.count(DBdata, peaks = cons.peaks, bUseSummarizeOverlaps = T,
      #                           summits = F,
      #                           bParallel = T,
      #                           bScaleControl = T,
      #                           bSubControl = T)
      #
      # }
      print("DBdata.count:")
      dba.show(DBdata.count) %>% print()

      print("if summits in use: (0: summits in use and the width is estimated by DiffBind; NULL: summits is not in use)")
      print(DBdata.count$summits)
      print(" ")


      write.table(cons.peaks %>% as.data.frame() %>% dplyr::select(seqnames, start, end, width, strand) %>% distinct(),
        "./consensus.peak.bed",
        quote = F, col.names = T, row.names = F, sep = "\t", eol = "\n"
      )
    }


    dba.save(DBdata.count,
      file = "DBA", dir = ".", pre = "dba_", ext = "RData",
      bRemoveAnalysis = T, bRemoveBackground = FALSE,
      bCompress = FALSE
    )
  }

  if (file.exists("./DBdata.count.score.global_binding_matrix.txt") == T) {
    # Load the DBdata.count.RData -----
    print("global_binding matrix already exist; return DBdata.count")
  } else {
    print("global_binding matrix counting...")
    print(DBdata.count$samples)
    # DBdata.count$peaks
    print(DBdata.count$config)

    info <- dba.show(DBdata.count)
    libsizes <- cbind(
      LibReads = info$Reads,
      FRiP = info$FRiP, PeakReads = round(info$Reads * info$FRiP)
    )
    rownames(libsizes) <- info$ID
    libsizes
    write.table(libsizes,
      file = paste0("./DBdata.count.libsizes.txt"),
      sep = "\t", quote = FALSE, row.names = TRUE
    )

    # dba.count.score plots -----
    # score<-c("DBA_SCORE_READS", "DBA_SCORE_CONTROL_READS",
    #          "DBA_SCORE_READS_FOLD", "DBA_SCORE_READS_MINUS", "DBA_SCORE_RPKM",
    #          "DBA_SCORE_RPKM_FOLD", "DBA_SCORE_RPKM_MINUS", "DBA_SCORE_SUMMIT",
    #          "DBA_SCORE_SUMMIT_ADJ", "DBA_SCORE_SUMMIT_POS"
    # ) %>%as.character()
    score <- c("DBA_SCORE_RPKM_FOLD") %>% as.character()

    # score<-c("DBA_SCORE_RPKM_FOLD"
    # ) %>%as.character()

    pdf(paste0("./DBdata.count.score.cluster.pdf"))
    par(mfrow = c(4, 3))
    global_matrix <- global_matrix.1 <- NULL
    a <- 0
    for (scor in score) {
      a <- a + 1
      DBdata.count.score <- dba.count(DBdata.count,
        peaks = peakset, score = get(scor)
      )
      method <- c("pearson", "spearman")
      for (m in method) {
        plot(
          DBdata.count.score,
          sub = paste0(score),
          attributes = c(DBA_CONDITION, DBA_REPLICATE),
          ColAttributes = DBA_CONDITION,
          # margin = 5,
          distMethod = m
        )
      }

      # hcluster heatmap ---
      title <- paste("DBdata.count.score", scor)
      heatmap.method <- list(
        hclust = c("ward.D2", "complete", "average", "mcquitty"),
        dist = c("euclidean", "manhattan", "correlation", "maximum")
      )
      for (hclust in 1:length(heatmap.method$hclust)) {
        for (dist in 1:length(heatmap.method$dist)) {
          hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
          dba.plotHeatmap(
            DBdata.count.score,
            attributes = c(DBA_CONDITION, DBA_REPLICATE),
            ColAttributes = DBA_CONDITION,
            correlations = FALSE,
            scale = "row",
            margin = 5,
            colScheme = hmap,
            main = paste("heatmap", heatmap.method$dist[dist], heatmap.method$hclust[hclust], sep = "_"),
            distMethod = heatmap.method$dist[dist],
            hclustfun = function(x) {
              hclust(x, method = heatmap.method$hclust[hclust])
            },
            sub = title
          )
        }
      }

      DBdata.count.score %>% dba.plotPCA(sub = title)

      if (a == 1) {
        global_matrix <- dba.peakset(DBdata.count.score,
          bRetrieve = T,
          DataType = DBA_DATA_FRAME
        ) %>%
          pivot_longer(
            cols = 4:(length(DBdata.count.score$samples$SampleID) + 3),
            names_to = "SampleID", values_to = scor
          )
        print(global_matrix %>% head())
      } else {
        global_matrix.1 <- dba.peakset(DBdata.count.score,
          bRetrieve = T,
          DataType = DBA_DATA_FRAME
        ) %>%
          pivot_longer(
            cols = 4:(length(DBdata.count.score$samples$SampleID) + 3),
            names_to = "SampleID", values_to = scor
          )
        print(global_matrix.1 %>% head())
        global_matrix <- global_matrix %>% left_join(global_matrix.1)
      }
    }

    dev.off()

    # box plot
    value.1 <- global_matrix %>%
      mutate(ID = paste(CHR, START, END, sep = "_"), .before = CHR) %>%
      dplyr::select(-c(CHR, START, END)) %>%
      mutate(
        stage = SampleID %>% str_extract("^\\w+"), .after = SampleID,
        SampleID = SampleID %>% str_remove("\\.H3.*")
      ) %>%
      mutate(stage = factor(
        stage,
        levels = level
      )) %>%
      arrange(stage)

    Plots <- list()
    a <- 0
    for (p in 4:length(colnames(value.1))) {
      a <- a + 1
      input <- value.1[, c(1:3, p)] %>%
        mutate(SampleID = SampleID %>%
          factor(levels = unique(SampleID)))
      colnames(input)[4] <- "value"

      outlier <- quantile(input$value, c(0.25, 0.75))
      if (outlier[1] - IQR(input$value) * 1.5 >= 0) {
        box.length <- c(
          outlier[1] - IQR(input$value) * 1.5,
          outlier[2] + IQR(input$value) * 1.5
        )
      } else {
        box.length <- c(
          0,
          outlier[2] + IQR(input$value) * 1.5
        )
      }
      print(head(input))
      Plots[[a]] <- input %>%
        ggplot() +
        geom_boxplot(aes(
          x = SampleID,
          y = value
        ), outlier.shape = NA) +
        xlab("Sample") +
        ylab("value") +
        scale_y_continuous(limits = box.length) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(title = colnames(value.1)[p])
      a <- a + 1
      outlier <- quantile((input$value + 1) %>% log2(), c(0.25, 0.75))
      if (outlier[1] - IQR((input$value + 1) %>% log2()) * 1.5 >= 0) {
        box.length <- c(
          outlier[1] - IQR((input$value + 1) %>% log2()) * 1.5,
          outlier[2] + IQR((input$value + 1) %>% log2()) * 1.5
        )
      } else {
        box.length <- c(
          0,
          outlier[2] + IQR((input$value + 1) %>% log2()) * 1.5
        )
      }
      Plots[[a]] <- input %>%
        ggplot() +
        geom_boxplot(aes(
          x = SampleID,
          y = (value + 1) %>% log2()
        ), outlier.shape = NA) +
        xlab("Sample") +
        ylab("log2(value)") +
        scale_y_continuous(limits = box.length) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(title = colnames(value.1)[p])
    }

    plot <-
      ggarrange(
        plotlist = Plots,
        labels = LETTERS[1:length(Plots)],
        ncol = 2,
        nrow = (length(Plots) / 2),
        align = "v",
        font.label = list(size = 12, color = "black", face = "bold", family = NULL, position = "top")
      )
    sacle.factor <- length(Plots) / 2 / 11
    ggsave(
      paste0("dba.count.score.value.boxplot.pdf"),
      plot,
      width = 12,
      height = 24 * sacle.factor,
      units = "in",
      device = "pdf"
    )


    write.table(
      global_matrix,
      file = paste0("./DBdata.count.score.global_binding_matrix.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = F
    )
  }
  return(DBdata.count)
}
pairwise.matrix <- function(res.matrix, tool) {
  level <- c(
    "cot",
    "em",
    "mm",
    "lm",
    "ed",
    "md",
    "ld",
    "dry",
    "sdlg3h",
    "sdlg16h",
    "sdlg27h",
    "s8d"
  )

  res.tool <-
    res.matrix %>%
    filter(tools %>% str_detect(tool)) %>%
    mutate(
      positive = ifelse(tools %>% str_detect("up"),
        pos.,
        neg.
      ),
      negative = ifelse(tools %>% str_detect("down"),
        pos.,
        neg.
      )
    ) %>%
    dplyr::select(positive, negative, DBS_number) %>%
    mutate(
      positive = positive %>% factor(
        levels = level
      ),
      negative = negative %>% factor(
        levels = level
      )
    ) %>%
    arrange(negative) %>%
    arrange(positive)
  res.tool.1 <- res.tool %>%
    pivot_wider(names_from = negative, values_from = DBS_number) %>%
    # relocate(cot, .before=em)%>%
    dplyr::select(c(positive, paste(res.tool$positive %>% as.character(), sep = ","))) %>%
    dplyr::rename("row.vs.col" = positive)
}
matrix.plot <- function(data, title, color) {
  level <- c(
    "cot",
    "em",
    "mm",
    "lm",
    "ed",
    "md",
    "ld",
    "dry",
    "sdlg3h",
    "sdlg16h",
    "sdlg27h",
    "s8d"
  )
  level <- level

  matrix <- data[, -1]
  matrix$positive <- data$row.vs.col
  x.melt <- reshape2::melt(matrix, id.vars = "positive") %>% dplyr::rename(negative = variable)
  plot <- x.melt %>%
    mutate(
      compaire.1 = positive %>% factor(levels = level),
      compaire.2 = negative %>% factor(levels = level)
    ) %>%
    ggplot(aes(x = compaire.1, y = compaire.2, fill = value)) +
    geom_tile(colour = "white", size = 0.25) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    coord_fixed() +
    # coord_flip() +
    scale_fill_gradientn(
      colours = color
      # limits=c(0, totalMergedpeak)
    ) +
    geom_text(aes(label = comma(value)),
      color = "black",
      size = 4, angle = 90, fontface = "bold"
    ) +
    ggtitle(title %>% as.character() %>% basename()) +
    theme(
      legend.key.height = unit(0.6, "cm"),
      plot.title = element_text(vjust = -20, size = 13),
      legend.text = element_text(angle = 90),
      legend.title = element_text(angle = 90, hjust = 1),
      legend.position = "top",
      axis.ticks = element_line(size = 0.5),
      plot.background = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(
        size = 10,
        face = "bold",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 10,
        angle = 90,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.y = element_text(vjust = 5)
    )
  print(paste("output:", file.path(outpath, paste0(title, "_pair_wise_comparison_matrix.pdf"))))
  pdf(file.path(outpath, paste0(title, "_pair_wise_comparison_matrix.pdf")), width = 9, height = 9)
  grid.newpage()
  print(plot,
    vp = viewport(
      width = 1,
      height = 1,
      angle = -90
    )
  )
  dev.off()

  # ggsave(file.path(outpath, paste0(title,"_pair_wise_comparison_matrix.pdf")), plot,
  #        width = 12, height = 12, units = "in",
  #        device = "pdf")
}
DiffBind.nor.DE.fun <- function(DBdata.count, peakset = NULL, color = color, log2FC = 1, q = 0.05) {
  print(paste("fold change cut-off:", 2^log2FC, "; q cut-off:", q))

  if (!is.null(peakset)) {
    print("specific peakset in use")
  } else {
    print("NULL peakset in use")
    peakset <- NULL
  }
  nor.parameter <- data.frame(
    nor = c("DBA_NORM_LIB"),
    lib = c("DBA_LIBSIZE_FULL"),
    BG = c(F),
    offset = c(F),
    title = c("LIB-FULL")
  )

  Plots <- list()
  a <- 0
  for (nor in 1:length(nor.parameter$title)) {
    a <- a + 1
    norm <- NULL
    norm <- dba.normalize(
      DBdata.count,
      method = DBA_ALL_METHODS,
      normalize = get(nor.parameter$nor[nor]),
      library = get(nor.parameter$lib[nor]),
      background = nor.parameter$BG[nor],
      offset = nor.parameter$offset[nor]
    )
    # norm<- dba.count(norm, peaks=NULL, score = DBA_SCORE_NORMALIZED)
    # norm.table <- dba.peakset(norm, peaks=NULL, bRetrieve = T, DataType = DBA_DATA_FRAME)

    norm <- dba.contrast(norm, categories = DBA_CONDITION, minMembers = 2)
    Plots[[a]] <- value.box.plot(norm, nor.parameter$title[nor], peakset = peakset)
    pdf(paste0("./DBdata.normalized.cluster.", nor.parameter$title[nor], ".pdf"))
    dba.plot(norm, nor.parameter$title[nor])
    dev.off()
    contrast.dba(norm, nor.parameter$title[nor])
  }


  plot.nor <-
    ggarrange(
      plotlist = Plots,
      labels = LETTERS[1:length(Plots) * 2],
      ncol = 1,
      nrow = length(Plots),
      align = "v",
      font.label = list(size = 12, color = "black", face = "bold", family = NULL, position = "top")
    )

  ggsave(
    paste0("dba.normalized.value.boxplot.pdf"),
    plot.nor,
    width = 12,
    height = 24 / 6,
    units = "in",
    device = "pdf"
  )

  ##### ---------
  # count pair-wise matrix: DEseq2, edgeR, and inBoth -------

  # totalMergedpeak<-DBAobject$totalMerged
  report.files <- data.frame(
    edgeR.report = list.dirs(file.path(outpath, "DBS_table")) %>%
      list.files("dba.report_edgeR", full.names = T)
  ) %>%
    mutate(
      pair_wise = edgeR.report %>% as.character() %>% basename() %>%
        str_extract("\\w+(?=(\\.dba))"),
      method = edgeR.report %>% as.character() %>% dirname() %>%
        str_extract("(?<=(DBS_table/)).*$")
    ) %>%
    left_join(
      data.frame(DEseq2.report = list.dirs(file.path(outpath, "DBS_table")) %>%
        list.files("dba.report_DEseq2", full.names = T)) %>%
        mutate(
          pair_wise = DEseq2.report %>% as.character() %>% basename() %>%
            str_extract("\\w+(?=(\\.dba))"),
          method = DEseq2.report %>% as.character() %>% dirname() %>%
            str_extract("(?<=(DBS_table/)).*$")
        )
    )

  report.files$method %>% unique()

  # creat folders ----
  for (m in report.files$method %>% unique()) {
    dir.create(file.path(outpath, "DBS_bed", m, "edgeR"), recursive = T)
    dir.create(file.path(outpath, "DBS_bed", m, "DEseq2"), recursive = T)
    dir.create(file.path(outpath, "DBS_bed", m, "inBoth"), recursive = T)
  }

  # pair-wise comparison peak number statistic ----
  res.allmethod <- NULL
  m <- "LIB-FULL"
  for (m in unique(report.files$method)) {
    res <- NULL
    report.files.1 <- report.files %>% filter(method == m)
    for (i in 1:length(report.files.1$pair_wise)) {
      # cmd<-paste0( "wc -l ", report.files.1$edgeR.report[i]%>%as.character(), " | awk '{print $1}'")
      # file.line=system(cmd, intern = TRUE)

      edgeR.d1 <- read.delim(
        report.files.1$edgeR.report[i] %>% as.character(),
        stringsAsFactors = F
      )
      edgeR.d1 <- edgeR.d1 %>%
        mutate(
          peakID = paste(seqnames, start, end, sep = "_"),
          .before = seqnames
        )
      # dplyr::count(edgeR.d1, Called2, Called1) %>% print()
      head(edgeR.d1)
      sig.edgeR.up <- edgeR.d1 %>%
        filter(FDR < q) %>%
        filter(Fold %>% as.numeric() > log2FC)
      # dplyr::count(sig.edgeR.up, Called2, Called1) %>% print()
      print(paste(min(sig.edgeR.up$Fold), max(sig.edgeR.up$Fold)))

      sig.edgeR.down <- edgeR.d1 %>%
        filter(FDR < q) %>%
        filter(Fold %>% as.numeric() < -log2FC)
      # dplyr::count(sig.edgeR.down, Called2, Called1)
      print(paste(min(sig.edgeR.down$Fold), max(sig.edgeR.down$Fold)))


      DEseq2.d1 <- read.delim(
        report.files.1$DEseq2.report[i] %>% as.character(),
        stringsAsFactors = F
      )
      DEseq2.d1 <- DEseq2.d1 %>%
        mutate(
          peakID = paste(seqnames, start, end, sep = "_"),
          .before = seqnames
        )
      # dplyr::count(DEseq2.d1, Called2, Called1)
      sig.DEseq2.up <- DEseq2.d1 %>%
        filter(FDR < q) %>%
        filter(Fold %>% as.numeric() > log2FC)
      # dplyr::count(sig.DEseq2.up, Called2, Called1)
      print(paste(min(sig.DEseq2.up$Fold), max(sig.DEseq2.up$Fold)))

      sig.DEseq2.down <- DEseq2.d1 %>%
        filter(FDR < q) %>%
        filter(Fold %>% as.numeric() < -log2FC)
      # dplyr::count(sig.DEseq2.down, Called2, Called1)
      print(paste(min(sig.DEseq2.down$Fold), max(sig.DEseq2.down$Fold)))

      sig.inBoth.up <- dplyr::intersect(sig.edgeR.up %>% dplyr::select(seqnames, start, end, peakID), sig.DEseq2.up %>% dplyr::select(seqnames, start, end, peakID))
      # print(paste(min(sig.inBoth.up$Fold), max(sig.inBoth.up$Fold)))

      sig.inBoth.down <- dplyr::intersect(sig.edgeR.down %>% dplyr::select(seqnames, start, end, peakID), sig.DEseq2.down %>% dplyr::select(seqnames, start, end, peakID))
      # print(paste(min(sig.inBoth.down$Fold), max(sig.inBoth.down$Fold)))

      write.table(sig.edgeR.up %>% dplyr::select(seqnames, start, end, peakID),
        file.path(
          outpath, "DBS_bed", m, "edgeR",
          paste0(report.files.1$pair_wise[i], "_q0pt05.FC2.edgeR.up.bed")
        ),
        sep = "\t", row.names = F, col.names = F, eol = "\n", quote = F
      )
      write.table(sig.edgeR.down %>% dplyr::select(seqnames, start, end, peakID),
        file.path(
          outpath, "DBS_bed", m, "edgeR",
          paste0(report.files.1$pair_wise[i], "_q0pt05.FC2.edgeR.down.bed")
        ),
        sep = "\t", row.names = F, col.names = F, eol = "\n", quote = F
      )

      write.table(sig.DEseq2.up %>% dplyr::select(seqnames, start, end, peakID),
        file.path(
          outpath, "DBS_bed", m, "DEseq2",
          paste0(report.files.1$pair_wise[i], "_q0pt05.FC2.DEseq2.up.bed")
        ),
        sep = "\t", row.names = F, col.names = F, eol = "\n", quote = F
      )
      write.table(sig.DEseq2.down %>% dplyr::select(seqnames, start, end, peakID),
        file.path(
          outpath, "DBS_bed", m, "DEseq2",
          paste0(report.files.1$pair_wise[i], "_q0pt05.FC2.DEseq2.down.bed")
        ),
        sep = "\t", row.names = F, col.names = F, eol = "\n", quote = F
      )


      write.table(sig.inBoth.up %>% dplyr::select(seqnames, start, end, peakID),
        file.path(
          outpath, "DBS_bed", m, "inBoth",
          paste0(report.files.1$pair_wise[i], "_q0pt05.FC2.inBoth.up.bed")
        ),
        sep = "\t", row.names = F, col.names = F, eol = "\n", quote = F
      )
      write.table(sig.inBoth.down %>% dplyr::select(seqnames, start, end, peakID),
        file.path(
          outpath, "DBS_bed", m, "inBoth",
          paste0(report.files.1$pair_wise[i], "_q0pt05.FC2.inBoth.down.bed")
        ),
        sep = "\t", row.names = F, col.names = F, eol = "\n", quote = F
      )

      res.1 <-
        data.frame(
          pair_wise = report.files.1$pair_wise[i],
          method = m,
          sig.edgeR.up = sig.edgeR.up$peakID %>% length(),
          sig.edgeR.down = sig.edgeR.down$peakID %>% length(),
          sig.DEseq2.up = sig.DEseq2.up$peakID %>% length(),
          sig.DEseq2.down = sig.DEseq2.down$peakID %>% length(),
          sig.inBoth.up = sig.inBoth.up$peakID %>% length(),
          sig.inBoth.down = sig.inBoth.down$peakID %>% length()
        )
      res <- res %>% bind_rows(res.1)
    }
    res.allmethod <- res.allmethod %>% bind_rows(res)
  }
  write.table(res.allmethod, file.path(outpath, "pair-wise_peak.number.txt"),
    quote = F, sep = "\t", row.names = F, eol = "\n"
  )
  # test
  # table.txt<-"/nas/qnap4_1/Project/Gm/ChIPseq_analysis/Gmax189/Analysis/Methylation/DiffBind/output/H3K27me3/all/IDR/DBS_table/LIB-FULL/LIB-FULL.cot_dry.dba.report_DEseq2.txt"
  # cmd<-paste0("awk 'NR == 1 || ($11 < 0.05 && ($9 > 2))' ",table.txt, " | grep -v '^seqnames' | wc -l")
  # cmd<-paste0("awk 'NR == 1 || ($11 < 0.05 && ($9 > 2 || $9 < -2))' ",table.txt, " | grep -v '^seqnames' | wc -l")
  # system(cmd)


  # generate pair-wise comparison peak number matrix -----

  # m="LIB-FULL"

  for (m in report.files$method %>% unique()) {
    m <- m %>% as.character()
    # if(c("res.edgeR", "res.DEseq2", "res.inBoth") %in% ls()!=T){
    res.matrix <- res.allmethod %>%
      filter(method == m) %>%
      dplyr::select(-method) %>%
      distinct()

    res.matrix <- res.matrix %>%
      pivot_longer(
        cols = 2:length(colnames(res.matrix)),
        names_to = "tools",
        values_to = "DBS_number"
      ) %>%
      separate(
        col = pair_wise,
        into = c("pos.", "neg."),
        sep = "_"
      )
    res.edgeR <- pairwise.matrix(res.matrix, "edgeR")
    res.DEseq2 <- pairwise.matrix(res.matrix, "DEseq2")
    res.inBoth <- pairwise.matrix(res.matrix, "inBoth")

    dir.create(file.path(outpath, "DBS_matrix", m), recursive = T)
    # output pair-wise comparison matrix --------
    write.table(res.edgeR,
      file.path(outpath, "DBS_matrix", m, "pair_wise_matrix_q0pt05.FC2.edgeR.txt"),
      quote = F, row.names = F, sep = "\t", eol = "\n"
    )

    write.table(res.DEseq2,
      file.path(outpath, "DBS_matrix", m, "pair_wise_matrix_q0pt05.FC2.DEseq2.txt"),
      quote = F, row.names = F, sep = "\t", eol = "\n"
    )

    write.table(res.inBoth,
      file.path(outpath, "DBS_matrix", m, "pair_wise_matrix_q0pt05.FC2.inBoth.txt"),
      quote = F, row.names = F, sep = "\t", eol = "\n"
    )

    # }else{
    #   res.edgeR<-read.delim(file.path(outpath, "DBS_matrix", m, "pair_wise_matrix_q0pt05.FC2.edgeR.txt"), stringsAsFactors = F)
    #   res.DEseq2<-read.delim(file.path(outpath, "DBS_matrix", m, "pair_wise_matrix_q0pt05.FC2.DEseq2.txt"), stringsAsFactors = F)
    #   res.inBoth<-read.delim(file.path(outpath, "DBS_matrix", m, "pair_wise_matrix_q0pt05.FC2.inBoth.txt"), stringsAsFactors = F)
    # }

    # pair-wise comparison heatmap -------
    matrix.plot(res.edgeR, file.path("DBS_matrix", m, paste0(m, "_edgeR_q0pt05_FC2")), color)
    matrix.plot(res.DEseq2, file.path("DBS_matrix", m, paste0(m, "_DEseq2_q0pt05_FC2")), color)
    matrix.plot(res.inBoth, file.path("DBS_matrix", m, paste0(m, "_inBoth_q0pt05_FC2")), color)
  }
}
