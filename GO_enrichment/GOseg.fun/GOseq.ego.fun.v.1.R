#########################################
# GOseq.ego.fun.v.1.R
# revised date: 2024.05.06
# reviser: YCW
###
## be called by function:
# GOseq.hyper.FDR.DEG.fun.v.2.R
###############

# source.path <- getCurrentFileLocation()
if (!require("clusterProfiler", quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
}

if (!require("topGO", quietly = TRUE)) {
    BiocManager::install("topGO")
}

if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
if (!require("jamba", quietly = TRUE)) {
    remotes::install_github("jmw86069/jamba")
}

if (!require("Rgraphviz", quietly = TRUE)) {
    BiocManager::install("Rgraphviz")
}



require(goseq)
require(dplyr)
require(stringr)
require(tidyr)
require(GO.db)
require(ggplot2)
require(clusterProfiler)

# source.path <- "/bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Toolbox/GO_enrichment"
# print(source.path)
# print(file.path(source.path %>% list.dirs() %>% list.files("jamenrich-enrichdf2er.R", full.names = T)))
# source(file.path(source.path %>% list.dirs() %>% list.files("jamenrich-enrichdf2er.R", full.names = T)))
# print(file.path(source.path %>% list.dirs() %>% list.files("jamenrich-import.R", full.names = T)))
# source(file.path(source.path %>% list.dirs() %>% list.files("jamenrich-import.R", full.names = T)))

GO.seq.ego <- function(pwf, GO_DB, termfile = termfile, q.cut.off) {
    termfile <- termfile
    head(termfile) %>% print()
    result <- list()
    head(pwf)
    GOI.list <- rownames(pwf %>% filter(DEgenes == 1)) %>% unique()
    head(GOI.list)
    colnames(GO_DB) <- c("geneID", "ID", "Ontology")
    head(GO_DB)

    print("genes in GO term list:")
    genes <- GO_DB$geneID %>%
        unique() %>%
        length()
    print(genes)

    print("input GOI Number:")
    print(length(GOI.list))

    print("GOIs in GO term list:")
    GOIs <- length(which(GOI.list %in% GO_DB$geneID))
    print(GOIs)

    ego <- goseq(pwf, gene2cat = GO_DB, method = "Hypergeometric", use_genes_without_cat = T) # Wallenius
    head(ego, n = 10)
    ego$qval <- p.adjust(ego$over_represented_pvalue, method = "BH")
    head(ego)
    nrow(ego)
    result[["ego"]] <- ego

    # # use_genes_without_cat test ----
    #   go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric",use_genes_without_cat=T)
    #   go.wall.BP.pwf <- goseq(pwf,gene2cat=GO_BP, method="Hypergeometric", use_genes_without_cat=F)
    #   go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
    #   go.wall.BP.pwf$qval <- p.adjust(go.wall.BP.pwf$over_represented_pvalue,method="BH")
    #   plot(
    #     log10(go.wall.BP.pwf[, 8]),
    #     log10(go.wall.BP[match(go.wall.BP.pwf[, 1], go.wall.BP[, 1]), 8]),
    #     xlab = "use_genes_without_cat=F, log10(qval)",
    #     ylab = "use_genes_without_cat=T, log10(qval)",
    #     xlim = c(-3, 0),
    #     ylim = c(-3, 0)
    #   )
    #   abline(0,1,col=3,lty=2)
    #   abline(h=-2, col="red",lty=4)
    #   abline(v=-2, col="red",lty=4)
    #   #x<-log10(0.01)


    #   # different method test ----
    #   go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric",use_genes_without_cat=T)
    #   go.wall.BP.pwf <- goseq(pwf,gene2cat=GO_BP, method="Wallenius", use_genes_without_cat=T)
    #   go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
    #   go.wall.BP.pwf$qval <- p.adjust(go.wall.BP.pwf$over_represented_pvalue,method="BH")
    #   plot(
    #     log10(go.wall.BP.pwf[, 8]),
    #     log10(go.wall.BP[match(go.wall.BP.pwf[, 1], go.wall.BP[, 1]), 8]),
    #     xlab = "Wallenius, log10(qval)",
    #     ylab = "Hypergeometric, log10(qval)",
    #     xlim = c(-3, 0),
    #     ylim = c(-3, 0)
    #   )
    #   abline(0,1,col=3,lty=2)
    #   abline(h=-2, col="red",lty=4)
    #   abline(v=-2, col="red",lty=4)
    #   #x<-log10(0.01)

    #   # different method, different bg test ----
    #   go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric",use_genes_without_cat=T)
    #   go.wall.BP.pwf <- goseq(pwf,gene2cat=GO_BP, method="Wallenius", use_genes_without_cat=F)
    #   go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
    #   go.wall.BP.pwf$qval <- p.adjust(go.wall.BP.pwf$over_represented_pvalue,method="BH")
    #   plot(
    #     log10(go.wall.BP.pwf[, 8]),
    #     log10(go.wall.BP[match(go.wall.BP.pwf[, 1], go.wall.BP[, 1]), 8]),
    #     xlab = "Wallenius, use_genes_without_cat=F, log10(qval)",
    #     ylab = "Hypergeometric, use_genes_without_cat=T, log10(qval)",
    #     xlim = c(-3, 0),
    #     ylim = c(-3, 0)
    #   )
    #   abline(0,1,col=3,lty=2)
    #   abline(h=-2, col="red",lty=4)
    #   abline(v=-2, col="red",lty=4)
    #   #x<-log10(0.01)
    ## -----

    enriched.ego <- subset(ego, ego$qval < q.cut.off)
    head(enriched.ego)
    # if have Go terms qval lower than 0.05, draw a tree.
    if (nrow(enriched.ego) > 0) {
        goIDs <- enriched.ego[, 1]
        pvalue <- enriched.ego[, 8]

        enriched.ego.annot <- merge(enriched.ego, termfile, by.x = "category", by.y = "ID", all.x = TRUE)
        head(enriched.ego.annot)
        ### for table ------
        enriched.ego.annot.1 <- enriched.ego.annot %>%
            dplyr::rename(ID = category, FDR = qval, geneHits = numDEInCat, pathGenes = numInCat, pvalue = over_represented_pvalue) %>%
            left_join(GO_DB) %>%
            filter(geneID %in% GOI.list) %>%
            group_by(ID) %>%
            mutate(
                entries = paste(geneID, collapse = ",")
            ) %>%
            ungroup() %>%
            dplyr::select(-geneID) %>%
            distinct() %>%
            mutate(
                GeneRatio = paste(geneHits, GOIs, sep = "/"),
                BgRatio = paste(pathGenes, genes, sep = "/")
            ) %>%
            as.data.frame()
        # head(enriched.ego.annot.1)

        enriched.ego.annot.1 <- enriched.ego.annot.1[order(enriched.ego.annot.1$FDR), ]
        head(enriched.ego.annot.1)
        result[["fin.result"]] <- enriched.ego.annot.1
    }
    if (nrow(enriched.ego) == 0) {
        enriched.go.wall.annot.BP <- NULL
        result[["fin.result"]] <- NULL
    }
    return(result)
}
DAG.GOseq.fun <- function(GOseq.result, GO_DB, GO.ontology, termfile, GOI.list, q.cut.off, outpath = outpath) {
    print(GO.ontology)
    GO_DB$Ontology %>% unique()
    GOseq.result$ontology %>% unique()

    tmp <- GOseq.result[which(GOseq.result$ontology %>% is.na()), ]

    if (length(which(!tmp$category %in% termfile$ID)) != 0) {
        print("Please check the GO term annotation files to ensure that each term includes information about its ontology classification.")
    }
    # GO_DB[which(GO_DB$ID%in% tmp$category),]
    # termfile[which(termfile$ID%in% tmp$category),]

    head(GOseq.result)
    colnames(GO_DB) <- c("geneID", "ID", "Ontology")
    head(GO_DB)
    genes <- GO_DB[, 1] %>%
        unique()
    print(genes %>%
        length())

    GOIs <- GOI.list[which(GOI.list %in% genes)]
    print(GOIs %>%
        length())

    # GOseq.result %>% head()
    GOseq.result.1 <- GOseq.result %>%
        # na.omit() %>%
        dplyr::rename(ID = category, FDR = qval, geneHits = numDEInCat, pathGenes = numInCat, pvalue = over_represented_pvalue) %>%
        left_join(GO_DB) %>%
        filter(geneID %in% GOIs) %>%
        group_by(ID) %>%
        mutate(
            Genes = paste(geneID %>% unique(), collapse = ",")
        ) %>%
        ungroup() %>%
        dplyr::select(-geneID) %>%
        distinct() %>%
        mutate(
            GeneRatio = paste(geneHits, GOIs %>%
                length(), sep = "/"),
            BgRatio = paste(pathGenes, genes %>%
                length(), sep = "/")
        ) %>%
        as.data.frame() %>%
        dplyr::select("ID", "pvalue", "geneHits", "term", "FDR", "Genes", "GeneRatio", "BgRatio")

    nrow(GOseq.result.1)
    GOseq.result.1$Genes %>% head()
    head(GOseq.result.1)
    colnames(GOseq.result.1)

    enrich_res <- enrichDF2enrichResult(GOseq.result.1,
        qvalueCutoff = q.cut.off,
        pAdjustMethod = "BH",
        keyColname = "ID",
        # pathGenes = "pathGenes",
        geneColname = "Genes",
        geneHits = "geneHits",
        geneRatioColname = "GeneRatio",
        geneDelim = "[,/ ]+",
        geneSep = ",",
        pvalueColname = "pvalue",
        qvalueColname = "FDR",
        descriptionColname = "ID",
        msigdbGmtT = NULL,
        verbose = FALSE
    )
    # enrich_res@gene
    enrich_res@ontology <- GO.ontology
    enrich_res@universe <- GO_DB$geneID %>% unique()
    eGO_DB <- GO_DB %>%
        filter(ID %in% GOseq.result.1$ID) %>%
        distinct()

    enrich_res@geneSets <- split(as.character(eGO_DB$geneID), as.character(eGO_DB$ID))
    print(outpath)

    pdf(outpath)
    tryCatch(
        {
            enrich_res %>% plotGOgraph(useFullNames = TRUE, .NO.CHAR = 35)
        },
        error = function(e) {
            message("An error occurred: ", e$message)
        }
    )
    dev.off() # Ensure the PDF device is turned off in case of an error
}
