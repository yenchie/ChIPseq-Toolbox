annoGene <- function(peak, promoter = c(-1000, 0)) {
    print(promoter)
    nrow(gene.df)
    unique(gene.df$geneID) %>% length()

    ref.annt <- gene.df %>%
        mutate(
            promoter.start = ifelse(strand == ">", start + promoter[1], end - promoter[2]),
            promoter.end = ifelse(strand == ">", start + promoter[2], end - promoter[1])
        ) %>%
        dplyr::select(seqnames, promoter.start, promoter.end, geneID, Name) %>%
        distinct() %>%
        dplyr::rename(start = promoter.start, end = promoter.end) %>%
        mutate(feature = "promoter") %>%
        bind_rows(
            gene.df %>%
                dplyr::select(seqnames, start, end, geneID, Name) %>%
                distinct() %>%
                mutate(feature = "genic")
        ) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = T)

    head(ref.annt)

    hits <- findOverlaps(ref.annt, peak) %>%
        as.data.frame() %>%
        mutate(across(everything(), as.character))
    nrow(hits)

    res <- hits %>%
        left_join(peak %>% as.data.frame() %>% rownames_to_column(var = "subjectHits")) %>%
        left_join(ref.annt %>% as.data.frame() %>% rownames_to_column(var = "queryHits") %>%
            dplyr::rename(gene.start = start, gene.end = end, Strand = strand, gene.width = width))

    head(res)
    nrow(res)
    genes <- unique(res$geneID)
    length(genes)
    return(genes)
}
