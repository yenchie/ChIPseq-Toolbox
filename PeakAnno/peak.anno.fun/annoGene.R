annoGene <- function(peak.gr, promoter = c(-1000, 0)) {
    ## note:
    # per transcript
    # ls(envir = .GlobalEnv) %>% print()

    if (!("transcript.id.df" %in% c(ls(envir = .GlobalEnv) %>% as.character()))) {
        print("missing annotation table in the environment")
        return(NULL)
    } else {
        unique(transcript.id.df$geneID) %>% length()
        ref <- transcript.id.df

        ref.annt <- ref %>%
            filter(feature == "mRNA") %>%
            mutate(
                promoter.start = ifelse(strand == "+", start + promoter[1], end - promoter[2]),
                promoter.end = ifelse(strand == "+", start + promoter[2], end - promoter[1])
            ) %>%
            dplyr::select(seqnames, promoter.start, promoter.end, strand, geneID, transcript.ID) %>%
            distinct() %>%
            dplyr::rename(start = promoter.start, end = promoter.end) %>%
            mutate(feature = "promoter") %>%
            bind_rows(
                ref %>% # filter(feature != "mRNA") %>%
                    dplyr::select(seqnames, start, end, strand, geneID, transcript.ID, feature) %>%
                    distinct()
            ) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = T)

        head(ref.annt)
        unique(ref.annt$feature)

        hits <- findOverlaps(ref.annt, peak.gr)
        hits.df <- hits %>%
            as.data.frame() %>%
            mutate(across(everything(), as.character))

        head(hits.df)
        # p <- Pairs(ref.annt, peak.gr, hits = hits)


        if (nrow(hits.df) != 0) {
            nrow(hits.df)
            priority <- c("promoter", "five_prime_UTR", "CDS", "three_prime_UTR", "mRNA", "intergenic")
            res <- hits.df %>%
                left_join(peak.gr %>% as.data.frame() %>% rownames_to_column(var = "subjectHits")) %>%
                left_join(ref.annt %>% as.data.frame() %>% rownames_to_column(var = "queryHits") %>%
                    dplyr::rename(gene.start = start, gene.end = end, Strand = strand, gene.width = width)) %>%
                dplyr::select(-c(queryHits, subjectHits)) %>%
                mutate(peakID = paste(seqnames, start, end, sep = "__")) %>%
                group_by(peakID) %>%
                arrange(factor(feature, levels = priority)) %>%
                mutate(priority_mark = ifelse(row_number() == 1, 1, seq_along(feature))) %>%
                mutate(feature = ifelse(all(feature == "mRNA"), "intron", feature)) %>% ## only apear in intron across all transcripts annoted
                ungroup() %>%
                mutate(distance_to_TSS = ifelse(Strand == "+",
                    min(c(abs(start - gene.start), abs(end - gene.start))),
                    min(c(abs(start - gene.end), abs(end - gene.end)))
                ))

            head(res %>% as.data.frame())
            nrow(res)

            unalign <- peak.gr %>%
                as.data.frame() %>%
                rownames_to_column(var = "subjectHits") %>%
                filter(!subjectHits %in% hits.df$subjectHits) %>%
                distinct() %>%
                mutate(peakID = paste(seqnames, start, end, sep = "__"))
            head(unalign)
            nrow(unalign)

            priority <- c("promoter", "five_prime_UTR", "CDS", "three_prime_UTR", "intron", "intergenic")

            feature.count <- data.frame(
                feature = "intergenic",
                n = length(peak.gr) - length(hits.df$subjectHits %>% unique())
            ) %>%
                bind_rows(res %>%
                    filter(priority_mark == 1) %>%
                    dplyr::count(feature)) %>%
                mutate(proportion = ((n / sum(n)) %>% round(3)) * 100) %>%
                dplyr::select(-n) %>%
                mutate(feature = feature %>% factor(levels = priority)) %>%
                arrange(feature)
            print(feature.count)
            genes <- unique(res$geneID)

            annoPeak <- new("annoPeak",
                annoGene = genes,
                gr.ref = ref.annt,
                gr.peak = peak.gr,
                promoterRegion = promoter,
                Annotation = res,
                annoStat = feature.count,
                peakNum = length(peak.gr)
            )
        } else {
            res <- NULL
            feature.count <- data.frame(
                feature = "intergenic",
                n = length(peak.gr) - length(hits.df$subjectHits %>% unique())
            ) %>%
                mutate(proportion = ((n / sum(n)) %>% round(3) * 100)) %>%
                dplyr::select(-n)
            print(feature.count)
            annoPeak <- new("annoPeak",
                gr.ref = ref.annt,
                gr.peak = peak.gr,
                promoterRegion = promoter,
                annoStat = feature.count,
                peakNum = length(peak.gr)
            )
        }

        return(annoPeak)
        # return(res)
    }
}
