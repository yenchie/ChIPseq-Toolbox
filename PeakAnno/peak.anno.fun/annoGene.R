## update: 2024.11.01


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
        head(ref)

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
                    distinct() %>%
                    mutate(feature = ifelse(feature %>% str_detect("exon.1"), "first_exon",
                        feature %>% str_remove(".\\d+$")
                    ))
            ) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = T)

        tss.annt <- ref %>%
            filter(feature == "mRNA") %>%
            mutate(
                tss = ifelse(strand == "+", start, end)
            ) %>%
            dplyr::select(seqnames, tss, strand, geneID, transcript.ID) %>%
            distinct() %>%
            mutate(start = tss, end = tss + 1) %>%
            dplyr::select(-tss) %>%
            mutate(feature = "TSS") %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = T)
        head(tss.annt)


        head(ref.annt)
        unique(ref.annt$feature) %>% print()

        #ref.annt<- c(ref.annt, tss.annt)
        hits <- findOverlaps(ref.annt, peak.gr)
        hits.df <- hits %>%
            as.data.frame() %>%
            mutate(across(everything(), as.character))

        head(hits.df)
        nrow(hits.df)
        # p <- Pairs(ref.annt, peak.gr, hits = hits)

        if (nrow(hits.df) != 0) {
            nrow(hits.df)
            priority <- c("promoter", "five_prime_UTR", "first_exon", "exon", "three_prime_UTR", "mRNA", "intergenic")
            res <- hits.df %>%
                left_join(peak.gr %>% as.data.frame() %>% mutate(start = start - 1) %>% rownames_to_column(var = "subjectHits")) %>%
                left_join(ref.annt %>% as.data.frame() %>% mutate(start = start - 1) %>% rownames_to_column(var = "queryHits") %>%
                    dplyr::rename(gene.start = start, gene.end = end, Strand = strand, gene.width = width)) %>%
                dplyr::select(-c(queryHits, subjectHits)) %>%
                mutate(peakID = paste(seqnames, start, end, sep = "__")) %>% 
                group_by(peakID, feature, transcript.ID) %>%
                mutate(
                    start_distance = ifelse(
                        feature == "promoter", 
                            ifelse(
                                Strand == "+", start - (gene.end+1), (start - (gene.start-1))
                            ),
                            ifelse( Strand == "+", start - (gene.start), (start - (gene.end)))),
                    end_distance = ifelse(
                        feature == "promoter", 
                            ifelse(
                                Strand == "+", end - (gene.end+1), (end - (gene.start-1))
                            ),
                            ifelse(Strand == "+", end - (gene.start), (end - (gene.end)))),
                    distance_to_TSS = ifelse(feature == "mRNA",
                        ifelse((start_distance < 0 & end_distance > 0) |
                            (start_distance > 0 & end_distance < 0),
                        0, min(abs(start_distance), abs(end_distance))
                        ), 
                        ifelse(feature == "promoter", 
                            ifelse((start_distance < 0 & end_distance > 0) |
                            (start_distance > 0 & end_distance < 0),
                            0, -min(abs(start_distance), abs(end_distance))
                        ), NA)
                    )
                ) %>%
                ungroup() %>%
                group_by(peakID) %>%
                arrange(factor(feature, levels = priority))%>%
                mutate(priority_mark = ifelse(row_number() == 1, 1, seq_along(feature))) %>%
                mutate(feature.1 = ifelse(all(feature == "mRNA"), "intron", feature)) %>% ## only apear in intron across all transcripts annoted
                ungroup()  %>%
                dplyr::select(-c(start_distance, end_distance))

            head(res %>% as.data.frame())
            tail(res %>% as.data.frame())
            head(res %>% filter(feature== "mRNA")%>% as.data.frame())
            tail(res %>% filter(feature== "mRNA")%>% as.data.frame())
            nrow(res)

            unique(res$feature) %>% print()
            unique(res$feature.1) %>% print()

            unalign <- peak.gr %>%
                as.data.frame() %>%
                mutate(start = start - 1) %>%
                rownames_to_column(var = "subjectHits") %>%
                filter(!subjectHits %in% hits.df$subjectHits) %>%
                distinct() %>%
                mutate(peakID = paste(seqnames, start, end, sep = "__"))
            head(unalign)
            nrow(unalign)

            priority <- c("promoter", "five_prime_UTR", "first_exon", "exon", "three_prime_UTR", "intron", "intergenic")

            feature.count <- data.frame(
                feature.1 = "intergenic",
                n = length(peak.gr) - length(hits.df$subjectHits %>% unique())
            ) %>%
                bind_rows(res %>%
                    filter(priority_mark == 1) %>%
                    dplyr::count(feature.1)) %>%
                mutate(proportion = ((n / sum(n)) %>% round(3)) * 100) %>%
                dplyr::select(-n) %>%
                mutate(feature.1 = feature.1 %>% factor(levels = priority)) %>%
                arrange(feature.1)
            print(feature.count)
            genes <- unique(res$geneID)

            annoPeak <- new("annoPeak",
                annoGene = genes,
                gr.ref = ref.annt,
                gr.peak = peak.gr,
                gr.ref.tss = tss.annt,
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
                gr.ref.tss = tss.annt,
                promoterRegion = promoter,
                annoStat = feature.count,
                peakNum = length(peak.gr)
            )
        }

        return(annoPeak)
        # return(res)
    }
}
