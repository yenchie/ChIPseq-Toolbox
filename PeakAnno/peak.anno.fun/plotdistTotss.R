# plotdistTotss.R
# author: YCW
# date: 2024.06.05
###
# for annoPeak Object


plotdistTotss <- function(annoPeakList, tssRegion = c(0, 1), limits = c(-1, 1)) {
    print(tssRegion)
    res.peak.annt <- NULL
    for (set in names(annoPeakList)) {
        print(set)
        annoPeak.1 <- annoPeakList[[set]]

        peak.gr <- annoPeak.1@gr.peak
        gr.ref <- annoPeak.1@gr.ref
        gr.ref <- gr.ref[(elementMetadata(gr.ref)$feature == "mRNA")]
        # elementMetadata(gr.ref)$geneID %>%
        #     unique() %>%
        #     length()

        gr.tss <- flank(gr.ref, width = -1, start = TRUE)

        peaks_df <- as.data.frame(peak.gr) %>% rowid_to_column(var = "queryHits")
        genes_df <- as.data.frame(gr.tss) %>%
            rowid_to_column(var = "subjectHits")
        head(genes_df)
        head(peaks_df)

        distances <- distanceToNearest(peak.gr, gr.tss, select = "arbitrary", ignore.strand = F) %>%
            as.data.frame() %>%
            left_join(peaks_df, by = c("queryHits")) %>%
            left_join(genes_df, by = c("subjectHits")) %>%
            mutate(distance = ifelse(distance == 0, distance,
                ifelse(strand.y == "+",
                    ifelse(start.y < start.x, distance, -distance),
                    ifelse(end.y > end.x, distance, -distance)
                )
            )) %>%
            dplyr::select(queryHits, distance) %>%
            mutate(peaksetID = set) %>%
            distinct()

        head(distances)
        nrow(distances)

        res.peak.annt <- res.peak.annt %>% bind_rows(distances)
    }

    head(res.peak.annt)
    nrow(res.peak.annt)

    bin_width <- 0.01
    res.peak.annt.1 <- res.peak.annt %>%
        mutate(
            distance_kb = distance / 10^3,
            bin = floor(distance_kb / bin_width) * bin_width
        ) %>%
        group_by(bin, peaksetID) %>%
        mutate(bin_count = n()) %>%
        ungroup() %>%
        dplyr::select(bin, bin_count, peaksetID) %>%
        distinct() %>%
        group_by(peaksetID) %>%
        mutate(proportion = bin_count / sum(bin_count)) %>%
        ungroup()

    res.peak.annt.1 %>%
        ggplot() +
        geom_col(aes(x = bin, y = proportion, fill = peaksetID), position = "identity", alpha = 0.5) +
        scale_x_continuous(limits = limits) +
        ggtitle("distance to TSS") +
        xlab("distance to nearest TSS(kb) (bin: 10b)") +
        ylab("percentage (per all)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12, face = "bold")
        ) +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_y_sqrt()

    # res.peak.annt %>%
    #     mutate(peaksetID = peaksetID %>% factor(levels = c(names(annoPeakList) %>% rev()))) %>%
    #     arrange(peaksetID) %>%
    #     ggplot() +
    #     # geom_density(aes(distance / 10^3, ..scaled.., color = peaksetID), position = "identity") +
    #     geom_histogram(aes(distance / 10^3, y = (..count..) / sum(..count..), fill = peaksetID), bins = 100, alpha = 0.8, position = "identity") +
    #     ggtitle("distance to TSS") +
    #     xlab("distance to nearest TSS (Kb)") +
    #     ylab("density") +
    #     scale_x_continuous(limits = c(-2.5, 2.5)) +
    #     # scale_y_continuous(trans = "log10") +
    #     theme_bw() +
    #     theme(
    #         axis.text = element_text(size = 12, face = "bold"),
    #         axis.title = element_text(size = 12, face = "bold")
    #     ) +
    #     guides(fill = guide_legend(reverse = TRUE)) +
    #     scale_y_sqrt()
}
# fin -----
