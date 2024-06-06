# plotanntBar.R
# author: YCW
# date: 2024.06.05
###
# for annoPeak Object


plotanntBar <- function(annoPeakList) {
    res.peak.annt <- NULL
    for (set in names(annoPeakList)) {
        annot.1 <- as.data.frame(annoPeakList[[set]]@annoStat) %>%
            mutate(peaksetID = set)

        res.peak.annt <- res.peak.annt %>% bind_rows(annot.1)
    }

    level <- c("promoter", "five_prime_UTR", "CDS", "intron", "three_prime_UTR", "intergenic")
    col <- c(
        "#A6CEE3",
        "#1F78B4", "#33A02C", "#E31A1C",
        "#FF7F00", "#CAB2D6"
    )

    col <- setNames(col, level)

    res.peak.annt %>%
        mutate(peaksetID = peaksetID %>% factor(levels = c(names(annoPeakList) %>% rev()))) %>%
        mutate(feature = feature %>% factor(levels = level %>% rev())) %>%
        arrange(feature) %>%
        arrange(peaksetID) %>%
        ggplot(aes(x = peaksetID, y = proportion, fill = feature, order = as.numeric(feature))) +
        geom_col() +
        ggtitle("Feature Distribution") +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12, face = "bold")
        ) +
        scale_fill_manual(values = col %>% rev()) +
        coord_flip() +
        guides(fill = guide_legend(reverse = TRUE))
}
# fin -----
