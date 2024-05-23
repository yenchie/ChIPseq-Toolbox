# plotAnntBar.fun.v.1.R
# author: YCW
# date: 2024.03.08
###
# revised date:
# - plot bar plot for percentage of peak feature Distribution


plotAnntBar <- function(peakAnnoList) {
  res.peak.annt <- NULL
  for (set in names(peakAnnoList)) {
    annot.1 <- as.data.frame(peakAnnoList[[set]]@anno)
    unique(annot.1$annotation)
    annot.2 <- annot.1 %>%
      mutate(
        annt =
          ifelse(
            annotation %>% str_detect("Exon"),
            ifelse(annotation %>% str_detect("exon 1 of"),
              "1st Exon",
              "Other Exon"
            ),
            ifelse(
              annotation %>% str_detect("Intron"),
              ifelse(
                annotation %>% str_detect("intron 1 of"),
                "1st Intron",
                "Other Intron"
              ),
              annotation
            )
          )
      ) %>%
      dplyr::count(annt) %>%
      mutate(
        percentage = (n / sum(n)) * 100,
        peaksetID = set
      )
    res.peak.annt <- res.peak.annt %>% bind_rows(annot.2)
  }

  promoter.annt <- c(grep("Promoter", unique(res.peak.annt$annt), value = T))
  if (promoter.annt %>% str_detect("\\d")) {
    sorted_promoter.annt <- promoter.annt[order(as.numeric(gsub("\\D", "", promoter.annt)), decreasing = T)]
    col <- data.frame(
      Feature = c(
        sorted_promoter.annt,
        "5' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "3' UTR",
        "Downstream (<=300bp)", "Distal Intergenic"
      ),
      colors = c(
        brewer.pal(n = length(sorted_promoter.annt) - 1, name = "Blues") %>% rev(), "#A6CEE3",
        "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
        "#FF7F00", "#CAB2D6"
      )
    )
  } else {
    sorted_promoter.annt <- promoter.annt
    col <- data.frame(
      Feature = c(
        sorted_promoter.annt,
        "5' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "3' UTR",
        "Downstream (<=300bp)", "Distal Intergenic"
      ),
      colors = c(
        "#A6CEE3",
        "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
        "#FF7F00", "#CAB2D6"
      )
    )
  }

  col.1 <- col %>% filter(Feature %in% unique(res.peak.annt$annt))
  res.peak.annt %>%
    mutate(peaksetID = peaksetID %>% factor(levels = c(names(peakAnnoList) %>% rev()))) %>%
    mutate(annt = annt %>% factor(levels = col$Feature %>% rev())) %>%
    arrange(annt) %>%
    arrange(peaksetID) %>%
    ggplot(aes(x = peaksetID, y = percentage, fill = annt, order = as.numeric(annt))) +
    geom_col() +
    ggtitle("Feature Distribution") +
    xlab("") +
    ylab("Percentage (%)") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    scale_fill_manual(values = col.1$colors %>% rev()) +
    coord_flip() +
    guides(fill = guide_legend(reverse = TRUE))
}
# fin -----
