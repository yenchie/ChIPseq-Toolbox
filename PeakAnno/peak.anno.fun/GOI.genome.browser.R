# author: YCW
# date:2024.05.23

# dependency:
# beds2grl.fun.R
# Gmax189.GB.Rdata

# USAGE:
# GOI.genome.browser(grl = grl, geneID = NULL, seqname = NULL, range = NULL)


## features display ----
# gff.v1.1 <- c("./Gmax_189_gene_exons.gff3.gz") %>%
#     read.delim(header = F, stringsAsFactors = F, skip = 1)
#  gene.df <- gff.v1.1 %>%
#      filter(V3 == "mRNA") %>%
#      mutate(
#          geneID = V9 %>% str_extract("(?<=Parent=)\\w+\\.*\\d*"),
#          Name = V9 %>% str_extract("(?<=Name=)\\w+\\.*\\d*(?=;)"),
#          transcript.ID = V9 %>% str_extract("(?<=pacid=)\\w+"),
#          seqnames = V1, start = (V4 - 1), end = V5,
#          strand = ifelse(V7 == "+", ">", "<"),
#          feature = "gene", name = "gene"
#      ) %>%
#      dplyr::select(seqnames, start, end, geneID, Name, transcript.ID, strand, feature, name) %>%
#      distinct()

# transcript.id.df <- gff.v1.1 %>%
#     filter(V3 != "gene") %>%
#     filter(V3 != "exon") %>%
#     mutate(
#         transcript.ID = V9 %>% str_extract("(?<=ID=PAC:)\\w+"),
#         seqnames = V1, start = (V4 - 1), end = V5,
#         strand = ifelse(V7 == "+", ">", "<"),
#         feature = V3, name = "transcript"
#     ) %>%
#     dplyr::select(seqnames, start, end, transcript.ID, strand, feature, name) %>%
#     distinct() %>%
#     left_join(gene.df %>% dplyr::select(geneID, Name, transcript.ID)) %>%
#     group_by(geneID) %>%
#     mutate(alter = Name %>% str_extract("\\d+$") %>% as.numeric()) %>%
#     ungroup()
# save(gene.df, transcript.id.df, file = "./000/Toolbox/PeakAnno/Gmax189.GB.Rdata")


GOI.genome.browser <- function(grl = grl, geneID = NULL, seqname = NULL, range = NULL) {
    res <- NULL
    for (x in names(grl)) {
        d1 <- grl[[x]] %>%
            as.data.frame() %>%
            mutate(
                name = x, feature = "peak"
            )
        res <- res %>% bind_rows(d1)
    }
    head(res) %>% print()
    unique(res$name) %>% print()

    if (!geneID %>% is.null()) {
        seqname <- gene.df$seqnames[which(gene.df$geneID == geneID)] %>% unique()
        range <- c(gene.df$start[which(gene.df$geneID == geneID)], gene.df$end[which(gene.df$geneID == geneID)])
        minx <- min(range) / 10^6 - 0.003
        maxx <- max(range) / 10^6 + 0.003
    } else {
        seqname <- seqname
        range <- range
        minx <- min(range) / 10^6
        maxx <- max(range) / 10^6
    }
    print(paste("display region:", seqname, paste(range, collapse = "-")))


    p.d1 <- res %>%
        bind_rows(gene.df) %>%
        bind_rows(transcript.id.df) %>%
        filter(seqnames %in% seqname) %>%
        filter((minx * 10^6 < start) & (end < maxx * 10^6))

    if (all(unique(res$name) %in% unique(p.d1$name))) {
        print("all peakset display in the region")
    } else {
        p.d1 <- data.frame(
            name =
                unique(res$name)[which(!unique(res$name) %in% unique(p.d1$name))]
        ) %>%
            mutate(feature = "peak") %>%
            bind_rows(p.d1)
    }

    peakset.plot.df <- p.d1 %>%
        rownames_to_column("order") %>%
        mutate(o.ID = ifelse(name == "transcript",
            paste(transcript.ID, seqnames, start, end, sep = "__"),
            ifelse(name == "gene", paste(seqnames, start, end, sep = "__"),
                paste(order, seqnames, start, end, sep = "__")
            )
        )) %>%
        dplyr::select(seqnames, start, end, strand, geneID, transcript.ID, feature, name, Name, alter, o.ID) %>%
        distinct() %>%
        mutate(label.coord = (start + end + 1) / 2) %>%
        as.data.frame() %>%
        mutate(
            geneID = ifelse(feature == "gene", geneID,
                ifelse(feature == "mRNA", Name, NA)
            ),
            alter = ifelse(name == "transcript",
                paste("transcript", alter, sep = "."), name
            ),
            strand = ifelse(feature == "peak",
                NA, strand
            ),
            feature = ifelse(feature == "peak",
                name, ifelse(feature == "mRNA", "intron", feature)
            )
        ) %>%
        pivot_longer(
            cols = c(start, end),
            values_to = "coord.",
            names_to = "type"
        )
    tail(peakset.plot.df)
    nrow(peakset.plot.df)

    line.width <- data.frame(
        feature = c(names(grl), "gene", "five_prime_UTR", "CDS", "three_prime_UTR", "intron"),
        line.width = c(rep(10, length(grl)), 10, 5, 10, 5, 2),
        colors = c(rainbow(n = length(grl)), "blue", "#417041", "#0031c3", "#10a410", "gray")
    )

    line.width <- line.width %>% filter(feature %in% unique(peakset.plot.df$feature))

    p <- peakset.plot.df %>%
        mutate(
            alter = alter %>%
                factor(levels = c(unique(peakset.plot.df$alter) %>% rev())),
            name = name %>% factor(levels = unique(peakset.plot.df$name)),
            feature = feature %>% factor(levels = c(line.width$feature) %>% unique())
        ) %>%
        mutate(
            facet = ifelse(name == "transcript", "transcript", "peak")
        ) %>%
        {
            ggplot(.) +
                geom_line(
                    aes(
                        x = (coord. / 10^6), y = alter, group = o.ID, color = feature,
                        linewidth = feature
                    ),
                    alpha = 1, na.rm = F
                ) +
                geom_text(aes(group = o.ID, x = (label.coord / 10^6), y = alter, label = geneID),
                    vjust = 2.5, position = "identity", check_overlap = T, size = 3
                ) +
                geom_text(aes(group = o.ID, x = (label.coord / 10^6), y = alter, label = strand), color = "white", size = 3) +
                ylab("peakset") +
                xlab("position (Mbase)") +
                scale_x_continuous(
                    labels = scales::number_format(accuracy = 0.01),
                    limits = c(minx, maxx + (10^-3))
                ) +
                scale_discrete_manual(
                    aesthetic = "linewidth",
                    values = line.width$line.width
                ) +
                scale_color_manual(
                    values = line.width$colors
                ) +
                theme(
                    axis.text.x = element_text(size = 20, face = "bold"),
                    axis.text.y = element_text(size = 16, face = "bold"),
                    axis.title = element_text(size = 20, face = "bold"),
                    legend.position = "bottom",
                    legend.text = element_text(size = 14, face = "bold"),
                    legend.title = element_text(size = 14, face = "bold"),
                    strip.background = element_rect(colour = "white", fill = "white"),
                    strip.text = element_blank(),
                    panel.background = element_rect(
                        fill = "white",
                        colour = "white",
                        size = 0.5, linetype = "solid"
                    ),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_line(
                        size = 0.25, linetype = "solid",
                        colour = "gray"
                    ),
                    plot.background = element_rect(fill = "#eceae1")
                ) +
                facet_grid(facet ~ ., scales = "free_y", space = "free_y")
        }
    return(p)
}
## fin---
