# peak.Venn.20240422.R
# YCW
# 2024.04.22


library(ChIPpeakAnno)
library(tidyverse)
path <- "/nas/qnap4_1/Project/Gm/ChIPseq_analysis/Gmax189/Analysis/Methylation/Peakcalling/output/PE/K4K27"
beds <- path %>%
    list.dirs() %>%
    list.files("sorted.narrowPeak", full.names = T)

peaks <- map_df(beds, ~ {
    d1 <- read.delim(.x, stringsAsFactors = F, header = F) %>%
        mutate(IP = paste0(.x %>% as.character() %>% basename() %>% str_extract(paste(c("CG", "CHH", "CHG"), collapse = "|"))))
    d1
}) %>%
    dplyr::rename(seqnames = V1, start = V2, end = V3)
grl <- GRangesList()

for (ip in unique(d1$IP)) {
    print(ip)
    gr <- d1 %>% filter(IP == ip)
    gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE, starts.in.df.are.0based = T)
    print(gr)
    ID <- paste(s, ip, sep = ".")
    grl[[ID]] <- gr
}
grl[["validated"]] <- validated.bi.peaks
names(grl)
promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 0) # signed your definition of promoter
options(ChIPseeker.downstreamDistance = 0)
# 0base
grl1 <- NULL
grl1 <- grl
head(grl1)
makeVennDiagram(grl1,
    NameOfPeaks = c(names(grl1)),
    scaled = FALSE, euler.d = FALSE, totalTest = lengths(grl1) %>% max(),
    connectedPeaks = "keepAll",
    fill = c("#009E73", "#F0E442", "#D55E00", "#0072B2"), # circle fill color
    col = c("#D55E00", "#0072B2", "#009E73", "#F0E442"), # circle border color
    cat.col = c("#D55E00", "#0072B2", "#009E73", "#F0E442")
)
