# # annotate.variables.R
# # features display ----
# # input: gff
# # output:
# ## gene.df
# ### "seqnames", "start", "end", "geneID", "transcript.ID", "transcript", "strand", "feature", "name"=c("transcript")
# ## transcript.id.df
# ### "seqnames", "start", "end", "transcript", "strand", "feature", "name"=c("transcript"), "geneID", "transcript.ID"
# # 0 base system


# gff.v1.1 <- c("/homes/yenching/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/assembly/Gmax189_gene_exons.gff3") %>%
#     read.delim(header = F, stringsAsFactors = F, skip = 1)
# head(gff.v1.1)
# gene.df <- gff.v1.1 %>%
#     filter(V3 == "mRNA") %>%
#     mutate(
#         geneID = V9 %>% str_extract("(?<=Parent=)\\w+\\.*\\d*"),
#         transcript.ID = V9 %>% str_extract("(?<=Name=)\\w+\\.*\\d*(?=;)"),
#         transcript = V9 %>% str_extract("(?<=pacid=)\\w+"),
#         seqnames = V1, start = (V4 - 1), end = V5,
#         strand = V7,
#         feature = "gene", name = "gene"
#     ) %>%
#     dplyr::select(seqnames, start, end, geneID, transcript.ID, transcript, strand, feature, name) %>%
#     distinct()
# head(gene.df)
# colnames(gene.df)
# transcript.id.df <- gff.v1.1 %>%
#     filter(V3 != "gene") %>%
#     filter(V3 != "exon") %>%
#     mutate(
#         transcript = V9 %>% str_extract("(?<=ID=PAC:)\\w+"),
#         seqnames = V1, start = (V4 - 1), end = V5,
#         strand = V7,
#         feature = V3, name = "transcript"
#     ) %>%
#     dplyr::select(seqnames, start, end, transcript, strand, feature, name) %>%
#     distinct() %>%
#     left_join(gene.df %>% dplyr::select(geneID, transcript.ID, transcript)) %>%
#     group_by(geneID) %>%
#     # mutate(alter = Name %>% str_extract("\\d+$") %>% as.numeric()) %>%
#     ungroup()
# head(transcript.id.df)
# colnames(transcript.id.df)
# save(gene.df, transcript.id.df, file = "./000/Toolbox/PeakAnno/Gmax189.GB.Rdata")


# library(tidyverse)
# gff.a6 <- c("/bcst/JYL/JYL_qnap/db/Gm/GmaxWm82.a6.v1_Gmax880/annotation/Gmax_880_Wm82.a6.v1.gene_exons.gff3") %>%
#     read.delim(header = F, stringsAsFactors = F, skip = 6)
# head(gff.a6, n=10)
# unique(gff.a6$V3)
# gene.df <- gff.a6 %>%
#     filter(V3 == "mRNA") %>%
#     mutate(
#         geneID = V9 %>% str_extract("(?<=Parent=).+(?=(.Wm82.a6.v1))"),
#         transcript.ID = V9 %>% str_extract("(?<=Name=).+(?=;pacid=)"),
#         transcript = V9 %>% str_extract("(?<=pacid=)\\w+"),
#         seqnames = V1, start = (V4 - 1), end = V5,
#         strand = V7,
#         feature = "gene", name = "gene"
#     ) %>%
#     dplyr::select(seqnames, start, end, geneID, transcript.ID, transcript, strand, feature, name) %>%
#     distinct()
# head(gene.df)
# colnames(gene.df)
# transcript.id.df <- gff.a6 %>%
#     filter(V3 != "gene") %>%
#     filter(V3 != "CDS") %>%
#     mutate(
#         transcript = V9 %>% str_extract("(?<=pacid=)\\w+"),
#         seqnames = V1, start = (V4 - 1), end = V5,
#         strand = V7,
#         # feature = V3, 
#         feature =  ifelse(V3 =="exon", V9 %>% str_extract("(?<=Wm82.a6.v1.).+(?=;Parent=)"), V3 ),
#         name = "transcript"
#     ) %>%
#     dplyr::select(seqnames, start, end, transcript, strand, feature, name) %>%
#     distinct() %>%
#     left_join(gene.df %>% dplyr::select(geneID, transcript.ID, transcript)) 

# head(transcript.id.df,  n=10)%>%as.data.frame()
# unique(transcript.id.df$feature)
# colnames(transcript.id.df)
# nrow(transcript.id.df)
# unique(transcript.id.df$geneID)%>% length() #48387
# unique(transcript.id.df$transcript.ID)%>% length() #80374
# # transcript.id.df%>% filter(geneID =="Glyma.04G182900")%>%as.data.frame()%>% tail()

# save(gene.df, transcript.id.df, file = "./000/Toolbox/PeakAnno/Gmax880.GB.Rdata")
