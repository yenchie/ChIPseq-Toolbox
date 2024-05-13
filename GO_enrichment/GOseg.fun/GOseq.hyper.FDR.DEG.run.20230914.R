source("/nas/cluster_data/homes/yenching/JYL/JYL_qnap_2/YCWang/0_Script/000/Toolbox/GO_enrichment/GOseq.hyper.FDR.DEG.fun.v.2.R")
# datapath<-"/bcst/JYL/JYL_qnap_2/YCWang/Project/Gm/RNAseq/Gm_seed_develop/GO_enrichment/AdjUp_mm_306TF_id"
datapath <- "/homes/yenching/testing/GO_enrichment/RO_test.Data"
# name="/bcst/JYL/JYL_qnap_2/YCWang/Project/Gm/RNAseq/Gm_seed_develop/GO_enrichment/adjUP/output/GOterm/adjUp_sdlg_GO0006355_regulation_of_transcription.TF_GOseq.enrichment.txt"


log_file <- file.path(datapath, "output.log")
sink(log_file, append = TRUE)
GOseq.hyper.FDR.DEG.fun(datapath)
sink()
