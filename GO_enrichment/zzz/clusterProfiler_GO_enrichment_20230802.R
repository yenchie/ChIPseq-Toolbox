# clusterProfiler_GO_enrichment_20230802.R

# Load required libraries
# packageurl <-
#   "https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz"
# install.packages(packageurl, repos = NULL, type = "source")
# remove.packages("clusterProfiler")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")

library(clusterProfiler); library(enrichplot); library(dplyr); library(stringr); library(tidyr); library(GO.db); library(ggplot2)

rm(list = ls())

#read all the table. 
featurepath<- "/bcst/JYL/JYL_qnap/Project/Gm/ChIPseq_analysis/Gmax189/Methylation/Analysis/Other_Analysis/Bivalent_H3K4me3_H3K27me3/GO_enrichment/input/feature"
all.genes<-read.table(file.path(featurepath, "Gmax_189_Gene_Model.lengths.txt"), sep="\t",stringsAsFactors=FALSE, header=FALSE)
GO_BP<-read.table(file.path(featurepath, "Gmax.soybase.GO_Biological_Process.txt"), sep="\t",stringsAsFactors=FALSE)%>% dplyr::rename(GO=V2, gene=V1) %>% mutate(Ontology= "GO_BP")
GO_MF<-read.table(file.path(featurepath, "Gmax.soybase.GO_Molecular_Function.txt"), sep="\t",stringsAsFactors=FALSE)%>% dplyr::rename(GO=V2, gene=V1) %>% mutate(Ontology= "GO_MF")
GO_CC<-read.table(file.path(featurepath, "Gmax.soybase.GO_Cellular_Component.txt"), sep="\t",stringsAsFactors=FALSE)%>% dplyr::rename(GO=V2, gene=V1) %>% mutate(Ontology= "GO_CC")
termfile<-read.table(file.path(featurepath, "GOterm.txt"), sep ="\t", stringsAsFactors=FALSE,header=FALSE,fill=TRUE,quote="")%>% dplyr::rename(GO=V1, descritpion=V2, taq=V3)

GO.all<- bind_rows(GO_BP, GO_MF) %>% bind_rows(GO_CC)
GO.all$Ontology<-GO.all$Ontology%>%as.character()
unique(GO.all$GO)%>%length()
unique(termfile$GO)%>%length()
# query GO term description


query.GO<-GO.all%>% left_join(termfile) %>% dplyr::select(-c(taq, V4)) %>% 
  filter(descritpion%>%is.na())%>% dplyr::select(GO)%>%distinct()
term=list()
for(go in query.GO$GO){
  
  if(GOTERM[[go]]%>%is.null()){
    print(paste("non description", go))
    cat("--------------------------------------\n")
  }else{
    term[[go]]<-GOTERM[[go]]
    print(GOTERM[[go]])
    cat("--------------------------------------\n")
  }
}

termfile.1<- data.frame(GO=sapply(term, function(x) x@GOID ),
                        descritpion = sapply(term, function(x) x@Term))
termfile<- termfile%>%bind_rows(termfile.1)
GO.all <- GO.all%>% left_join(termfile) %>% dplyr::select(-c(taq, V4))

# write.table(GO.all, file.path(featurepath, "Gmax.soybase.GO_all_20230802.txt"), sep="\t", quote = F, row.names = F, col.names = F)

# no description: GO:0090450
# provided by a data.frame of GO (column 1) and gene (column 2) direct annotation this function will building gene to GO and GO to gene mapping, with directly and undirectly (ancestor GO term) annotation.
# TERM2GENE<- buildGOmap((GO_BP%>%dplyr::select(GO, gene)) )
# class(TERM2GENE)

data.path<- "/bcst/JYL/JYL_qnap/Project/Gm/ChIPseq_analysis/Gmax189/Methylation/Analysis/DiffBind/Data/DEG/GO_enrichment"
filelist=list.files( file.path(data.path),pattern="_gene_id|_TF_id", full.names = T)
#gene.id<- read.delim(filelist[1], stringsAsFactors = F, header = F)$V1%>%as.vector()



for (file in 1:length(filelist)) {
  gene.id<- read.delim(filelist[file], stringsAsFactors = F, header = F)$V1%>%as.vector()
  name<-filelist[file]%>%as.character%>%basename()%>%str_remove("_id")
  
  for (go in unique(GO.all$Ontology)){
    print(paste("start:",  name, go))
    
    TERM2GENE<- GO.all%>%
      dplyr::filter(Ontology==go) 
    
    print(which(gene.id%in% TERM2GENE$gene)%>%length())
    count(TERM2GENE, GO)
    n_distinct(TERM2GENE$gene)
    ego <- enricher(
      gene.id ,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      # universe	background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
      universe = NULL, 
      qvalueCutoff = 0.05,
      minGSSize=5,
      TERM2GENE = TERM2GENE[c("GO", "gene")],
      TERM2NAME = TERM2GENE[c("GO", "descritpion")]
    )
    ego.t <- ego%>%as.data.frame() %>% mutate(Ontology = go)
    # GeneRatio = genes of interest in the gene set / total genes of interest.  
    #             the fraction of differentially expressed genes found in the gene set.
    # BgRatio = size of the geneset / size of all of the unique genes in the collection of genesets
    
    p1<-dotplot(ego, color = "qvalue", showCategory = length(ego.t$ID%>%unique()),  
                font.size = 12, title = paste( name, go, sep="_")) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=35)) 
    
    if(length(ego.t$ID%>%unique())>0){
      p2<-dotplot(ego, color = "qvalue", showCategory = 20,  
                  font.size = 12, title = paste( name, go, sep="_")) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=35)) 
    }else{
      print(paste("no GO term enriched in",  name, go))
      p2<-dotplot(ego, color = "qvalue", showCategory = length(ego.t$ID%>%unique()),
                  font.size = 12, title = paste( name, go, sep="_")) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=35)) 
    }
    
    
    
    if(length(ego.t$ID%>%unique())<20){
      zoom<- 1
    }else{
      zoom<- length(ego.t$ID%>%unique())/20
    }
   
    print(zoom)
    print(paste("output:", file.path(data.path, "output", "plot"), name, go))
    
    ggsave(
      filename = paste( name, go,"dotplot_all_terms.pdf", sep="_"),
      plot = p1,
      path = file.path(data.path, "output", "plot"),
      width = 10,
      height = 16*zoom, 
      limitsize = FALSE
      
    )
    ggsave(
      filename = paste( name, go,"dotplot_top20_qval_terms.pdf", sep="_"),
      plot = p2,
      path = file.path(data.path, "output", "plot"),
      width = 10,
      height = 12
      
    )
    write.table(ego.t, file.path(data.path, "output", "GOterm", paste( name, go,"clusterprofile.enrichment.txt", sep="_")), sep="\t", quote = F, row.names = F, col.names = T, eol="\n")

  }
 
}

source("./YCWang/Rscript/1_Function/CMD_GO_USE_THIS_ONE/clusterProfiler_GO_enrichment.function.v1.R")
datapath<-"/bcst/JYL/YCWang//testing/GO_enrichment/test_data"
GO_enricher(datapath, datapath, 0.05)



