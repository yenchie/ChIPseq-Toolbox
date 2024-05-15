# GOseg.fun (v.1)

## date: 20240514
## author:YCW

## Description:
run GO enrichment analysis via [goseq](https://bioconductor.org/packages/goseq/) and plotting.


## Dependency:
[goseq](https://bioconductor.org/packages/goseq/)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("goseq")
```

[GO.db](https://bioconductor.org/packages/release/data/annotation/html/GO.db.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GO.db")
```

[topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("topGO")
```

[Rgraphviz](https://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rgraphviz")
```

[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
```

install.packages("dplyr")
install.packages("tidyr")
install.packages("stringr")
install.packages("ggplot2")

[jamba](https://rdrr.io/github/jmw86069/jamba/)
```R
install.packages("remotes")
remotes::install_github("jmw86069/jamba")
```

## Input:
- Path of gene ID files (suffix of file names should be .id). 
 file format: no header single column list

## Output: 

|output
|--|GOterm
|--|--|[file.name]_GOseq.enrichment.txt
|--|plot
|--|--|[file.name]_dotplot_top10_terms.pdf
|--|--|BP_[file.name]_dotplot_all_terms.pdf
|--|--|CC_[file.name]_dotplot_all_terms.pdf
|--|--|MF_[file.name]_dotplot_all_terms.pdf
|--|--|[file.name]_DAG.BP0.05.pdf
|--|--|[file.name]_DAG.CC0.05.pdf
|--|--|[file.name]_DAG.MF0.05.pdf

## Reference:
### five files need to be put on the working space:
- GOterm.txt
- Gene_Model.lengths.txt
- GO_Molecular_Function.txt
- GO_Cellular_Component.txt
- GO_Biological_Process.txt

- Gmax189
- RO_v.3 (currently)

## Usage: 

```R
GOseq.hyper.FDR.DEG.fun(datapath)
```

## Example:

```R
files.sources <- list.files("/Toolbox/GO_enrichment/GOseg.fun", ".R$", ignore.case = T, full.names = T)
sapply(files.sources, source)
datapath <- "./testing/GO_enrichment/RO_test.Data"

log_file <- file.path(datapath, "output.log")
sink(log_file, append = TRUE)
GOseq.hyper.FDR.DEG.fun(datapath)
sink()
```

## Scripts:
- GOseq.hyper.FDR.DEG.fun.v.2.R
- GOseq.ego.fun.v.1.R
- jamenrich-import.R
- jamenrich-enrichdf2er.R

## Structure:
```R
GOseq.hyper.FDR.DEG.fun(datapath)
    GO.seq.ego(pwf, GO_DB, termfile = termfile, q.cut.off)
    DAG.GOseq.fun(GOseq.result, GO_DB, GO.ontology, termfile, GOI.list, q.cut.off, outpath = outpath)
        enrichDF2enrichResult()
```


