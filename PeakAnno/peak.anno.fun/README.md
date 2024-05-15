# peak.anno.fun (v.1)

## date: 20240515
## author: YCW

## Description:
Perform peak (region) annotation using Chipseeker to obtain annotated genes, features, TSS-relative distances, and a Venn diagram illustrating peak overlap and gene intersection.

## Dependency:

[GenomicFeatures](https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
```
[ChIPseeker](https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPseeker")
```
[DOSE](https://bioconductor.org/packages/release/bioc/html/DOSE.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DOSE")
```
```R
install.packages("dplyr") 
install.packages("stringr") 
install.packages("tibble") 
install.packages("forcats") 
install.packages("ggplot2") 
install.packages("ggpubr")
```

## Input:
- option1: Character vector containing file paths of peak files. Peak set ID will use file names as labels and output gene id files (if chose).
e.g. 
```R
beds<- c("path1", "path2", ..., "pathN")
```

- option2: Table storing file paths and peak set IDs: Column 1 represents paths, and column 2 contains peak set IDs as desired which will be used as labels and output gene id files (if chose)

e.g. 
|path|ID|
|---|---|
|./path1.bed|Apeakset|

## Output: 

```
..
|--[peak_set_ID1]ChIPseeker.geneAnno.[date].gene.id
|--[peak_set_ID2]ChIPseeker.geneAnno.[date].gene.id
|-- ...
|--ChIPseeker.gene.number.[date].txt
|--ChIPseeker.plot.[date].pdf

```

## Note:
### Reference:
#### 
- Gmax189

### Parameter:
```R
promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 0)
options(ChIPseeker.downstreamDistance = 0)
tssRegion = c(-1000, -1)
level = "gene", sameStrand = F

```

## Usage: 

```R
peakAnnoList <-peak.Anno.fun(files, outpath, date, plot = F, output.gene.table = T)

# plot not veen digram not supported more than 6 peak sets

```

## Example:

```R
files.sources <- list.files("./Toolbox/PeakAnno/peak.anno.fun", ".R$",
    ignore.case = T, full.names = T
)
sapply(files.sources, source)

ls()
beds.path <- "./DiffBind/output/H3K4me3/q.01/DBS_bed/LIB-FULL/inBoth"
beds <- beds.path %>% list.files(full.names = T)

outpath <- beds.path %>% str_replace("DBS_bed", "DBG")
dir.create(outpath, mode="0770", recursive = T)
print(outpath)

date = "20240515"
log_file <- file.path(outpath, "output.log")
sink(log_file, append = TRUE)
peakAnnoList <- peak.Anno.fun(files = beds, outpath, date, plot = F, output.gene.table = T)
sink()
```

## Scripts:
- peakannt.chipseeker.fun.v.1.R
- plotAnntBar.fun.v.1.R

## Structure:
```R
peak.Anno.fun(files, outpath, date, plot = F, output.gene.table = T)
    ChIPseeker::annotatePeak(
        TxDb = txdb,
        tssRegion = c(-1000, -1), # signed your definition for TSS
        verbose = TRUE,
        assignGenomicAnnotation = TRUE,
        genomicAnnotationPriority = c(
        "Promoter", "5UTR", "3UTR", "Exon", "Intron",
        "Downstream", "Intergenic"),
        level = "gene", sameStrand = F)

    plotAnntBar(peakAnnoList)

```


