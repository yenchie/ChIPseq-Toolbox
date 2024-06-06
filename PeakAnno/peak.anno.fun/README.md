# peakannt.annoGene.fun (v.1)

## date: 2024006
## author: YCW

## Description:
Perform peak (region) annotation using GenomicFeatures to obtain annotated genes, features, nearestTSS distances, and a Venn diagram illustrating peak overlap and gene intersection.

## Dependency:

[GenomicFeatures](https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
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

- option3: named grl object

## Output: 

```
..
|--[peak_set_ID1].annoGene.geneAnno.[date].gene.id
|--[peak_set_ID2].annoGene.geneAnno.[date].gene.id
|-- ...
|--annoGene.gene.number.[date].txt
|--annoGene.plot.[date].pdf

```

## Note:
### Reference:
#### 
depend on variables, gene.df and transcript.id.df in loaded RData
- Gmax189
'''R
load("./000/Toolbox/PeakAnno/peak.anno.fun/Gmax189.GB.Rdata",
    ignore.case = T, full.names = T
)
'''
### Parameter:
```R
annoGene(peak.gr, promoter = c(-1000, 0))
tssRegion = c(0, 1)
sameStrand = F
priority = c("promoter", "five_prime_UTR", "CDS", "three_prime_UTR", "intron", "intergenic")
```

## Usage: 

```R
annoPeakList <- peak.anno.fun(files=NULL, grl = grl, outpath = outpath, date = date, plot = F, output.gene.table = F)
# plot not veen digram not supported more than 6 peak sets

```

## Example:

```R
files.sources <- "~/000/Toolbox/PeakAnno/peak.anno.fun"
sapply(
    list.files(files.sources,
        ".R$",
        ignore.case = T, full.names = T
    ),
    source
)
load(list.files(files.sources,
    "\\.Rdata$",
    ignore.case = T, full.names = T
), envir = .GlobalEnv)
ls()

beds.path <- "./DiffBind/output/H3K4me3/q.01/DBS_bed/LIB-FULL/inBoth"
beds <- beds.path %>% list.files(full.names = T)

outpath <- beds.path %>% str_replace("DBS_bed", "DBG")
dir.create(outpath, mode="0770", recursive = T)
print(outpath)

date = "20240515"
log_file <- file.path(outpath, "output.log")
sink(log_file, append = TRUE)
annoPeakList <- peak.anno.fun(files=beds, outpath = outpath, date = date, plot = F, output.gene.table = F)
sink()
```

## Scripts:
- peakannt.annoGene.fun.v.1.R
- annoPeak.R
- annoGene.R
- plotanntBar.fun.v.1.R
- plotdistTotss.R
- beds2grl.fun.R
- GOI.genome.browser.R

## Structure:
```R
peak.anno.fun(files, outpath, date, plot = F, output.gene.table = T, Txdb=Txdb)
    beds2grl(files)
    annoGene(peak.gr, promoter = c(-1000, 0))
    plotanntBar(peakAnnoList)
    plotdistTotss(annoPeakList, tssRegion = c(0, 1), limits = c(-1, 1))


GOI.genome.browser.R(grl = grl, geneID = NULL, seqname = NULL, range = NULL)
    beds2grl(files)
    
```


-----------------------
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
files.sources <- "~/000/Toolbox/PeakAnno/peak.anno.fun"
sapply(
    list.files(files.sources,
        ".R$",
        ignore.case = T, full.names = T
    ),
    source
)
load(list.files(files.sources,
    "\\.Rdata$",
    ignore.case = T, full.names = T
))
ls()

beds.path <- "./DiffBind/output/H3K4me3/q.01/DBS_bed/LIB-FULL/inBoth"
beds <- beds.path %>% list.files(full.names = T)

outpath <- beds.path %>% str_replace("DBS_bed", "DBG")
dir.create(outpath, mode="0770", recursive = T)
print(outpath)

date = "20240515"
log_file <- file.path(outpath, "output.log")
sink(log_file, append = TRUE)
peakAnnoList <- peak.Anno.fun(files = beds, outpath, date, plot = F, output.gene.table = T, Txdb=Txdb)
sink()
```

## Scripts:
- peakannt.chipseeker.fun.v.1.R
- plotAnntBar.fun.v.1.R
- beds2grl.fun.R
- GOI.genome.browser.R

## Structure:
```R
peak.Anno.fun(files, outpath, date, plot = F, output.gene.table = T, Txdb=Txdb)
    beds2grl(files)
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


GOI.genome.browser.R(grl = grl, geneID = NULL, seqname = NULL, range = NULL)
    beds2grl(files)
    
```




