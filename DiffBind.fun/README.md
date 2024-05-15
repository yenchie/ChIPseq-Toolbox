# DiffBind.fun (v.1)

## date: 20240515
## author: YCW

## Description:
Do DiffBind and output results

## Dependency:

[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DiffBind")
```

install.packages("dplyr")
install.packages("stringr")
install.packages("tidyr")
install.packages("purrr")
install.packages("tibble")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("scales")
install.packages("grid")
install.packages("reshape2")

## Input:
- path of sample sheet
formate:

|SampleID|bamReads|bamControl|Factor|Condition|ControlID|Tissue|Treatment|Replicate|Peaks|PeakCaller|
|--------|--------|----------|------|---------|---------|------|---------|---------|-----|----------|
|ABBCC   |bam.path|ctr.bam.path|NA  |NA       |ABBCC.ctr|NA    |NA       |1        |peak.path|macs/narrow/bed|


## Output: 

```
..
|--output.log
|--consensus.peak.bed
|--dba_DBA.RData
|--dba.count.score.value.boxplot.pdf
|--DBdata.count.score.cluster.pdf
|--DBdata.count.libsizes.txt
|--DBdata.count.score.global_binding_matrix.txt
|--dba.normalized.value.boxplot.pdf
|--DBdata.normalized.cluster.LIB-FULL.pdf
|--pair-wise_peak.number.txt
|--DBS_matrix
    |--LIB-FULL
        |--LIB-FULL_DEseq2_q0pt05_FC2_pair_wise_comparison_matrix.pdf
        |--LIB-FULL_edgeR_q0pt05_FC2_pair_wise_comparison_matrix.pdf
        |--LIB-FULL_inBoth_q0pt05_FC2_pair_wise_comparison_matrix.pdf
        |--pair_wise_matrix_q0pt05.FC2.DEseq2.txt
        |--pair_wise_matrix_q0pt05.FC2.edgeR.txt
        |--pair_wise_matrix_q0pt05.FC2.inBoth.txt
|--DBS_table
    |--LIB-FULL
        |--LIB-FULL.[condition1]_[condition2].dba.report_DEseq2.txt
        |--LIB-FULL.[condition1]_[condition2].dba.report_edgeR.txt
        |--LIB-FULL.[condition1]_[condition3].dba.report_DEseq2.txt
         ...
|--DBS_bed
    |--LIB-FULL
        |--DEseq2
            |--[condition1]_[condition2]_q0pt05.FC2.DEseq2.down.bed
            |--[condition1]_[condition2]_q0pt05.FC2.DEseq2.up.bed
            |--[condition1]_[condition3]_q0pt05.FC2.DEseq2.down.bed
             ...
        |--edgeR
            |--[condition1]_[condition2]_q0pt05.FC2.edgeR.down.bed
            |--[condition1]_[condition2]_q0pt05.FC2.edgeR.up.bed
            |--[condition1]_[condition3]_q0pt05.FC2.edgeR.down.bed
             ...
        |--inBoth
            |--[condition1]_[condition2]_q0pt05.FC2.inBoth.down.bed
            |--[condition1]_[condition2]_q0pt05.FC2.inBoth.up.bed
            |--[condition1]_[condition3]_q0pt05.FC2.inBoth.down.bed
             ...
```

## Note:

### Parameter:
#### Config:
Th: 0.05
doBlacklist = F
doGreylist = F
bParallel = T
bUseSummarizeOverlaps = T
bScaleControl = T
bSubControl = T
bRemoveDuplicates = F 
minOverlap = 1
summits = T
filter = 0
method = DBA_ALL_METHODS
consensus = DBA_CONDITION (stages)

#### 
```R
# score to visualization  and in output table
    score <- c("DBA_SCORE_RPKM_FOLD")

    # other options
    #score<-c("DBA_SCORE_READS", "DBA_SCORE_CONTROL_READS",
    #          "DBA_SCORE_READS_FOLD", "DBA_SCORE_READS_MINUS", "DBA_SCORE_RPKM",
    #          "DBA_SCORE_RPKM_FOLD", "DBA_SCORE_RPKM_MINUS", "DBA_SCORE_SUMMIT",
    #          "DBA_SCORE_SUMMIT_ADJ", "DBA_SCORE_SUMMIT_POS"
    # )
```
```R
# method for correlation display in correlation plot
method <- c("pearson", "spearman")

# method for heatmap
heatmap.method <- list(
        hclust = c("ward.D2", "complete", "average", "mcquitty"),
        dist = c("euclidean", "manhattan", "correlation", "maximum")
      )

# method for normalized
 nor.parameter <- data.frame(
    nor = c("DBA_NORM_LIB"),
    lib = c("DBA_LIBSIZE_FULL"),
    BG = c(F),
    offset = c(F),
    title = c("LIB-FULL")
  )

```

## Usage: 

```R
DiffBind.build.fun(samplesheet_path = samplesheet_path, outpath = outpath, peakset = NULL) 
## build DBA item and return DBdata.count (also save as RData (dba_DBA.RData) in outpath)

DiffBind.nor.DE.fun(DBdata.count, peakset = NULL, color = color, log2FC = 1, q = 0.05)
## do dba.compare and output result in outpath

```

## Example:

```R
files.sources <- list.files("./Toolbox/DiffBind", ".R$", ignore.case = T, full.names = T)
sapply(files.sources, source)

outpath <- "./DiffBind/output"
samplesheet_path <- "./samplesheet.csv"
samples <- read.delim(samplesheet_path, stringsAsFactors = F, sep = ",")

log_file <- file.path(outpath, "output.log")
sink(log_file, append = TRUE)
print(samples)
DBdata.count <- DiffBind.build.fun(samplesheet_path, outpath)
dba.show(DBdata.count)
color <- c("white", "darkgoldenrod1", "orangered2", "coral3")
DiffBind.nor.DE.fun(DBdata.count, color = color, log2FC = 1, q = 0.05)
sink()

```

## Scripts:
- diffbind.fun.v.1.0.R

## Structure:
```R
DiffBind.build.fun(samplesheet_path = samplesheet_path, outpath = outpath, peakset = NULL)
    DiffBind::dba(
                sampleSheet = samples, 
                attributes = c(DBA_ID, DBA_REPLICATE, DBA_FACTOR, DBA_CONDITION),
                config = list(
                th = th,
                AnalysisMethod = DBA_ALL_METHODS,
                doBlacklist = F,
                doGreylist = F
                )
                )

    DiffBind::dba.count(DBdata,
                peaks = cons.peaks, bUseSummarizeOverlaps = T,
                filter = 0,
                summits = T,
                bParallel = T,
                bScaleControl = T,
                bSubControl = T
                )


DiffBind.nor.DE.fun(DBdata.count, peakset = NULL, color = color, log2FC = 1, q = 0.05)  
    DiffBind::dba.normalize(
                DBdata.count,
                method = DBA_ALL_METHODS,
                normalize = get(nor.parameter$nor[nor]),
                library = get(nor.parameter$lib[nor]),
                background = nor.parameter$BG[nor],
                offset = nor.parameter$offset[nor]
                )

    value.box.plot(data, title, peakset = NULL)

    DiffBind::dba.contrast(norm, categories = DBA_CONDITION, minMembers = 2)   

    dba.plot(data, title)   

    contrast.dba(norm, title)
        DiffBind::dba.analyze(norm, method = DBA_ALL_METHODS)
        DiffBind::dba.show(dba.analyz, bContrasts = T, th = th)

    pairwise.matrix(res.matrix, tool)

    matrix.plot(data, title, color)

```


