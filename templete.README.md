# FUN.name (v.1)

## date: 
## author:

## Description:
[fun](https://bioconductor.org/packages/fun/)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("fun")
```

## Dependency:

## Input:
- 

## Output: 

```
|output
|--
|--|--
```

## Note:
### Reference:
#### 
- Gmax189


## Usage: 

```R
fun(x)
```

## Example:

```R
files.sources <- list.files("./", ".R$", ignore.case = T, full.names = T)
sapply(files.sources, source)
datapath <- "./output"

log_file <- file.path(datapath, "output.log")
sink(log_file, append = TRUE)
fun(x)
sink()
```

## Scripts:
- 

## Structure:
```R
fun(x)
    fun.1(x)
    fun.2(x)
        fun.2.1(x)
```


