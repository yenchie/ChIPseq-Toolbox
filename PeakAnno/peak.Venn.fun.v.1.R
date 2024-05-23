# peak.Venn.fun.v.1.R
# YCW
# 2024.05.15
# 2024.04.22

# USAGE
# source("/bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Toolbox/PeakAnno/peak.Venn.fun.v.1.R")
# beds.path <- "./beds"
# beds <- data.frame(path = beds.path %>% list.files(full.names = T, "bed"), ID = LETTERS[seq( from = 1, to = length(beds.path %>% list.files(full.names = T, "bed")))] )
# beds
# outpath <- "./"

# date = "20240515"
# peak.venn(files = beds, outpath, date)


library(ChIPpeakAnno)
library(tidyverse)

peak.venn <- function(files, outpath, date, colors=NULL) {
    print(outpath)
    ## input -----
    if (!is.data.frame(files)) {
        if (file.exists(files) %>% which() %>% length() == length(files)) {
            files <- data.frame(
                path = files
            ) %>%
                mutate(ID = path %>% as.character() %>% basename() %>% str_remove(".\\w+$"))
            nrow(files) %>% print()
            head(files) %>% print()
            print("0-base coordinate system in use \n please check your input files format") # names for peakset
        } else {
            print(" please recheck the input table. \n files formate: col1: path; col2: ID # names of peaksets")
            print("0-base coordinate system in use \n please check your input files format")
            return()
        }
    } else if (is.data.frame(files)) {
        colnames(files) <- c("path", "ID")
        if (file.exists(files$path) %>% which() %>% length() == nrow(files)) {
            nrow(files) %>% print()
            head(files) %>% print()
        } else {
            print(" please recheck the input table. \n files formate: col1: path; col2: ID # names of peaksets")
            print("0-base coordinate system in use \n please check your input files format")
        }
    } else {
        print(" please recheck the input table. \n files formate: col1: path; col2: ID # names of peaksets")
        print("0-base coordinate system in use \n please check your input files format")
        return()
    }

    grl1 <- NULL
    grl <- GRangesList()
    for (id in unique(files$ID)) {
        if (files$path[files$ID == id] %>% file.size() != 0) {
            gr <- read.delim(files$path[files$ID == id], header = FALSE)
            colnames(gr)[1:3] <- c("seqnames", "start", "end")
            gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE, starts.in.df.are.0based = T)
            grl[[id]] <- gr
        } else {
            print(paste(id, "show none region in the providing file"))
            grl[[id]] <- GRanges()
        }
    }

    head(grl)
    length(grl)
    names(grl)

    grl1 <- grl
    if(colors=NULL){
        colors = rainbow(length(grl1))
    }

    print(file.path(outpath, paste0("peak.venn.", date, ".pdf")))
    pdf(file.path(outpath, paste0("peak.venn.", date, ".pdf")))
    tryCatch(
        {
            makeVennDiagram(grl1,
                NameOfPeaks = c(names(grl1)),
                scaled = FALSE, euler.d = FALSE, totalTest = lengths(grl1) %>% max(),
                connectedPeaks = "keepAll",
                fill = colors, # circle fill color
                col = colors, # circle border color
                cat.col = colors
            )
        },
        error = function(e) {
            cat("An error occurred:", conditionMessage(e), "\n")
        }
    )

    dev.off()
}
