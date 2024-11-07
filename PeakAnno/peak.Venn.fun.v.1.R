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

peak.venn <- function(files, outpath, date, colors = NULL, connectedPeaks = "merge", outname= "outname") {
    print(outpath)
    print(outname)
    ## input -----
    if (!is.data.frame(files)) {
        if (file.exists(files) %>% which() %>% length() == length(files)) {
            files <- data.frame(
                path = files
            ) %>%
                mutate(ID = path %>% as.character() %>% basename() %>% str_remove(".\\w+$"))
            nrow(files) %>% print()
            head(files) %>% print()
            print("0-base coordinate system in use \n please check your input files format") 
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
    if (colors %>% is.null()) {
        colors <- rainbow(length(grl1))
    }

    print(file.path(outpath, paste0("peak.venn.", date, ".pdf")))
    pdf(file.path(outpath, paste0(outname, ".peak.venn.", date, ".pdf")), width = 18, height = 18)

    if(length(grl1)==2){
        Cat.dist<- c( 0.01, 0.01)
        Cat.pos <- c(5, -5)
        Cex=3
        Cat.cex=3
    }else if(length(grl1)==3){
        Cat.dist<- c( 0.02, 0.02, 0.02)
        Cat.pos <- c(-10, 10, 180)
        Cex=3
        Cat.cex=2

    }else if(length(grl1)==4){
        Cat.dist<- c( 0.2, 0.2, 0.1, 0.1)
        Cat.pos <- c(0, 0, 0, 0)
        Cex=3
        Cat.cex=2
        grl1<- grl1[c(1,4,2,3)]

    }else if(length(grl1)==5){
        Cat.dist<- c( 0.2, 0.2, 0.21, 0.18, 0.2)
        Cat.pos <- c(5, -30, -120, 150, 30)
        Cex=2
        Cat.cex=2
        grl1<-grl1[c(1, 5 ,4 ,3 , 2)]
    }

    tryCatch(
            {
            makeVennDiagram(grl1,
                NameOfPeaks = c(paste(names(lengths(grl1)), lengths(grl1), sep="\n")),
                cex = Cex,  cat.cex=Cat.cex, 
                fontface = 2,
                scaled = T, euler.d = T, totalTest = lengths(grl1) %>% max(), 
                cat.fontface = 3,
                cat.dist = Cat.dist, # Modified
                cat.pos =  Cat.pos, # Modified
                connectedPeaks = connectedPeaks,
                fill = colors, # circle fill color
                col = "black", # circle border color
                cat.col = "black"
            )
            },
            error = function(e) {
                cat("An error occurred:", conditionMessage(e), "\n")
            }
        )
        
   
    dev.off()
}
