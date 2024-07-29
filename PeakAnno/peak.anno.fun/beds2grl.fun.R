# author: YCW
# date:2024.05.23


beds2grl <- function(files) {
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

    return(grl)
}
# fin-----
