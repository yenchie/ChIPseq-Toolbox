## peakannt.annoGene.fun.v.1.R
# author: YCW
# date:2024.06.05
# revised:
# reviser:

# dependency:
# beds2grl.fun.R
# annoGene.R
# annoPeak.R

## function ----
anno.gene.fun <- function(x) {
    return(x@annoGene)
}


## core fun -----

print("files formate: col1: path; col2: ID # names of peaksets")
print("0-base coordinate system in use \n please check your input files format")

peak.anno.fun <- function(files = NULL, grl = NULL, outpath = NULL, date = date, plot = F, output.gene.table = T) {
    if (outpath %>% is.null()) {
        outpath <- "./"
    }
    print(outpath)

    if (is.null(grl)) {
        grl <- beds2grl(files)
    } else {
        print("grl as input")
    }

    head(grl)
    length(grl)
    names(grl)

    grl1 <- grl[which(lengths(grl) != 0)]
    head(grl1)
    length(grl1)
    names(grl1)

    annoPeakList <- lapply(grl, annoGene)
    genes <- lapply(annoPeakList, anno.gene.fun)

    head(genes)
    length(genes)
    names(genes)

    if (output.gene.table == T) {
        gene.table <- gene.table.0 <- NULL
        for (l in names(grl)) {
            d1 <- NULL
            x <- NULL
            if (l %in% names(genes)) {
                x <- genes[[l]]
                head(x)
                d1 <- data.frame(geneID = x %>% unlist()) %>%
                    mutate(group = l)
                print(nrow(d1))
                gene.table <- gene.table %>% bind_rows(d1)

                write.table(d1$geneID, file.path(outpath, paste0(l, "annoGene.geneAnno.", date, ".gene.id")), quote = F, row.names = F, sep = "\t", col.names = F)
            } else {
                d1 <- NULL
                print(paste(id, "none annotated genes, output empty file"))
                write.table(d1, file.path(outpath, paste0(l, ".annoGene.geneAnno.", date, ".gene.id")), quote = F, row.names = F, sep = "\t", col.names = F)
                gene.table.0 <- gene.table.0 %>% bind_rows(data.frame(group = l, n = 0))
            }
        }
        print(file.path(outpath, paste0(l, ".annoGene.geneAnno.", date, ".gene.id")))

        head(gene.table)
        gene.count <- gene.table %>%
            dplyr::count(group) %>%
            bind_rows(gene.table.0) %>%
            dplyr::rename(count = n)

        print(file.path(outpath, paste0("annoGene.gene.number.", date, ".txt")))
        write.table(gene.count, file.path(outpath, paste0("annoGene.gene.number.", date, ".txt")), quote = F, row.names = F, sep = "\t")
    }

    if (plot == T) {
        # plotting -----
        print(file.path(outpath, paste0("annoGene.plot.", date, ".pdf")))
        pdf(file.path(outpath, paste0("annoGene.plot.", date, ".pdf")))
        if(length(grl1)==2){
                cat.dist<- c( 0.01, 0.01)
                cat.pos <- c(5, -5)

            }else if(length(grl1)==3){
                cat.dist<- c( 0.02, 0.02, 0.02)
                cat.pos <- c(-10, 10, 180)

            }else if(length(grl1)==4){
                cat.dist<- c( 0.2, 0.2, 0.1, 0.1)
                cat.pos <- c(5, 0, 0, 0)

            }else if(length(grl1)==5){
                cat.dist<- c( 0.2, 0.2, 0.21, 0.18, 0.2)
                cat.pos <- c(5, -30, -120, 150, 30)

            }else{

            }

        colors <- rainbow(length(grl1))  
        connectedPeaks="keepAll"    
        tryCatch(
            {
                makeVennDiagram(grl1,
                    NameOfPeaks = c(paste(names(lengths(grl1)), lengths(grl1), sep="\n")),
                cex = 3,  cat.cex=2, 
                fontface = 2,
                scaled = FALSE, euler.d = FALSE, totalTest = lengths(grl1) %>% max(), 
                cat.fontface = 4,
                cat.dist = cat.dist, # Modified
                cat.pos =  cat.pos, # Modified
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


        print(plotanntBar(annoPeakList))

        print(plotdistTotss(annoPeakList))

        # gene annotation intersection ------------
        tryCatch(
            {
                print(vennplot(genes, by = "Vennerable"))
            },
            error = function(e) {
                cat("An error occurred:", conditionMessage(e), "\n")
            }
        )
        tryCatch(
            {
                print(vennplot(genes, by = "gplots"))
            },
            error = function(e) {
                cat("An error occurred:", conditionMessage(e), "\n")
            }
        )

        Vgene <- Venn(genes)
        print(plot(Vgene, doWeights = F))
        V <- Vgene
        tryCatch(
            {
                if (length(genes) <= 3) {
                    VennList <- compute.Venn(V, doWeights = F, type = "circles")
                } else if (length(genes) == 4) {
                    VennList <- compute.Venn(V, doWeights = F, type = "ellipses")
                } else {
                    VennList <- compute.Venn(V, doWeights = F, type = "ChowRuskey")
                }
            },
            error = function(e) {
                VennList <- compute.Venn(V, doWeights = F, type = "AWFE")
                # cat("An error occurred:", conditionMessage(e), "\n")
            }
        )

        pcentFun <- function(x) {
            100 * (x / sum(x))
        }
        tryCatch(
            {
                Weight <- V@IndicatorWeight %>%
                    as.data.frame() %>%
                    rownames_to_column()
                Weight_sorted <- Weight[order(match(Weight$rowname, VennList@FaceLabels$FaceName)), ]
                areas <- Weight_sorted$.Weight
                names(areas) <- Weight_sorted$rowname
                areasPcent <- round(pcentFun(areas), digits = 2)
                VennList@FaceLabels$Signature <- paste0(areas, "\n", areasPcent %>% round(1), "%")

                # plot(VennList, show = list(FaceText = c("signature"), DarkMatter = F))
                gp <- VennThemes(VennList, increasingLineWidth = F, colourAlgorithm = "signature")
                modify_vector <- function(element, option, value) {
                    if (option %in% names(element)) {
                        element[[option]] <- value
                    }
                    return(element)
                }

                gp[["SetText"]] <- lapply(gp[["SetText"]], function(x) modify_vector(x, "fontsize", 12))
                gp[["FaceText"]] <- lapply(gp[["FaceText"]], function(x) modify_vector(x, "fontsize", 12))

                gridExtra::grid.arrange(
                    grid::grid.grabExpr(height = 4, width = 6, plot(
                        VennList,
                        show = list(
                            FaceText = c("signature"),
                            DarkMatter = F
                        ),
                        gp = gp
                    )),
                    top = textGrob(
                        paste("Intersection of Gene Sets in Venn"),
                        gp = gpar(fontsize = 15)
                    )
                )
            },
            error = function(e) {
                cat("An error occurred:", conditionMessage(e), "\n")
            }
        )
        dev.off()
    }

    return(annoPeakList)
}
## fin---
