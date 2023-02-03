

#' @title This function is a wrapper to \code{heatmap.2} that displays
#' quantitative data in the \code{Biobase::exprs()} table of an object of
#' class \code{MSnSet}
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param distance The distance used by the clustering algorithm to compute
#' the dendrogram. See \code{help(heatmap.2)}.
#'
#' @param cluster the clustering algorithm used to build the dendrogram.
#' See \code{help(heatmap.2)}
#'
#' @param dendro A boolean to indicate fi the dendrogram has to be displayed
#'
#' @return A heatmap
#'
#' @author Alexia Dorffer
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- 'peptide'
#' metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
#' indices <- GetIndices_WholeLine(metacell.mask)
#' wrapper.heatmapD(obj)
#'
#' @export
#'
wrapper.heatmapD <- function(obj,
    distance = "euclidean",
    cluster = "complete",
    dendro = FALSE
    ) {
    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)[["Condition"]]
    for (j in seq_len(length(colnames(qData)))) {
        .ncol <- ncol(Biobase::pData(obj))
        colnames(qData)[j] <- paste(
            as.character(
                Biobase::pData(obj)[j, seq.int(from = 2, to = .ncol)]),
            collapse = " "
        )
    }
    heatmapD(qData, conds, distance, cluster, dendro)
}




#' @title This function is a wrapper to \code{heatmap.2} that displays
#' quantitative data in the \code{Biobase::exprs()} table of an object of
#' class \code{MSnSet}
#'
#' @param qData A dataframe that contains quantitative data.
#'
#' @param conds A vector containing the conditions
#'
#' @param distance The distance used by the clustering algorithm to compute
#' the dendrogram. See \code{help(heatmap.2)}
#'
#' @param cluster the clustering algorithm used to build the dendrogram.
#' See \code{help(heatmap.2)}
#'
#' @param dendro A boolean to indicate fi the dendrogram has to be displayed
#'
#' @return A heatmap
#'
#' @author Florence Combes, Samuel Wieczorek, Enor Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10), ]
#' level <- 'peptide'
#' metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
#' indices <- GetIndices_WholeLine(metacell.mask)
#' qData <- Biobase::exprs(obj)
#' conds <- Biobase::pData(obj)[["Condition"]]
#' heatmapD(qData, conds)
#'
#'
#' @export
#'
heatmapD <- function(qData,
    conds,
    distance = "euclidean",
    cluster = "complete",
    dendro = FALSE) {
    
    pkgs.require(c('stats', 'dendextend', "gplots", 'grDevices', 'RColorBrewer'))
    
    
    .data <- matrix(qData,
        ncol = ncol(qData),
        byrow = FALSE,
        dimnames = list(rownames(qData), colnames(qData))
    )
    colors <- c(
        seq(-3, -2, length = 100),
        seq(-2, 0.5, length = 100),
        seq(0.5, 6, length = 100)
    )
    heatmap.color <- grDevices::colorRampPalette(c("green", "red"))(n = 1000)

    # samples label color
    x <- t(.data)
    x[is.na(x)] <- -1e5
    dist <- dist(x, method = distance)
    hcluster <- stats::hclust(dist, method = cluster)
    pal <- GetColorsForConditions(conds)
    
    cols_branches <- pal
    dend1 <- stats::as.dendrogram(hcluster)
    dend1 <- dendextend::color_branches(
        dend1, 
        k = length(conds), 
        col = cols_branches
        )
    col_labels <- dendextend::get_leaves_branches_col(dend1)

    if (dendro) {
        .dendro <- "row"
    } else {
        .dendro <- "none"
    }
    p <- gplots::heatmap.2(
        x = t(.data),
        distfun = function(x) {
            x[is.na(x)] <- -1e5
            dist(x, method = distance)
        },
        hclustfun = function(x) {
            x[is.na(x)] <- -1e5
            stats::hclust(x, method = cluster)
        },
        dendrogram = .dendro,
        Rowv = TRUE,
        col = heatmap.color,
        density.info = "none",
        key = TRUE,
        trace = "none",
        scale = "none",
        # srtCol=45,
        labCol = "",
        margins = c(4, 12),
        cexRow = 1.5 + ncol(.data) * -0.011,
        keysize = 1.5,
        # lhei = c(1.5, 9),
        # lwid = c(1.5, 4),
        # lmat = rbind(4:3, 2:1),
        colRow = col_labels
    )
}


#' @title  xxx
#' 
#' @description 
#' This function is inspired from the function \code{heatmap.2}
#' that displays quantitative data in the \code{Biobase::exprs()} table of an 
#' object of
#' class \code{MSnSet}. For more information, please refer to the help
#' of the heatmap.2 function.
#'
#' @param x A dataframe that contains quantitative data.
#'
#' @param col colors used for the image. Defaults to heat colors (heat.colors).
#'
#' @param srtCol angle of column conds, in degrees from horizontal
#'
#' @param labCol character vectors with column conds to use.
#'
#' @param labRow character vectors with row conds to use.
#'
#' @param key logical indicating whether a color-key should be shown.
#'
#' @param key.title main title of the color key. If set to NA no title will
#' be plotted.
#'
#' @param main main title; default to none.
#'
#' @param ylab y-axis title; default to none.
#'
#' @return A heatmap
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(100)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
#' indices <- GetIndices_WholeLine(metacell.mask)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' heatmapForMissingValues(qData)
#'
#' @export
#'
heatmapForMissingValues <- function(x,
    col = NULL,
    srtCol = NULL,
    labCol = NULL,
    labRow = NULL,
    key = TRUE,
    key.title = NULL,
    main = NULL,
    ylab = NULL) {
        
    pkgs.require(c('grDevices', 'graphics'))

        if (is.null(col))
            col <- grDevices::heat.colors(100)
        
        scale01 <- function(x, low = min(x), high = max(x)) {
            x <- (x - low) / (high - low)
            x
        }

        offsetCol <- 0.5
        offsetRow <- 0.5
        srtRow <- NULL
        colRow <- NULL
        colCol <- NULL
        xlab <- NULL
        key.par <- list()
        margins <- c(5, 5)
        sepcolor <- "white"
        na.color <- "white"
        keysize <- 1.5
        breaks <- NULL
        na.rm <- TRUE

        if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
            stop("`x' must be a numeric matrix")
        }
        nr <- di[1]
        nc <- di[2]
        if (nr <= 1 || nc <= 1) {
            stop("`x' must have at least 2 rows and 2 columns")
        }
        x <- x[nr:1, ]
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
        cexCol <- 0.2 + 1 / log10(nc)
        cexRow <- 0.2 + 1 / log10(nr)
        iy <- seq_len(nr)
        breaks <- length(col) + 1
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
            length = breaks
        )

        nbr <- length(breaks)
        ncol <- length(breaks) - 1

        min.breaks <- min(breaks)
        max.breaks <- max(breaks)
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks
        lhei <- c(keysize, 4)
        lwid <- c(keysize, 4)
        lmat <- rbind(4:3, 2:1)
        lmat[is.na(lmat)] <- 0

        op <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(op))
        graphics::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

        graphics::par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)


        graphics::image(seq_len(nc), seq_len(nr), x,
            xlim = 0.5 + c(0, nc), ylim = 0.5 +
                c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
            breaks = breaks
        )


        if (!is.null(labCol)) {
            graphics::axis(1, 
                seq_len(nc),
                label = labCol, 
                las = 2, 
                line = -0.5 + offsetCol, 
                tick = 0, 
                cex.axis = cexCol, 
                hadj = NA,
                padj = 0
            )
        } else {
            adjCol <- c(1, NA)
            xpd.orig <- graphics::par("xpd")
            graphics::par(xpd = NA)
            xpos <- graphics::axis(1, 
                seq_len(nc),
                label = rep("", nc), 
                las = 2,
                tick = 0
            )
            graphics::text(
                x = xpos, 
                y = graphics::par("usr")[3] - (1 + offsetCol) *
                    graphics::strheight("M"), 
                label = labCol, 
                adj = adjCol,
                cex = cexCol, 
                srt = srtCol, 
                col = colCol
            )
            graphics::par(xpd = xpd.orig)
        }


        if (!is.null(labRow)) {
            graphics::axis(4, 
                iy,
                label = labRow, 
                las = 5, 
                line = -0.5 + offsetRow,
                tick = 0, 
                cex.axis = cexRow, 
                hadj = 0, 
                padj = NA
            )
        } else {
            xpd.orig <- graphics::par("xpd")
            graphics::par(xpd = NA)
            ypos <- graphics::axis(4, 
                iy,
                label = rep("", nr), 
                las = 2,
                line = -0.5, 
                tick = 0
            )
            
            .strw <- graphics::strwidth("M")
            graphics::text(
                x = graphics::par("usr")[2] + (1 + offsetRow) * .strw,
                y = ypos, 
                label = labRow, 
                adj = c(0, NA), 
                cex = cexRow,
                srt = srtRow, 
                col = colRow
            )
            graphics::par(xpd = xpd.orig)
        }

        graphics::par(mar = c(margins[1], 0, 0, 0))
        graphics::plot.new()
        graphics::par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))

        graphics::plot.new()
        if (!is.null(main)) {
            graphics::title(main, cex.main = 1.5 * op[["cex.main"]])
        }


        if (key) {
            mar <- c(5, 4, 2, 1)
            graphics::par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
            if (length(key.par) > 0) {
                do.call(par, key.par)
            }

            tmpbreaks <- breaks
            min.raw <- min.breaks
            max.raw <- max.breaks

            z <- seq(min.raw, max.raw, by = min(diff(breaks) / 100))
            graphics::image(
                z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                xaxt = "n", yaxt = "n"
            )
            graphics::par(usr = c(0, 1, 0, 1))
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, label = lv)

            xargs$side <- 1
            do.call(graphics::axis, xargs)
            key.xlab <- "Intensity value"

            graphics::mtext(
                side = 1, 
                key.xlab, 
                line = graphics::par("mgp")[1], 
                padj = 0.5,
                cex = graphics::par("cex") * graphics::par("cex.lab")
            )

            if (is.null(key.title)) {
                graphics::title("Color Key")
            }
        }
    }

