
#' @title Builds a violinplot from a dataframe
#'
#' @param obj xxx
#'
#' @param conds xxx
#'
#' @param keyId xxx
#'
#' @param legend A vector of the conditions (one condition per sample).
#'
#' @param pal xxx
#'
#' @param subset.view xxx
#'
#' @return A violinplot
#'
#' @author Samuel Wieczorek, Anais Courtier
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot
#' legend <- conds <- Biobase::pData(obj)$Condition
#' key <- "Protein_IDs"
#' violinPlotD(obj, conds, key, legend, subset.view = seq_len(10))
#'
#'
#' @export
#'
violinPlotD <- function(obj,
                        conds,
                        keyId,
                        legend = NULL,
                        pal = NULL,
                        subset.view = NULL) {
    pkgs.require(c('stats', 'vioplot', 'graphics'))
    
    graphics::plot.new()

    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    } else {
        qData <- Biobase::exprs(obj)
    }

    if (missing(conds)) {
        stop("'conds' is missing.")
    }

    if (length(subset.view) == 0) {
        subset.view <- NULL
    }

    if (is.null(legend)) {
        legend <- conds
        for (i in unique(conds)) {
            legend[which(conds == i)] <- paste0(i, "_", 
                seq_len(length(which(conds == i))))
        }
    }

    if (!is.null(subset.view)) {
        if (is.null(keyId) || missing(keyId)) {
            stop("'keyId' is missing.")
        } else {
            if (!grep(keyId, colnames(Biobase::fData(obj)))) {
                stop("'keyId' does not belong to metadata")
            }
        }
    }

    myColors <- NULL
    if (is.null(pal)) {
        warning("Color palette set to default.")
        myColors <- GetColorsForConditions(conds,
            ExtendPalette(length(unique(conds))))
    } else {
        if (length(pal) != length(unique(conds))) {
            warning("The color palette has not the same dimension as
                the number of samples")
            myColors <- GetColorsForConditions(conds,
                ExtendPalette(length(unique(conds))))
        } else {
            myColors <- GetColorsForConditions(conds, pal)
        }
    }

    graphics::plot.window(
        xlim = c(0, ncol(qData) + 1),
        ylim = c(
            min(stats::na.omit(qData)), 
            max(stats::na.omit(qData))
            )
    )
    graphics::title(
        ylab = "Log (intensity)",
        xlab = "Samples"
    )

    for (i in seq_len(ncol(qData))) {
        vioplot::vioplot(
            stats::na.omit(
                qData[, i]), 
                col = myColors[i], 
                add = TRUE, 
                at = i)
    }


    graphics::axis(2, yaxp = c(
        floor(min(stats::na.omit(qData))),
        floor(max(stats::na.omit(qData))), 5
    ), las = 1)

    if (!is.null(legend)) {
        if (is.vector(legend)) {
            N <- 1
        } else {
            N <- ncol(legend)
        }

        for (i in seq_len(N)) {
            graphics::axis(
                side = 1,
                at = seq_len(ncol(qData)),
                labels = if (is.vector(legend)) {
                    legend
                } else {
                    legend[, i]
                },
                line = 2 * i - 1
            )
        }
    }

    # Display of rows to highlight (index of row in subset.view)
    if (!is.null(subset.view)) {
        idVector <- keyId
        pal <- ExtendPalette(length(subset.view), "Dark2")

        n <- 0
        for (i in subset.view) {
            n <- n + 1
            for (c in seq_len(ncol(qData) - 1)) {
                graphics::segments(
                    y0 = qData[i, c], 
                    y1 = qData[i, c + 1], 
                    x0 = c, 
                    x1 = c + 1, 
                    pch = 16, 
                    col = pal[n], 
                    lwd = 2)
                graphics::points(
                    y = qData[i, c], 
                    x = c, 
                    pch = 16, 
                    col = pal[n])
            }
            graphics::points(
                y = qData[i, ncol(qData)], 
                x = ncol(qData), 
                pch = 16, 
                col = pal[n])
        }
        graphics::legend("topleft",
            legend = Biobase::fData(obj)[subset.view, keyId],
            lty = 1,
            lwd = 2,
            col = pal,
            pch = 16,
            bg = "transparent",
            bty = "n"
        )
    }
}
