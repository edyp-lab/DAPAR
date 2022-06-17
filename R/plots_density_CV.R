
#' Builds a densityplot of the CV of entities in the Biobase::exprs() table.
#' of an object \code{MSnSet}. The variance is calculated for each
#' condition present
#' in the dataset (see the slot \code{'Condition'} in the \code{Biobase::pData()} table).
#'
#' @title Distribution of CV of entities
#'
#' @param obj An object of class \code{MSnSet}
#'
#' @param ... arguments for palette.
#'
#' @return A density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' wrapper.CVDistD_HC(Exp1_R25_pept)
#'
#' @export
#'
wrapper.CVDistD_HC <- function(obj, ...) {
    if (nrow(obj) == 0) {
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)[, "Condition"]
    CVDistD_HC(qData, conds, ...)
}



#' Builds a densityplot of the CV of entities in the Biobase::exprs() table
#' of a object. The CV is calculated for each condition present
#' in the dataset (see the slot \code{'Condition'} in the \code{Biobase::pData()} table)
#'
#' @title Distribution of CV of entities
#'
#' @param qData A dataframe that contains quantitative data.
#'
#' @param conds A vector of the conditions (one condition per sample).
#'
#' @param pal xxx
#'
#' @return A density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' conds <- Biobase::pData(Exp1_R25_pept)[, "Condition"]
#' CVDistD_HC(Biobase::exprs(Exp1_R25_pept), conds)
#' pal <- ExtendPalette(2, "Dark2")
#' CVDistD_HC(Biobase::exprs(Exp1_R25_pept), conds, pal)
#'
#' @import highcharter
#' @importFrom stats density var
#'
#' @export
#'
CVDistD_HC <- function(qData,
                       conds = NULL,
                       pal = NULL) {
    if (is.null(conds)) {
        warning("The vector of conditions is empty. The plot cannot be drawn.")
        return(NULL)
    }

    conditions <- unique(conds)
    n <- length(conditions)


    if (is.null(pal)) {
        pal <- ExtendPalette(n)
    } else {
        if (length(pal) != n) {
            warning("The color palette has not the same dimension as the number of samples. Set to default.")
            pal <- ExtendPalette(n)
        }
    }


    h1 <- highchart() %>%
        my_hc_chart(chartType = "spline", zoomType = "x") %>%
        hc_colors(pal) %>%
        hc_legend(enabled = TRUE) %>%
        hc_xAxis(title = list(text = "CV(log(Intensity))")) %>%
        hc_yAxis(title = list(text = "Density")) %>%
        hc_tooltip(
            headerFormat = "",
            pointFormat = "<b>{series.name}</b>: {point.y} ",
            valueDecimals = 2
        ) %>%
        my_hc_ExportMenu(filename = "logIntensity") %>%
        hc_plotOptions(
            series = list(
                connectNulls = TRUE,
                marker = list(
                    enabled = FALSE
                )
            )
        )

    minX <- maxX <- 0
    maxY <- 0
    for (i in 1:n) {
        if (length(which(conds == conditions[i])) > 1) {
            t <- apply(
                qData[, which(conds == conditions[i])], 1,
                function(x) 100 * var(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
            )
            tmp <- data.frame(
                x = density(t, na.rm = TRUE)$x,
                y = density(t, na.rm = TRUE)$y
            )

            ymaxY <- max(maxY, tmp$y)
            xmaxY <- tmp$x[which(tmp$y == max(tmp$y))]
            minX <- min(minX, tmp$x)
            maxX <- max(maxX, 10 * (xmaxY - minX))


            h1 <- h1 %>% hc_add_series(data = tmp, name = conditions[i])
        }
    }

    h1 <- h1 %>%
        hc_chart(
            events = list(
                load = JS(paste0("function(){
                         var chart = this;
                         this.xAxis[0].setExtremes(", minX, ",", maxX, ");
                         this.showResetZoom();}"))
            )
        )

    return(h1)
}
