

#' Densityplot of quantitative proteomics data over samples.
#'
#' @title Builds a densityplot from a dataframe
#'
#' @param obj xxx
#'
#' @param legend A vector of the conditions (one condition
#' per sample).
#'
#' @param pal xxx
#'
#' @return A density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' densityPlotD_HC(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' pal <- ExtendPalette(2, "Dark2")
#' densityPlotD_HC(Exp1_R25_pept, pal = pal)
#'
#' @import highcharter
#'
#' @export
#'
densityPlotD_HC <- function(obj,
                            legend = NULL,
                            pal = NULL) {
    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)$Condition

    if (is.null(legend)) {
        legend <- Biobase::pData(obj)[, "Condition"]
    }

    myColors <- NULL
    if (is.null(pal)) {
        myColors <- GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
    } else {
        if (length(pal) != length(unique(conds))) {
            warning("The color palette has not the same dimension as the number of samples. Set to default.")
            myColors <- GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
        } else {
            myColors <- GetColorsForConditions(conds, pal)
        }
    }

    h1 <- highchart() %>%
        hc_title(text = "Density plot") %>%
        my_hc_chart(chartType = "spline", zoomType = "x") %>%
        hc_colors(myColors) %>%
        hc_legend(enabled = TRUE) %>%
        hc_xAxis(title = list(text = "log(Intensity)")) %>%
        hc_yAxis(title = list(text = "Density")) %>%
        hc_tooltip(
            headerFormat = "",
            pointFormat = "<b> {series.name} </b>: {point.y} ",
            valueDecimals = 2
        ) %>%
        my_hc_ExportMenu(filename = "densityplot") %>%
        hc_plotOptions(
            series = list(
                animation = list(
                    duration = 100
                ),
                connectNulls = TRUE,
                marker = list(
                    enabled = FALSE
                )
            )
        )

    if (is.null(legend)) {
        legend <- paste0("series", 1:ncol(qData))
    }

    for (i in 1:ncol(qData)) {
        tmp <- data.frame(
            x = density(qData[, i], na.rm = TRUE)$x,
            y = density(qData[, i], na.rm = TRUE)$y
        )

        h1 <- h1 %>% hc_add_series(data = list_parse(tmp), name = legend[i])
    }

    return(h1)
}
