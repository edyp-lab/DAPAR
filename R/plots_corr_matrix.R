#' @title Displays a correlation matrix of the quantitative data of the
#' \code{Biobase::exprs()} table
#' 
#' @description Builds a correlation matrix based on a \code{MSnSet} object.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param rate A float that defines the gradient of colors.
#'
#' @param showValues xxx
#'
#' @return A colored correlation matrix
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept)
#' wrapper.corrMatrixD_HC(Exp1_R25_pept)
#'
#'
#' @export
#'
wrapper.corrMatrixD_HC <- function(obj, rate = 0.5, showValues = TRUE) {
    
    if (!requireNamespace("stats", quietly = TRUE)) {
        stop("Please install stats: BiocManager::install('stats')")
    }
    
    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    data <- stats::cor(qData, use = "pairwise.complete.obs")
    corrMatrixD_HC(data, samplesData, rate, showValues)
}




#' @title Displays a correlation matrix of the quantitative data of the
#' \code{Biobase::exprs()} table.
#'
#' @param object The result of the \code{cor} function.
#'
#' @param samplesData A dataframe in which lines correspond to samples and
#' columns to the meta-data for those samples.
#'
#' @param rate The rate parameter to control the exponential law for
#' the gradient of colors
#'
#' @param showValues xxx
#'
#' @return A colored correlation matrix
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept)
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' res <- cor(qData, use = "pairwise.complete.obs")
#' corrMatrixD_HC(res, samplesData)
#'
#' @import highcharter
#'
#' @export
#'
corrMatrixD_HC <- function(object, 
    samplesData = NULL, 
    rate = 0.5, 
    showValues = TRUE) {
    if (!requireNamespace("stats", quietly = TRUE)) {
        stop("Please install stats: BiocManager::install('stats')")
    }
    
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Please install dplyr: BiocManager::install('dplyr')")
    }
    
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Please install tidyr: BiocManager::install('tidyr')")
    }
    
    if (!requireNamespace("tibble", quietly = TRUE)) {
        stop("Please install tibble: BiocManager::install('tibble')")
    }
    df <- as.data.frame(object)
    .sData <- samplesData
    if (!is.null(.sData)) {
        for (j in seq_len(ncol(df))) {
            names(df)[j] <- paste(as.character(.sData[j, 2:ncol(.sData)]),
                collapse = " "
            )
        }
    }
    is.num <- vapply(df, is.numeric, FUN.VALUE = NA)
    df[is.num] <- lapply(df[is.num], round, 2)
    dist <- NULL

    x <- y <- names(df)

    df <- dplyr::tbl_df(cbind(x = y, df)) %>%
        tidyr::gather(y, dist, -x) %>%
        dplyr::mutate(
            x = as.character(x),
            y = as.character(y)
        ) %>%
        dplyr::left_join(
            tibble::tibble(
                x = y,
                xid = seq(length(y)) - 1
                ), 
            by = "x") %>%
        dplyr::left_join(
            tibble::tibble(
                y = y,
                yid = seq(length(y)) - 1
                ), 
            by = "y")

    ds <- df %>%
        dplyr::select("xid", "yid", "dist") %>%
        list_parse2()

    fntltp <- JS("function(){
    return this.series.xAxis.categories[this.point.x] + ' ~ ' +
    this.series.yAxis.categories[this.point.y] + ': <b>' +
    Highcharts.numberFormat(this.point.value, 2)+'</b>';
        ; }")
    cor_colr <- list(
        list(0, "#FF5733"),
        list(0.5, "#F8F5F5"),
        list(1, "#2E86C1")
    )
    highchart() %>%
        my_hc_chart(chartType = "heatmap") %>%
        hc_xAxis(categories = y, title = NULL) %>%
        hc_yAxis(categories = y, title = NULL) %>%
        hc_add_series(data = ds) %>%
        hc_plotOptions(
            series = list(
                boderWidth = 0,
                dataConditions = list(enabled = TRUE),
                dataLabels = list(enabled = showValues)
            )
        ) %>%
        hc_tooltip(formatter = fntltp) %>%
        hc_legend(
            align = "right", layout = "vertical",
            verticalAlign = "middle"
        ) %>%
        hc_colorAxis(stops = cor_colr, min = rate, max = 1) %>%
        my_hc_ExportMenu(filename = "corrMatrix")
}
