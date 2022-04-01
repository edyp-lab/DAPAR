#' 
#' 
#' #' Builds a correlation matrix based on a \code{MSnSet} object.
#' #' 
#' #' @title Displays a correlation matrix of the quantitative data of the
#' #' \code{Biobase::exprs()} table
#' #' 
#' #' @param obj An object of class \code{MSnSet}.
#' #' 
#' #' @param rate A float that defines the gradient of colors.
#' #' 
#' #' @return A colored correlation matrix
#' #' 
#' #' @author Alexia Dorffer
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' wrapper.corrMatrixD(Exp1_R25_pept)
#' #' 
#' #' 
#' #' @export
#' #' 
#' wrapper.corrMatrixD <- function(obj, rate=5){
#'   qData <- Biobase::exprs(obj)
#'   samplesData <- Biobase::pData(obj)
#'   corrMatrixD(qData, samplesData, rate)
#' }

#' Builds a correlation matrix based on a \code{MSnSet} object. 
#' 
#' @title Displays a correlation matrix of the quantitative data of the
#' \code{Biobase::exprs()} table
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.corrMatrixD_HC(Exp1_R25_pept)
#'  
#' @importFrom stats cor
#' 
#' @export
#' 
wrapper.corrMatrixD_HC <- function(obj, rate=0.5, showValues=TRUE){
  if (is.null(obj)) {
    warning('The dataset is NULL and cannot be shown')
    return(NULL)
  } else if (nrow(obj) == 0) {
    warning('The dataset is empty and cannot be shown')
    return(NULL)
  }
  
  qData <- Biobase::exprs(obj)
  samplesData <- Biobase::pData(obj)
  data <- cor(qData,use = 'pairwise.complete.obs')
  corrMatrixD_HC(data,samplesData, rate, showValues)
}




#' Correlation matrix based on a \code{MSnSet} object.
#' 
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' res <- cor(qData,use = 'pairwise.complete.obs')
#' corrMatrixD_HC(res, samplesData)
#' 
#' @import highcharter
#' @importFrom dplyr tbl_df mutate left_join select
#' @importFrom tidyr gather
#' @importFrom tibble tibble
#' @importFrom stats cor
#' 
#' @export
#' 
corrMatrixD_HC <- function(object,samplesData = NULL, rate = 0.5, showValues=TRUE) {
  
  df <- as.data.frame(object)
  
  if (!is.null(samplesData)){
    for (j in 1:ncol(df)){
      names(df)[j] <- paste(as.character(samplesData[j,2:ncol(samplesData)]), 
                            collapse =" ")
    }
  }
  is.num <- sapply(df, is.numeric)
  df[is.num] <- lapply(df[is.num], round, 2)
  dist <- NULL
  
  x <- y <- names(df)
  
  df <- dplyr::tbl_df(cbind(x = y, df)) %>% 
    tidyr::gather(y, dist, -x) %>% 
    dplyr::mutate(x = as.character(x),
                  y = as.character(y)) %>% 
    dplyr::left_join(tibble(x = y,
                            xid = seq(length(y)) - 1), by = "x") %>% 
    dplyr::left_join(tibble(y = y,
                            yid = seq(length(y)) - 1), by = "y")
  
  ds <- df %>% 
    dplyr::select("xid", "yid", "dist") %>% 
    list_parse2()
  
  fntltp <- JS("function(){
                  return this.series.xAxis.categories[this.point.x] + ' ~ ' +
                         this.series.yAxis.categories[this.point.y] + ': <b>' +
                         Highcharts.numberFormat(this.point.value, 2)+'</b>';
               ; }")
  cor_colr <- list( list(0, '#FF5733'),
                    list(0.5, '#F8F5F5'),
                    list(1, '#2E86C1')
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
      )) %>% 
    hc_tooltip(formatter = fntltp) %>% 
    hc_legend(align = "right", layout = "vertical",
              verticalAlign="middle") %>% 
    hc_colorAxis(  stops= cor_colr,min=rate,max=1) %>%
    my_hc_ExportMenu(filename = "corrMatrix")
}