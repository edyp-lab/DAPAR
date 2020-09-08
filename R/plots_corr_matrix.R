#' 
#' 
#' #' Builds a correlation matrix based on a \code{MSnSet} object.
#' #' 
#' #' @title Displays a correlation matrix of the quantitative data of the
#' #' \code{exprs()} table
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
#' #' @importFrom Biobase exprs pData
#' #' 
#' #' @export
#' #' 
#' wrapper.corrMatrixD <- function(obj, rate=5){
#'   qData <- Biobase::exprs(obj)
#'   samplesData <- Biobase::pData(obj)
#'   corrMatrixD(qData, samplesData, rate)
#' }

#' Builds a correlation matrix based on a \code{MSnSet} object. 
#' Same as the function \code{\link{wrapper.corrMatrixD}} but uses the package \code{highcharter}
#' 
#' @title Displays a correlation matrix of the quantitative data of the
#' \code{exprs()} table
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param rate A float that defines the gradient of colors.
#' 
#' @param showDataLabels xxx
#' 
#' @return A colored correlation matrix
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' wrapper.corrMatrixD_HC(Exp1_R25_pept)
#'  
#' @importFrom Biobase exprs pData
#' @importFrom stats cor
#' 
#' @export
#' 
wrapper.corrMatrixD_HC <- function(obj, rate=0.5, showDataLabels=TRUE){
  qData <- Biobase::exprs(obj)
  samplesData <- Biobase::pData(obj)
  data <- cor(qData,use = 'pairwise.complete.obs')
  corrMatrixD_HC(data,samplesData, rate, showDataLabels)
}

#' 
#' 
#' #' Correlation matrix based on a \code{MSnSet} object
#' #' 
#' #' @title Displays a correlation matrix of the quantitative data of the
#' #' \code{exprs()} table.
#' #' 
#' #' @param qData A dataframe of quantitative data.
#' #' 
#' #' @param samplesData A dataframe where lines correspond to samples and 
#' #' columns to the meta-data for those samples.
#' #' 
#' #' @param gradientRate The rate parameter to control the exponential law for 
#' #' the gradient of colors
#' #' 
#' #' @return A colored correlation matrix
#' #' 
#' #' @author Florence Combes, Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' qData <- Biobase::exprs(Exp1_R25_pept)
#' #' samplesData <- Biobase::pData(Exp1_R25_pept)
#' #' corrMatrixD(qData, samplesData)
#' #' 
#' #' @export
#' #' 
#' #' @importFrom stats pexp cor
#' #' @importFrom ggplot2 element_text qplot theme labs scale_fill_gradientn
#' #' @importFrom reshape2 melt
#' #' 
#' corrMatrixD <- function(qData, samplesData, gradientRate = 5){
#'   Var1 <- Var2 <- value <- NULL
#'   
#'   for (j in 1:length(colnames(qData))){
#'     colnames(qData)[j] <- paste(as.character(samplesData[j,2:ncol(samplesData)]), 
#'                                 collapse =" ")
#'   }
#'   
#'   z <- cor(qData,use = 'pairwise.complete.obs')
#'   text <- element_text(colour="black", size = 16, face = "bold")
#'   d <- qplot(x = Var1, 
#'              y = Var2, 
#'              data = melt(z), 
#'              fill = value, 
#'              geom = "tile") +
#'     theme(axis.text = element_text(size=16),
#'           axis.title = element_text(size=20, face="bold"),
#'           axis.text.x = element_text(angle=30, vjust=1, hjust=1),
#'           legend.text = text,
#'           legend.title = text) +
#'     labs(x = "", y = "") +
#'     
#'     scale_fill_gradientn (
#'       colours = colorRampPalette (c ("white", "lightblue","darkblue")) (101),
#'       values = c(pexp(seq(0,1,0.01), rate=gradientRate),1), limits=c(0,1))
#'   
#'   plot(d)
#' }
#' 



#' Correlation matrix based on a \code{MSnSet} object. Same as the 
#' function \link{corrMatrixD} but uses the package \code{highcharter}
#' 
#' @title Displays a correlation matrix of the quantitative data of the
#' \code{exprs()} table.
#' 
#' @param object The result of the \code{cor} function.
#' 
#' @param samplesData A dataframe in which lines correspond to samples and 
#' columns to the meta-data for those samples.
#' 
#' @param rate The rate parameter to control the exponential law for 
#' the gradient of colors
#' 
#' @param showDataLabels xxx
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
corrMatrixD_HC <- function(object,samplesData = NULL, rate = 0.5, showDataLabels=TRUE) {
  
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
        dataLabels = list(enabled = showDataLabels)
      )) %>% 
    hc_tooltip(formatter = fntltp) %>% 
    hc_legend(align = "right", layout = "vertical",
              verticalAlign="middle") %>% 
    hc_colorAxis(  stops= cor_colr,min=rate,max=1) %>%
    my_hc_ExportMenu(filename = "corrMatrix")
}