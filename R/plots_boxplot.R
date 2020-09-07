#' Boxplot for quantitative proteomics data using the library \code{highcharter}
#' 
#' @title Builds a boxplot from a dataframe using the library \code{highcharter}
#' 
#' @param qData Numeric matrix 
#' 
#' @param conds xxx
#' 
#' @param keyId xxxx
#' 
#' @param legend A vector of the conditions (one condition per sample).
#' 
#' @param palette xxx
#' 
#' @param subset.view A vector of index indicating rows to highlight
#' 
#' @return A boxplot
#' 
#' @author Samuel Wieczorek, Anais Courtier, Enora Fremy
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs((Exp1_R25_pept))
#' conds <- legend <- Biobase::pData(Exp1_R25_pept)[["Condition"]]
#' key <- "Protein_group_IDs"
#' boxPlotD_HC(qData, conds, key, legend, c(rep('blue',3), rep('green',3)), 1:10)
#' 
#' @import highcharter
#' 
#' @importFrom grDevices colorRampPalette boxplot.stats
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @importFrom stats na.exclude
#' 
#' @export
#' 
boxPlotD_HC <- function(qData, conds, keyId=NULL, legend=NULL, palette = NULL, subset.view=NULL){
  
  if (is.null(qData)){
    warning('The dataset in NULL and cannot be shown')
    return(NULL)
  }
  
  if(missing(conds))
    stop("'conds' is missing.")
  
  if (is.null(legend)) {
    legend <- conds
    for (i in unique(conds))
      legend[which(conds==i)] <- paste0(i, '_', 1:length(which(conds==i)))
  }
  
  if (!is.null(subset.view)) {
    if (is.null(keyId)|| missing(keyId))
      stop("'keyId' is missing.")
  }
  
  
  add_variable_to_series_list <- function(x, series_list, key_vector, value_vector){
    base::stopifnot(length(key_vector) == length(value_vector))
    base::stopifnot(length(series_list) == length(key_vector))
    series_list[[x]][length(series_list[[x]])+1]<- value_vector[x]
    names(series_list[[x]])[length(series_list[[x]])]<- key_vector[x]
    return(series_list[[x]])
  }
  
  
  # From highcharter github pages:
  hc_add_series_bwpout <- function(hc, value, by, ...) {
    z = lapply(levels(by), function(x) {
      bpstats = boxplot.stats(value[by == x])$stats
      outliers = c()
      for (y in na.exclude(value[by == x])) {
        if ((y < bpstats[1]) | (y > bpstats[5]))
          outliers = c(outliers, list(which(levels(by)==x)-1, y))
      }
      outliers
    })
    hc %>%
      hc_add_series(data = z, type="scatter", ...)
  }
  
  
  gen_key_vector <- function(variable, num_times){
    return(rep(variable, num_times))
  } 
  
  
  gen_boxplot_series_from_df <- function(value, by,...){
    value<- base::as.numeric(value)
    by<- base::as.factor(by)
    box_names<- levels(by)
    z=lapply(box_names, function(x) {
      boxplot.stats(value[by==x])$stats
    })
    tmp<- lapply(seq_along(z), function(x){
      var_name_list<- list(box_names[x])
      #tmp0<- list(names(df)[x])
      names(var_name_list)<- "name"
      index<- x-1
      tmp<- list(c(index,  z[[x]]))
      tmp<- list(tmp)
      names(tmp)<- "data"
      tmp_out<- c(var_name_list, tmp)
      #tmp<- list(tmp)
      return(tmp_out)
      
    })
    return(tmp)
  }
  
  
  # Usage: 
  #series<- gen_boxplot_series_from_df(value = df$total_value, by=df$asset_class)
  
  
  ## Boxplot function:
  make_highchart_boxplot_with_colored_factors <- function(value, by, chart_title="Boxplots",
                                                         chart_x_axis_label="Values", show_outliers=FALSE,
                                                         boxcolors=NULL, box_line_colors=NULL){
    by <- as.factor(by)
    box_names_to_use <- levels(by)
    series <- gen_boxplot_series_from_df(value = value, by=by)
    
    if(is.null(boxcolors)){
      cols <- rep("#FFFFFF", ncol(qData))
      #cols<- viridisLite::viridis(n= length(series), alpha = 0.5) # Keeping alpha in here! (COLORS FOR BOXES ARE SET HERE)
    } else {
      cols <- boxcolors
    }
    
    
    if(is.null(box_line_colors)){
      if(base::nchar(cols[[1]])==9){
        cols2<- substr(cols, 0,7) # no alpha, pure hex truth, for box lines 
      } else {
        cols2<- cols
      }
      
    } else {
      cols2<- box_line_colors
    }
    
    # Injecting value 'fillColor' into series list
    kv<- gen_key_vector(variable = "fillColor", length(series)) 
    series2<- lapply(seq_along(series), function(x){ add_variable_to_series_list(x = x, series_list = series, key_vector = kv, value_vector = cols) })
    
    if(show_outliers == TRUE){
      hc<- highcharter::highchart() %>%
        highcharter::hc_chart(type="boxplot", inverted=FALSE) %>%
        highcharter::hc_title(text=chart_title) %>%
        highcharter::hc_legend(enabled=FALSE) %>%
        highcharter::hc_xAxis(type="category", categories=box_names_to_use, title=list(text=chart_x_axis_label)) %>%
        highcharter::hc_yAxis(title = list(text = "Log (intensity)")) %>%
        highcharter::hc_add_series_list(series2) %>%
        hc_add_series_bwpout(value = value, by=by, name="Outliers", colorByPoint = TRUE) %>%
        hc_plotOptions(series = list(
          marker = list(
            symbol = "circle"
          ),
          grouping=FALSE
        )) %>%
        highcharter::hc_colors(cols2) %>%
        highcharter::hc_exporting(enabled=TRUE)
      
    } else{
      hc<- highcharter::highchart() %>%
        highcharter::hc_chart(type="boxplot", inverted=FALSE) %>%
        highcharter::hc_title(text=chart_title) %>%
        highcharter::hc_legend(enabled=FALSE) %>%
        highcharter::hc_xAxis(type="category", categories=box_names_to_use, title=list(text=chart_x_axis_label)) %>%
        highcharter::hc_yAxis(title = list(text = "Log (intensity)")) %>%
        highcharter::hc_add_series_list(series2) %>%
        hc_plotOptions(series = list(
          marker = list(
            symbol = "circle"
          ),
          grouping=FALSE
        )) %>%
        highcharter::hc_colors(cols2) %>%
        highcharter::hc_exporting(enabled=TRUE)
    }
    hc
  }
  
  #palette <- BuildPalette(conds, palette)
  if (is.null(palette)){palette <- rep("#FFFFFF", ncol(qData))
  } else {
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  df <- data.frame(cbind(categ =rep(colnames(qData),nrow(qData)),
                         value = as.vector(apply(qData, 1, function(x) as.vector(x)))))
  df$value<- base::as.numeric(df$value)
  hc <- make_highchart_boxplot_with_colored_factors(value = df$value, 
                                                    by=df$categ, 
                                                    chart_title = "",
                                                    chart_x_axis_label = "Samples",
                                                    show_outliers = TRUE, 
                                                    boxcolors = palette,
                                                    box_line_colors = "black")
  
    # Display of rows to highlight (index of row in subset.view) 
  if(!is.null(subset.view)){
    idVector <- keyId
    pal=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(subset.view))    
    n=0
    for(i in subset.view){
      n=n+1
      dfSubset <- data.frame(y = as.vector(qData[i,],mode='numeric'), x = as.numeric(factor(names(qData[i,])))-1, stringsAsFactors = FALSE)
      hc<-hc %>%
        highcharter::hc_add_series(type= "line",data=dfSubset,color=pal[n], dashStyle = "shortdot",name=idVector[i],
                                   tooltip=list(enabled=T,headerFormat ="",pointFormat="{point.series.name} : {point.y: .2f} ") )
      
    }
  }  
  
  hc
  
}

