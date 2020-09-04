



#' @title Builds a boxplot from a dataframe
#' 
#' @param obj xxx
#' 
#' @param conds xxx
#' 
#' @param legend A vector of the conditions (one string per sample).
#' 
#' @param palette xxx
#' 
#' @return A boxplot
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @seealso \code{\link{densityPlotD}}
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' boxPlotD(Exp1_R25_pept, conds)
#' 
#' @importFrom Biobase exprs pData
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#' 
boxPlotD <- function(obj,conds, legend=NULL,palette=NULL){
  qData <- Biobase::exprs(obj)
  if (is.null(palette)){
    pal <- RColorBrewer::brewer.pal(length(unique(conds)),"Dark2")[1:length(unique(conds))]
    
    for (i in 1:ncol(qData)){
      palette[i] <- pal[ which(conds[i] == unique(conds))]
    }
    
  }else{
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  
  boxplot(qData
          ,las = 1
          , col = palette
          , cex = 2
          , axes=TRUE
          , xaxt = "n"
          , ylab = "Log (intensity)"
          , pt.cex = 4
          , horizontal = FALSE
  )
  
  
  if( is.null(legend)){legend <- Biobase::pData(obj)$Condition}
  axis(side=1,at = 1:ncol(qData), label = legend)
  #mtext("Samples", side=1,  line=(6+ncol(legend)), cex.lab=1, las=1)
  
  abline(h=0) 
  
}





#' @title Builds a boxplot from a dataframe using the library \code{highcharter}
#' 
#' @param obj xxx
#' 
#' @param legend A vector of the conditions (one condition per sample).
#' 
#' @param palette xxx
#' 
#' @return A boxplot
#' 
#' @author Samuel Wieczorek
#' 
#' @seealso \code{\link{densityPlotD_HC}}
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' legend <- Biobase::pData(Exp1_R25_pept)[,"Sample.name"]
#' boxPlotD_HC(Exp1_R25_pept, legend)
#' 
#' @importFrom Biobase exprs
#' @import highcharter
#' 
#' @export
#' 
boxPlotD_HC <- function(obj, legend=NULL, palette = NULL){
  
  qData <- Biobase::exprs(obj)
  if( is.null(legend)){legend <- Biobase::pData(obj)[,"Sample.name"]}
  if (is.null(palette)){palette <- rep("#FFFFFF", ncol(qData))
  } else {
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  bx <-boxplot(qData, na.rm=TRUE)
  df_outlier <- data.frame(x=bx$group-1,y = bx$out)
  
  tmp <- NULL
  for (i in 1:ncol(qData)){
    tmp <- c(tmp, rep(paste(paste0(rep("A", i), collapse=""),legend[i], sep='_'),nrow(qData)))
  }
  
  df <- data.frame(values = as.vector(qData,mode='numeric'),samples = tmp, stringsAsFactors = FALSE)
  
  hc <- highcharter::hcboxplot(x=df$values, var = df$samples, colorByPoint = TRUE, outliers = TRUE) %>%
    hc_chart(type="column") %>%
    hc_yAxis(title = list(text = "Log (intensity)")) %>%
    hc_xAxis(title = list(text = "Samples"), categories=legend) %>%
    hc_colors(palette) %>%
    hc_add_series(type= "scatter",df_outlier) %>%
    hc_tooltip(enabled = FALSE) %>%
    hc_plotOptions(
      
      boxplot= list(
        
        fillColor= c('lightgrey'),
        lineWidth= 3,
        medianColor= 'grey',
        medianWidth= 3,
        stemColor= '#A63400',
        stemDashStyle= 'dot',
        stemWidth= 1,
        whiskerColor= '#3D9200',
        whiskerLength= '20%',
        whiskerWidth= 3
      ),
      scatter = list(
        marker=list(
          fillColor = 'white',
          lineWidth = 0.5,
          lineColor = 'grey'
        )
      )
    )
  
  
  hc
  
}

















#' Builds a heatmap of the quantitative proteomic data of a 
#' \code{MSnSet} object.
#' 
#' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
#' quantitative data in the \code{exprs()} table of an object of
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
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- mvFilter(Exp1_R25_pept[1:1000], "wholeMatrix", 6)
#' wrapper.heatmapD(obj)
#' }
#' 
#' @importFrom Biobase exprs pData
#' 
#' @export
#' 
wrapper.heatmapD  <- function(obj, distance="euclidean", cluster="complete", 
                              dendro = FALSE){
  qData <- Biobase::exprs(obj)
  conds <- Biobase::pData(obj)[['Condition']]
  for (j in 1:length(colnames(qData))){
    colnames(qData)[j] <- paste(as.character(Biobase::pData(obj)[j,2:ncol(Biobase::pData(obj))]), 
                                collapse =" ")
  }
  heatmapD(qData, conds, distance, cluster, dendro)
}




#' Heatmap of the quantitative proteomic data of a \code{MSnSet} object
#' 
#' @title This function is a wrapper to \code{\link{heatmap.2}} that displays 
#' quantitative data in the \code{exprs()} table of an object of
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
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- mvFilter(Exp1_R25_pept[1:1000], "wholeMatrix", 6)
#' qData <- Biobase::exprs(obj)
#' conds <- pData(obj)[['Condition']]
#' heatmapD(qData, conds)
#' }
#' 
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#' 
heatmapD <- function(qData, conds, distance="euclidean", cluster="complete", dendro = FALSE){
  ##Check parameters
  # paramdist <- c("euclidean", "manhattan") 
  # if (!(distance %in% paramdist)){
  #     stop("Param distance is not correct.")
  #     return (NULL)
  # }
  # 
  # paramcluster <- c("ward.D", "average")
  # if (!(cluster %in%  paramcluster)){
  #     stop("Param clustering is not correct.")
  #     return (NULL)
  # }
  
  
  # if (isTRUE(dendro) && getNumberOfEmptyLines(qData) != 0)  {
  #     stop("Your dataset contains empty lines: the dendrogram cannot 
  # be computed.
  #         Please filter or impute missing values before.")
  #     return (NULL)
  # }
  # else {
  .data <- matrix(qData, 
                  ncol = ncol(qData), 
                  byrow = FALSE,
                  dimnames = list(rownames(qData), colnames(qData))
  )
  colors = c(seq(-3, -2, length=100),
             seq(-2, 0.5, length=100),
             seq(0.5, 6, length=100))
  heatmap.color <- colorRampPalette(c("green", "red"))(n = 1000)
  
  # samples label color
  x=t(.data)
  x[is.na(x)] <- -1e5
  dist= dist(x, method=distance)
  hcluster = hclust(dist, method=cluster)
  palette<-NULL
  palette.init <- RColorBrewer::brewer.pal(8,"Dark2")[1:length(unique(conds))]
  for (i in 1:length(conds)){
    palette[i] <- palette.init[which(conds[i] == unique(conds))]
  }
  cols_branches<-palette
  dend1 <- as.dendrogram(hcluster)
  dend1 <- dendextend::color_branches(dend1, k = length(conds), col = cols_branches)
  col_labels <- dendextend::get_leaves_branches_col(dend1)
  
  if (dendro){ .dendro = "row"} else {.dendro = "none"}
  p <- gplots::heatmap.2(
    x=t(.data),
    distfun = function(x) {
      x[is.na(x)] <- -1e5
      dist(x, method=distance)
    },
    hclustfun = function(x) {
      x[is.na(x)] <- -1e5
      hclust(x, method=cluster)
    },
    dendrogram =.dendro,
    Rowv=TRUE,
    col=heatmap.color ,
    density.info='none',
    key=TRUE,
    trace="none",
    scale="none",
    #srtCol=45,
    labCol="",
    margins=c(4,12),
    cexRow= 1.5 + ncol(.data)*-0.011,
    keysize = 1.5,
    lhei = c(1.5, 9),
    lwid = c(1.5, 4),
    lmat = rbind(4:3, 2:1),
    colRow = col_labels
    
  )
  #    }
}


#' Heatmap inspired by the heatmap.2 function.
#' 
#' @title This function is inspired from the function \code{\link{heatmap.2}} 
#' that displays quantitative data in the \code{exprs()} table of an object of
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- mvFilter(Exp1_R25_pept, "wholeMatrix", 6)
#' qData <- Biobase::exprs(obj)
#' heatmap.DAPAR(qData)
#' 
#' @export
#' 
heatmap.DAPAR <- 
  function (x, 
            col = heat.colors(100),
            srtCol=NULL,
            labCol = NULL,
            labRow = NULL,
            key = TRUE, 
            key.title = NULL,
            main = NULL,  
            ylab = NULL) 
  {
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    
    offsetCol <- 0.5
    offsetRow = 0.5
    srtRow = NULL
    colRow = NULL
    colCol = NULL 
    xlab = NULL
    key.par = list()
    margins = c(5, 5)
    sepcolor = "white"
    na.color = "white"
    keysize = 1.5
    breaks <- NULL
    na.rm = TRUE
    
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
      stop("`x' must have at least 2 rows and 2 columns")
    x <- x[nr:1,]
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    cexCol = 0.2 + 1/log10(nc)
    cexRow = 0.2 + 1/log10(nr)
    iy <- 1:nr
    breaks <- length(col) + 1
    breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                  length = breaks)
    
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
    
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    
    
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
            c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
          breaks = breaks)
    
    
    if (!is.null(labCol)) 
    {
      axis(1, 1:nc, label = labCol, las = 2, line = -0.5 + 
             offsetCol, tick = 0, cex.axis = cexCol, hadj = NA, 
           padj = 0)
    } else {
      adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, label = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
             strheight("M"), label = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, col = colCol)
      par(xpd = xpd.orig)
    }
    
    
    if (!is.null(labRow) ) {
      axis(4, iy, label = labRow, las = 5, line = -0.5 + offsetRow, 
           tick = 0, cex.axis = cexRow, hadj = 0, padj = NA)
    }
    else {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, label = rep("", nr), las = 2, 
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
           y = ypos, label = labRow, adj = c(0,NA), cex = cexRow, 
           srt = srtRow, col = colRow)
      par(xpd = xpd.orig)
    }
    
    par(mar = c(margins[1], 0, 0, 0))
    plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    
    plot.new()
    if (!is.null(main)) 
      title(main, cex.main = 1.5 * op[["cex.main"]])
    
    
    if (key) {
      mar <- c(5, 4, 2, 1)
      par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
      if (length(key.par) > 0) 
        do.call(par, key.par)
      
      tmpbreaks <- breaks
      min.raw <- min.breaks
      max.raw <- max.breaks
      
      z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at = xv, label = lv)
      
      xargs$side <- 1
      do.call(axis, xargs)
      key.xlab <- "Intensity value"
      
      mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
            cex = par("cex") * par("cex.lab"))
      
      if (is.null(key.title)) 
        title("Color Key")
    }
    
  }



# 
# rep.row <-function(x,n){
#   matrix(rep(x,each=n),nrow=n)
# }
# 
# rep.col<-function(x,n){
#   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
# }

###--------------------------------------------------------------------
# heatmap_HC <- function(qData, col=heat.colors(100),labCol)
# {
#   conds_v <- c(rep.row(colnames(qData), nrow(qData)))
#   lines_v <- c(rep.col(1:nrow(qData), ncol(qData)))
#   data <- tibble(id=lines_v, condition=conds_v, value=round(c(qData), digits=2 ))
#   
# #   fntltp <- JS("function(){
# #                return this.point.x + ' ' +  this.series.yAxis.categories[this.point.y] + ':<br>' +
# #                Highcharts.numberFormat(this.point.value, 2);
# # }")
# 
#   # plotline <- list(
#   #   color = "#fde725", value = 3, width = 2, zIndex = 5,
#   #   label = list(
#   #     text = "", verticalAlign = "top",
#   #     style = list(color = "black"), textAlign = "left",
#   #     rotation = 0, y = -5)
#   # )
#   
#  # highchart2() %>%
#     hchart(data, "heatmap", hcaes(x = condition, y = id, value = value)) %>% 
#     hc_colorAxis(stops = color_stops(100, col),type = "linear") %>% 
#     # hc_yAxis(reversed = FALSE, offset = -20, tickLength = 0,
#     #          gridLineWidth = 0, minorGridLineWidth = 0,
#     #          labels = list(style = list(fontSize = "8px"))) %>% 
#     hc_tooltip(enabled = FALSE) %>% 
#     #hc_xAxis(plotLines = list(plotline)) %>%
#     hc_title(text = "MEC repartition") %>% 
#     hc_legend(layout = "vertical", verticalAlign = "top",
#               align = "left", valueDecimals = 0)
#   
# }