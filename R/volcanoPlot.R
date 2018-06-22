
##' Plots a volcanoplot after the differential analysis.
##' Typically, the log of Fold Change is represented on the X-axis and the
##' log10 of the p-value is drawn on the Y-axis. When the \code{threshold_pVal}
##' and the \code{threshold_logFC} are set, two lines are drawn respectively on
##' the y-axis and the X-axis to visually distinguish between differential and
##' non differential data.
##' 
##' @title Volcanoplot of the differential analysis
##' @param FC A vector of the log(fold change) values of the differential
##' analysis.
##' @param pVal A vector of the p-value values returned by the differential
##' analysis.
##' @param threshold_pVal A floating number which represents the p-value that
##' separates differential and non-differential data.
##' @param threshold_logFC A floating number which represents the log of the
##' Fold Change that separates differential and non-differential data.
##' @param conditions A list of the names of condition 1 and 2 used for the
##' differential analysis.
##' @return A volcanoplot
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' keepThat <- mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' qData <- Biobase::exprs(obj)
##' sTab <- Biobase::pData(obj)
##' limma <- limmaCompleteTest(qData,sTab)
##' diffAnaVolcanoplot(limma$FC[,1], limma$P_Value[,1])
diffAnaVolcanoplot <- function(FC=NULL, 
                                pVal=NULL, 
                                threshold_pVal=1e-60, 
                                threshold_logFC=0, 
                                conditions=NULL){

xtitle <- paste("log2 ( mean(",conditions[2],") / mean(",conditions[1],") )",
                sep="")


if (is.null(FC)||is.null(pVal)) {

    p <- plot(-1,-1
            , xlab = xtitle
            , ylab="- log10 ( p-value )"
            ,xlim = range(0,1)
            , xaxt='n'
            ,yaxt='n'
    )
    return (NULL)
}

x <- FC
y <- -log10(pVal)

colorCode <- c("gray", "orange")
color <- rep(colorCode[1], length(y))

for (i in 1:length(y)){
    if ( (y[i] >= threshold_pVal) && (abs(x[i]) >= threshold_logFC) ){
    color[i] <- colorCode[2]
    }
}

p <- plot(x
            , y
            , xlab = xtitle
            , ylab=  "- log10 ( p-value )"
            , xlim = c(-max(abs(x)), max(abs(x)))
            , ylim = c(0, max(y))
            , col = color
            , pch = 16
            , las = 1
            , cex = 1
            , cex.lab = 1.5
            , cex.axis = 1.5
            , cex.main = 3
)

abline(h = threshold_pVal, col = "gray")
abline(v = threshold_logFC, col = "gray")
abline(v = -threshold_logFC, col = "gray")

return(p)
}


##' Plots an interactive volcanoplot after the differential analysis.
##' Typically, the log of Fold Change is represented on the X-axis and the
##' log10 of the p-value is drawn on the Y-axis. When the \code{threshold_pVal}
##' and the \code{threshold_logFC} are set, two lines are drawn respectively on
##' the y-axis and the X-axis to visually distinguish between differential and
##' non differential data. With the use of the package Highcharter, a 
##' customizable tooltip appears when the user put the mouse's pointer over 
##' a point of the scatter plot.
##' 
##' @title Volcanoplot of the differential analysis
##' @param df A dataframe which contains the following slots :
##' x : a vector of the log(fold change) values of the differential analysis,
##' y : a vector of the p-value values returned by the differential analysis.
##' index : a vector of the rowanmes of the data.
##' This dataframe must has been built with the option stringsAsFactors set 
##' to FALSE. There may be additional slots which will be used to show 
##' informations in the tooltip. The name of these slots must begin with the 
##' prefix "tooltip_". It will be automatically removed in the plot.
##' @param threshold_pVal A floating number which represents the p-value that
##' separates differential and non-differential data.
##' @param threshold_logFC A floating number which represents the log of the
##' Fold Change that separates differential and non-differential data.
##' @param conditions A list of the names of condition 1 and 2 used for the
##' differential analysis.
##' @param clickFunction A string that contains a JavaScript function used to 
##' show info from slots in df. The variable this.index refers to the slot 
##' named index and allows to retrieve the right row to show in the tooltip
##' @return An interactive volcanoplot
##' @author Samuel Wieczorek
##' @examples
##' library(highcharter) 
##' library(tidyverse)
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' keepThat <- mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' qData <- Biobase::exprs(obj)
##' sTab <- Biobase::pData(obj)
##' data <- limmaCompleteTest(qData,sTab)
##' df <- data.frame(x=data$FC, y = -log10(data$P_Value),index = as.character(rownames(obj)))
##' colnames(df) <- c("x", "y", "index")
##' tooltipSlot <- c("Sequence", "Score")
##' df <- cbind(df,Biobase::fData(obj)[tooltipSlot])
##' colnames(df) <- gsub(".", "_", colnames(df), fixed=TRUE)
##' if (ncol(df) > 3){
##'     colnames(df)[4:ncol(df)] <- 
##'     paste("tooltip_", colnames(df)[4:ncol(df)], sep="")}
##' hc_clickFunction <- JS("function(event) {
##' Shiny.onInputChange('eventPointClicked', [this.index]);}")
##' cond <- c("25fmol", "10fmol")
##' diffAnaVolcanoplot_rCharts(df, 2.5, 1, cond,hc_clickFunction) 
diffAnaVolcanoplot_rCharts <- function(df, 
                                        threshold_pVal=1e-60, 
                                        threshold_logFC=0, 
                                        conditions=NULL, 
                                        clickFunction=NULL){
    
    xtitle <- paste("log2 ( mean(",
                    conditions[2],
                    ") / mean(",
                    conditions[1],
                    ") )",
                    sep="")
    
    if (is.null(clickFunction)){
        clickFunction <- 
            JS("function(event) {Shiny.onInputChange('eventPointClicked', [this.index]+'_'+ [this.series.name]);}")
    }
    
    
    df <- cbind(df, 
                g=ifelse(df$y >= threshold_pVal & abs(df$x) >= threshold_logFC, "g1", "g2")
    )
    
    
    i_tooltip <- which(startsWith(colnames(df),"tooltip"))
    txt_tooltip <- NULL
    for (i in i_tooltip){
        t <- txt_tooltip <- paste(txt_tooltip,"<b>",gsub("tooltip_", "", 
                                                        colnames(df)[i], 
                                                        fixed=TRUE), 
                                 " </b>: {point.", colnames(df)[i],"} <br> ", 
                                 sep="")
    }
    
    h1 <-  hchart(df, "scatter", hcaes(x,y,group=g)) %>%
        hc_colors(c("orange", "grey")) %>%
        my_hc_chart(zoomType = "xy",chartType="scatter") %>%
        hc_legend(enabled = FALSE) %>%
        hc_yAxis(title = list(text="-log10(pValue)"),
                 plotBands = list(list(from= 0, to = threshold_pVal, color = "lightgrey")),
                 plotLines=list(list(color= "grey" , width = 2, value = 0, zIndex = 5))) %>%
        hc_xAxis(title = list(text = "FC"),
                 plotBands = list(list(from= -threshold_logFC, to = threshold_logFC, color = "lightgrey")),
                 plotLines=list(list(color= "grey" , width = 2, value = 0, zIndex = 5))) %>%
        hc_tooltip(headerFormat= '',pointFormat = txt_tooltip) %>%
        hc_plotOptions( series = list( animation=list(duration = 100),
                                       cursor = "pointer", 
                                       point = list( events = list( 
                                           click = clickFunction ) ) ) ) %>%
        my_hc_ExportMenu(filename = "volcanoplot")
    
    return(h1)
}
