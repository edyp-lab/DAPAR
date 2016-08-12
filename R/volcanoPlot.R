
##' Plots a volcanoplot after the differential analysis.
##' Typically, the log of Fold Change is represented on the X-axis and the
##' log10 of the p-value is drawn on the Y-axis. When the \code{threshold_pVal}
##' and the \code{threshold_logFC} are set, two lines are drawn respectively on
##' the y-axis and the X-axis to visually distinguish between differential and
##' non differential data.
##' 
##' @title Volcanoplot of the differential analysis
##' @param logFC A vector of the log(fold change) values of the differential
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
##' data(UPSpep25)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' data <- wrapper.diffAnaLimma(UPSpep25, condition1, condition2)
##' diffAnaVolcanoplot(data$logFC, data$P.Value)
diffAnaVolcanoplot <- function(logFC=NULL, 
                                pVal=NULL, 
                                threshold_pVal=1e-60, 
                                threshold_logFC=0, 
                                conditions=NULL){

xtitle <- paste("log2 ( mean(",conditions[2],") / mean(",conditions[1],") )",
                sep="")


if (is.null(logFC)||is.null(pVal)) {

    p <- plot(-1,-1
            , xlab = xtitle
            , ylab="- log10 ( p-value )"
            ,xlim = range(0,1)
            , xaxt='n'
            ,yaxt='n'
    )
    return (NULL)
}

x <- logFC
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
