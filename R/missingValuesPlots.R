##' This method plots a histogram which represents the distribution of the 
##' number of missing values (NA) per lines (ie proteins).
##' 
##' @title Histogram of missing values per lines
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend The indice of the column name's in \code{pData()} tab 
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A histogram
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' mvPerLinesHisto(UPSprotx2)
mvPerLinesHisto <- function(obj, indLegend="auto", showValues=TRUE){
    .data <- exprs(obj)
if (indLegend == "auto") { indLegend <- c(2:length(colnames(pData(obj))))}
    for (j in 1:length(colnames(.data))){
        noms <- NULL
        for (i in 1:length(indLegend)){
            noms <- paste(noms, pData(obj)[j,indLegend[i]], sep=" ")
            }
        colnames(.data)[j] <- noms
    }

    coeffMax <- .1

    NbNAPerCol <- colSums(is.na(.data))
    NbNAPerRow <- rowSums(is.na(.data))
    #par(mar = c(10,3, 3, 3))

    nb.col <- dim(.data)[2] 
    nb.na <- NbNAPerRow
    temp <- table(NbNAPerRow)
    nb.na2barplot <- c(temp,rep(0,1+ncol(.data)-length(temp)))

    if (sum(NbNAPerRow) == 0){
        nb.na2barplot <- rep(0,1+ncol(.data))
    }
        x <- barplot(nb.na2barplot[-1], 
                main = "# lines by # of NA",
                xlab = "# NA per lines",
                names.arg = as.character(c(1:(ncol(.data)))), 
                col = c(rep("lightgrey",nb.col-1), "red"),
                ylim = c(0, max(1,nb.na2barplot[-1])), 
                las=1,
                cex.names=1.5,
                cex.axis=1.5
                )

        if (showValues) {
                text(x, 
#                 nb.na2barplot[-1] + coeffMax*(max(nb.na2barplot[-1])), 
                0,
                nb.na2barplot[-1], 
                cex= 1.5,
                pos = 3,
                col = rep("black",nb.col)) 
    }
}


##' This method plots a histogram of missing values.
##' 
##' @title Histogram of missing values
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param indLegend The indices of the column name's in \code{pData()} tab
##' @param showValues A logical that indicates wether numeric values should be
##' drawn above the bars.
##' @return A histogram
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' mvHisto(UPSprotx2,showValues=TRUE)
mvHisto <- function(obj, indLegend="auto", showValues=FALSE){

        if (indLegend == "auto") { 
            indLegend <- c(2:length(colnames(pData(obj))))
        }

    .data <- exprs(obj)
    colnames(.data) <- pData(obj)$Label

    coeffMax <- .1
    pal <- getPaletteForLabels(obj)

    NbNAPerCol <- colSums(is.na(.data))
    NbNAPerRow <- rowSums(is.na(.data))

    if (sum(NbNAPerCol) == 0) {if (sum(NbNAPerCol) == 0){
        NbNAPerCol <- rep(0,1+ncol(.data))
    }} 
    x <- barplot(NbNAPerCol, 
                main = "# NA per columns",
                col=pal,
                las=1,
                ylim = c(0, max(1,NbNAPerCol)),
                names.arg = c(1:18), 
                cex.names=1.5,
                cex.axis=1.5,
                axisnames = FALSE
                )

        par(xpd = TRUE)
        text(x, -5,
            labels = colnames(.data),
            srt = 45,
            adj=1,
            cex=1.4)

        if (showValues) {
            text(x, 
                0, 
                NbNAPerCol, 
                cex= 1.5,
                pos = 3) 
        }

}


##' Plots a heatmap of the quantitative data. Each column represent one of
##' the conditions in the object of class \code{\link{MSnSet}} and 
##' the color is proportional to the mean of intensity for each line of
##' the dataset.
##' The lines have been sorted in order to vizualize easily the different
##' number of missing values. A white square is plotted for missing values.
##' 
##' @title Heatmap of missing values
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A heatmap
##' @author Samuel Wieczorek, Thomas Burger
##' @examples data(UPSprotx2)
##' mvImage(UPSprotx2)
mvImage <- function(obj){
    ### build indices of conditions
    indCond <- list()
    LabelNames <- unique(pData(obj)$Label)
    for (i in LabelNames) {
        indCond <- append(indCond, list(which(i == pData(obj)$Label)))
    }
    indCond <- setNames(indCond, as.list(c("cond1", "cond2")))

    fData(obj) <- cbind(fData(obj),
                        nNA1 = apply(as.matrix(exprs(obj)[,indCond$cond1]), 1, 
                                    function(x) sum(is.na(x))))
    fData(obj) <- cbind(fData(obj),
                        nNA2 = apply(as.matrix(exprs(obj)[,indCond$cond2]), 1,
                                    function(x) sum(is.na(x))))
    o <- order(((fData(obj)$nNA1 +1)^2) / (fData(obj)$nNA2 +1))
    yy <- obj[o,]

    for (i in 1:nrow(yy)){
        k <- order(exprs(yy)[i,indCond$cond1])
        exprs(yy)[i,rev(indCond$cond1)] <- exprs(yy)[i, k]
        .temp <- mean(exprs(yy)[i,rev(indCond$cond1)], na.rm = TRUE)
        exprs(yy)[i,which(!is.na(exprs(yy)[i,indCond$cond1]))] <- .temp

        k <- order(exprs(yy)[i,indCond$cond2])
        exprs(yy)[i,indCond$cond2] <- exprs(yy)[i, k+length(indCond$cond1)]
        .temp <- mean(exprs(yy)[i,indCond$cond2], na.rm = TRUE)
        exprs(yy)[i,length(indCond$cond1) + 
                    which(!is.na(exprs(yy)[i,indCond$cond2]))] <- .temp
    }
    
    colfunc <- colorRampPalette(c("yellow", "red"))

    q <- heatmap.2(exprs(yy),
                dendrogram = "none",
                Rowv=FALSE,
                Colv = FALSE,
                col = colfunc(100),
                density.info='none',
                key=TRUE,
                trace="none",
                scale="none",
                srtCol= 0,
                labCol=pData(obj)$Label,
                labRow = "",
                ylab = "Peptides / proteins",
                main = "Missing values heatmap"
    ) 
}



##' This method plots a scatter plot which represents the distribution of
##' missing values.
##' The colors correspond to the different conditions (slot Label in in the
##' dataset of class \code{\link{MSnSet}}).
##' The x-axis represent the mean of intensity for one condition and one
##' entity in the dataset (i. e. a protein) 
##' whereas the y-axis count the number of missing values for this entity
##' and the considered condition.
##' The data have been jittered for an easier vizualisation.
##' 
##' @title Distribution of missing values with respect to intensity values
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param threshold An integer for the intensity that delimits MNAR and 
##' MCAR missing values.
##' @return A scatter plot
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' mvTypePlot(UPSprotx2)
mvTypePlot <- function(obj,threshold=0){

    pal <- unique(getPaletteForLabels(obj))
    color <- NULL
    col.legend <- c(1:length(pal))
    data <- exprs(obj)
    
    conditions <- pData(obj)[,"Label"]
    i <- 1
    d <- data[,which(conditions==unique(conditions)[i])]
    color <- rep(pal[i],nrow(d))
    
    for (i in 2:length(unique(conditions))){
        d <- rbind(d, data[,which(conditions==unique(conditions)[i])])
        color <- rbind(color, rep(pal[i],
                        nrow(data[,which(conditions==unique(conditions)[i])])))
    }
    
    total <- apply(d,1,mean ,na.rm=TRUE)
    nbNA <- apply(d,1,function(x) length(which(is.na(x) == TRUE)))
    plot(total,
            rep(-1,length(total)),
            xlim = range(total, na.rm = TRUE),
            ylim = c(0, ncol(data)/length(unique(conditions))),
            xlab = "Mean of quantity values", 
            ylab = "# of missing values",
            main =  "Missing values repartition")
    
    if (sum(nbNA) > 0)
    {
        points(total,
                jitter(nbNA, 0.9), 
                col = color,
                pch = 16,
                cex=0.5)
        
        lines(lowess(nbNA ~ total,
                        delta = 0.1 * diff(range(total, na.rm = TRUE))), 
                col="black", lwd=2)
        abline(v=threshold, col="blue", lwd=3)
    }
    
    legend("topright"         
            , legend = unique(pData(obj)$Label)
            , col = col.legend
            , pch = 15 
            , bty = "n"
            , pt.cex = 2
            , cex = 1
            , horiz = FALSE
            , inset=c(0,0)
    )
}
