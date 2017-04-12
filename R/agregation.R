##' This function computes the number of proteins that are only defined by 
##' specific peptides, shared peptides or a mixture of two. 
##' 
##' @title computes the number of proteins that are only defined by 
##' specific peptides, shared peptides or a mixture of two.
##' @param matUnique The adjacency matrix with only specific peptides.
##' @param matShared The adjacency matrix with both specific and 
##' shared peptides.
##' @return A list
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' MShared <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' MUnique <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, TRUE)
##' getProteinsStats(MUnique,MShared)
getProteinsStats <- function(matUnique, matShared){
    if (is.null(matUnique) || is.null(matShared) ||
        !is.matrix(matUnique) || !is.matrix(matShared)){return(NULL)}
    
    t <- setdiff(union(rownames(matUnique), rownames(matShared)), 
                 intersect(rownames(matUnique), rownames(matShared)))
    sharedPeptides <- matShared[t,]
    sharedPeptides <- sharedPeptides[,-which(colSums(sharedPeptides)==0)]
    protOnlyUnique <- setdiff(union(colnames(sharedPeptides), 
                                    colnames(matShared)), 
                              intersect(colnames(sharedPeptides), 
                                        colnames(matShared)))
    
    protOnlyShared <- setdiff(union(colnames(matUnique), 
                                    colnames(matShared)), 
                              intersect(colnames(matUnique), 
                                        colnames(matShared)))
    a <- union(protOnlyUnique, protOnlyShared)
    b <- colnames(matShared)
    protMix <- setdiff(union(union(protOnlyUnique, protOnlyShared),
                             colnames(matShared)),
                       intersect(union(protOnlyUnique, protOnlyShared),
                                 colnames(matShared)))

    return (list(protOnlyUniquePep =protOnlyUnique,
                  protOnlySharedPep =protOnlyShared,
                  protMixPep = protMix))
}




##' This function creates a column for the protein dataset after agregation 
##' by using the previous peptide dataset.
##' 
##' @title creates a column for the protein dataset after agregation by
##'  using the previous peptide dataset.
##' @param peptideData A data.frame of meta data of peptides. It is the fData 
##' of the MSnset object.
##' @param matAdj The adjacency matrix used to agregate the peptides data.
##' @param columnName The name of the column in fData(peptides_MSnset) that 
##' the user wants to keep in the new protein data.frame.
##' @return A vector
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' data <- Biobase::fData(Exp1_R25_pept[1:1000])
##' name <- "organism"
##' BuildColumnToProteinDataset(data, M, name )
BuildColumnToProteinDataset <- function(peptideData, matAdj, columnName){
nbProt <- ncol(matAdj)
newCol <- rep("", nbProt)

for (p in 1:nbProt){
    listeIndicePeptides <- which(matAdj[,p] == 1)
    listeData <- unique(peptideData[listeIndicePeptides,columnName])
    newCol[p] <- paste(listeData, collapse = ", ")
}
return(newCol)
}





##' This function computes the number of peptides used to aggregate proteins.
##' 
##' @title Compute the number of peptides used to aggregate proteins
##' @param M A "valued" adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @return A vector of boolean which is the adjacency matrix 
##' but with NA values if they exist in the intensity matrix.
##' @author Alexia Dorffer
##' @examples
##' library(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' CountPep(M)
CountPep <- function (M) {
    z <- M
    z[z!=0] <- 1
    return(z)
}



##' Method to create a plot with proteins and peptides on
##' a MSnSet object (peptides)
##' 
##' @title Function to create a histogram that shows the repartition of
##' peptides w.r.t. the proteins
##' @param mat An adjacency matrix.
##' @return A histogram  
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' mat <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein.group.IDs")
##' GraphPepProt(mat)
GraphPepProt <- function(mat){
    if (is.null(mat)){return (NULL)} 

    #mat <- as.matrix(mat)
    t <- t(mat)
    t <- apply(mat, 2, sum, na.rm=TRUE)
    tab <- table(t)
    position <- seq(1, length(tab),by=3)
    label <- names(tab)

    #par(mar=c(6,4,4,8) + 0.1)#, mgp=c(3,0.5,0)
    barplot(tab, 
            xlim=c(1, length(tab)),
            xlab="Nb of peptides", ylab="Nb of proteins",
            names.arg=label, xaxp=c(1, length(tab), 3), las=1
            , col = "orange")

}




##' Method to create a binary matrix with proteins in columns and peptides 
##' in lines on a MSnSet object (peptides)
##' 
##' @title Function matrix of appartenance group
##' @param obj.pep An object (peptides) of class \code{\link{MSnbase}}.
##' @param protID The name of proteins ID column 
##' @param unique A boolean to indicate whether only the unique peptides must 
##' be considered (TRUE) or if the shared peptides have to 
##' be integrated (FALSE).
##' @return A binary matrix  
##' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept) 
##' BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein.group.IDs", TRUE)
BuildAdjacencyMatrix <- function(obj.pep, protID, unique=TRUE){
    
    data <- Biobase::exprs(obj.pep)
    PG <- Biobase::fData(obj.pep)[,protID]
    PG.l <- strsplit(as.character(PG), split=";", fixed=TRUE)
    
    Un1 <- unlist(PG.l)
    X<- sparseMatrix(i = rep(seq_along(PG.l), lengths(PG.l)),
                     j=as.integer(factor(Un1, levels = unique(Un1))),
                     x=1, 
                     dimnames=list(rownames(data),as.character(unique(Un1))))
    
    if (unique == TRUE){
        X <- X[which(rowSums(as.matrix(X))==1),]
        X <- X[,which(colSums(as.matrix(X))>0)]
    }
    
    return(X)
}



##' This function computes the intensity of proteins based on the sum of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins with the sum of the intensities
##' of their peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' SumPeptides(M, Biobase::exprs(Exp1_R25_pept[1:1000]))
SumPeptides <- function(matAdj, expr){
    expr <- expr[rownames(matAdj),]
    expr[is.na(expr)] <- 0
    Mp <- t(matAdj) %*% expr
    .temp <- expr
    .temp[!is.na(.temp)] <- 1
    .temp[is.na(.temp)] <- 0
    pep <- t(matAdj) %*% .temp
    
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)
    
    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)
    return(res)
    
}

##' This function computes the intensity of proteins as the mean of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins as the mean of the intensities
##' of their peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' matAdj <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' MeanPeptides(matAdj, Biobase::exprs(Exp1_R25_pept[1:1000]))
MeanPeptides <- function(matAdj, expr){
    expr <- expr[rownames(matAdj),]
    expr[is.na(expr)] <- 0
    Mp <- t(matAdj) %*% expr
    .temp <- expr
    .temp[!is.na(.temp)] <- 1
    .temp[is.na(.temp)] <- 0
    pep <- t(matAdj) %*% .temp
    Mp <- Mp / pep
    
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)
    
    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)
    return(res)
    
}


##' This function computes the intensity of proteins as the sum of the 
##' intensities of their n best peptides.
##' 
##' @title Compute the intensity of proteins as the sum of the 
##' intensities of their n best peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @param n The maximum number of peptides used to aggregate a protein.
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' matAdj <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' TopnPeptides(matAdj, Biobase::exprs(Exp1_R25_pept[1:1000]), 3)
TopnPeptides <-function(matAdj,expr,n){
    
    #Get the indices of the n peptides with the best median
    med <- apply(expr[rownames(matAdj),], 1, median)
    xmed <- as(matAdj * med, "dgCMatrix")
    


   for (c in 1:ncol(matAdj)){
        v <- order(xmed[,c],decreasing=TRUE)[1:n]
       l <- v[which((xmed[,c])[v] != 0)]
       
       if (length(l) > 0){
            diff <- setdiff( which(matAdj[,c] == 1), l)
            if (length(diff)) {matAdj[diff,c] <- 0}
       }
    }

    
#     test <- function(A,B,n){
#         v <- order(B,decreasing=TRUE)[1:n]
#         l <- v[which((B)[v] != 0)]
# 
#         if (length(l) > 0){
#             diff <- setdiff( which(A == 1), l)
#             if (length(diff)) {A[diff] <- 0}
#         }
#     return(A)
#     }
# 
#     time1 <- system.time(t <- mapply(test, 
    #split(matAdj, col(matAdj)), split(xmed, col(xmed)),n))
# 
# print(time0)
# print(time1)
#     
    
    
    res <- SumPeptides(matAdj, expr)

    return(res)
}



##' Method to agregate with a method peptides to proteins on
##' a MSnSet object (peptides)
##' 
##' @title Function agregate peptides to proteins 
##' @param obj.pep An object (peptides) of class \code{\link{MSnbase}}.
##' @param protID The name of proteins ID column 
##' @param method The method used to aggregate the peptides into proteins.
##' Values are "sum", "mean" or "sum on top n" : do the sum / mean of intensity
##' on all peptides belonging to proteins. Default is "sum"
##' @param matAdj An adjacency matrix
##' @param n The number of peptides considered for the aggregation.
##' @return An object of class \code{\link{MSnbase}} with proteins  
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples 
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' mat <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, TRUE)
##' pepAgregate(Exp1_R25_pept[1:1000], protID, "sum overall", mat)
pepAgregate <- function (obj.pep, protID, method="sum overall", 
                         matAdj=NULL, 
                         n=NULL){
    #Check the validity of parameters
    parammethod <- c("sum overall", "mean", "sum on top n") 
    if (sum(is.na(match(method, parammethod) == TRUE)) > 0){return (NULL)}
    if (is.null(matAdj)){warning("Adjacency matrix is missing.")
        return (NULL)}
    
    if (!is.na(match(method, "sum on top n")) && is.null(n)){
        warning("With the top n method, the parameter n must not be NULL.")
        return (NULL)}
    
    condname <- Biobase::pData(obj.pep)$Experiment
    condition <- Biobase::pData(obj.pep)
    expr <- 2^(Biobase::exprs(obj.pep))
    
    if(method == "sum overall"){ res <- SumPeptides(matAdj, expr)}
    else if  (method == "mean"){ res <- MeanPeptides(matAdj, expr)}
    else if (method == "sum on top n"){ res <- TopnPeptides(matAdj, expr, n) }
    
    Mp <- as.matrix(res$matfin)
    Mp[Mp == 0] <- NA
    Mp[is.nan(Mp)] <- NA
    Mp[is.infinite(Mp)] <-NA
    
    
    
    pep <- as.matrix(res$nbpep)
    protId <- res$idprot
    fd <- data.frame(protId, pep)
    
    obj <- MSnSet(exprs = log2(Mp), 
                  fData = fd, 
                  pData = Biobase::pData(obj.pep))
    obj@experimentData@other  <- list(obj@experimentData@other,
                                      typeOfData ="protein")
    
    return(obj)
}

