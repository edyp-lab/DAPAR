# peptide-level t-test
##########

## Build design matrix X
X <- as.matrix(BuildAdjacencyMatrix(data, 'Protein.group.IDs', unique = FALSE))
X.spec <- as.matrix(BuildAdjacencyMatrix(data, 'Protein.group.IDs', unique = TRUE))

y <- exprs(data)
y <- y[rownames(X), ]
n <- n1+n2 # number of samples in c1 and c2

## Keep track of proteins that are lost in aggregation step
unsup.prot <- !(colnames(X) %in% colnames(X.spec))
q <- nrow(X) # Number of peptides
p <- ncol(X) # Number of proteins


# Two ways to compute it (function groupttest below)
peptide.spec.based.tmp <- apply(X.spec, 2, FUN=function(mask) t.test(x=as.vector(y[mask == 1, 1:n1]), y=as.vector(y[mask == 1, -(1:n1)]), var.equal=TRUE)$p.value)
peptide.spec.based.tmp <- groupttest(X.spec,y)

# then:
peptide.spec.based.pv <- rep(1, ncol(X))
peptide.spec.based.pv[!unsup.prot] <- peptide.spec.based.tmp





#######################
### attention, à rendre générique: pour l'instant ne marche qu'avec n1=3 (car codé en dur)
##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param qData A matrix of quantitative data, without any missing values.
##' @param X A matrix of quantitative data, without any missing values.
##' @return xxxxx
##' @author Thomas Burger, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' keepThat <- mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' sTab <- Biobase::pData(obj)
##' qData <- Biobase::exprs(obj)
##' gttest <- groupttest(X, qData)
groupttest <- function(MatAdj, expr){
  nProt <- dim(MatAdj)[2]
  res <- rep(0,nProt)
  
  for(i in 1:nProt){
    #print(i)
    index <- names(which(MatAdj[,i]==1))
    if(length(index)== 0){
      res[i] <- 1
    } else{
      peptidestotest <- expr[which(rownames(expr)%in% index),]
      if(length(index)== 1){
        res[i] <- t.test(x=peptidestotest[1:3], y=peptidestotest[-(1:3)], var.equal=TRUE)$p.value
      } else{
        res[i] <- t.test(x=peptidestotest[,1:3], y=peptidestotest[,-(1:3)], var.equal=TRUE)$p.value
      }
    }
  }
  return(res)
}  