# peptide-level t-test
##########


## Build design matrix X
X <- as.matrix(BuildAdjacencyMatrix(obj, 'Protein_group_IDs', unique = FALSE))
X.spec <- as.matrix(BuildAdjacencyMatrix(obj, 'Protein_group_IDs', unique = TRUE))

y <- exprs(obj)
y <- y[rownames(X), ]
n1 <- n2 <- 3
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
##' obj <- Exp1_R25_pept
##' protID <- "Protein_group_IDs"
##' keepThat <-  mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' X.spec <- BuildAdjacencyMatrix(obj, protID,  unique = TRUE)
##' qData <- Biobase::exprs(obj)
##' gttest <- groupttest(X.spec, qData)
groupttest <- function(MatAdj, expr){
  nProt <- dim(MatAdj)[2]
  res <- rep(0,nProt)
  
  for(i in 1:nProt){
    print(i)
    index <- names(which(MatAdj[,i]==1))
    if(length(index)== 0){
      res[i] <- 1
    } else{
      peptidestotest <- expr[index,]
      if(length(index)== 1){
        res[i] <- t.test(x=peptidestotest[1:3], y=peptidestotest[-(1:3)], var.equal=TRUE)$p.value
      } else{
        res[i] <- t.test(x=peptidestotest[,1:3], y=peptidestotest[,-(1:3)], var.equal=TRUE)$p.value
      }
    }
  }
  return(res)
}  




##'
##' @title xxxxxx
##' @param qData A matrix of quantitative data, without any missing values.
##' @param sTab xxxx
##' @param X.spec A matrix of quantitative data, without any missing values.
##' @param contrast Indicates if the test consists of the comparison of each 
##' biological condition versus 
##' each of the other ones (contrast=1; 
##' for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
##' or each condition versus all others (contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
##' H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
##' @param type xxxxx
##' @return xxxxx
##' @author Thomas Burger, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept
##' protID <- "Protein_group_IDs"
##' keepThat <-  mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' X.spec <- BuildAdjacencyMatrix(obj, protID,  unique = TRUE)
##' qData <- Biobase::exprs(obj)
##' compute.group.t.tests <- groupttest(qData, sTab, X.spec)
compute.group.t.tests <- function(qData,sTab, X.spec, contrast="OnevsOne", type="Student"){
  
  
  switch(type,
         Student=.type <- TRUE,
         Welch=.type <- FALSE)
  
  
  
  res<-list()
  logFC <- list()
  P_Value <- list()
  
  nbComp <- NULL
  
  sTab.old <- sTab
  Conditions.f <- factor(sTab$Condition, levels=unique(sTab$Condition))
  sTab <- sTab[unlist(lapply(split(sTab, Conditions.f), function(x) {x['Sample.name']})),]
  qData <- qData[,unlist(lapply(split(sTab.old, Conditions.f), function(x) {x['Sample.name']}))]
  Conditions <- sTab$Condition
  
  Cond.Nb<-length(levels(Conditions.f))
  
  
  if(contrast=="OnevsOne"){
    nbComp <- Cond.Nb*(Cond.Nb-1)/2
    
    for(i in 1:(Cond.Nb-1)){
      for (j in (i+1):Cond.Nb){
        
        c1Indice <- which(Conditions==levels(Conditions.f)[i])
        c2Indice <- which(Conditions==levels(Conditions.f)[j])
        
        res.tmp <- apply(qData[,c(c1Indice,c2Indice)], 1, 
                         function(x) {
                           t.test(x~Conditions[c(c1Indice,c2Indice)],  var.equal=.type)
                         })
        p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
        m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
        m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
        m1.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[1])))[1]
        m2.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[2])))[1]
        logFC.tmp <- m1.tmp - m2.tmp
        if (grepl(levels(Conditions.f)[1], m2.name)){logFC.tmp <- -logFC.tmp}
        
        txt <- paste(levels(Conditions.f)[i],"_vs_",levels(Conditions.f)[j], sep="")
        
        logFC[[paste(txt, "logFC", sep="_")]] <- logFC.tmp
        P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
      }
    }
  } ##end Contrast==1
  
  if(contrast=="OnevsAll"){
    nbComp <- Cond.Nb
    
    for(i in 1:nbComp){
      
      c1 <- which(Conditions==levels(Conditions.f)[i])
      
      Cond.t.all<-c(1:length(Conditions))
      Cond.t.all[c1]<-levels(Conditions.f)[i]
      Cond.t.all[-c1]<-"all"
      
      res.tmp <- apply(qData, 1, 
                       function(x) {
                         t.test(x~Cond.t.all, var.equal=.type)
                       })
      
      p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
      m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
      m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
      m1.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[1])))[1]
      m2.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[2])))[1]
      logFC.tmp <- m1.tmp - m2.tmp
      if (grepl(levels(Conditions.f)[i], m2.name)){logFC.tmp <- -logFC.tmp}
      
      txt <- paste(levels(Conditions.f)[i],"_vs_(all-",levels(Conditions.f)[i],")", sep="")
      logFC[[paste(txt, "logFC", sep="_")]] <- logFC.tmp
      P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
    }
  } # End Contrast=2
  
  
  res.l <- list(
    logFC = as.data.frame(logFC),
    P_Value = as.data.frame(P_Value)
  )
  colnames(res.l$logFC) <- names(logFC)
  colnames(res.l$P_Value) <- names(P_Value)
  
  return(res.l) 
  
}