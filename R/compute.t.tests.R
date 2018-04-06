

##' This function is a wrapper xxxxx
##'
##' @title xxxxx
##' @param obj An object of class \code{MSnSet} with no missing values
##' @param design An integer that reflects the type of comparisons. Available values are 1 (One vs One) or (One vs All)
##' @return xxxxxxx
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' ttest <- wrapper.t_test_Complete(obj, 1)
wrapper.t_test_Complete <- function(obj, design){
    
    qData <- Biobase::exprs(obj)
    conds <- pData(obj)[,"Label"]
    switch(design,
           OnevsOne=contrast <- 1,
           OnevsAll=contrast <- 2)
    
    ttest <- compute.t.tests(qData,conds, contrast)
    
    return (ttest)
}




##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param qData A matrix of quantitative data, without any missing values.
##' @param Conditions A vector of factor which indicates the name of the 
##' biological condition for each replicate. 
##' @param Contrast Indicates if the test consists of the comparison of each 
##' biological condition versus 
##' each of the other ones (Contrast=1; 
##' for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
##' or each condition versus all others (Contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
##'  H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
##' @return A list of two items : FC and P_Value; both are dataframe. The first one contains
##' the logFC values of all the comparisons (one column for one comparison), the second one contains
##' the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes
##' are identical and correspond to the description of the comparison. 
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' ttest <- compute.t.tests(obj, 1)
compute.t.tests <- function(qData,Conditions, Contrast=1){

res<-list()
FC <- list()
P_Value <- list()

nbComp <- NULL
suffixes <- c("pval", "m1", "m2", "logFC")

Conditions.f <- factor(Conditions)
#Cond<-levels(Conditions.f)
Cond.Nb<-length(levels(Conditions.f))


    if(Contrast==1){
        nbComp <- Cond.Nb*(Cond.Nb-1)/2

        for(i in 1:(Cond.Nb-1)){
            for (j in (i+1):Cond.Nb){
    
                c1Indice <- which(Conditions==levels(Conditions.f)[i])
                c2Indice <- which(Conditions==levels(Conditions.f)[j])
    
                res.tmp <- apply(qData[,c(c1Indice,c2Indice)], 1, 
                                 function(x) {
                   t.test(x~Conditions[c(c1Indice,c2Indice)], data=qData, var.equal=F)
                })
                p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
                m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
                m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
                FC.tmp <- m1.tmp - m2.tmp
                
                txt <- paste(unique(Conditions[c1Indice]),"-vs-",unique(Conditions[c2Indice]), sep="")
                
                FC[[paste(txt, "logFC", sep="_")]] <- FC.tmp
                P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
            }
        }
    } ##end Contrast==1

    if(Contrast==2){
        nbComp <- Cond.Nb
        
        for(i in 1:nbComp){
            
            c1<-which(Conditions==levels(Conditions.f)[i])
           
            Cond.t.all<-c(1:length(Conditions))
            Cond.t.all[c1]<-levels(Conditions.f)[i]
            Cond.t.all[-c1]<-"all"
            
            res.tmp <- apply(qData, 1, 
                             function(x) {
                                 t.test(x~Cond.t.all, data=qData, var.equal=FALSE)
                             })
            
            p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
            m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
            m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
            FC.tmp <- m1.tmp - m2.tmp
            
            txt <- paste(unique(Conditions[c1]),"-vs-(all-",unique(Conditions[c1]),")", sep="")
            
            FC[[paste(txt, "logFC", sep="_")]] <- FC.tmp
            P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
        }
    } # End Contrast=2
    
    
    res.l <- list(
              FC = as.data.frame(FC),
              P_Value = as.data.frame(P_Value)
    )
    
    return(res.l) 
    
}
