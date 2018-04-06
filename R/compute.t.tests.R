

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
##' @param qData xxxx
##' @param Conditions xxxxx
##' @param Contrast xxxxx
##' @return xxxxxx
##' @author Florence Combes
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findLapalaBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceLapala(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' ttest <- wrapper.limmaCompleteTest(obj, 1)
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

       # nbCompa.1v1 <- Cond.Nb*(Cond.Nb-1)/2
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
                
                FC[[paste(txt, "FC", sep="_")]] <- FC.tmp
                P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
            }
        }
 
    } ##end Contrast==1

    if(Contrast==2){
        
        #nbCompa.1vsAll<-Cond.Nb
        nbComp <- nbCompa.1vsAll
        
        for(i in 1:nbCompa.1vsAll){
            
            c1<-which(Conditions==levels(Conditions.f)[i])
           
            Cond.t.all<-c(1:length(Conditions))
            Cond.t.all[c1]<-levels(Conditions.f)[i]
            Cond.t.all[-c1]<-"all"
            
            p.tmp <- m1.tmp <- m2.tmp <- FC.tmp <- c()
            for(k in 1:dim(qData)[1]){
                t<-t.test(as.numeric(qData[k,]) ~ Cond.t.all, data=qData, var.equal=FALSE)
                p.tmp[k]<-t$p.value
                m1.tmp[k]<-as.numeric(t$estimate[1])
                m2.tmp[k]<-as.numeric(t$estimate[2])
                FC.tmp[k] <- as.numeric(t$estimate[1])-as.numeric(t$estimate[2])
            }
            
            res.tmp<-cbind(p.tmp, m1.tmp, m2.tmp, FC.tmp)
            txt <- paste(unique(Conditions[c1Indice]),"-vs-(all-",unique(Conditions[c1Indice]),")", sep="")
            colnames(res.tmp)<-paste(txt, suffixes, sep="_")
            
            res[[txt]]<-res.tmp
        }
    } # End Contrast=2
    
    
    res.l <- list(
              FC = as.data.frame(FC),
              P_Value = as.data.frame(P_Value)
    )
    
    return(res.l) 
    
}
