


#dat<-readRDS("Wilson.2.4.6.norm.MSnset")


#####################DESIGN#######################

Conditions<-c("KO12","KO12","KO36","KO36","KO36","KO36","KO44","KO44","KO44","KO44","KO44","KO44")
Conditions.f <- factor(Conditions)
Cond<-levels(Conditions.f)
Cond.Nb<-length(levels(Conditions.f))



wrapper.t_test_Complete <- function(obj, design){
    
    qData <- Biobase::exprs(obj)
    conds <- pData(obj)[,"Label"]
    switch(design,
           OnevsOne=contrast <- 1,
           OnevsAll=contrast <- 2)
    
    ttest <- t_test_Complete(qData,conds, contrast)
    
    return (ttest)
}





compute.t.tests <- function(qData,Conditions, Contrast=1){

    res<-list()
    if(Contrast==1){

        nbCompa.1v1<-Cond.Nb*(Cond.Nb-1)/2
        c <- 1
        
        for(i in 1:(Cond.Nb-1)){
            for (j in (i+1):Cond.Nb){
    
                c1Indice <- which(Conditions==levels(Conditions.f)[i])
                c2Indice <- which(Conditions==levels(Conditions.f)[j])
    
                p.tmp <- m1.tmp <- m2.tmp <- c()
                for(k in 1:dim(qData)[1]){
                    t<-t.test(as.numeric(qData[k,c(c1Indice,c2Indice)]) ~ Conditions[c(c1Indice,c2Indice)], data=qData, var.equal=F)
                    p.tmp[k]<-t$p.value
                    m1.tmp[k]<-as.numeric(t$estimate[1])
                    m2.tmp[k]<-as.numeric(t$estimate[2])
                }
        
                res.tmp <- cbind(p.tmp, m1.tmp, m2.tmp)
                colnames(res.tmp)<-c(paste("p", unique(Conditions[c1Indice]),"VS",unique(Conditions[c2Indice]), sep="."),
                                 paste("m1", unique(Conditions[c1Indice]),"VS",unique(Conditions[c2Indice]), sep="."),
                                 paste("m2", unique(Conditions[c1Indice]),"VS",unique(Conditions[c2Indice]), sep=".")
                )
                res[[c]] <- res.tmp
                c=c+1
            }
        }
 
    } ##en Contrast==1

    if(Contrast==2){
        
        nbCompa.1vsAll<-Cond.Nb

        
        for(i in 1:nbCompa.1vsAll){
            
            c1<-which(Conditions==levels(Conditions.f)[i])
           
            Cond.t.all<-c(1:length(Conditions))
            Cond.t.all[c1]<-levels(Conditions.f)[i]
            Cond.t.all[-c1]<-"all"
            
            p.tmp<-m1.tmp<-m2.tmp<-c()
            for(k in 1:dim(qData)[1]){
                t<-t.test(as.numeric(qData[k,]) ~ Cond.t.all, data=qData, var.equal=F)
                p.tmp[k]<-t$p.value
                m1.tmp[k]<-as.numeric(t$estimate[1])
                m2.tmp[k]<-as.numeric(t$estimate[2])
            }
            
            res.tmp<-cbind(p.tmp, m1.tmp, m2.tmp)
            colnames(res.tmp)<-c(paste("p.", unique(Conditions[c1]),"VS(all-",unique(Conditions[c1]), ")",sep=""),
                                 paste("m1.", unique(Conditions[c1]),"VS(all-",unique(Conditions[c1]), ")", sep=""),
                                 paste("m2.", unique(Conditions[c1]),"VS(all-",unique(Conditions[c1]), ")", sep="")
            )
            
            res[[i]]<-res.tmp
        }
    }
    
    
    
    res.df <- as.data.frame(res)
    indPVal <- seq(from=1,by=3, length.out=3)
    indm1 <- seq(from=2,by=3, length.out=3)
    indm2 <- seq(from=3,by=3, length.out=3)
    
    res.l <- list(
        FC = res.df[,indm1]/res.df[,indm2],
        P_Value = res.df[,indPVal]
    )
    
    colnames(res.l$FC) <- paste(cn, "logFC",sep="_")
    colnames(res.l$P_Value) <- paste(cn, "pval",sep="_")
    
    FC <- P_Value <- data.fra
    for (i in 1:length(res)){
        #FC <- cbind(FC, res[[i]][,2]/res[[i]][,3])
        P_Value <- cbind(P_Value, as.data.frame(res[[i]])[,1])
    }
    return(res) 
    
}