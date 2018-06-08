
##############################################################
#Fonction realisant un test de contrastes entre conditions a l'aide du 
#package LIMMA.
#
#En entree:
#qData: tableau de donn?es sans valeurs manquantes avec chaque replicat en 
#colonne
#Conditions: indique le numero/la lettre de la condition biologique auquel 
#appartient chaque replicat

#
#Contrast: indique si l'on souhaite tester chaque condition biologique 
#contre chacune (Contrast=1; par exemple H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#ou chaque condition contre toutes les autres (Contrast=2; par exemple 
#H0:"C1=(C2+C3)/2" vs H1:"C1!=(C2+C3)/2", etc. si on a 3 conditions ).
#
#En sortie :
#Objet fit renvoye par la fonction eBayes de LIMMA.
##QGG, Aout 2015
##############################################################


##' This function is a wrapper of \code{limmaCompleteTest} for objects of class \code{MSnSet}
##'
##' @title Computation of the hierarchical design matrix : 2-level case : Bio-Tech or 
##' Bio-Analytical or Tech-Analytical.
##' @param obj An object of class \code{MSnSet} with no missing values
##' @param design An integer that reflects the type of comparisons. Available values are 1 (One vs One) or 2 (One vs All)
##' @return xxxxxxxxxx
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' limma <- wrapper.limmaCompleteTest(obj, 1)
wrapper.limmaCompleteTest <- function(obj, design){
    
    qData <- Biobase::exprs(obj)
    RepBio <- RepTech <- factor(1:nrow(pData(obj)))
    #RepBio <- RepTech <- factor(1:nrow(pData(obj)))
    conds <- factor(pData(obj)[,"Label"])
    switch(design,
           OnevsOne=contrast <- 1,
           OnevsAll=contrast <- 2)
    
    limma <- limmaCompleteTest(qData,conds,RepBio, RepTech, contrast)
    
    return (limma)
}







##' This function is a limmaCompleteTest
##' 
##' @title Computes a hierarchical differential analysis
##' @param qData A matrix of quantitative data, without any missing values.
##' @param Conditions A vector of factor which indicates the name of the 
##' biological condition for each replicate. 
##' @param RepBio A vector of factor which indicates the number of the bio rep 
##' for each replicate. 
##' @param RepTech A vector of factor which indicates the number of the tech 
##' rep for each replicate.
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
##' @author Quentin Giai-Gianetto
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' lapala <- findMECBlock(obj)
##' obj <- wrapper.impute.detQuant(obj)
##' obj <- reIntroduceMEC(obj, lapala)
##' obj <- wrapper.impute.detQuant(obj)
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(obj)
##' RepBio <- RepTech <- factor(1:6)
##' conds <- factor(Biobase::pData(obj)[,"Label"])
##' limma <- limmaCompleteTest(qData,conds,RepBio, RepTech)
limmaCompleteTest <- function(qData,Conditions, RepBio, RepTech, Contrast=1){
    
     make.design.2=function(qData, Condition, RepBio){
        #Renome the levels of factor
        levels(Condition)=c(1:length(levels(Condition)))
        levels(RepBio)=c(1:length(levels(RepBio)))
        
        #Initial design matrix
        design=model.matrix(qData[1,]~0+Condition:RepBio)
        
        #Remove empty columns in the design matrix
        design=design[,(apply(design,2,sum)>0)]
        #Remove identical columns in the design matrix
        coldel=-1
        for (i in 1:(length(design[1,])-1)){
            d2=as.matrix(design[,(i+1):length(design[1,])]);
            for (j in 1:length(d2[1,])){
                d2[,j]=d2[,j]-design[,i];
            }
            e=as.matrix(rnorm(length(design[,1]),10,1));
            sd2=t(e)%*%d2
            liste=which(sd2==0)
            coldel=c(coldel,liste+i)
        }
        design=design[,(1:length(design[1,]))!=coldel]
        colnames(design)=make.names(colnames(design))
        return(design)
    }
    
    
    
    #######################
    make.design.3=function(qData,Condition,RepBio,RepTech){
        #Rename the levels of factor
        levels(Condition)=c(1:length(levels(Condition)))
        levels(RepBio)=c(1:length(levels(RepBio)))
        levels(RepTech)=c(1:length(levels(RepTech)))
        
        
        #Initial design matrix
        design=model.matrix(qData[1,]~0+Condition:RepBio:RepTech)
        
        #Remove empty columns in the design matrix
        design=design[,(apply(design,2,sum)>0)]
        
        #Remove identical columns in the design matrix
        coldel=-1
        for (i in 1:(length(design[1,])-1)){
            d2=as.matrix(design[,(i+1):length(design[1,])]);
            for (j in 1:length(d2[1,])){
                d2[,j]=d2[,j]-design[,i];
            }
            e=as.matrix(rnorm(length(design[,1]),10,1));
            sd2=t(e)%*%d2
            liste=which(sd2==0)
            coldel=c(coldel,liste+i)
        }
        design=design[,(1:length(design[1,]))!=coldel]
        colnames(design)=make.names(colnames(design))
        return(design)
    }
    
    
    
    aggreg.column.design=function(design,Condition){
        nb.cond=length(levels(Condition))
        name.col=colnames(design)
        name.cond=NULL
        nb.col=NULL
        for (i in 1:nb.cond){
            col.select=NULL
            col.name.begin=paste("Condition",i, sep = "")
            nc=nchar(col.name.begin)
            for (j in 1:length(design[1,])){
                if (substr(name.col[j], 1, nc)==col.name.begin){
                    col.select=c(col.select,j)
                }
            }
            name.aggreg=NULL
            for (j in 1:length(col.select)){
                name.aggreg=paste(name.aggreg,name.col[col.select[j]],sep="+")
            }
            name.aggreg=substr(name.aggreg, 2, nchar(name.aggreg))
            name.cond=c(name.cond,name.aggreg)
            nb.col=c(nb.col,length(col.select))
        }
        return(list(name.cond,nb.col))
    }
    
    
    make.contraste.1.1=function(design,Condition){
        nb.cond=length(levels(Condition))
        r=aggreg.column.design(design,Condition)
        label.agg=r[[1]]
        nb.agg=r[[2]]
        k=1
        contra=rep(0,sum(1:(nb.cond-1)))
        for (i in 1:(nb.cond-1)){
            for (j in (i+1):nb.cond){
                contra[k]=c(paste("(",label.agg[i],")/",
                                  nb.agg[i],"-(",label.agg[j],")/",
                                  nb.agg[j]))
                k=k+1
            }
        }
        return(contra)
    }
    
    
    
    make.contraste.2.1=function(design,Condition){
        nb.cond=length(levels(Condition))
        r=aggreg.column.design(design,Condition)
        label.agg=r[[1]]
        nb.agg=r[[2]]
        k=1
        #contra=rep(0,sum(1:(nb.cond-1)))
        contra=rep(0,nb.cond)
        for (i in 1:(nb.cond)){
            contra[k]=c(paste("(",label.agg[i],")/",nb.agg[i]))
            nb=sum(nb.agg[(1:nb.cond)[(1:nb.cond)!=i]])
            for (j in (1:nb.cond)[(1:nb.cond)!=i]){
                contra[k]=c(paste(contra[k],"-(",label.agg[j],")/",nb))
            }
            k=k+1
        }
        return(contra)
    }
    
    
    
    
    ### Begin of the main function
    
    Conditions <- factor(Conditions)
    RepBio <- factor(RepBio)
    RepTech <- factor(RepTech)
    
    lt=length(qData[1,])
    lc=length(Conditions)
    lb=length(RepBio)
    lte=length(RepTech)
    
    #Check the correct length of vectors
    if (length(levels(Conditions))!=lc){
        if(lt==lc){
            if (lb==lc){
                if (lte==lc){  
                    #CHeck if the number of factors in Conditions is less than 
                    #the one in Bio.Rep, itself must
                    #be less than tech.Rep
                    if (length(levels(Conditions))<=length(levels(RepBio))){
                        if (length(levels(RepBio))<=length(levels(RepTech))){
                            #Get the number of hierarchical levels
                            if (length(levels(RepBio))==lt){niveau=1;}
                            else{ 
                              if (length(levels(RepTech))==lt){niveau=2;}
                                else{
                                  if (length(levels(RepTech))==length(levels(RepBio))){niveau=2;}
                                 else { niveau=3;}
                                }
                            }
                            #Non hierarchical case (only one type of replicate)
                            if (niveau==1){
                                nb_cond=length(levels(Conditions))
                                #CGet the number of replicates per condition
                                nb_Rep=rep(0,nb_cond)
                                for (i in 1:nb_cond){
                                    nb_Rep[i]=sum((Conditions==levels(Conditions)[i]))
                                }
                                #####################
                                #Compute the design matrix
                                design=matrix(0,lt,nb_cond)
                                n0=1
                                coln=NULL
                                for (j in 1:nb_cond){
                                    coln=c(coln,paste("Condition",j,collapse=NULL,sep=""))
                                    design[(n0:(n0+nb_Rep[j]-1)),j]=rep(1,length((n0:(n0+nb_Rep[j]-1))))
                                    n0=n0+nb_Rep[j]
                                }
                                colnames(design)=coln
                                
                                #####################
                                #Compute the contrast matrices
                                if (Contrast==1){contra=make.contraste.1.1(design,Condition=Conditions);}
                                if (Contrast==2){contra=make.contraste.2.1(design,Condition=Conditions);}
                                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                                fit=eBayes(contrasts.fit(lmFit(qData, design), cmtx))
                                # return(fit)  
                            }
                            if (niveau==2){
                                #####################
                                #Compute the design matrix
                                design=make.design.2(qData,Condition=Conditions,RepBio=RepBio)
                                #####################
                                #Compute the contrast matrices
                                if (Contrast==1){
                                    contra=make.contraste.1.1(design,Condition=Conditions);}
                                if (Contrast==2){
                                    contra=make.contraste.2.1(design,Condition=Conditions);}
                                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                                fit <- eBayes(contrasts.fit(lmFit(qData, design), cmtx))
                                #return(fit)
                            }
                            if (niveau==3){
                                #####################
                                #Compute the design matrix
                                design=make.design.3(qData,Condition=Conditions,RepBio=RepBio,RepTech=RepTech)
                                #####################
                                #Compute the contrast matrices
                                if (Contrast==1){
                                    contra=make.contraste.1.1(design,Condition=Conditions);}
                                if (Contrast==2){
                                    contra=make.contraste.2.1(design,Condition=Conditions);}
                                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                                fit <- eBayes(contrasts.fit(lmFit(qData, design), cmtx))
                                #return(fit)
                            }
                        }else{cat("Problem: the number of factors in Bio.Rep has to be less or equal to the number of factor in tech.Rep.!\n")}
                        
                    }else{cat("Problem: the number of factors in Conditions has to be less or equal to the number of factor in bio.Rep.!\n")}
                    
                }else{cat("Problem: the length of the vector tech.Rep must be equal to the length of the vector bio.Rep!\n")}
                
            }else{cat("Problem: the length of the vector bio.Rep must be equal to the length of the vector Conditions!\n")}
            
        }else{cat("Problem: the length of the vector Conditions must be equal to the length of the input matrix!\n")}
        
    }else{cat("Problem: the factor-vector Conditions must contain several replicates in each condition  !\n")}
    
    
    
    ##colnames (string treatment in R ...)
    #how many comparisons have been done (and thus how many columns of pval)
    Compa.Nb <- dim(fit$p.value)[2]
    
    #pas topTable, ou adapt?. 
    #ramener ttes les pvalues et colnames corrects. 
    
    res.tmp <- topTable(fit,number=Inf, sort.by="none")
    res <- cbind(res.tmp[,1:Compa.Nb], fit$p.value)
    #names(res) <- gsub(".", "_", names(res), fixed=TRUE)
    
    
    #empty colnames vector
    cn<-c()
    for (i in 1:Compa.Nb){
        
        #not the same syntax to pars if Contast=1 or Contrast=2
        if(Contrast==1){
            compa <- stringr::str_match_all(colnames(fit$p.value)[i],"[[:space:]]Condition([[:digit:]]+)")[[1]]
            cn[i] <- paste(levels(Conditions)[as.numeric(compa[1,2])], "_vs_",levels(Conditions)[as.numeric(compa[2,2])], sep="")
        }
        if(Contrast==2){
            #hierarchic only
            #compa<-str_match_all(colnames(fit$p.value)[i], "[[:space:]]Condition([[:digit:]]+)[[:space:]]")[[1]]
            #cn[i]<-paste(levels(Conditions)[as.numeric(compa[1,2])], "vs(all-",levels(Conditions)[as.numeric(compa[1,2])], ")", sep="")
            
            #hier and non hier
            compa<-str_match_all(colnames(fit$p.value)[i], "[[:space:]]Condition([[:digit:]]+)")[[1]]
            cn[i]<-paste(levels(Conditions)[as.numeric(compa[1,2])], "_vs_(all-",levels(Conditions)[as.numeric(compa[1,2])], ")", sep="")
        }
    }
    
     
    res.l <- list(
      FC = as.data.frame(res[,1:Compa.Nb]),
      P_Value = as.data.frame(res[,-(1:Compa.Nb)] )
    )
    
    colnames(res.l$FC) <- paste(cn, "FC",sep="_")
    colnames(res.l$P_Value) <- paste(cn, "pval",sep="_")
    ## end colnames
    
    return(res.l)
}







