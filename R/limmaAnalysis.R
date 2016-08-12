
##############################################################
#Fonction r?alisant un test de contrastes entre conditions ? l'aide du package LIMMA.
#
#En entr?e:
#qData: tableau de donn?es sans valeurs manquantes avec chaque r?plicat en colonne
#Conditions: indique le num?ro/la lettre de la condition biologique auquel appartient chaque r?plicat

#
#Contrast: indique si l'on souhaite tester chaque condition biologique contre chacune (Contrast=1; par exemple H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#ou chaque condition contre toutes les autres (Contrast=2; par exemple H0:"C1=(C2+C3)/2" vs H1:"C1!=(C2+C3)/2", etc. si on a 3 conditions ).
#
#En sortie :
#Objet fit renvoy? par la fonction eBayes de LIMMA.
##QGG, Aout 2015
##############################################################





##' This function is a limmaCompleteTest
##' 
##' @title Computes a hierarchical differential analysis
##' @param qData A matrix of quantitative data, without any missing values.
##' @param Conditions A vector of factor which indicates the name of the biological condition for each replicate. 
##' @param RepBio A vector of factor which indicates the number of the bio rep for each replicate. 
##' @param RepTech A vector of factor which indicates the number of the tech rep for each replicate.
##' @param Contrast Indicates if the test concists of the comparison of each biological condition versus 
##' each of the other ones (Contrast=1; for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
##' or each condition versus all others (Contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
##'  H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
##' @return fdsfdgfdg
##' @author Quentin Giai-Gianetto
##' @examples
##' data(UPSpep25)
##' obj <- wrapper.mvImputation(UPSpep25, "QRILC")
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(obj)
##' RepBio <- RepTech <- factor(1:6)
##' conds <- factor(c(rep(condition1, 3), (rep(condition2, 3))))
##' limma <- limmaCompleteTest(qData,conds,RepBio, RepTech)
limmaCompleteTest <- function(qData,Conditions, RepBio, RepTech, Contrast=1){
    #####################################################################
    #Fonctions permettant de cr?er des design et contrastes
    #######################
    
    
    
    #Compute of the hierarchical design matrix : 2-level case : Bio-Tech or Bio-Analytical or Tech-Analytical .
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
    #Compute of the hierarchical design matrix : 3-levels case : Bio, Tech and Analytical.
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
    
    
    
    
    
    #######################
    #Aggregation of columns from the same condition to build the contrast matrix.
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
                if (substr(name.col[j], 1, nc)==col.name.begin){col.select=c(col.select,j)}
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
    
    
    #######################
    #Compute the contrast matrix: case a condition vs one condition.
    make.contraste.1.1=function(design,Condition){
        nb.cond=length(levels(Condition))
        r=aggreg.column.design(design,Condition)
        label.agg=r[[1]]
        nb.agg=r[[2]]
        k=1
        contra=rep(0,sum(1:(nb.cond-1)))
        for (i in 1:(nb.cond-1)){
            for (j in (i+1):nb.cond){
                contra[k]=c(paste("(",label.agg[i],")/",nb.agg[i],"-(",label.agg[j],")/",nb.agg[j]))
                k=k+1
            }
        }
        return(contra)
    }
    
    
    #######################
    #Compute the contrast matrix: case a condition vs all other conditions.
    make.contraste.2.1=function(design,Condition){
        nb.cond=length(levels(Condition))
        r=aggreg.column.design(design,Condition)
        label.agg=r[[1]]
        nb.agg=r[[2]]
        k=1
        contra=rep(0,sum(1:(nb.cond-1)))
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
    
    
    
    #####################################################################
    #Main function
    
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
                    #CHeck if the number of factors in Conditions is less than the one in Bio.Rep, itself must
                    #be less than tech.Rep
                    if (length(levels(Conditions))<=length(levels(RepBio))){
                        if (length(levels(RepBio))<=length(levels(RepTech))){
                            #Get the number of hierarchical levels
                            if (length(levels(RepBio))==lt){niveau=1;}
                            else{ if (length(levels(RepTech))==lt){niveau=2;}
                                else{niveau=3;}
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
                                if (Contrast==1){contra=make.contraste.1.1(design,Condition=Conditions);}
                                if (Contrast==2){contra=make.contraste.2.1(design,Condition=Conditions);}
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
                                if (Contrast==1){contra=make.contraste.1.1(design,Condition=Conditions);}
                                if (Contrast==2){contra=make.contraste.2.1(design,Condition=Conditions);}
                                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                                fit <- eBayes(contrasts.fit(lmFit(qData, design), cmtx))
                                #return(fit)
                            }
                        }
                        else{cat("Problem: the number of factors in Bio.Rep has to be less or equal to the number of factor in tech.Rep.!\n")}
                        
                    }
                    else{cat("Problem: the number of factors in Conditions has to be less or equal to the number of factor in bio.Rep.!\n")}
                    
                }
                else{cat("Problem: the length of the vector tech.Rep must be equal to the length of the vector bio.Rep!\n")}
                
            }
            else{cat("Problem: the length of the vector bio.Rep must be equal to the length of the vector Conditions!\n")}
            
        }
        else{cat("Problem: the length of the vector Conditions must be equal to the length of the input matrix!\n")}
        
    }
    else{cat("Problem: the factor-vector Conditions must contain several replicates in each condition  !\n")}
    
    res <- topTable(fit,number=Inf)
   
        return(res)
}


