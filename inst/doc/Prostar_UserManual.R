### R code from vignette source 'Prostar_UserManual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: installProstaRBiocond (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite("Prostar")


###################################################
### code chunk number 3: runProstarStandalone (eval = FALSE)
###################################################
## library(Prostar)
## Prostar()


###################################################
### code chunk number 4: installDAPARBiocond (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite("DAPAR")


###################################################
### code chunk number 5: sessionLog (eval = FALSE)
###################################################
## getProcessingInfo(obj)


###################################################
### code chunk number 6: diffAnalysis (eval = FALSE)
###################################################
## res <- diffAnaLimma(imputed_dataset, condition1, condition2)
## obj <- diffAnaSave(imputed_dataset, res, "limma", condition1, condition2)


###################################################
### code chunk number 7: installDAPARBiocondJAVAHOME (eval = FALSE)
###################################################
## Sys.setenv("JAVA_HOME", "")


###################################################
### code chunk number 8: sessioninfo
###################################################
toLatex(sessionInfo())


