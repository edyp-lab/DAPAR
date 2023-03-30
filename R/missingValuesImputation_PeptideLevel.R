


#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @description 
#' This method is a wrapper to the function \code{impute.mle()} of the 
#' package \code{imp4p} adapted to an object of class \code{MSnSet}.
#' It does not impute MEC missing values.
#'
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @return The \code{Biobase::exprs(obj)} matrix with imputed values instead 
#' of missing
#' values.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10), ]
#' level <- 'peptide'
#' metacell.mask <- match.metacell(GetMetacell(obj), c("Missing POV", "Missing MEC"), level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj.imp.na <- wrapper.impute.mle(obj)
#'
#' @export
#'
#'
wrapper.impute.mle <- function(obj) {
    
    pkgs.require('imp4p')
    
    if (missing(obj))
        stop("'obj' is required.")

    tmp_cond <- Biobase::pData(obj)$Condition
    cond <- factor(tmp_cond, levels = unique(tmp_cond))
    res <- imp4p::impute.mle(Biobase::exprs(obj), conditions = cond)

    Biobase::exprs(obj) <- res
    obj <- UpdateMetacellAfterImputation(obj)

    return(obj)
}




#' @title Missing values imputation using the LSimpute algorithm.
#' 
#' @description 
#' This method is a wrapper to the function \code{impute.mi()} of the package
#' \code{imp4p} adapted to an object of class \code{MSnSet}.
#'
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param nb.iter Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param nknn Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param selec Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param siz Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param weight Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param ind.comp Same as the function \code{mi.mix} in the package
#' \code{imp4p}
#'
#' @param progress.bar Same as the function \code{mi.mix} in the package
#' \code{imp4p}
#'
#' @param x.step.mod Same as the function \code{estim.mix} in the package
#' \code{imp4p}
#'
#' @param x.step.pi Same as the function \code{estim.mix} in the package
#' \code{imp4p}
#'
#' @param nb.rei Same as the function \code{estim.mix} in the package
#' \code{imp4p}
#'
#' @param method Same as the function \code{estim.mix} in the package
#' \code{imp4p}
#'
#' @param gridsize Same as the function \code{estim.mix} in the package
#' \code{imp4p}
#'
#' @param q Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param q.min Same as the function \code{impute.pa} in the package
#' \code{imp4p}
#'
#' @param q.norm Same as the function \code{impute.pa} in the package
#' \code{imp4p}
#'
#' @param eps Same as the function \code{impute.pa} in the package \code{imp4p}
#'
#' @param methodi Same as the function \code{mi.mix} in the package \code{imp4p}
#'
#' @param lapala xxxxxxxxxxx
#'
#' @param distribution The type of distribution used. Values are \code{unif}
#' (default) or \code{beta}.
#'
#' @return The \code{Biobase::exprs(obj)} matrix with imputed values 
#' instead of missing values.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(100)]
#' level <- 'peptide'
#' metacell.mask <- match.metacell(GetMetacell(obj), c("Missing POV", "Missing MEC"), level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj.imp.na <- wrapper.dapar.impute.mi(obj, nb.iter = 1, lapala = TRUE)
#' obj.imp.pov <- wrapper.dapar.impute.mi(obj, nb.iter = 1, lapala = FALSE)
#'
#' @export
#'
#'
wrapper.dapar.impute.mi <- function(obj,
                                    nb.iter = 3,
                                    nknn = 15,
                                    selec = 600,
                                    siz = 500, 
                                    weight = 1,
                                    ind.comp = 1,
                                    progress.bar = FALSE,
                                    x.step.mod = 300,
                                    x.step.pi = 300,
                                    nb.rei = 100,
                                    method = 4,
                                    gridsize = 300,
                                    q = 0.95,
                                    q.min = 0,
                                    q.norm = 3,
                                    eps = 0,
                                    methodi = "slsa",
                                    lapala = TRUE,
                                    distribution = "unif") {

    pkgs.require('imp4p')

    if (missing(obj))
        stop("'obj' is required.")



    ## order exp and pData table before using imp4p functions
    tmp <- Biobase::pData(obj)$Condition
    conds <- factor(tmp, levels = unique(tmp))
    sample.names.old <- Biobase::pData(obj)$Sample.name
    sTab <- Biobase::pData(obj)
    qData <- Biobase::exprs(obj)
    new.order <- unlist(lapply(split(sTab, conds), function(x) {
        x["Sample.name"]
    }))
    qData <- Biobase::exprs(obj)[, new.order]
    sTab <- Biobase::pData(obj)[new.order, ]

    conditions <- as.factor(sTab$Condition)
    repbio <- sTab$Bio.Rep
    reptech <- sTab$Tech.Rep


    dat.slsa <- imp4p::impute.rand(tab = qData, conditions = conditions)


    res <- imp4p::estim.mix(tab = qData,
                            tab.imp = dat.slsa,
                            conditions = conditions
                            )


    born <- imp4p::estim.bound(tab = qData, conditions = conditions)
    proba <- imp4p::prob.mcar.tab(born$tab.upper, res)


    data.mi <- imp4p::mi.mix(tab = qData,
                             tab.imp = dat.slsa,
                             prob.MCAR = proba,
                             conditions = conditions,
                             repbio = repbio,
                             reptech = reptech,
                             nb.iter = nb.iter
                             )

    if (lapala == TRUE) {
        data.final <- impute.pa2(tab = data.mi,
                                 conditions = conditions,
                                 q.min = q.min,
                                 q.norm = q.norm,
                                 eps = eps,
                                 distribution = distribution
                                 )
        } else {
          data.final <- data.mi
          }


    # restore previous order
    colnames(data.final) <- new.order
    data.final <- data.final[, sample.names.old]

    Biobase::exprs(obj) <- data.final

    msg <- paste("Missing values imputation using imp4p")
    obj@processingData@processing <- c(obj@processingData@processing, msg)

    obj@experimentData@other$imputation.method <- "imp4p"
    # na.type <- "Missing POV"
    # if (isTRUE(lapala)) {
    #     na.type <- "Missing"
    # }

    obj <- UpdateMetacellAfterImputation(obj)

    return(obj)
}



#' @title Generator of simulated values
#'
#' @param n An integer which is the number of simulation (same as in rbeta)
#'
#' @param min An integer that corresponds to the lower bound of the interval
#'
#' @param max An integer that corresponds to the upper bound of the interval
#'
#' @param param1 An integer that is the first parameter of rbeta function.
#'
#' @param param2 An integer that is second parameter of rbeta function.
#'
#' @return A vector of n simulated values
#'
#' @author Thomas Burger
#'
#' @examples
#' translatedRandomBeta(1000, 5, 10, 1, 1)
#'
#' @export
#'
#'
#'
translatedRandomBeta <- function(n, min, max, param1 = 3, param2 = 1) {
    pkgs.require('stats')

    scale <- max - min
    simu <- stats::rbeta(n, param1, param2)
    res <- (simu * scale) + min
    return(res)
}




#'
#' @title Missing values imputation from a \code{MSnSet} object
#' 
#' @description 
#' This method is a wrapper to the function \code{impute.pa2()} adapted to
#' objects of class \code{MSnSet}.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param q.min A quantile value of the observed values allowing defining the
#' maximal value which can be generated. This maximal value is defined by the
#' quantile q.min of the observed values distribution minus eps.
#' Default is 0 (the maximal value is the minimum of observed values minus eps).
#'
#' @param q.norm A quantile value of a normal distribution allowing defining
#' the minimal value which can be generated. Default is 3 (the minimal value
#' is the maximal value minus qn*median(sd(observed values)) where sd is the
#' standard deviation of a row in a condition).
#'
#' @param eps A value allowing defining the maximal value which can be
#' generated. This maximal value is defined by the quantile q.min of the
#' observed values distribution minus eps. Default is 0.
#'
#' @param distribution The type of distribution used. Values are \code{unif}
#' (default) or \code{beta}.
#'
#' @return The object \code{obj} which has been imputed
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' obj.imp.pa2 <- wrapper.impute.pa2(Exp1_R25_pept[seq_len(100)], 
#' distribution = "beta")
#'
#' @export
#'
#'
wrapper.impute.pa2 <- function(obj,
                               q.min = 0,
                               q.norm = 3,
                               eps = 0,
                               distribution = "unif") {

    if (missing(obj))
        stop("'obj' is required.")


    ## order exp and pData table before using imp4p functions
    tmp <- Biobase::pData(obj)$Condition
    conds <- factor(tmp, levels = unique(tmp))
    sample.names.old <- Biobase::pData(obj)$Sample.name
    sTab <- Biobase::pData(obj)
    new.order <- unlist(lapply(split(sTab, conds), function(x) {
        x["Sample.name"]
    }))
    qData <- Biobase::exprs(obj)[, new.order]
    sTab <- Biobase::pData(obj)[new.order, ]


    tab <- qData
    conditions <- as.factor(sTab$Condition)
    tab_imp <- impute.pa2(tab, conditions, q.min, q.norm, eps, distribution)

    # restore previous order
    colnames(tab_imp) <- new.order
    tab_imp <- tab_imp[, sample.names.old]

    Biobase::exprs(obj) <- tab_imp
    obj <- UpdateMetacellAfterImputation(obj)

    return(obj)
}






#'
#' @title Missing values imputation from a \code{MSnSet} object
#' 
#' @description 
#' 
#' This method is a variation to the function \code{impute.pa()} from the 
#' package \code{imp4p}.
#'
#' @param tab An object of class \code{MSnSet}.
#'
#' @param conditions A vector of conditions in the dataset
#'
#' @param q.min A quantile value of the observed values allowing defining the
#' maximal value which can be generated. This maximal value is defined by the
#' quantile q.min of the observed values distribution minus eps.
#' Default is 0 (the maximal value is the minimum of observed values minus eps).
#'
#' @param q.norm A quantile value of a normal distribution allowing defining
#' the minimal value which can be generated. Default is 3 (the minimal value
#' is the maximal value minus qn*median(sd(observed values)) where sd is the
#' standard deviation of a row in a condition).
#'
#' @param eps A value allowing defining the maximal value which can be
#' generated. This maximal value is defined by the quantile q.min of the
#' observed values distribution minus eps. Default is 0.
#'
#' @param distribution The type of distribution used. Values are unif or beta.
#'
#' @return The object \code{obj} which has been imputed
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' obj.imp <- wrapper.impute.pa2(Exp1_R25_pept[seq_len(100)], 
#' distribution = "beta")
#'
#' @export
#'
impute.pa2 <- function(tab,
                       conditions, 
                       q.min = 0, 
                       q.norm = 3, 
                       eps = 0, 
                       distribution = "unif"
                       ) {

    pkgs.require('stats')

    tab_imp <- tab
    qu <- apply(tab_imp, 2, stats::quantile, na.rm = TRUE, q.min)
    nb_cond <- length(levels(conditions))
    nb_rep <- rep(0, nb_cond)
    k <- 1
    j <- 1
    for (i in seq_len(nb_cond)) {
        nb_rep[i] <- sum((conditions == levels(conditions)[i]))
        sde <- apply(
            tab_imp[, seq.int(from = k, to = (k + nb_rep[i] - 1))], 
            1, 
            stats::sd,
            na.rm = TRUE
        )
        while (j < (k + nb_rep[i])) {
            if (distribution == "unif") {
                tab_imp[which(is.na(tab_imp[, j])), j] <- stats::runif(
                    n = sum(is.na(tab_imp[, j])),
                    min = qu[j] - eps - q.norm * stats::median(sde, 
                        na.rm = TRUE),
                    max = qu[j] - eps
                )
            } else if (distribution == "beta") {
                tab_imp[which(is.na(tab_imp[, j])), j] <- translatedRandomBeta(
                    n = sum(is.na(tab_imp[, j])),
                    min = qu[j] - eps - q.norm * stats::median(sde, na.rm = TRUE),
                    max = qu[j] - eps,
                    param1 = 3,
                    param2 = 1
                )
            }
            j <- j + 1
        }
        k <- k + nb_rep[i]
    }

    return(tab_imp)
}
