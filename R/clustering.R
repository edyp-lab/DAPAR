#' Average protein/peptide abundances for each condition studied
#'
#' @description Calculate the average of the abundances for each protein in
#' each condition for an ExpressionSet or MSnSet. Needs to have the array
#' expression data ordered in the same way as the phenotype data (columns of
#' the array data in the same order than the condition column in the phenotype
#' data).
#'
#' @param ESet_obj ExpressionSet object containing all the data
#' 
#' @return a dataframe in wide format providing (in the case of 3 or more
#' conditions) the means of intensities for each protein/peptide
#' in each condition. If there are less than 3 conditions, an error message is
#' returned.
#' 
#' @author Helene Borges
#' 
#' @examples
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj <- MetaCellFiltering(obj, indices, cmd='delete')
#' averageIntensities(obj$new)
#' 
#' @export
#' 
#' @importFrom dplyr bind_cols as_tibble group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom Mfuzz standardise
#' 
averageIntensities <- function(ESet_obj){
    intensities <- exprs(ESet_obj)
    sTab <- pData(ESet_obj)
    sTab$Condition <- as.factor(sTab$Condition)
    intensities_t <- as.data.frame(t(intensities))
    intensities_t <- dplyr::bind_cols(intensities_t, condition = sTab$Condition, sample = rownames(intensities_t))
    tbl_intensities <- dplyr::as_tibble(intensities_t, rownames = NA)
    longer_intensities <- tbl_intensities %>%
        tidyr::pivot_longer(-c(condition,sample), names_to = "feature", values_to = "intensity")
    mean_intensities <- longer_intensities %>%
        dplyr::group_by(condition, feature) %>%
        dplyr::summarise(mean = mean(intensity))
    mean_intensities_wide <- tidyr::pivot_wider(mean_intensities,
                                                names_from = condition,
                                                values_from = mean)
    averaged_data <- as.data.frame(mean_intensities_wide)
    rownames(averaged_data) <- averaged_data[,1]
    
    return(averaged_data)
    
}



#' Prepare the data for subsequent clustering
#'
#' @description The first step is to standardize the data (with the \code{Mfuzz}
#' package). Then the function checks that these data are clusterizable or not
#'  (use of [diptest::dip.test()] to determine whether the distribution is 
#'  unimodal or #'  multimodal). Finally, it determines the "optimal" k by 
#'  the Gap statistic approach.
#'
#' @param standards a matrix or dataframe containing only the standardized 
#' mean intensities returned by the function [standardiseMeanIntensities()]
#' 
#' @param b Parameter B of the function [gap_cluster()] 
#' 
#' @return a list of 2 elements:
#' * dip_test: the result of the clusterability of the data
#' * gap_cluster: the gap statistic obtained with the function 
#' [cluster::clusGap()].
#' 
#' @author Helene Borges
#' 
#' @examples
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot[1:100]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj <- MetaCellFiltering(obj, indices, cmd='delete')
#' averaged_means <- averageIntensities(obj$new)
#' only_means <- dplyr::select_if(averaged_means, is.numeric)
#' only_features <- dplyr::select_if(averaged_means, is.character)
#' means <- purrr::map(purrr::array_branch(as.matrix(only_means), 1),mean)
#' centered <- only_means - unlist(means)
#' centered_means <- dplyr::bind_cols(feature = dplyr::as_tibble(only_features), 
#' dplyr::as_tibble(centered))
#' checkClust <- checkClusterability(centered_means, b=100)
#' 
#' @export
#' 
#' @importFrom cluster clusGap
#' @importFrom diptest dip.test
#' 
checkClusterability <- function(standards, b = 500){
    
    if(methods::is(standards, "data.frame")){
        # if there are columns with something other than numeric values
        if(!is.numeric(standards)){
            # then we only keep the columns with numeric values
            standards <- dplyr::select_if(standards, is.numeric)
        }
        # we transform into a matrix to be usable by dip.test
        standards <- as.matrix(standards)
    }else if(!methods::is(standards, "matrix") || !methods::is(standards, "data.frame")){
        stop("Input must be a matrix or a dataframe")
    }
    # check the clusterability of the data
    dip_res <- diptest::dip.test(x = standards)
    print("dip test done")
    # d.power = 2 corresponds to the Tibshirani criteria. B = 500 is to have 
    # stable results from one simulation to another.
    gap_cluster <- cluster::clusGap(standards, FUNcluster = kmeans, nstart = 20, K.max = 10, d.power = 2, B = b)
    
    return(list(
        "dip_test" = dip_res,
        "gap_cluster" = gap_cluster
    ))
}

#' Visualize the clusters according to pvalue thresholds
#' 
#' @param dat the standardize data returned by the function [checkClusterability()]
#' 
#' @param clust_model the clustering model obtained with dat.
#' 
#' @param adjusted_pValues vector of the adjusted pvalues obtained for each 
#' protein with a 1-way ANOVA (for example obtained with the function 
#' [wrapperClassic1wayAnova()]).
#'  
#' @param FDR_th the thresholds of FDR pvalues for the coloring of the profiles.
#' The default (NULL) creates 4 thresholds: 0.001, 0.005, 0.01, 0.05
#' For the sake of readability, a maximum of 4 values can be specified.
#' 
#' @param ttl title for the plot.
#' 
#' @param subttl subtitle for the plot.
#' 
#' @return a ggplot object
#' 
#' @author Helene Borges
#' 
#' @examples
#' library(dplyr)
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj <- MetaCellFiltering(obj, indices, cmd='delete')
#' expR25_ttest <- compute_t_tests(obj$new)
#' averaged_means <- averageIntensities(obj$new)
#' only_means <- dplyr::select_if(averaged_means, is.numeric)
#' only_features <- dplyr::select_if(averaged_means, is.character)
#' means <- purrr::map(purrr::array_branch(as.matrix(only_means), 1),mean)
#' centered <- only_means - unlist(means)
#' centered_means <- dplyr::bind_cols(feature = dplyr::as_tibble(only_features), 
#' dplyr::as_tibble(centered))
#' difference <- only_means[,1] - only_means[,2]
#' clusters <- as.data.frame(difference) %>%
#'     dplyr::mutate(cluster = dplyr::if_else(difference > 0, 1,2))
#' vizu <- visualizeClusters(dat = centered_means,
#'                           clust_model = as.factor(clusters$cluster),
#'          adjusted_pValues = expR25_ttest$P_Value$`25fmol_vs_10fmol_pval`,
#'                           FDR_th = c(0.001,0.005,0.01,0.05),
#'                           ttl = "Clustering of protein profiles")
#'                           
#' @export
#' 
#' @importFrom stringr str_glue
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom methods is
#' 
visualizeClusters <- function(dat, 
                              clust_model, 
                              adjusted_pValues, 
                              FDR_th = NULL, 
                              ttl = "", 
                              subttl = ""){
    if(is.null(FDR_th)){
        FDR_th <- c(0.001,0.005,0.01,0.05)
    }else if(length(FDR_th) > 4){
        message("Too many FDR thresholds provided. Please do not exceed 4 values")
        return(NULL)
    }
    
    str_try <- stringr::str_glue("<{FDR_th}")
    str_max <- stringr::str_glue(">{max(FDR_th)}")
    
    dat$FDR_threshold <- cut(adjusted_pValues, breaks = c(-Inf,FDR_th,Inf), labels = c(str_try, str_max))
    desc_th <- FDR_th[order(FDR_th, decreasing = TRUE)]
    str_desc <- stringr::str_glue("<{desc_th}")
    
    dat$FDR_threshold <- factor(dat$FDR_threshold,
                                levels = c(str_max, str_desc))
    
    dat$adjusted_pvalues <- adjusted_pValues
    if(methods::is(clust_model, "factor")){
        dat$cluster <- clust_model
    }else if(methods::is(clust_model, "kmeans")){
        dat$cluster <- as.factor(clust_model$cluster)
        
    }else if(methods::is(clust_model, "APResult")){
        dat$cluster <- NA
        for(k in seq_len(length(clust_model@clusters))){
            dat$cluster[clust_model@clusters[[k]]] <- k
        }
        
    }else{
        message("Wrong model. Currently, only the Kmeans and the Affinity
                Propagation models are supported.")
        return(NULL)
    }
    
    melted <- reshape2::melt(dat,
                             id.vars = c("feature", "cluster", "adjusted_pvalues", "FDR_threshold"),
                             value.name = "intensity",
                             variable.name = "Condition")
    
    palette2use <- c("#DEDEDE","#FFC125", "#FF8C00", "#CD3700", "#8B0000")
    names(palette2use) <- levels(melted$FDR_threshold)
    
    colScale <- ggplot2::scale_colour_manual(name = "FDR threshold",
                                             values = palette2use)
    
    
    melted <- dplyr::arrange(melted, desc(adjusted_pvalues))
    return(ggplot2::ggplot(data = melted,
                           ggplot2::aes(x = Condition, y = intensity, col = FDR_threshold)) +
               ggplot2::geom_line(ggplot2::aes(group = feature)) +
               ggplot2::facet_wrap(~ cluster) +
               colScale +
               ggplot2::xlab("Condition") +
               ggplot2::ylab("Standardized averaged intensity") +
               ggplot2::labs(title = ttl, subtitle = subttl) +
               ggplot2::theme(
                   plot.title = ggplot2::element_text(size = 15, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 12),
                   panel.spacing = ggplot2::unit(1.5, "lines"),
                   panel.background = ggplot2::element_rect(fill = NA),
                   panel.grid.major = ggplot2::element_line(linetype = 'solid',
                                                            size = 0.25,
                                                            colour = "#f6f6f6"),
                   panel.grid.minor = ggplot2::element_blank()
               )
    )
}



#' Run a clustering pipeline of protein/peptide abundance profiles.
#'
#' @description This function does all of the steps necessary to obtain a
#' clustering model and its graph from average abundances of proteins/peptides.
#'  It is possible to carry out either a kmeans model or an affinity
#'  propagation model. See details for exact steps.
#'  
#' @param obj ExpressionSet or MSnSet object.
#' 
#' @param clustering_method character string. Three possible values are 
#' "kmeans", "affinityProp" and "affinityPropReduced. See the details section 
#' for more explanation.
#' 
#' @param conditions_order vector specifying the order of the Condition factor
#' levels in the phenotype data. Default value is NULL, which means that it is 
#' the order of the condition present in the phenotype data of "obj" which is
#' taken to create the profiles.
#' 
#' @param k_clusters integer or NULL. Number of clusters to run the kmeans
#' algorithm. If `clustering_method` is set to "kmeans" and this parameter is
#' set to NULL, then a kmeans model will be realized with an optimal number of
#' clusters `k` estimated by the Gap statistic method. Ignored for the Affinity
#'  propagation model.
#'  
#' @param adjusted_pvals vector of adjusted pvalues returned by the 
#' [wrapperClassic1wayAnova()]
#' 
#' @param ttl the title for the final plot
#' 
#' @param subttl the subtitle for the final plot
#' 
#' @param FDR_thresholds vector containing the different threshold
#' values to be used to color the profiles according to their adjusted pvalue.
#' The default value (NULL) generates 4 thresholds: [0.001, 0.005, 0.01, 0.05].
#'  Thus, there will be 5 intervals therefore 5 colors: the pvalues <0.001,
#'  those between 0.001 and 0.005, those between 0.005 and 0.01, those between
#'  0.01 and 0.05, and those> 0.05. The highest given value will be considered
#'  as the threshold of insignificance, the profiles having a pvalue> this
#'  threshold value will then be colored in gray.
#'  
#' @details The first step consists in averaging the abundances of
#' proteins/peptides according to the different conditions defined in the
#' phenotype data of the expressionSet / MSnSet. Then we standardize the data 
#' if there are more than 2 conditions. If the user asks to realize a kmeans 
#' model without specifying the desired number of
#' clusters (`clustering_method =" kmeans "` and `k_clusters = NULL`), the 
#' function checks data's clusterability and estimates a number of clusters k 
#' using the gap statistic method. It is advise however to specify a k for the 
#' kmeans, because the gap stat gives the smallest possible k, whereas in 
#' biology a small number of clusters can turn out to be uninformative. 
#' If you want to run a kmeans but you don't know what number of clusters to 
#' give, you can let the pipeline run the first time without specifying 
#' `k_clusters`, in order to view the profiles the first time and choose by the
#' following is a more appropriate value of k.
#' If it is assumed that the data can be structured with a large number of
#' clusters, it is recommended to use the affinity propagation model instead.
#' This method simultaneously considers all the data as exemplary potentials,
#' unlike hard clustering (kmeans) which initializes with a number k of points
#' taken at random. The "affinityProp" model will use a q parameter set to NA,
#' meaning that exemplar preferences are set to the median of non-Inf values
#' in the similarity matrix (set q to 0.5 will be the same). The
#' "affinityPropReduced" model will use a q set to 0, meaning that exemplar
#' preferences are set to the sample quantile with threshold 0 of non-Inf
#' values. This should lead to a smaller number of final clusters.
#' 
#' @author Helene Borges
#' 
#' @return a list of 2 elements: "model" is the clustering model, "ggplot" is 
#' the ggplot of profiles clustering.
#' 
#' @references
#' Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number of 
#' data clusters via the Gap statistic. 
#' *Journal of the Royal Statistical Society* B, 63, 411–423.
#'
#' Frey, B. J. and Dueck, D. (2007) Clustering by passing messages between 
#' data points. *Science* 315, 972-976. 
#' DOI: \href{https://science.sciencemag.org/content/315/5814/972}{10.1126/science.1136800}
#'
#' @examples
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj <- MetaCellFiltering(obj, indices, cmd='delete')
#' expR25_ttest <- compute_t_tests(obj$new)
#' wrapperRunClustering(obj = obj$new, 
#' adjusted_pvals = expR25_ttest$P_Value$`25fmol_vs_10fmol_pval`)
#' 
#' @export
#' 
#' @importFrom apcluster apcluster
#' @importFrom forcats as_factor fct_relevel
#' @importFrom stats kmeans
#' 
wrapperRunClustering <- function(obj, clustering_method, conditions_order = NULL, k_clusters = NULL, adjusted_pvals, ttl = "", subttl = "", FDR_thresholds = NULL){
    res <- list("model" = NULL, "ggplot" = NULL)
    # reorder conditions if requested
    if(!is.null(conditions_order)){
        # check that given levels are correct in number...
        if(length(conditions_order) == length(levels(forcats::as_factor(obj@phenoData@data$Condition)))){
            # ...and have same labels
            if(all(conditions_order %in% levels(forcats::as_factor(obj@phenoData@data$Condition)))){
                obj@phenoData@data$Condition <- forcats::fct_relevel(obj@phenoData@data$Condition, conditions_order)
            }else{
                valid_labels <- str_c(levels(forcats::as_factor(obj@phenoData@data$Condition)), collapse = ", ")
                message(stringr::str_glue("Wrong labels given. The valid labels are {valid_labels}"))
                return(NULL)
            }
        }else{
            y <- length(conditions_order)
            n <- length(levels(forcats::as_factor(obj@phenoData@data$Condition)))
            message(stringr::str_glue("Wrong number of given levels. There are currently {n} levels for Condition. You provided {y} levels."))
            return(NULL)
        }
    }
    
    averaged_means <- averageIntensities(obj)
    averaged_means <- na.omit(averaged_means)
    
    sTab <- pData(obj)
    only_means <- dplyr::select_if(averaged_means, is.numeric)
    only_features <- dplyr::select_if(averaged_means, is.character)
    if(length(levels(as.factor(sTab$Condition))) == 2){
        print("there are only two conditions")
        means <- purrr::map(purrr::array_branch(as.matrix(only_means), 1),mean)
        centered <- only_means - unlist(means)
        centered_means <- dplyr::bind_cols(feature = dplyr::as_tibble(only_features), dplyr::as_tibble(centered))
        
        difference <- only_means[,1] - only_means[,2]
        clusters <- as.data.frame(difference) %>%
            dplyr::mutate(cluster = dplyr::if_else(difference > 0, 1,2))
        
        res$model <- as.factor(clusters$cluster)
        
        res$ggplot <- visualizeClusters(dat = centered_means,
                                        clust_model = res$model,
                                        adjusted_pValues = adjusted_pvals,
                                        FDR_th = FDR_thresholds,
                                        ttl = ttl,
                                        subttl = subttl
        )
    }else if(length(levels(as.factor(sTab$Condition))) > 2){
        eset <- ExpressionSet(assayData = as.matrix(only_means))
        eset_s <- Mfuzz::standardise(eset)
        standards <- eset_s@assayData$exprs
        
        standardized_means <- dplyr::bind_cols(feature = dplyr::as_tibble(only_features), dplyr::as_tibble(standards))
        names(standardized_means) <- c("feature", names(dplyr::as_tibble(standards)))
        if(clustering_method == "affinityProp"){
            res$model <- apcluster::apcluster(apcluster::negDistMat(r=2), standards)
            
        }else if(clustering_method == "affinityPropReduced"){
            print("running affinity propagation reduced algorithm")
            res$model <- apcluster::apcluster(apcluster::negDistMat(r=2), standards, q=0)
        }else if(clustering_method == "kmeans"){
            if(is.null(k_clusters)){
                res <- list("model" = NULL, "dip" = NULL, "ggplot" = NULL)
                checked_means <- checkClusterability(standards)
                res$dip <- checked_means$dip_test
                best_k <- cluster::maxSE(checked_means$gap_cluster$Tab[, "gap"],
                                         checked_means$gap_cluster$Tab[, "SE.sim"],
                                         method="Tibs2001SEmax")
                res$model <- kmeans(standards, centers = best_k, nstart = 25)
            }else if(k_clusters > 1){
                best_k <- k_clusters
                res$model <- kmeans(standards, centers = best_k, nstart = 25)
            }else{ # corresponds to the case k = 1 so no need for clustering
                one_cluster <- as.factor(rep_len(1, nrow(standardized_means)))
                res$ggplot <- visualizeClusters(dat = standardized_means,
                                                clust_model = one_cluster,
                                                adjusted_pValues = adjusted_pvals,
                                                FDR_th = FDR_thresholds,
                                                ttl = ttl,
                                                subttl = subttl
                )
                return(res)
            }
        }else{
            message("Wrong method given. Valid names are affinityProp, affinityPropReduced and kmeans.")
            return(NULL)
        }
        print("generating gglot...")
        res$ggplot <- visualizeClusters(dat = standardized_means,
                                        clust_model = res$model,
                                        adjusted_pValues = adjusted_pvals,
                                        FDR_th = FDR_thresholds,
                                        ttl = ttl,
                                        subttl = subttl
        )
    }
    
    return(res)
}
