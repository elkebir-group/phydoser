
#' Generate the distinguishing feature family for a set of trees
#' @param treeList a list of trees as binary matrices
#' @return a list of lists of featurettes.
#' @description   For each tree in the input list, generateDistFeat finds the corresponding distinguishing feature family which
#' is a set of distinguishing features made up of root to mutation paths (featurettes) that unqiuely distinguishes that tree from
#' all other trees in the input list.
#' @export
#'
#' @examples
#' generateDistFeat(phydata$trees)
#' \dontrun{generateDistFeat(AML38$trees)}
#'
generateDistFeat <- function(treeList){
    dff <- tryCatch({
        enumerateDF(treeList)
    },
    error = function(err){
        print(paste("MY_ERROR:  ",err))
        return(list())}
    )


    for(i in 1:length( dff)){
        dff[i] <- lapply(dff[i], stringr::str_split, ",")
    }

    return(dff)
}

#' Determine the number of cells to sequence in a follow-up single cell experiment
#' @param trees a list of binary matrices representing trees with clones as the rows and mutations as the columns
#' @param fmatries a list of frequency matrices with samples as the rows and mutation frequencies (VAFs) as the columns
#' @param distFeat a optional list of distinguishing features in the same order as the trees
#' @param gamma a number that determines the confidence level of the users of observing the tree during sequencing, defaults to 0.95
#' @param fnr a number between 0 and 1 that represents the false negative rate of the sequencing technology, defaults to 0
#' @param kstar_quant the quantile that determines
#'
#' @return a list that contains the value of k* and a dataframe that contains the value of k for for each tree and sample
#' @export
#' @examples
#' phydose(trees = phydata$trees, fmatrices = phydata$fmatrices)
#' phydose(trees = phydata$trees umatrices = phdata$u_list)
#' \dontrun{
#' phydose(trees = AML38$trees, fmatrices = AML38$fmatrices, fnr= 0.05)
#' }
phydose <- function(trees, fmatrices=NULL, umatrices = NULL,
                    distFeat = NULL,
                    gamma = 0.95, fnr = 0, kstar_quant = 1){


    if(!is.list(trees)){
        trees <- list(trees)
    }




    # if(!is.list(distFeat)){
    #     distFeat <- list(distFeat)
    # }

    phyDOSE <- data.frame(stringsAsFactors = F)
    if(is.null(distFeat)){

        distFeat <- generateDistFeat(trees)

    }

    if(length(distFeat) == 0){
        return(list(kstar=NA, kT.df=NULL))
    }


    for(i in 1:length(trees)){

        if(is.null(umatrices)){

            if(is.matrix(fmatrices) || length(fmatrices)==1){
                singleFmat <- TRUE
                if(is.list(fmatrices)){
                    fmat <- fmatrices[[1]]
                }else{
                    fmat <- fmatrices
                }
            }else{
                fmat <- fmatrices[[i]]
            }




            u_mat <- CalcU(trees[[i]],fmat)
        }
        else{
            u_mat <- umatrices[[i]]
        }

        # u_mat <- try({
        #         logger.trace("class(userInput) = %s", class(userInput))
        #
        #     }, silent=TRUE)
        #
        #
        # if ( "try-error" %in% class(u_mat) ) {
        #     err_msg <- geterrmessage()
        #     print("Cannot calculate clonal prevelance matrix U from the trees and fmatrics, check inputs and try again")
        #     # logging of error message
        #     # detection and handling of particular error strings
        #      stop()# if necessary with user friendly error strings
        # }

        if(is.matrix(u_mat)){

            samples <- nrow(u_mat)
        }else{
            u_mat <- as.matrix(u_mat)
            samples <- 1
        }
        distFeatFam <- distFeat[[i]]
        curr_dff <- .RenameDF(distFeatFam, u_mat)
        for(j in 1:samples){



            u <- u_mat[j,]
            k_t <- num_cells(curr_dff, u, fnr, gamma)
            tree_name <- paste("tree",i)
            tree_df <- data.frame(tree = tree_name,
                                  tree_id = i,
                                  sample = j,
                                  cells = k_t,
                                  stringsAsFactors = F)
            phyDOSE <- dplyr::bind_rows(tree_df, phyDOSE)

        }



    }
    kstar <- .CalcKStar(phyDOSE)
    return(list(kstar=kstar, kTdata = phyDOSE))
}
.CalcKStar <- function(df, kstar_quant=1){
    phyDOSE <- df %>%dplyr:: group_by(sample) %>%
        dplyr::summarize(cells = quantile(cells, kstar_quant)) %>%
        dplyr::ungroup() %>%
        dplyr::summarize(kstar = min(cells)) %>%
        dplyr::ungroup()

    return(phyDOSE$kstar)
}


.FMultiplier <- function(fmat, mult){
    return(fmat*mult)
}

.RemoveWS <- function(x, ws= ""){
    return(x[x!=ws])
}

    .RenameDF <- function(df_fam, u){
        clones <- stringr::str_split(colnames(u), pattern = " ")
        df_rename <- list()

        for(i in 1:length(df_fam)){
            df <- df_fam[[i]][[1]]

            dfclones <- stringr::str_split(df, " ")
            dfclones <- lapply(dfclones, .RemoveWS)

           df_new <- character()
                for(d in dfclones){
                    for(c in clones){
                        if(identical(sort(d), sort(c))){
                               df_new  <- c(  paste(c, collapse = " "), df_new)
                               break
                        }
                    }

                }

              df_rename[[i]] <- df_new
            }
            # all(dfclones %in%)
            #
            # new_df_vector <- vector()
            # df <- df_fam[[i]]
            #
            # for(featurette in df){
            #
            #     featurette <- stringr::str_split(featurette, pattern = " ")[[1]]
            #
            #
            #     for(c in clones){
            #
            #         m <- match(featurette, c)
            #         m <- m[!is.na(m)]
            #         if(length(m) == length(featurette)){
            #             #print("true")
            #             new_df_vector <-unique(c(new_df_vector,paste(c, collapse = " ")))
            #
            #         }
            #     }
            # }
            #
            # df_rename[[i]] <- new_df_vector
            #

        return(df_rename)
    }


#' Generates the clonal prevalence matrix U by mutlplying the frequency matrix by
#' the inverse of the tree matrix
#' @param treeMat  a binary matrix representing a tree with clones as rows and mutations as columns
#' @param f a frequency matrix with samples as rows and mutations as columns
#' @export
#' @return a matrix of clonal prevelance with samples as rows and clones as columns
CalcU <- function(treeMat, f){
    f <- f[, colnames(treeMat)]
    t_inv <- solve(treeMat)
    u <- f %*% t_inv
    return(u)
}

# PhyDOSE <- function(filename, conf_level=0.95, fn_rate=0, fmult=1){
#
#     enum_err = FALSE
#     dff <- tryCatch({
#         enumerateDF(filename)
#         },
#              error = function(err){
#                  print(paste("MY_ERROR:  ",err))
#                  return(list())}
#              )
#     if(length(dff) == 0){
#         return(list(k_star=-1, graphs=NULL))
#     }
#     phyDOSE <- data.frame(stringsAsFactors = F)
#     resu_all <- main(filename, fmult)
#     resu <- resu_all$main
#     graphs <- resu_all$graphs
#     #Rgraphviz::plot(graphs[[1]])
#     for(i in 1:length(resu)){
#         input <- resu[[i]]
#         curr_df <- dff[[i]]
#         trees <- input$trees
#         f_matrix <- input$f
#         if(is.matrix(f_matrix)){
#
#             samples <- nrow(f_matrix)
#         }else{
#             samples <- 1
#         }
#         u_list <- input$u
#         for(j in 1:samples){
#             k_t <- calc_cells(curr_df,
#                               sample = j,
#                               tree_num = i,
#                               u_list = u_list,
#                               fn_rate = fn_rate,
#                               gamma  = conf_level)
#             tree_name <- paste("tree",i)
#             tree_df <- data.frame(tree = tree_name,
#                                   tree_num = i,
#                                   sample = j,
#                                   cells = k_t,
#                                   stringsAsFactors = F)
#             phyDOSE <- bind_rows(tree_df, phyDOSE)
#         }
#     }
#     phyDOSE.OUT <- phyDOSE %>% group_by(sample) %>%
#         summarize(cells = max(cells)) %>%
#         ungroup() %>%
#         summarize(k_star = min(cells)) %>%
#         ungroup()
#
#     k_star <- phyDOSE.OUT$k_star
#
#     phyDOSE.FILTER <- (phyDOSE[phyDOSE$cells == k_star][1,])
#     best_sample <- phyDOSE.FILTER$sample
#     #print(phyDOSE.FILTER)
#     #phyDOSE.TREE <- (phyDOSE[phyDOSE$cells == k_star][1,])$tree_num
#     best_tree <- phyDOSE.FILTER$tree_num
#     print(best_tree)
#
#     # generate data for plot
#     gamma_list <- seq(0.0, 0.99, by=0.001)
#     cells <- numeric(length(gamma_list))
#     index <- 1
#     curr_df <- dff[[1]]
#     input <- resu[[1]]
#     u_list <- input$u
#     for(g in gamma_list){
#         cell_list <- calc_cells(curr_df,
#                                 sample = best_sample,
#                                 tree_num = best_tree,
#                                 u_list = u_list,
#                                 fn_rate = fn_rate,
#                                 gamma  = g)
#         cells[index] = cell_list
#
#         index = index + 1
#     }
#     res <- data.frame(k=cells, gamma = gamma_list)
#
#
#     # return all relevant information
#     ret_all <- list(k_star=k_star, graphs=graphs, plot_vals=res, dff=dff)
#     return (ret_all)
# }
