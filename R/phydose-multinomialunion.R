

##########################################################################################
#
# Multinomial Tail Probability Power Calculation for Multiple Distinguishing Features
#
##############################################################################################



calc_cells <- function(df_fam, sample, tree_num, u_list, fn_rate, gamma){
    
    u_mat <- u_list[[sample]]
    df_fam <- rename_df(df_fam, u_mat)
    u <- u_mat[sample,]
    sample_k <- num_cells(df_fam, u, fn_rate, gamma)
    
    return(sample_k$cells)
}


##################################################################
#
#   Num_Cells is the main function to conduct the power calculation
#   on the number of cells required to see each distinguishing feature
#   at least once with probability gamma.
#
#  @param dfFamily - a distinguishing feature family where each 
#                      item in the list is a distinguishing feature comprised 
#                      of featurettes
#
#  @param u          - a vector indicating the probability of each featurette (aka clone)
#  
#  @param gamma       - the probability with which you want to see each distinguishing feature
#                     - at least once 
##################################################################
num_cells <- function(dfFamily, u, fn, gamma){
    
    #if dfFamily is passed a
    if(!is.list(dfFamily)){
        dfFamily <- list(dfFamily)
    }
    
    dfFamily <- filter_df(dfFamily, u)
    num_sets <- length(dfFamily)
    set_indices = 1:num_sets
    cells <- 0
    prob_not_met <- TRUE
    probs <- NA
    
    
    
    
    #if the list of distinguishing features with nonzero featurettes is empty then return Inf cells
    if(length(dfFamily) == 0){
        return(list(cells = Inf))
        
    }else if(length(dfFamily) == 1){
        path_length_vec <-.GetPathLengths(dfFamily[[1]])
        prob_vec <- u[dfFamily[[1]]]
        prob_vec <- .AdjustFN(prob_vec, path_length_vec, fn)
        prob_vec <- c(prob_vec, 1- sum(prob_vec))
        lower_bound <- c(rep(0, length(prob_vec)-1), -Inf)
        cells <- pmultinom::invert.pmultinom(lower = lower_bound, probs = prob_vec, target.prob = gamma, method="exact")
    }else{
        
        
        while(prob_not_met){
            cells = cells + 1
            #print(paste0("cells:", cells))
            p <- 0
            j <- 1
            while(j <= num_sets){
                combos <- combn(set_indices, j, simplify = FALSE)
                for(c in combos){
                    #print(c)
                    num_df_sets = length(c)
                    
                    if(num_df_sets == 1){
                        prob_vec <- u[dfFamily[[c]]]
                        path_length_vec <-.GetPathLengths(dfFamily[[c]])
                        
                        prob_vec <- .AdjustFN(prob_vec, path_length_vec, fn)
                        
                        item_prob <- calc_prob(prob_vec, cells)
                        
                        
                    }else{ #take the union of the sets
                        intersection <- NA
                        for(i in c){
                            
                            intersection <- base::union(dfFamily[[i]], intersection)
                        }
                        
                        intersection <- intersection[!is.na(intersection)]
                        prob_vec <- u[intersection]
                        path_length_vec <-.GetPathLengths(intersection)
                        
                        prob_vec <- .AdjustFN(prob_vec, path_length_vec, fn)
                        
                        #print(intersection)
                        item_prob <- calc_prob(prob_vec, cells)
                        
                    }
                    
                    if(num_df_sets %% 2 == 0){
                        p = p - item_prob
                    }else{
                        p = p + item_prob
                    }
                    
                    
                    
                }
                
                j= j + 1
                
            }
            
            probs[cells] = p
            
            if(probs[cells] >= gamma){
                prob_not_met = FALSE
            }
        }
    }
    
    
    
    return(list(cells = cells))
    
    
}


###########################################################
#
#  Helper function to calcule the tail probability of a
#   distinguishing feature for a particular number of cells
#  
#  @param prob_vec - the vector of probabilities for a distinguishing feature
# 
#  @param cells - the desired number of cells for which the multinomial tail
#                 probability should be calculated
# 
#
###########################################################



calc_prob <- function(prob_vec, cells){
    
    
    
    if(sum(prob_vec) == 1){
        lower_bounds <- c(rep(0, length(prob_vec)))
        
    }else{
        prob_vec <- c(prob_vec, "x"=1- sum(prob_vec))
        lower_bounds <- c(rep(0, length(prob_vec) -1), -Inf)
    }
    
    #threshold and renormalize
    prob_vec <- ifelse(prob_vec < 1e-11, 0, prob_vec)
    prob_vec <- prob_vec/sum(prob_vec)
    
    
    if(any(prob_vec==0)){
        return(Inf)
    }
    
    
    
    
    
    
    
    
    pmult <- pmultinom::pmultinom(lower = lower_bounds, size = cells, probs = prob_vec, method = "exact")
    
    return(pmult)
}

.GetPathLengths <- function(distFeat){
    paths <- numeric(length(distFeat))
    
    for(i in 1:length(distFeat)){
        feat <- distFeat[i]
        feat <-stringr::str_split(feat, " ")[[1]]
        paths[i] <- length(feat)
    }
    return(paths)
}

.AdjustFN <- function(probs, path_vec, fn){
    adjustvec <- (1-fn)^path_vec
    return(adjustvec * probs)
    
}


####################################################################################


#####################################################################################
#
#             Helper Function to remove any distinguishing features from the family
#                  that there is zero probability of seeing
#         
#
#   @param dfFamily - a list of distinguishing feautures
# 
#   @param u - - a vector indicating the probability of each featurette (aka clone)

#  returns a filtered dfFamily where all featurettes of all distinguishing features
#          have non-zero probabilities

###########################################################################

#removes and distinguishing features that have a featurette with probability 0
filter_df <- function(dfFamily, u){
    #scan each distinguishing feature and check if any featurette has a prob of zero
    #remove the entire distinguishing feature if it does
    fil_list <- list()
    index <- 1
    for(df in dfFamily){
        if(length(df) > 0){
            if(all(u[names(u) %in% df] > 1e-5)){
                fil_list[[index]] <- df
                index <- index + 1
            }
            
        }
        
    }
    
    return(fil_list)
}



