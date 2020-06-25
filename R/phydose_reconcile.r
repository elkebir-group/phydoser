#takes in dataframe wth all the trees  scs_file_name, number of cells to draw)
#returns a boolean vector with whether or not that tree was identified by that sample

#' Calculate the support for each distinguishing feature in a given set of single-cells
#' @param dffs a list of distinguishing feature families for trees in the set
#' @param a dataframe of sampled single cells
#'
#' @return a dataframe that returns the support for each tree in the set of input cells
Reconcile <- function(dffs, cells){
  tree_support <- numeric(0)

  for(i in 1:length(dffs)){

    dfFamily <- dffs[[i]]

    list.df <- convert_dfFamily_to_scs(dfFamily, colnames(cells))
    pres <- unlist(lapply(list.df, is_present, cells))
    supported_df <- data.frame()
    for(k in 1:length(pres)){
      if(pres[k] > 0){
        supported_df <- dplyr::bind_rows(list.df[[k]], supported_df)
      }
    }
    supported_df <- supported_df %>% dplyr::distinct()
    if(nrow(supported_df) > 0){
      support <- count_cells(supported_df,cells )
    }else{
      support <- 0
    }

    tree_support[i] <- support
  }
  support.df <- data.frame(tree = 1:length(dffs), support = tree_support)
  return(support.df)
}

# scs_sim <- function(tree_df, cells, scs.dat,df_folder_path,  trials, stat, fn, write_file){
#
#
#
#   bool_replace <- ifelse(cells > nrow(scs.dat), TRUE, FALSE)
#   cells <- ceiling(cells)
#   #num_success <- rep(0, nrow(tree_df))
#   scs.all <- data.frame()
#   tree_support <- data.frame()
#   for(i in 1:trials){
#
#     scs.samples <- sample_n(scs.dat, cells, replace = bool_replace)
#
#     for(j in 1:nrow(tree_df)){
#       tree <- tree_df[j,]
#
#       dfFamily <- df_family(file.path(df_folder_path, tree$tree_file_name))
#       list.df <- convert_dfFamily_to_scs(dfFamily, colnames(scs.samples))
#       pres <- unlist(lapply(list.df, is_present, scs.samples))
#       supported_df <- data.frame()
#       for(k in 1:length(pres)){
#         if(pres[k] > 0){
#           supported_df <- bind_rows(list.df[[k]], supported_df)
#         }
#       }
#       supported_df <- supported_df %>% distinct()
#       if(nrow(supported_df) > 0){
#         support <- count_cells(supported_df,scs.samples )
#       }else{
#         support <- 0
#       }
#
#
#
#       df <-unlist(lapply(dfFamily, paste, collapse = ":"))
#       df_res <- as.data.frame(df)
#       df_res$df_id =1:length(df)
#       df_res$success <- pres
#       df_res$support <- support
#       df_res$trial <- i
#       df_res$tree_file_name <- unique(tree$tree_file_name)
#       df_res$tree_num <- unique(tree$tree_num)
#       df_res$cells <- cells
#
#       #print(tail(df_res))
#       tree_support <- bind_rows(df_res, tree_support)
#
#
#     }
#
#
#
#     scs.all <- bind_rows(scs.samples %>% mutate(trial = i), scs.all)
#   }
#
#
#
#
#
#
#
#   scs.all$stat <- stat
#   scs.all$fn <- fn
#
#   #write.csv(scs.all,write_file, row.names=F)
#   return(tree_support)
#
#
# }


convert_dfFamily_to_scs <- function(dfList, cnames){
  data.list <- lapply(dfList, convert_to_vector, cnames)
  return(data.list)
}

convert_to_vector <- function(vec, cnames){
  if(substr(cnames[1],0,1) == "X"){
    muts <- as.numeric(str_extract(cnames, '[0-9]+'))
  }
  else{
    muts <- cnames
  }
  #
  names.df <- data.frame(n =cnames)
  vec2 <- stringr::str_split(vec, pattern = " ", simplify = FALSE)
  for(v in vec2){
    x <- ifelse(muts %in% v, 1, 0)
    names.df <- cbind(names.df,x)
  }
  names.df <- as.data.frame(t(names.df[,-1]))
  colnames(names.df) <- cnames
  rownames(names.df) <- NULL

  return(names.df)


}


count_cells <- function(df, scs.dat){
  df_support <- nrow(dplyr::semi_join(scs.dat, df, by= colnames(scs.dat)))

  return(df_support)

}

is_present <- function(df, scs.dat){
  temp <- dplyr::semi_join(df, scs.dat, by = colnames(scs.dat))
  pres <- ifelse(nrow(temp) == nrow(df), 1, 0)

  return(pres)

}


