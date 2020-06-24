
#' Read in set of candidate trees in the given directory and
#' convert to a binary matrix with clones as the rows and mutations as columns
#'
#' @param fname the path to directory where the trees are located
#' @return A list of binary matrices representing all trees located in the directory
#' @export
ReadTrees <- function(fname){

    fil <- readLines(fname)
    trees_graphs <- .GetAllTrees(fil)
    return(trees_graphs)

}


#' Creates a list of of the parents for each node
#'
#' @param vecs the portion of the file that contains the edge list of the tree
#' @return a list of vectors for each node that contains the names of all ancestors of that node
.CreateParentList <- function(vecs){

  parents_all <- list()
  nodes <- character(0)
  branches <- stringr::str_split(vecs, pattern = " ")
  for(b in branches){

    for(i in 1:length(b)){
      if(!(b[i] %in% nodes)){
        nodes <- c(b[i], nodes)
      }
      #it's a child so add the immediate parent
      if( i > 1 ){

        parents_all[b[i]] <- b[1]
      }
    }
  }

  root_node <- nodes[!(nodes %in% names(parents_all))]



  continue <- T
  while(continue){
    cur_len <- lapply(parents_all, length)
    for(p in names(parents_all)){

      if( p != root_node){
        for(b in parents_all[[p]]){

          if(b != root_node){
            parents_all[[p]] <- unique(c(parents_all[[p]], parents_all[[b]]))
          }
        }


      }


    }

    new_len <- lapply(parents_all, length)

    if(identical(new_len, cur_len)){
      break;
    }

  }

  parents_all[[root_node]] <- NA

  return(parents_all)

}


#' Create a key that maps the columns of the trees and frequencies matrices to
#' the mutations
#' @param mat an input matrix where the column names are used as the key
#' @return a dataframe creating the mapping of column ids to mutations
.CreateKey <- function(mat){
  muts <- colnames(mat)
  ids <- 0:(length(muts)-1)
  key <- data.frame(id = ids, mutation = muts, stringsAsFactors = F)
  return(key)
}


#' Read in data from a directory
#' @param dir the directory where the files are stored
#' @return a list
#' @export

ReadJointDir <- function(dir, fmultiplier=1){


    files <- list.files(dir , include.dirs = TRUE)



  trees <- list()
  fmatrices <- list()
  u_list <- list()
  graphs <- list()
  for(f in files){
   #print(f)
   file_out <- ReadJoint(file.path(dir,f), fmultiplier)
   trees <- c(trees, file_out$trees)
   fmatrices <- c(fmatrices, file_out$fmatrices)
   u_list <- c(u_list, file_out$u_list)
   graphs <- c(graphs, file_out$graph)





  }
  results <- list(trees = trees, fmatrices = fmatrices, u_list = u_list, graph = graphs )

  return(results)
}
#' Read a file containing a candidate set of trees and an accompnaying F matrix for each tree.
#' @param fname the path to the file that is to be read
#' @param fmultiplier a multiplier for the frequency values, defaults to 1
#' @param fvalues a string specifying whether the "fplus" or "fminus" values should be read, defaults to "fplus"
#' @return a list
#' @export
#' @importFrom magrittr %>%
ReadJoint <- function(fname, fmultiplier=1, fvalues = "fplus"){



    fil <- readLines(fname)


    trees_dots <- .GetAllTrees(fil)


    trees <- trees_dots$trees

    dots <- trees_dots$graph

    fmats <- .GetAllF(fil, fvalues)


    #key <- .CreateKey(trees[[1]])
    if(fmultiplier != 1 && fmultiplier > 0){
      fmats <- lapply(fmats, .FMultiplier, fmultiplier)
    }

    u_list <- list()
    for (i in 1:length(fmats)){
      tree <- trees[[i]]
      fmat <- fmats[[i]]



      u_list[[i]] <- CalcU(tree,fmat)

      # if(!is.matrix(fmat)){
      #   fmat<- t(as.matrix(fmat))
      # }


      # temp.df <- data.frame(mutation = colnames(tree), stringsAsFactors = F)
      # temp.df <- dplyr::left_join(temp.df, key, by="mutation")
      # colnames(tree) <- temp.df$id
      # colnames(fmat) <- temp.df$id
      #
      # fmats[[i]] <- fmat
      # trees[[i]] <- tree
    }
    result <- list(trees = trees, fmatrices = fmats, u_list = u_list, graph = dots )

    return(result)
}



#' Helper function to read in all trees from a file and convert to binary matrices
#'
#' @param charVec a character vector with all lines of the file containing the trees
#' and frequency matrices alternating
#'
#' @return a list of trees in the format of a binary matrix with clones as rows and mutations as columns
.GetAllTrees <- function(charVec){
  tree_count <- 0
  s <- NULL
  e <- NULL
  trees <- list()
  dots <- list()
  for(i in 1:(length(charVec) + 1)){

    # if(!is.null(s) && !is.null(e)){
    #   vecs <- charVec[s:e]
    #   tree_count <- tree_count + 1
    #   print(s)
    #   print(e)
    #   trees[[tree_count]] <- .CreateTree(vecs)
    #
    #   dots[[tree_count]] <- .CreateGraphDot(vecs)
    #
    #   s <- NULL
    #   e <- NULL
    # }

    if(i > length(charVec)){
      break
    }




    if(stringr::str_detect(charVec[i], "sites")|| (stringr::str_detect(charVec[i], "tree [0-9]+") && !is.null(s))){
      e <- i -1
      vecs <- charVec[s:e]
      tree_count <- tree_count + 1

      trees[[tree_count]] <- .CreateTree(vecs)

      dots[[tree_count]] <- .CreateGraphDot(vecs)

     if((stringr::str_detect(charVec[i], "tree [0-9]+"))){
       s <- i + 1
     }else{
       s <- NULL
     }
      e <- NULL
    }

    if(stringr::str_detect(charVec[i], "tree [0-9]+") && is.null(s) ){
      s <- i + 1
    }

  }

  return(list(trees = trees, graph = dots))
}

#' Converts all frequency matrices in the character vector to a list of frequency matrices
#'
#' @param charVec a character vector that contains all of the text from the input file
#' that has trees and frequency matrices alternating
#'
#' @return a list of frequency matrices
#'
.GetAllF <- function(charVec, fvalues){
  f_count <- 0
  s <- NULL
  e <- NULL
  fmat <- list()
  for(i in 1:(length(charVec) + 1)){

    if(!is.null(s) && !is.null(e) ){

      vecs <- charVec[s:e]

      f_count <- f_count + 1
      fmat[[f_count]] <- .CreateF(vecs, fvalues)
      s <- NULL
      e <- NULL
    }

    if(i > length(charVec)){
      break
    }
    if(stringr::str_detect(charVec[i], "sample_index")){
      s <- i
    }

    if(stringr::str_detect(charVec[i], "tree [0-9]+") && !is.null(s)){
      e <- i-1
    }

    if(i == length(charVec)){
      e <- i
    }
  }

  return(fmat)
}


#' Cleans the raw dataframe frequency matrix and converts it into a matrix with samples as rows and mutations as columns
#' @param df the frequency dataframe read in from the input file
#'
#' @return a matrix representing the frequency matrix with samples as rows and columns as mutations
#' @importFrom magrittr %>%
.CleanFmatrix <-function(df){
  fmatrix <- df %>% dplyr::select(sample_index, character_label, fplus) %>%
    dplyr::mutate(fplus = as.numeric(fplus)) %>%
    tidyr::pivot_wider(id_cols = sample_index, names_from = character_label, values_from = fplus)  %>%
    dplyr::select(-sample_index)
  fmatrix <- as.matrix(fmatrix)

  return(fmatrix)
}



.CleanFmatrixFminus <-function(df){
  fmatrix <- df %>% dplyr::select(sample_index, character_label, fminus) %>%
    dplyr::mutate(fplus = as.numeric(fminus)) %>%
    tidyr::pivot_wider(id_cols = sample_index, names_from = character_label, values_from = fminus)  %>%
    dplyr::select(-sample_index)
  fmatrix <- as.matrix(fmatrix)

  return(fmatrix)
}

#' Creates an F matrix from character vector
#' @param vecs the character vector from which the frequency matrix is created
#' @param fvaluse a string either "fplus" or "fminus" specifying with frequencies should be read
#' @return a matrix of the frequencies
.CreateF<- function(vecs, fvalues){


      vecs <- stringr::str_replace(vecs,"#", "")

      header <- unlist(stringr::str_split(vecs[1], "\t"))

      header <- c(header[1:6], c("fminus", "fplus"))

      vecs <- stringr::str_split(vecs[-1], "\t")
      fmatrix <-do.call("rbind", vecs)
      colnames(fmatrix) <- header

      fmatrix <- as.data.frame(fmatrix, stringsAsFactors=F)
      if(fvalues == "fplus"){
        fmatrix <- .CleanFmatrix(fmatrix)
      }else{
        fmatrix <- .CleanFmatrixFminus(fmatrix)
      }


  return(fmatrix)
}

#' Create a tree in dot format
#' @param vec character vector with the edge list
#' @return a character vector representing the tree in dot format
.CreateGraphDot<- function(vec){
  # input "0 3" "1 2" "2 4" "3 1"
  # output nodes (vector)
  # nodes <- c("0", "1", "2", "3", "4")

  # edges
  # edgeList <- list(0=list(edges=c("3")),
  #                 1=list(edges=c("2")),
  #                 2=list(edges=c("4")),
  #                 3=list(edges=c("1")))

  # graph <- new("graphNEL", nodes=nodes, edgeL=edgeList, edgemode="undirected")
  dot <- "digraph T{"
  nodes <- c()
  # edgeList <- list()
  for(i in 1:length(vec)){
    edgestr <- strsplit(vec[[i]], ' ')
    if (!(edgestr[[1]][1] %in% nodes)){
      dot <- paste(dot,edgestr[[1]][1], ";")
      nodes <- c(nodes, edgestr[[1]][1])
    }
    if (!(edgestr[[1]][2] %in% nodes)){
      dot <- paste(dot,edgestr[[1]][2], ";")
      nodes <- c(nodes, edgestr[[1]][2])
    }
  }
  for(i in 1:length(vec)){
    edgestr <- strsplit(vec[[i]], ' ')
    mychar <- edgestr[[1]][1]
    newchar <- edgestr[[1]][2]
    dot <- paste(dot, mychar, "->", newchar, ";")
  }
  dot <- paste(dot,"}")

  return(dot)
}


#' sort a matrix that represents a tree
#' @param mat a binary matrix
#' @return a matrix where the columns are sorted by descending number of 1's an the rows
#' and rows are sorted ascending order of total 1's
.SortMat <- function(mat){
  colOrder <- sort(-colSums(mat))

  rowOrder <- sort(rowSums(mat))
  return(mat[names(rowOrder), names(colOrder)])

}

#' Creates a binary matrix tree with clones as rows and mutations as columns from a vector of string edge lists
#' @param vecs a vector of strings with a space separator where the left string represents the parent and the
#' right string represents the child
#' @return a binary matrix with clones as rows and mutations as columns
.CreateTree <- function(vecs){


  parent_list <- .CreateParentList(vecs)
  numMuts <- length(parent_list)
  treeMat <-  matrix(rep(0, numMuts**2),
                     nrow=numMuts,
                     ncol=numMuts,
                     dimnames=list(clones=as.character(0:(numMuts-1)),
                                   mutations=names(parent_list))
  )

  root_mut <- unlist(lapply(lapply(parent_list, is.na), any))
  root_node <- names(root_mut)[root_mut]
  treeMat[,root_node] <- 1

  clones <- list()
  edge <- 1
  clones[[edge]] <- c(root_node)

  for(v in names(parent_list)){
    if(v == root_node){
      next
    }
    edge <- edge + 1
    new_clone <- v

    treeMat[edge,v] <-1
    for(p in parent_list[[v]]){
      new_clone <- c(new_clone,p)
      treeMat[edge,p] <-1
    }

    clones[[edge]] <- new_clone
  }

  treeMat <- .SortMat(treeMat)




  rownames(treeMat) <-apply(treeMat, MARGIN=1, .NameClones)
  return(treeMat)
}




.NameClones <- function(clone, sep=" "){


    clone <- paste(names(clone[clone==1]), collapse = sep)

     return(clone)
}
