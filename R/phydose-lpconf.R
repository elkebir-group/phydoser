# tree <- phydata$trees[[1]]
#
# dff <- phydataDFF[[1]]
#
#
#
# fmat <- phydata$fmatrices[[1]]
#
# depth <- 10000
#
# fmatMinus <- list()
# fmatPlus <- list()
# for(i in 1:length(phydata$trees)){
#   tree <- phydata$trees[[i]]
#   fmat <- phydata$fmatrices[[i]]
#   fmat_list <- genConf(fmat)
#
#   dff <- DF(phydataDFF[[i]], phydata$u_list[[i]])
#
#   fmax <- solveF(tree, fmat_list$fminus, fmat_list$fplus,dff, maxObj=T)
#   fmin <- solveF(tree, fmat_list$fminus, fmat_list$fplus,dff, maxObj=F)
#
#   fmatMinus[[i]] <- fmax
#   fmatPlus[[i]] <- fmin
# }
#
#
#
#
# genConf <- function(fmat, sample=1, conf.level = 0.95, depth=1000){
#   reads <- floor(fmat[sample,]*depth)
#   bconf <- binom.confint(reads, rep(depth, length(reads)), conf.level=conf.level, methods="wilson")
#
#
#   fminus <- bconf$lower
#   names(fminus) <- colnames(fmat)
#   fplus <- bconf$upper
#   names(fplus) <- colnames(fmat)
#
#   return(list(fminus= fminus, fplus=fplus))
# }


#' Find the frequencies to use for the k star confidence interval
#' @param tree a binary tree matrix
#' @param fminus the lower bound of the frequency confidence intervals
#' @param fplus the upper bound of the frequency confidence intervals
#' @param dff the distinguishing feature of the tree
#' @param maxObj a bolean determining if the objective should be maximized for minimized
#' @return a list of lists of featurettes.
#' @description   For each tree in the input list, generateDistFeat finds the corresponding distinguishing feature family which
#' is a set of distinguishing features made up of root to mutation paths (featurettes) that unqiuely distinguishes that tree from
#' all other trees in the input list.
#' @export
#'


solveF <- function(tree, fminus, fplus, dff, maxObj= T){

  dff <- .MatchCloneNames(dff, rownames(tree))
  if( requireNamespace("lpSolve", quietly = TRUE)){
    objval <- -1

    tree_inv <- t(solve(tree))
    for(i in 1:length(dff)){

      const <- .generateConstraints(tree_inv, fminus, fplus, dff[[i]], maxObj)

      obj_constants <- c(rep(0, length(fminus)), 1)

      obj_dir <- ifelse(maxObj, "max", "min")



      opt <- lpSolve::lp(direction = obj_dir,
                         objective.in = obj_constants,
                         const.mat <- const$constraints,
                         const.dir =  const$direction,
                         const.rhs = const$rhs)

      if(opt$objval > objval){
        fmat <- matrix(opt$solution[1:length(fminus)], nrow=1, ncol=length(fminus))
        objval <- opt$objval
      }


    }



    colnames(fmat) <- names(fminus)
    return(fmat)
  }else{
    return(NULL)
  }
}

.generateConstraints <- function(tree_inv, fminus, fplus, dff, maxObj=T){
  tree_inv_ord <- tree_inv[, names(fminus)]
  diagMat <- diag(x=1, length(fminus), length(fminus))
  colnames(diagMat) <- names(fminus)

  #add left hand side of constraints to ensure f is non-negative and less than 1
  constraints <- rbind(diagMat, diagMat)

  #add constraints for the sum condition
  constraints <- rbind(constraints, tree_inv_ord)



  dff_clones <- tree_inv_ord[rownames(tree_inv_ord) %in% trimws(dff),]

  z_col <- c(rep(0, nrow(constraints)), rep(-1, length(dff)))

  constraints <- rbind(constraints, dff_clones)
  constraints <- cbind(constraints, z_col)

  rhs <- c(fminus, fplus, rep(0, nrow(tree_inv_ord)), rep(0, length(dff)))

  dir <- c( rep(">=", nrow(diagMat)), rep("<=", nrow(diagMat)), rep(">=", (nrow(tree_inv_ord))))

  if(maxObj){
    dir_z <- rep(">=", length(dff))
  }else{
    dir_z <- rep("<=", length(dff))
  }

  dir <- c(dir, dir_z)
  return(list(constraints = constraints, rhs = rhs, direction = dir))
}


.MatchCloneNames <- function(df_fam, clone_names){
  clones <- stringr::str_split(clone_names, pattern = " ")
  df_rename <- list()

  for(i in 1:length(df_fam)){
    df <- df_fam[[i]][[1]]

    dfclones <- stringr::str_split(df, " ")
    dfclones <- lapply(dfclones, function(x) x[x != ""])

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


  return(df_rename)
}


# phydose(phydata$trees, fmatrices = fmatMinus)
# phydose(phydata$trees, fmatrices = fmatPlus)
# phydose(phydata$trees, fmatrices = phydata$fmatrices)
