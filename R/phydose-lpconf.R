




#' Find the frequencies to use for the k star confidence interval
#' @param tree a binary tree matrix
#' @param fminus the lower bound of the frequency confidence intervals
#' @param fplus the upper bound of the frequency confidence intervals
#' @param dff the distinguishing feature of the tree
#' @return a list of lists of featurettes.
#' @description   For each tree in the input list, generateDistFeat finds the corresponding distinguishing feature family which
#' is a set of distinguishing features made up of root to mutation paths (featurettes) that unqiuely distinguishes that tree from
#' all other trees in the input list.
#' @export
#'


solveFLower <- function(tree, fminus, fplus, dff){

  dff <- .MatchCloneNames(dff, rownames(tree))

  if( requireNamespace("lpSolve", quietly = TRUE)){
    objval <- -1

    tree_inv <- t(solve(tree))
    for(i in 1:length(dff)){

      const <- .generateConstraintsLower(tree_inv, fminus, fplus, dff[[i]])

      obj_constants <- c(rep(0, length(fminus)), 1)

      obj_dir <- "max"



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


solveFUpper2 <- function(tree, fminus, fplus, dff){

  dff <- .MatchCloneNames(dff, rownames(tree))

  if( requireNamespace("lpSolve", quietly = TRUE)){
    objval <- -1

    tree_inv <- t(solve(tree))
    for(i in 1:length(dff)){

      const <- .generateConstraintsUpper2(tree_inv, fminus, fplus, dff[[i]])

      obj_constants <- c(rep(0, length(fminus)), 1)

      obj_dir <- "min"



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


#' Find the frequencies to use for the k star confidence interval
#' @param tree a binary tree matrix
#' @param fminus the lower bound of the frequency confidence intervals
#' @param fplus the upper bound of the frequency confidence intervals
#' @param dff the distinguishing feature of the tree
#' @return a list of lists of featurettes.
#' @description   For each tree in the input list, generateDistFeat finds the corresponding distinguishing feature family which
#' is a set of distinguishing features made up of root to mutation paths (featurettes) that unqiuely distinguishes that tree from
#' all other trees in the input list.
#' @export
#'

solveFUpper <- function(tree, fminus, fplus, dff){

  dff <- .MatchCloneNames(dff, rownames(tree))

  if( requireNamespace("lpSolve", quietly = TRUE)){
    objval <- -Inf

    tree_inv <- t(solve(tree))
    for(i in 1:length(dff)){

      const <- .generateConstraintsUpper(tree_inv, fminus, fplus, dff[[i]])

      obj_constants <- c(rep(0, length(fminus)), 1, rep(0, length(dff[[i]])))

      obj_dir <- "min"
      bin.vec <- (length(obj_constants) - length(dff[[i]])):length(obj_constants)


      opt <- lpSolve::lp(direction = obj_dir,
                         transpose.constraints = F,
                         objective.in = obj_constants,
                         const.mat <- const$constraints,
                         const.dir =  const$direction,
                         binary.vec <- bin.vec,
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

.generateConstraintsUpper2 <- function(tree_inv, fminus, fplus, dff){
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


  dir_z <- rep(">=", length(dff))


  dir <- c(dir, dir_z)
  return(list(constraints = constraints, rhs = rhs, direction = dir))
}


.generateConstraintsLower <- function(tree_inv, fminus, fplus, dff){
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


  dir_z <- rep(">=", length(dff))


  dir <- c(dir, dir_z)
  return(list(constraints = constraints, rhs = rhs, direction = dir))
}


.generateConstraintsUpper <- function(tree_inv, fminus, fplus, dff){
  tree_inv_ord <- tree_inv[, names(fminus)]
  diagMat <- diag(x=1, length(fminus), length(fminus))
  colnames(diagMat) <- names(fminus)
  bigM <- 10000

  #add left hand side of constraints to ensure f is non-negative and less than 1
  constraints <- rbind(diagMat, diagMat)

  #add constraints for the sum condition
  constraints <- rbind(constraints, tree_inv_ord)



  dff_clones <- tree_inv_ord[rownames(tree_inv_ord) %in% trimws(dff),]

  z_col <- c(rep(0, nrow(constraints)), rep(-1, length(dff)))
  binary_vars <- length(dff)

  for(i in 1:(binary_vars+1)){
    constraints <- cbind(constraints,rep(0, nrow(constraints)))
  }
  if(binary_vars == 1){
    uconstraints <- c(dff_clones, -1, -bigM)
  }else{
    diagBigM <-  diag(x=-1*bigM, binary_vars, binary_vars)
    uconstraints <- cbind(dff_clones, rep(-1, binary_vars))
    uconstraints <- cbind(uconstraints, diagBigM)

  }




  constraints <- rbind(constraints, uconstraints)
  constraints <-  rbind(constraints, c(rep(0, ncol(constraints)-binary_vars), rep(1, binary_vars)))
  const.mat <- t(constraints)
  const.count <- ncol(const.mat)

  rhs <- c(fminus, fplus, rep(0, nrow(tree_inv_ord)), rep(0, length(dff)))
  rhs <- c(rhs, binary_vars-1)

  dir <- c( rep(">=", nrow(diagMat)), rep("<=", nrow(diagMat)), rep(">=", (nrow(tree_inv_ord))))


  dir_z <- rep("<=", length(dff))


  dir <- c(dir, dir_z, "=")

  return(list(constraints = const.mat, rhs = rhs, direction = dir))
}



#A simple example for min/min: if you have: min z where z = min(x1,x2,x3),
#then write: min z where z >= x1-M b1, z >= x2-M b2, z >= x3-M b3, b1+b2+b3=2 and b1,b2,b3 binary.

.MatchCloneNames <- function(df_fam, clone_names){
  clones <- stringr::str_split(clone_names, pattern = " ")
  df_rename <- list()

  for(i in 1:length(df_fam)){
    df <- df_fam[[i]]

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
