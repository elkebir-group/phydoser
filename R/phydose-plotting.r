

# fooPlot <- function(){
#   test <- "digraph graphname {a -> b -> c;  b -> d;}"
#   return(DiagrammeR::grViz(test))
# }

#' Convert a binary tree to dot format
#' @param tree a binary matrix
#' @return a string representing the string in dot format
#' @export
CreateDOT <- function(tree){
  nodes <- colnames(tree)
  tree_edges <- .GetEdges(tree)
  tree_edges <- tree_edges %>% dplyr::mutate(edge = sprintf("%s -> %s", from, to))

  node_strings <- paste(nodes, collapse = " ; ")
  all_edges <-  paste(tree_edges$edge, collapse = " ; ")
  digraph <- paste0("digraph T{ ", node_strings, " ; ", all_edges, " ; }")
  return(digraph)
}
CreateTreeGraph <- function(tree, distFeat=NULL, featurette=NULL, colors=NULL){

  if(is.null(colors)){
    colors = c( "#e41a1c", "#377eb8" ,"#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#f781bf", "#999999" )
    if(length(distFeat) > length(colors)){
      black <- rep("#000000", length(distFeat)-length(colors))
      colors < c(colors, black)
    }
  }

  nodes <- colnames(tree)
  tree_edges <- .GetEdges(tree) %>% dplyr::mutate(color=NA)
  node.df <- data.frame(node = colnames(tree))
  if(is.null(colors)){
    colors = c( "#e41a1c", "#377eb8" ,"#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#f781bf", "#999999" )
    if(length(distFeat) > length(colors)){
      black <- rep("#000000", length(distFeat)-length(colors))
      colors < c(colors, black)
    }
  }



  if(!is.null(distFeat)){
    color_index <- 1
    for(df in distFeat){
      color <- colors[color_index]
      featurette <- unlist(stringr::str_split(df, " "))
      node_strings = ifelse(nodes %in% featurette, paste0(nodes, " [color = ", color, "]"), nodes )
      nodes <- nodes[!(nodes %in% featurette)]

      tree_edges <- tree_edges %>%
        dplyr::mutate(color = ifelse( (from %in% featurette && to %in% featurette), color, ""))
    }


  } else if(!is.null(featurette)){
    featurette <- unlist(stringr::str_split(featurette, " "))
    node_strings = ifelse(nodes %in% featurette, paste0(nodes, " [color = ", colors[1], "]"), nodes )

    tree_edges <- tree_edges %>%
      dplyr::mutate(color = ifelse(from %in% featurette || to %in% featurette, colors[1], ""))

  }else{
    node_strings <- nodes
    tree_edges$color <- "#000000"
  }


  node_strings <- paste(node_strings, collapse = " ; ")

  tree_edges <- tree_edges %>% dplyr::mutate(color = ifelse(is.na(color), "#000000",  color),
                                             edge = sprintf("%s %s %s", from, " -> ", to))
  edge_strings_colors <- tree_edges$color
  all_edges  <- sprintf("%s [color = %s%s%s]", tree_edges$edge, '"', edge_strings_colors, '"')
  all_edges <- paste(all_edges, collapse =" ;")
  digraph <- paste0("digraph T{ ", node_strings, " ; ", all_edges, " ; }")
  # digraph <- stringr::str_replace_all(digraph, "\", "")
  DiagrammeR::grViz(digraph)
  return(digraph)
  # for(i in 1:length(nodes)){
  #   #edL[[i]] <- list(edges=5-i, weights=runif(1))
  #   print(i)
  #   nodeEdges <- dplyr::filter(tree_edges, from== nodes[i])
  #   print(nodeEdges)
  #   edL[[i]] <- list(edges= nodeEdges$to)
  # }


}

#' Plot a tree in the set
#' @param tree an input binary tree
#' @param distFeat a list of featurettes
#' @param featurette a vector of mutations
#' @param colors a vector of colors the length of distFeats or of length 1
#' @param u a named vector of clonal prevalence rates
#' @param f a named vector of frequencies
#' @return a DiagrammeR graph
#' @export
plotTree <- function(tree, distfeat = NULL, featurette = NULL, colors = NULL, u=NULL, f=NULL){

  if(requireNamespace("DiagrammeR", quietly = TRUE)){

    numNodes <- ncol(tree)
    nodes <- colnames(tree)



    if(is.null(colors)){
      colors = c( "#e41a1c", "#377eb8" ,"#4daf4a", "#984ea3", "#
                  ff7f00", "#ffff33", "#a65628","#f781bf", "#999999" )
      if(length(distfeat) > length(colors)){
        black <- rep("#000000", length(distfeat)-length(colors))
        colors < c(colors, black)
      }
    }


    default_color <- "#000000"
    tree_edges_df <- .GetEdges(tree) %>% dplyr::mutate(color = default_color)


    if(!is.null(distfeat)){
      color_index <- 1
      node_colors <- rep(default_color, length(nodes))


      if(!is.list(distfeat)){
        distfeat <- list(distfeat)
      }
      for(df in distfeat){
        color <- colors[color_index]
        featurette <- unlist(stringr::str_split(df, " "))
        default_nodes <- which(node_colors == default_color)
        nodes_in_featurette <- which(nodes %in% featurette)
        node_colors[intersect(default_nodes, nodes_in_featurette)] <- color
        #node_colors = ifelse(nodes %in% featurette && node_colors==default_color, color, default_color )

        edges_in_feature <- which(tree_edges_df$to %in% featurette)
        default_edges <- which(tree_edges_df$color == default_color)
        tree_edges_df$color[intersect(default_edges, edges_in_feature)] <- color
        color_index <- color_index + 1
      }


    } else if(!is.null(featurette)){
      featurette <- unlist(stringr::str_split(featurette, " "))
      node_colors <- rep(default_color, numNodes)
      node_colors[nodes %in% featurette] <- colors[1]


      tree_edges_df <- tree_edges_df %>%
        dplyr::mutate(color = ifelse(to %in% featurette, colors[1], default_color))

    }else{
      node_colors <- rep(default_color, numNodes)
      tree_edges_df$color <- default_color
    }
    node_labels <- ifelse(stringr::str_length(nodes) > 15, stringr::str_sub(nodes, 1, 15), nodes)

    fontsize <- ifelse(stringr::str_length(node_labels)  > 7, 6, 12)
    nodes_df <- DiagrammeR::create_node_df(n=numNodes,label = nodes, value = nodes,
                                           shape="ellipse", color=node_colors, fixedsize=F,fontsize=fontsize)


    from <- nodes_df$id[match(tree_edges_df$from, nodes_df$label)]
    to <- nodes_df$id[match(tree_edges_df$to, nodes_df$label)]






    if(!is.null(f)){
      if(is.matrix(f)){
        f <- f[1,]
      }
      f.df <- data.frame(label = names(f), f, stringsAsFactors = F)
      nodes_df <- nodes_df %>% dplyr::left_join(f.df, by = "label")


    }

    if(!is.null(u)){
      if(is.matrix(u)){
        u <- u[1,]
      }

      clones <- stringr::str_split(names(u), " ")
      u.df <- data.frame(label =unlist(lapply(clones, .CloneEnd)), u = u, stringsAsFactors = F)
      nodes_df <- dplyr::left_join(nodes_df,u.df, by="label")


    }






    if(!is.null(u) && is.null(f)){


      nodes_df <- dplyr::mutate(nodes_df, label = sprintf("%s\nu:%0.2f", label,u ))

    }else if(!is.null(f) && is.null(u)){
      nodes_df <- dplyr::mutate(nodes_df, label = sprintf("%s\nf:%0.2f", label,f ))

    }else if(!is.null(f) && !is.null(u)){
      nodes_df <- dplyr::mutate(nodes_df, label = sprintf("%s\nf:%0.2f u:%0.2f", label,f,u ))
      nodes_df$fontsize <- 5

    }else{
      nodes_df$label <- nodes_df$label
    }






    edges <- DiagrammeR::create_edge_df(from =from, to=to, color = tree_edges_df$color )




      treeGraph <- DiagrammeR::create_graph(nodes_df, edges, directed=T)
      #DiagrammeR::render_graph(treeGraph, layout="tree")

      return(treeGraph)
  }

}

.GetEdges <- function(tree){
  from <- character()
  to <- character()
  tree <- .SortMat(tree)
  clones <- apply(tree, MARGIN = 1, .NameClones)
  clone_list <- stringr::str_split(clones, " ")
  for(clone in clone_list){
    if(length(clone) > 1){
      start <- 1
      end <- 2
       while(end <= length(clone)){
         from <- c(from, clone[start])
         to <- c(to, clone[end])
         start <- start + 1
         end <- end + 1
       }
    }
  }

  edges <- data.frame(from = from, to = to, stringsAsFactors = F) %>% dplyr::distinct()



  return(edges)

}

.CloneEnd <- function(myvec){
  return(myvec[length(myvec)])
}

