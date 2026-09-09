plotCommonCluster = function(clusterList) {
  
  if(length(clusterList)==0) stop("clusterList is empty")
  
  #########################
  # Create member table
  #########################
  
  memberTab = NULL
  for(clusIdx in seq_along(clusterList)) {
    clus = clusterList[[clusIdx]]
    
    tmp = data.frame(
      cluster = clusIdx,
      evidence = as.character(clus[,1]),
      contributor = as.character(clus[,2])
    )
    #Unique node name
    tmp$node = paste0(tmp$evidence, ":C", tmp$contributor)
    memberTab = rbind(memberTab,tmp)
  }
  
  #########################
  # Create cluster nodes
  #########################
  clusterNodes = paste0("Clus",seq_along(clusterList))
  
  #Edges simply indicate cluster membership
  edgeTab = data.frame(
    from = clusterNodes[memberTab$cluster],
    to = memberTab$node
  )
  
  #########################
  # Create graph
  #########################
  g = igraph::graph_from_data_frame(edgeTab,directed=FALSE)
  
  #########################
  # Vertex properties
  #########################
  isCluster = igraph::V(g)$name %in% clusterNodes
  
  #Determine cluster number for every node
  vertexCluster = integer(igraph::vcount(g))
  
  #Cluster hubs
  vertexCluster[isCluster] = match(igraph::V(g)$name[isCluster],clusterNodes)
  
  #Evidence/contributor nodes
  vertexCluster[!isCluster] = memberTab$cluster[match(igraph::V(g)$name[!isCluster],memberTab$node)]
  
  #One colour per inferred cluster
  clusterCols = grDevices::hcl.colors(length(clusterList), palette="Dark 3")
  vertexCols = clusterCols[vertexCluster]
  
  #########################
  # Layout
  #########################
  layout = igraph::layout_with_kk(g)
  
  #########################
  # Plot
  #########################
  
  igraph::plot.igraph(
    g,layout = layout,
    vertex.color = vertexCols,
    vertex.frame.color = NA,
    
    #Cluster centres are visually distinct
    vertex.shape = ifelse(isCluster,"square","circle"),
    vertex.size = ifelse(isCluster,10,10),
    vertex.label = igraph::V(g)$name,
    vertex.label.cex = 0.8,
    edge.color = clusterCols[vertexCluster[match(edgeTab$from,igraph::V(g)$name)]],
    edge.width = 2 
  )
  invisible(g)
}