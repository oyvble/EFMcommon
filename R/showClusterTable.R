showClusterTable = function(clusterList, maxWidth = 150, fontSize = 12, palette = "Pastel 1", showN = TRUE) {
  
  if(!is.list(clusterList) || length(clusterList)==0) {
    stop("clusterList must be a non-empty list")
  }
  
  #########################
  # Obtain cluster members
  #########################
  
  memberList = lapply(clusterList, function(x) {
    if(is.null(dim(x)) || ncol(x)<2) {
      stop("Each cluster must be a matrix/table with evidence and contributor")
    }
    paste0(x[,1], ":C", x[,2])
  })
  
  nMembers = lengths(memberList)
  maxMembers = max(nMembers)
  
  #########################
  # Make rectangular table
  #########################
  
  memberList = lapply(memberList, function(x) {
    c(x, rep(NA_character_, maxMembers-length(x)))
  })
  
  tab = as.data.frame(
    memberList,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  #Column names
  clusterNames = paste0("Cluster ", seq_along(clusterList))
  
  if(showN) {
    clusterNames = paste0(clusterNames," (n=", nMembers, ")")
  }
  names(tab) = clusterNames
  
  #########################
  # Create colours
  #########################
  
  clusterCols = grDevices::hcl.colors(
    length(clusterList),palette = palette)
  
  #Use paler version for table body
  bodyCols = grDevices::adjustcolor(
    clusterCols, alpha.f = 0.30)
  
  #########################
  # Create gt table
  #########################
  
  #Create column-width specification
  widthSpec = rlang::new_formula(
    lhs = rlang::expr(gt::everything()),
    rhs = gt::px(maxWidth)
  )
  
  gtTab = gt::gt(tab) |>
    gt::sub_missing(missing_text = "") |>
    gt::cols_align(
      align = "left",
      columns = gt::everything()
    ) |>
    gt::cols_width(.list = list(widthSpec)) |>
    gt::tab_options(
      data_row.padding = gt::px(3),
      column_labels.padding = gt::px(5),
      table.font.size = gt::px(fontSize)
    )
  
  #########################
  # Colour each cluster
  #########################
  
  for(i in seq_along(clusterList)) {
    colName = clusterNames[i]
    
    #Body
    gtTab = gt::tab_style(
      gtTab,
      style = gt::cell_fill(
        color = bodyCols[i]
      ),
      locations = gt::cells_body(
        columns = tidyselect::all_of(colName)
      )
    )
    
    #Header
    gtTab = gt::tab_style(
      gtTab,
      style = list(
        gt::cell_fill(color = clusterCols[i]),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_column_labels(
        columns = tidyselect::all_of(colName)
      )
    )
  }

  return(gtTab)
}