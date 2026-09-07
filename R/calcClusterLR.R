#Helpfunction to evaluate cluster evidence
calcClusterLR = function(likEvidPerMarker,priorGlobalPerMarker, clusterMembers ) {
  markers = names(priorGlobalPerMarker)
  
  LRperMarker = setNames(rep(NA_real_,length(markers)),markers)
  for(marker in markers) {
    likVecProd = 1 #traverse each cluster member
    for(ridx in seq_len(nrow(clusterMembers))) {
      #      ridx = 1
      clusMember = clusterMembers[ridx,]
      sample = clusMember[1]
      ucontr = clusMember[2]
      likVec = likEvidPerMarker[[marker]][[sample]][,ucontr]
      likVecProd = likVecProd*likVec
    }
    LRperMarker[marker] = sum(likVecProd*priorGlobalPerMarker[[marker]])
  }      
  return(list(log10LR=sum(log10(LRperMarker)),LRperMarker=LRperMarker))
}  