#Helpfunction to evaluate posterior prob DC (Deconvolution) 
calcClusterDC = function(likEvidPerMarker,priorGlobalPerMarker, clusterList ) {
  markers = names(priorGlobalPerMarker)
  
  DClist = list()
  for(clusIdx in seq_along(clusterList)) {
#    clusIdx = 1
    DClist[[clusIdx]] = list()
    clusterMembers = clusterList[[clusIdx]]
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
      postProb = likVecProd*priorGlobalPerMarker[[marker]]
      DClist[[clusIdx]][[marker]] = postProb/sum(postProb) #ensure normalization
    }      
  }
  
  return(DClist)
}  