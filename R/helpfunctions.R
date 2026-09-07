
#helpfunction to indicate which values are similar (below threshold)
getDuplicatedValues = function(vec,thresh) {
  n = length(vec)
  if(n<2) return(FALSE)
  areUnequal = rep(TRUE,n) #default is that all are unequal
  for(uidx in 1:(n-1)) {
    # uidx=1
    if(!areUnequal[uidx]) next #skip if it was flagged as equal
    mxU1 = vec[uidx] #obtain Mx
    compareRange = (uidx+1):n #indices that are compared
    compareInds = compareRange[areUnequal[compareRange]]
    mxUother = vec[compareInds]
    hasSimilarMx = abs(mxU1-mxUother)<thresh
    hasSimilarMxCompared = compareInds[hasSimilarMx]
    areUnequal[hasSimilarMxCompared] = FALSE #set as false
  }
  return(!areUnequal)
}


#Helpfunction to get prior probs
getPriorProbs = function(freq,fst,knownAlleles=NULL) {
  priorList = calcGenoProb(freq = freq, nU = 1, fst = fst, nTyped = knownAlleles)
  priorProbs = setNames(priorList$Gprob,paste0(priorList$G[,1],"/",priorList$G[,2]))
  return(priorProbs)  
}


#Helpfunction to obtain an updated genotype name vector where missing alleles are accounted for (imputed with Qallele)
getMissingGenos = function(genoMat,alleles,Qallele="99") {
  genoMatUpdated = genoMat #copy
  genoMatUpdated[ !genoMat[,1]%in% alleles ,1] = Qallele #insert Q-allele for missing
  genoMatUpdated[ !genoMat[,2]%in% alleles ,2] = Qallele #insert Q-allele for missing
  
  #Need to swap with respect to sorting (important!)
  swap = genoMatUpdated[,2]<genoMatUpdated[,1]
  if(any(swap)) genoMatUpdated[swap,] = genoMatUpdated[swap,2:1]
  return(paste0(genoMatUpdated[,1],"/",genoMatUpdated[,2]))
}

#Helpfunction to obtain vector name from a matrix
getGenoVecName = function(genoMat) {
  paste0(genoMat[,1],"/",genoMat[,2])
}

#Convert back from vector names to matrix
getGenosAsMatrix = function(X) {
  return(t(matrix(unlist(strsplit(X,"/")),nrow=2)))
}          