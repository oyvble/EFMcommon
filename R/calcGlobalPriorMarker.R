#helpfunction to calculate the prior taking frequencies etc (wrapper)
calcGlobalPriorMarker = function(globalFreqs, allelesUse, Fst = 0, refList=NULL, Qallele = "99") {
  
  #Calculate global prior
  globalFreqs = globalFreqs[setdiff(allelesUse,Qallele)]
  if(Qallele%in%allelesUse) {
    globalFreqs = c(globalFreqs,1-sum(globalFreqs))
    names(globalFreqs)[length(globalFreqs)] = Qallele
  }
  
  #Obtain typed alleles 
  alleles = names(globalFreqs)
  nTypedAlleles = setNames(rep(0,length(globalFreqs)),alleles)
  for(ref in names(refList)) {
    #      ref = names(refList)[1]
    refAlleles = unlist(refList[[ref]]) #obtain alleles
    for(a in refAlleles) {
      #        a=refAlleles[2]
      indAllele <- a==alleles
      if(any(indAllele)) {
        nTypedAlleles[indAllele] = nTypedAlleles[indAllele] + 1 #add one
      }  else if(Qallele%in%alleles) {
        nTypedAlleles[Qallele] = nTypedAlleles[Qallele] + 1 #add one
      }
    }
  }
  priorGlobal = calcGenoProb(globalFreqs,nU=1,fst=Fst,nTyped=nTypedAlleles)
  pGlobal = setNames(priorGlobal$Gprob,getGenoVecName(priorGlobal$G))
  return(pGlobal)
}