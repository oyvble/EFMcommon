#Pairwise common-unknown donor LR:
calcCommonLRpair = function(sampleNames,markersAll,uRangePerSample,likEvidPerMarker,priorGlobalPerMarker, LRthresh=6) {
  
  pairResults = list()
  for(s1 in 1:(length(sampleNames)-1)) {
    sample1 = sampleNames[s1]
    pairResults[[sample1]] = list()
    for(s2 in (s1+1):length(sampleNames)) {
      #      s1 = 2;s2 = 3
      sample2 = sampleNames[s2]
      pairResults[[sample1]][[sample2]] = list()
      Urng1 = paste0(uRangePerSample[[sample1]])
      Urng2 = paste0(uRangePerSample[[sample2]])
      
      commonLR = matrix(NA_real_,nrow=length(Urng1),ncol=length(Urng2),dimnames = list(Urng1,Urng2))  
      commonLR_marker = array(NA_real_,dim = c(length(Urng1),length(Urng2),length(markersAll)),dimnames = list(Urng1,Urng2,markersAll))  
      for(u1 in Urng1) {
        for(u2 in Urng2) {
          #  u1=Urng1[1];u2=Urng2[2]
          LR_marker = setNames(rep(NA,length(markersAll)),markersAll)
          for(marker in markersAll) {
            lik1 = likEvidPerMarker[[marker]][[sample1]][,u1]
            lik2 = likEvidPerMarker[[marker]][[sample2]][,u2]
            LR_marker[marker] = sum(lik1*lik2*priorGlobalPerMarker[[marker]])
          }
          commonLR_marker[u1,u2,] = LR_marker
          commonLR[u1,u2] = sum(log10(LR_marker)) #on log10 scale
        }
      }
      pairResults[[sample1]][[sample2]] = list(log10LR=commonLR,log10LR_marker=commonLR_marker) 
    }
  }
  
  ####################################
  #Obtain clusters based on threshold#
  ####################################
  clusterList = list()  
  for(sample1 in names(pairResults)) {
    for(sample2 in names(pairResults[[sample1]])) {
      #      sample1 = sampleNames[1];sample2 = sampleNames[3]
      commonLR = pairResults[[sample1]][[sample2]]$log10LR
      commonIndMat = which(commonLR>=LRthresh,arr.ind = TRUE)
      if(nrow(commonIndMat)==0) next
      for(row in seq_len(nrow(commonIndMat))) {   #Add each to separate clusters:
        commonIdx = commonIndMat[row,]
        contr_Sample1 = rownames(commonLR)[commonIdx[1]]
        contr_Sample2 = colnames(commonLR)[commonIdx[2]]
        contrs = c(contr_Sample1,contr_Sample2)
        clusTab = cbind(
          Evid=c(sample1,sample2),
          Contr=contrs
        )
        clusterList[[length(clusterList)+1]] = list(
          members=clusTab,
          log10LR=commonLR[commonIdx[1], commonIdx[2]]
        )
      }
    }
  }
  return(list(pairResults=pairResults,clusterList=clusterList))
}