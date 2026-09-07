#' @title calcCommonCluster
#' @author Oyvind Bleka
#' @description Provides clustering of unknown contributors in samples
#' @details
#' Currently supported EuroForMix Family packages: (euroformix, EFMmps)
#'
#' Implements the FAC-clustering approach (greedy-search) from 
#' Bhembe et al (2026): "Salience, legitimacy, and credibility of two end-to-end forensic single cell probabilistic systems"
#'
#' Implementation structure
#' 1. Validate input
#' 2. Extract sample/marker/component information
#' 3. Harmonize genotype spaces and construct global priors
#' 4. Convert deconvolution posteriors to likelihood weights
#' 5. Initialize singleton clusters
#' 6. Calculate initial pairwise merge likelihoods
#' 7. Greedily merge, reusing unaffected calculations
#' 8. Return clusters
#' 
#' @param fitList List with fitted MLE-objects returned from a EuroForMix family program: List with sample elements which are again list elements: fitList[[sample]]=list(MLE,DC)
#' @param mergeThresh Threshold used for terminating the clustering: Note that this is relative LR
#' @param Fst The Fst used for the prior can be modified (over-writes existing)
#' @param mxThresh Threshold of whether multiple unknown contributors in an evididence should be merged
#' @param Qallele Name of compound drop-out alleles
#' @param returnPrepared Whether returning the pre-processed objects

# mergeThresh = 1; Fst=NULL; mxThresh=0.01; Qallele="99"; returnOnlyPrepared=FALSE
calcCommonCluster = function(fitList, mergeThresh = 0, Fst=NULL, mxThresh=0.01, Qallele="99", returnPrepared=FALSE) {
  sampleNames = names(fitList)
  if(is.null(sampleNames)) stop("Sample name must be given to the fitList")
  if(length(sampleNames)<=1) stop("Must have at least two samples")
  
  #########################################################
  #PART 1: THIS FIRST PART IS PREPARING LIKELIHOOD WEIGHTS# 
  #########################################################
  
  #Obtaining marker list per fitted sample
  markerList = lapply(fitList,function(x) x$MLE$prepareC$markerNames)
 
  #Traverse each sample and markers to obtain allele outcomes and prior info 
  alleleOutcomePerMarkerSample <- priorListPerMarker <- uRangePerSample <- list()
  samplesIgnore = character() #indicate which samples to ignore (because of no unknowns)
  for(sample in sampleNames) {
#    sample=sampleNames[3]
    c = fitList[[sample]]$MLE$prepareC
    
    #CHECK NUMBER OF UNKNOWNS DEFINED IN HYPOTHESIS
    NOC = c$nC #obtain number of defined unknowns
    sample_nKnowns = max(c$nKnowns)
    sample_nUnknowns = NOC - sample_nKnowns #obtain number of unknowns
    if(sample_nUnknowns==0) {
      samplesIgnore = c(samplesIgnore,sample) #ignore sample
      next #skip sample if no unknowns
    }
    
    #CHECK THAT MODEL WAS PROPERLY CONVERGED (proper Mx values)
    uRange = sample_nKnowns + seq_len(sample_nUnknowns)
    mxHatUnknowns =  fitList[[sample]]$MLE$fit$thetahat2[uRange]
    if(any(is.na(mxHatUnknowns))) {
      samplesIgnore = c(samplesIgnore,sample) #ignore sample
      next #SKIP SAMPLE IF NOT PROPERLY FITTED
    } 
    
    #Components that have very similar Mx are merged (skipping duplicates)
    if(sample_nUnknowns>=2) {
      isDuplicatedMxVal = getDuplicatedValues(mxHatUnknowns,mxThresh)
      uRange = uRange[!isDuplicatedMxVal] #update the unknown range
    }
    
    markers_Sample = markerList[[sample]]
    uRangePerSample[[sample]] = uRange
    for(marker in markers_Sample) {
#     marker = markersAll[1];sample = sampleNames[1]
        
      #Obtaining allele details (used for prior)
      mind = which(markerList[[sample]]==marker) #obtain marker index
      mrng = c$startIndMarker_nAlleles[mind] + seq_len(c$nAlleles[mind])
      sample_alleles = c$alleleNames[mrng]
      sample_freqs = setNames(c$freqs[mrng],sample_alleles)

      #calculate the prior
      genoProbList = calcGenoProb(sample_freqs,nU=1,fst=c$fst[mind],nTyped=c$maTyped[mrng])
      genoProbs = setNames(genoProbList$Gprob,getGenoVecName(genoProbList$G))
      
      #Insert info:
      isFirstElem = is.null(priorListPerMarker[[marker]]) #indicate if this is first elem
      if(isFirstElem) {
        priorListPerMarker[[marker]] = list()
        alleleOutcomePerMarkerSample[[marker]] = list()
      }
      
      priorListPerMarker[[marker]][[sample]] = genoProbs #store prior
      alleleOutcomePerMarkerSample[[marker]][[sample]] = sample_alleles#store allele names
    }    
  }
  sampleNames = setdiff(sampleNames,samplesIgnore) #update samples to consider
  markersAll = names(priorListPerMarker)
  
  #Structuring data with respect to common alleles (also adapt global prior here)
  likEvidPerMarker <- genosAllPerMarker <- priorGlobalPerMarker <- list() 
  for(marker in markersAll) {
    allelesList = alleleOutcomePerMarkerSample[[marker]]
    allelesAll = sort(unique(unlist(allelesList)))
    nAlleles = length(allelesAll)
       
    #Obtain the full genotype oucome
    genosAll = numeric()
    for(i in 1:nAlleles) genosAll = rbind(genosAll, cbind( allelesAll[rep(i,nAlleles - i + 1)], allelesAll[i:nAlleles] ))
    genosAllPerMarker[[marker]] = genosAll
    genoNameVec = getGenoVecName(genosAll)

    #Init likelihood evidence matrix and obtain info for calculating the global prior
    likEvidPerMarker[[marker]] = list()
    globalFreqs <- NULL
    refList = list()
    for(sample in sampleNames) {
      likEvidPerMarker[[marker]][[sample]] = matrix(1,nrow=length(genoNameVec),ncol=length(uRangePerSample[[sample]]),dimnames=list(genoNameVec,uRangePerSample[[sample]]))
      
      #Calculate the global prior (for full genotype outcome)
      freqs = fitList[[sample]]$MLE$model$popFreq[[marker]]
      freqs = freqs[names(freqs)%in%allelesAll]
      globalFreqs = c(globalFreqs,freqs[!names(freqs)%in%names(globalFreqs)])
      
      refData = fitList[[sample]]$MLE$model$refData
      if(!is.null(refData)) {
        refDataMarker = refData[[marker]] #check first this format
        if(is.null(refDataMarker)) refDataMarker = refData[[sample]][[marker]] #then this format
        
        if(!is.null(refDataMarker)) {
          newRefs = !names(refDataMarker)%in%names(refList) #obtain name of new refs
          refList = c(refList,refDataMarker[newRefs]) #add new refs to list if found
        }
      }
    }
    #Obtain Fst from first sample only:
    FstMarker = fitList[[sampleNames[1]]]$MLE$model$fst[1]
    if(!is.null(Fst)) FstMarker = Fst #use provided
    pGlobal = calcGlobalPriorMarker(globalFreqs, allelesAll,FstMarker,refList,Qallele)
    priorGlobalPerMarker[[marker]] = pGlobal[genoNameVec] #insert
  }

  #Calculating the likelihood weights by combining posterior and prior information
  for(sampleName in sampleNames) { #traverse sample wise first
#    sampleName=sampleNames[1]
    fitListSample = fitList[[sampleName]] #get fitted list for specific sample
    DCobj = fitListSample$DC
    if(is.null(DCobj)) {
      DCobj = deconvolve(fitListSample$MLE,alpha = 1) #ensure to list all outcomes
    }
    DCtab = DCobj$table3 #obtain full DC-tab
    if(is.null(DCtab)) stop("No deconvolution information obtained, cannot progress further!")
    sample_Urange = uRangePerSample[[sampleName]]  #Obtain unknown range for sample
    
    #traversing per-marker
    sample_markers = markerList[[sampleName]]
    for(marker in sample_markers) {
      allelesSample = alleleOutcomePerMarkerSample[[marker]][[sampleName]]
      if(is.null(allelesSample)) next #skip if not found
      DCtab_marker = DCtab[DCtab[,2]==marker,,drop=FALSE ] #filter relevant
      
      toGenoNameMat = genosAllPerMarker[[marker]] #get matrix
      toGenoNameMissing = getMissingGenos(toGenoNameMat,allelesSample,Qallele) #get updated genos
      
      #Obtaining genotype names (vectorized format)      
      priorProbs = priorListPerMarker[[marker]][[sampleName]]
      fromGenoPrior = names(priorProbs)
      insInd = match(toGenoNameMissing,fromGenoPrior) #obtain where to insert values
      if(any(is.na(insInd))) stop("Could not do genotype mapping!")
    #  all(fromGenoPrior[insInd]==toGenoNameMissing) #check
      #Obtain the posterior probabilites and insert to weights (After adjusting for prior)
      for(Cidx in sample_Urange) {
#        Cidx = sample_Urange[1]
        CidxTxt = paste0(Cidx)
        DCtab_U = DCtab_marker[paste0("C",CidxTxt)==DCtab_marker[,1], ,drop=FALSE]
        fromGenoPost = as.character(DCtab_U[,3])
        postProbs = setNames(as.numeric(DCtab_U[,4]),fromGenoPost) #obtain posterior probs
        
        #Calculate the likelihood weight and insert to universal matrix:
        likWeight = postProbs[fromGenoPrior]
        likWeight[is.na(likWeight)] = 0 #unlikely combination
        likWeight = likWeight/priorProbs
        #likWeight = likWeight/sum(likWeight)
        if(!all(is.finite(likWeight))) stop("Likelihood weights contained non-numerical values!")
        likEvidPerMarker[[marker]][[sampleName]][,CidxTxt] = likWeight[insInd] #insert to array
      }
    }  #end per marker
  } #end per sample
  
  #DIAGNOSTIC CHECK: sum_g p(g)*Lik(g)=1
  if(0) {
    for(sampleName in sampleNames) { #traverse sample wise first
      for(marker in markersAll) {
        priorProb = priorGlobalPerMarker[[marker]]
        print(t(likEvidPerMarker[[marker]][[sampleName]])%*%priorProb)
      }
    }
  }
  

  
  ##############################################
  #Part 2: Calculate multi-evidence per cluster#
  ##############################################
  
  #Obtain singletons of components
  clusterList = list()
  for(sample in names(uRangePerSample)) {
    for(contr in uRangePerSample[[sample]]) {
      clusterList[[length(clusterList)+1]] = rbind(c(sample,contr)) #must be a matrix
    }
  }

  #Init LR per separate clusters  
  nClus = length(clusterList)
  if(nClus==1) return(clusterList) #return if no multiple clusters
  separate_clusterLR = rep(NA,nClus)
  for(clusIdx in seq_len(nClus)) {
    clusObj = calcClusterLR(likEvidPerMarker,priorGlobalPerMarker,clusterList[[clusIdx]])
    separate_clusterLR[clusIdx] = clusObj$log10LR
  }
  
  #Internal helpfunction to get merged LR (two clusters)
  getMergedLR = function(clusList, idx1,idx2) {
    mergedClus = rbind(clusList[[idx1]],clusList[[idx2]]) #merge cluster
    
    #Avoid that two merged clusters contain same evidence names
    hasRepeatedSample = any(duplicated(mergedClus[,1])) #check if any has shared evidence name
    lr = NA
    if(!hasRepeatedSample)  { #if not shared evidence name we continue
      mergedLR = calcClusterLR(likEvidPerMarker,priorGlobalPerMarker,mergedClus)
      lr = mergedLR$log10LR #obtain LR
    }
    return(lr)
  }
  
  #Second: Consider the pairwise merging of clusters
  mergedClusLR = matrix(nrow=nClus*(nClus-1)/2,ncol=3)
  iter = 1 #counter
  for(clusIdx1 in 1:(nClus-1)) {
    for(clusIdx2 in (clusIdx1+1):nClus) {
      #        clusIdx1 =1;clusIdx2 = 2
      #Insert calculations
      mergedClusLR[iter,1] = clusIdx1
      mergedClusLR[iter,2] = clusIdx2
      mergedClusLR[iter,3] = getMergedLR(clusterList,clusIdx1,clusIdx2)
      iter = iter + 1
    }
  }
    
  #Third: Keep merge clusters until no changes
  while(TRUE) {

    #Obtaining next candidate: But stop if no improvement
    mergeScore = mergedClusLR[,3] - separate_clusterLR[mergedClusLR[,1]] - separate_clusterLR[mergedClusLR[,2]]
    maxScoreIdx = which.max(mergeScore)
    maxMergeScore = mergeScore[maxScoreIdx]
    
    #Stopping criterion: No more candidates to cluster
    if(length(maxMergeScore)==0 || maxMergeScore<mergeThresh) break 
    
    #Update Clusters
    maxScorePair = mergedClusLR[maxScoreIdx,]
    maxScorePairInds = as.integer(maxScorePair[1:2])
    maxScorePairVal = as.numeric(maxScorePair[3])
    updatedClus = rbind(clusterList[[maxScorePairInds[1]]],
                        clusterList[[maxScorePairInds[2]]])
    
    #Need to adjust existing indices
    #First copy:
    separate_clusterLR_updated = separate_clusterLR
    mergedClusLR_updated = mergedClusLR
    clusterList_updated = clusterList

    #Then make adjustments
    clusRange = seq_len(nClus) #current range
    nClus = nClus - 1 #one less (updated number of clusters)
    clusRange_updated = seq_len(nClus) #updated range
    keepElems = clusRange[-maxScorePairInds] #elements to keep
    clusRange_insPrev = seq_along(keepElems) #range to insert previous values
    
    #Updating variables:
    #Important: LAST CLUSTER IS THE UPDATED ONE (nClus <- 'nClus -1')
    separate_clusterLR_updated[clusRange_insPrev] = separate_clusterLR_updated[keepElems]
    separate_clusterLR_updated[nClus] = maxScorePairVal #insert updated score
    separate_clusterLR_updated = separate_clusterLR_updated[clusRange_updated] #update vector
    clusterList_updated[clusRange_insPrev] = clusterList_updated[keepElems]
    clusterList_updated[[nClus]] = updatedClus #insert updated cluster
    clusterList_updated = clusterList_updated[clusRange_updated] #update list
    
    #Adjusting mergedClusLR_updated to include the range of updated cluster numbers
      # Use old pair-index table only as a template for all pairs in new 1:nClus range.
      # Existing LR values are subsequently mapped from old cluster indices via keepElems.
    keepRows = mergedClusLR_updated[,1]%in%clusRange_updated & 
                mergedClusLR_updated[,2]%in%clusRange_updated
    mergedClusLR_updated = mergedClusLR_updated[keepRows,,drop=FALSE]
    
    #Must map previous vs updated pairwise LR values to correct position (updated)
    keyOld = paste(mergedClusLR[,1],mergedClusLR[,2])
    newInd1 = keepElems[mergedClusLR_updated[,1]]
    newInd2 = keepElems[mergedClusLR_updated[,2]]
    #hasNoNA = !(is.na(newInd1) | is.na(newInd2)) #MUST INCLUDE ALL
    keyOldKeep = paste(newInd1,newInd2)
    indMatchKey = match(keyOldKeep,keyOld)
    mergedClusLR_updated[,3] = mergedClusLR[indMatchKey,3] 
    
    #Need to update LR values for all pairs including new cluster
    rowsNewClusterInd = which(mergedClusLR_updated[,1]%in%nClus | mergedClusLR_updated[,2]%in%nClus)
    for(rowIdx in rowsNewClusterInd) {
#      rowIdx = rowsNewClusterInd[1]
      clusIdx1 = mergedClusLR_updated[rowIdx,1]
      clusIdx2 = mergedClusLR_updated[rowIdx,2]
      mergedClusLR_updated[rowIdx,3] = getMergedLR(clusterList_updated,clusIdx1,clusIdx2)
    }
    
    #Everything is updated (now update variables)
    separate_clusterLR = separate_clusterLR_updated
    clusterList = clusterList_updated
    mergedClusLR = mergedClusLR_updated
    
    if(nClus <= 1) break #no further clusters to combine
  }

  #Post-step: Obtaining DC output based on final clustering result:
  DClist = calcClusterDC(likEvidPerMarker,priorGlobalPerMarker,clusterList)

  #Decide what to return from function
  retList = list(clusterList=clusterList,DClist=DClist)
  if(returnPrepared) {
    retList$prepared = list(likEvidPerMarker=likEvidPerMarker,priorGlobalPerMarker=priorGlobalPerMarker,uRangePerSample=uRangePerSample)
  }
  return(retList)
}

