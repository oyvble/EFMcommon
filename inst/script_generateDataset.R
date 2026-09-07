#rm(list=ls())
#setwd("C:\\Users\\oyvbl\\Dropbox\\Forensic\\euroformix0\\calcCommon")
source("helpfiles_commonDC.R")
library(euroformix);sessionInfo() #load package
pkg = path.package("euroformix") #get package install folder

Qallele = "99"
fst = 0.03
kit = "ESX17" #defining kit to use (must be defined in getKit())
AT = 100 #analytical threshold used (global for all markers)
LRthresh = 6
includePrior = TRUE #whether prior should be considered
knownRefs = NULL #Indicate which references that are known
refData = NULL

knownAllelesMarkerList = list() #indicate known alleles per markers 

#Importing allele frequencies
freqFile = paste0(pkg,"/FreqDatabases/",kit,"_Norway.csv") #frequency file to use
popFreq =  freqImport(freqFile)[[1]] #need to select 1st population

#Importing evidence profiles:
evidfn = "ESX17_evidsCommon.csv"
evidfn_path = evidfn
#evidfn_path = paste0(pkg,"/examples/",evidfn)
evidData = sample_tableToList(tableReader(evidfn_path)) #read data
sampleNames = names(evidData) #obtain sample names
nSamples = length(sampleNames)

verbose = TRUE
#Performing model fit and DC for each sample separately
fitList = list()
for(i in seq_along(sampleNames)) {
#  i=1
  sampleName = sampleNames[i]
  NOC = length(strsplit(sampleName,"_")[[1]])
  mlefit = calcMLE(NOC,evidData[i],popFreq, kit=kit, BWS=FALSE,FWS=FALSE,AT = AT,fst = fst) 
  DC = deconvolve(mlefit,alpha=1)
  fitList[[i]] = list(MLE=mlefit,DC=DC) #attach
}

#Obtaining marker list
markerList = lapply(fitList,function(x) x$MLE$prepareC$markerNames)

#Perform pairwise comparison of common unknowns
commonUnknowns = NULL
for(i in 1:(nSamples-1)) {
  for(j in (i+1):nSamples) {
#    i=2;j=4
    commonLR = calcCommonLR(fitList[[i]]$MLE,fitList[[j]]$MLE,fitList[[i]]$DC,fitList[[j]]$DC)
    commonContrInds = which(commonLR>=LRthresh,arr.ind = TRUE)
    
    if(nrow(commonContrInds)==0) next
    newrows = cbind(Evid1=i,Evid2=j,commonContrInds)
    commonUnknowns = rbind(commonUnknowns,newrows)
  }
}
rownames(commonUnknowns) = NULL

#obtain clusters of common unknowns
commonList = list() #commonUnknowns #obtain common list
for(i in seq_len(nrow(commonUnknowns))) {
#  i=3
  match = commonUnknowns[i,]
  rows = rbind(match[c(1,3)],match[c(2,4)]) #put a evid/Uid row-wise
  colnames(rows) = c("Evid","Unknown")
  insIdx = as.integer() #default is first list element
  if(length(commonList)>0) {
    found = FALSE
    hasSame = rep(NA,length(commonList))
    for(j in seq_along(commonList)) {
      #j=1
      #CHECKING IF MATCH IS CONNECTED TO EXISTING LIST
      key1 = paste0(rows[,1],rows[,2]) #creating a key to compare from
      key2 = paste0(commonList[[j]][,1],commonList[[j]][,2]) #creating another key (compare to)
      hasSame[j] = any(key1%in%key2) #any already existing?
    }
    insIdx = which(hasSame)
  }
  
  #Add to list
  if(length(insIdx)==0) { #in case of now earlier found
    commonList[[length(commonList)+1]] = rows  #add as new element
  } else if(length(insIdx)==1) { #in case of previously found
    commonList[[insIdx]] = rbind(commonList[[insIdx]],rows)  #attach to existing
    commonList[[insIdx]] = unique(commonList[[insIdx]]) #keep only unique
  } else {
    cat("Note: A match was found to fit multiple clusters and hence ignored.")
  }
}
if(length(commonList)==0) return(NULL)

#We provide improved DC for the unknown in each cluster
DCgenoListClusters = list()
for(clusId in seq_along(commonList)) {
#  clusId = 1
  clusTab = commonList[[clusId]] #Obtain cluster of matches
  
  #Pre-step: Deside which component if same mixture
  nTable = table(clusTab[,1]) #count number of components per sample
  keepRows = rep(TRUE,nrow(clusTab)) #which rows to keep
  hasMultiple = as.integer(names(nTable)[nTable>1])
  for(i in hasMultiple) {
  #  i=hasMultiple
    evidRowInds = which(clusTab[,1]%in%i)
    checkCtrs = clusTab[evidRowInds,2] #extract contributors
    MxCtrs = fitList[[i]]$MLE$fit$thetahat2[checkCtrs] #obtain Mx estimates
    CtrKeepInd = which.max(MxCtrs) #indicate which contributor to use
    removeRows = evidRowInds[-CtrKeepInd] #indicate which rows to remove
    keepRows[removeRows] = FALSE #turn off
  }
  clusTab = clusTab[keepRows,]
  
  #Obtain full set of markers analysed
  markersCluster = unique(unlist(markerList[clusTab[,1]]))
  
  #obtaining DC-result per marker (across all cluster elements)
  DCgenoList = list()
  for(marker in markersCluster) {
#    marker = markersCluster[2]
    hasMarker = sapply(clusTab[,1],function(x) marker%in%markerList[[x]])
    clusTab_marker = clusTab[hasMarker,,drop=FALSE]
    
    #Extracting each DC object in the cluster
    DCtabList = list()
    alleleList = list() #store observed alleles per object
    for(oidx in seq_len(nrow(clusTab_marker))) {
  #    oidx=1
      sampleRow = clusTab_marker[oidx,]
      
      #Obtain posterior probabilities from per-object:
      fitListSample = fitList[[sampleRow[1]]] #get fitted list for specific sample
      DCtab = fitListSample$DC$table3 #obtain full DC-tab
      DCtab = DCtab[grepl(sampleRow[2],DCtab[,1]) & DCtab[,2]==marker, ] #filter relevant
      postProbs = setNames(as.numeric(DCtab[,4]),DCtab[,3])
      
      #Obtain prior prob to correct for;
      c = fitListSample$MLE$prepareC #obtain prepareC object
      mind = which(c$markerNames==marker) #get marker index
      mrng = c$startIndMarker_nAlleles[mind] + seq_len(c$nAlleles[mind]) #get range
      allelesSample = c$alleleNames[mrng]
      freqSample = setNames(c$freqs[mrng],allelesSample) #frequencies
      knowncounts = setNames(c$maTyped[mrng],allelesSample)
      knownAlleles = rep(names(knowncounts),knowncounts)
      fstMarker = c$fst[mind]
      
      #Obtain prior probability from frequencies only
      priorProbs = getPriorProbs(freqSample,fstMarker, knownAlleles)

      #Obtaining weights by dividing post/prio
      priorProbsSorted = priorProbs[match(names(postProbs),names(priorProbs))]
      #all(names(priorProbsSorted)==names(postProbs))
      likProbs = postProbs/priorProbsSorted
      likProbs[is.nan(likProbs)] = 0
      
      DCtabList[[oidx]] = likProbs/sum(likProbs) #normalize
      alleleList[[oidx]] = setdiff(allelesSample,Qallele)
    }
    
    #Expanding the genotype possibilities
    genosUnique = sort(unique(unlist(lapply(DCtabList,names))))
    genosTab = getGenoMatHelper(genosUnique)
    genosProbMat = matrix(NA,ncol=length(DCtabList),nrow=length(genosUnique),dimnames = list(genosUnique,NULL))
    
    #Assigning probabilities for each outcome (per object)
    for(oidx in seq_along(DCtabList)) {
    #    oidx=1
      DCprobs = DCtabList[[oidx]]
      genosProbMat[match(names(DCprobs),genosUnique),oidx] = DCprobs
      
      #insert the missing cells:
      isNArowIdx = which(is.na(genosProbMat[,oidx])) #get rows missing probs
      if(length(isNArowIdx)>0) {
        genosIsNA = genosTab[isNArowIdx,,drop=FALSE] #Obtain genotypes
        
        genosUpdated = getMissingGenos(genosIsNA,alleleList[[oidx]])
        matchInd = match(genosUpdated,names(DCprobs))
        genosProbMat[isNArowIdx,oidx] = DCprobs[matchInd]
        genosProbMat[is.na(genosProbMat[,oidx])] = 0 #force as zero if not found in list
      }
    }
    
    #calculate the product across all samples
    genoProd = apply(genosProbMat, 1, prod)  #NB: speed can improve here!
    
    #Apply prior (given known (unique) typed individuals for that marker and samples used in cluster)
    #Assumes same frequencies used
    allelesUnique = sort(unique(unlist(alleleList)))
    freqs = popFreq[[marker]][allelesUnique]
    if(any(is.na(freqs))) stop("Allele frequency not found. A global frequency set must be applied!")
    freqSUM = sum(freqs)
    if(freqSUM<1) { 
      freqs = c(freqs,1-freqSUM)
      names(freqs)[length(freqs)] = Qallele
    }
    if(!is.null(refData)) knownAlleles = NULL

    #Calculating a "common prior adjustment"
    #knownAlleles = knownAllelesMarkerList[[marker]]
    knownAlleles = NULL
    priorProbs = getPriorProbs(freqs,fstMarker, knownAlleles)
    priorProbsSorted = priorProbs[match(names(genoProd),names(priorProbs))]
    priorProbsSorted[is.na(priorProbsSorted)] = 0
#    all(names(genoProd)==names(priorProbsSorted))
    genoProd = genoProd*priorProbsSorted #adjust
    DCgenoList[[marker]] = genoProd/sum(genoProd) #normalize again and insert
  } #end for each marker
  DCgenoListClusters[[clusId]] = DCgenoList #insert cluster result
} #for each cluster

return(DCgenoListClusters)

#pdf(paste0("DCtest_clus",clusId,".pdf"),width=10,height=6)
#for(i in seq_along(DCgenoList)) barplot(DCgenoList[[i]],las=2,main=names(DCgenoList)[i])
#dev.off()

  