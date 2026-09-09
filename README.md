
# EFMcommon

Provides methods for evaluating common donor relationships between DNA
profiles using probabilistic deconvolution results from EuroForMix and
EFMmps. The package implements common-donor likelihood ratio and
deconvolution-based compatibility calculations for pairwise profile
comparisons, together with probabilistic methods for clustering multiple
profiles according to their inferred donor origin. The methods are
applicable to DNA mixture comparison and single-cell forensic genetic
analysis.

## Installation

Installation directly from source through GitHub (does not require
compilation):

``` r
install.packages("remotes")
remotes::install_github("oyvble/EFMcommon")
```

## Example of usage

### Part 1: Using euroformix

#### Step 1: Load software and data

``` r
library(euroformix)
library(EFMcommon)
pkg = path.package("EFMcommon") #get package install folder

#Settings:
kit = "ESX17" 
AT = 100 #analytical threshold used (global for all markers)

#Import allele frequencies
freqFile = paste0(pkg,"/examples/",kit,"_Norway.csv") #frequency file to use
popFreq =  euroformix::freqImport(freqFile)[[1]] #need to select 1st population

#Import evidence profiles:
evidfn = paste0(pkg,"/examples/ESX17_evidsCommon.csv")
evidData = euroformix::sample_tableToList(euroformix::tableReader(evidfn)) #read data
```

#### Step 2: Perform model fit and DC for each sample separately

``` r
fitList = list()
for(sampleName in names(evidData)) {
#  i=1
  NOC = length(strsplit(sampleName,"_")[[1]])
  mlefit = euroformix::calcMLE(NOC,evidData[sampleName],popFreq, kit=kit, BWS=FALSE,FWS=FALSE,AT = AT) 
  DC = euroformix::deconvolve(mlefit,alpha=1)
  fitList[[sampleName]] = list(MLE=mlefit,DC=DC) #attach
}
```

#### Step 3: Excecute full clustering

``` r
commonCluster = calcCommonCluster(fitList,mergeThresh = 1)
finalCluster = commonCluster$clusterList #Obtain final clustering
```
