#' @title calcGenoProb
#' @author Oyvind Bleka
#' @description Returns a list of joint genotypes with corresponding joing probabilities for unknowns contributors for a given marker.
#' @details A simplified version of calcGjoint function from euroformix package
#' The function returns the list of all possible genotypes, with corresponding probabilities. The allele-names in popFreq needs to be numbers.
#' @param freq A vector of allele frequencies for a given population.  
#' @param nU Number of unknowns
#' @param fst Assumed theta/fst-correction
#' @param nTyped contains a vector indicating number of typed allele
#' @return Glist A list with genotypes and genotype probabilities 
#' @export 

calcGenoProb = function(freq,nU=1,fst=0,nTyped=NULL) {
  if(length(fst)!=1) stop("Wrong input length for fst")
  if(length(nU)!=1) stop("Wrong input length for number of unknowns (nU)")
  if(nU==0) stop("You must specify at least one unknown to use this function!")
  if(!is.numeric(freq))  stop("freq argument must be numeric!") 
  
  sumsToOne = all.equal(1,sum(freq)) #checking if freqs sums to one (this must always be the case for EFM)
  if(is.character(sumsToOne)) warning("freq argument must sum to one!")
  
  #Function calculates genotypes for all possib
  nn = length(freq) #number of alleles
  nG <- nn*(1+nn)/2 #number of allele outcome
  #nG^nU #NUMBER OF ITERATIONS
  #print(nG^nU)
  
  #PRESTEP: GET GENOTYPE OUTCOME: SIMILAR TO getGlist function
  av <- names(freq)   
  
  #G-matrix is the vectorized upper triangular (1,1),(1,2),...,(1,n),(2,2),(2,3),...,(2,n),....,(n,n)
  G = numeric()
  for(i in 1:nn) G = rbind(G, cbind( av[rep(i,nn - i + 1)], av[i:nn] ))
  if(nrow(G)!=nG) warning("The size of the genotype outcome is not as expected!") #Must be true!
  ishomG <-  G[,1]==G[,2] #find G variants which are homozygous
  
  #COUNT NUMBER OF ALLELES IN refK:
  mkvec <- nTyped #setNames(nTyped,names(av))
  if(is.null(nTyped)) mkvec = rep(0,length(av))
  if(length(nTyped)!=length(av)) stop("Length of number of typed alleles did not match the number of alleles!")

  #####################################
  #GENERAL FORMULA WITH KINSHIP MODULE#
  #####################################
  
  P = function(Pi,ni,n) (ni*fst + (1-fst)*Pi)/(1+(n-1)*fst) #helpfunction for allele prob
  
  genoProb = function(mkvec2) { #helpfunction to get genotype prob for spec
    Gprob1 <- Gprob0 <- rep(NA,nG) #used to store genotype prob|typed alleles (all outcome)
    
    #prepare Gprobs for K0:
    for(i in 1:nG) { #for each genotype
      refG <- G[i,] #get proposed genotype
      indf <- match(refG,names(freq)) #get (allele) index of frequency
      tmpval = P(freq[indf[1]],ni=mkvec2[indf[1]],n=sum(mkvec2)) #prob. of 1st allele
      if(ishomG[i]) { 
       Gprob0[i] <- tmpval*P(freq[indf[2]],ni=mkvec2[indf[2]]+1,n=sum(mkvec2)+1) #get hom freq
      } else {
       Gprob0[i] <- 2*tmpval*P(freq[indf[2]],ni=mkvec2[indf[2]],n=sum(mkvec2)+1) #get het freq   
      }
    } #end for each type
    #sum(Gprob0) #MUST BE 1
    return(Gprob0) #return if no relatedness
  }
  
  #CREATE RECURSION FUNCTION TO CALCULATE pG for a given counted alleles and prev. given pG
  Gprob <- array(data = NA, dim = rep(nG,nU)) #create structure. Each index correspond to genotype in G
  
  #function returning genoprob for given start mkvec and pGeno
  calcGprob = function(mkvecTMP,track=NULL,Uk=1,pGenoTMP=1) { #Uk gives depth in recursion. Stops after finishing loop on Uk=1 (i.e. returning from all)
    #track=NULL for Uk=1
    pGenotmp <- pGenoTMP*genoProb(mkvecTMP) #update probabilities conditioned on prev typed alleles
    
    for(gind in seq_len(nG) )  { #for each genotype
      mkvectmp <- mkvecTMP #obtain a new copy here
      
      for(l in 1:2) { #UPDATE ALLELE COUNTER:
       indf <- which(G[gind,l]==av) #get index of frequency	 
       mkvectmp[indf] <-  mkvectmp[indf] + 1  #update counter of sampled allele
      }
      
      if(Uk<nU) { #if still more unknowns to calculate for
       calcGprob(mkvecTMP=mkvectmp,track=c(track,gind),Uk=Uk+1,pGenoTMP=pGenotmp[gind] )  #recurse deeper
      } else { #final unknown has been achievied: Need to store genotype and return back
       #print(paste0(paste0(track,collapse="-"),"-",gind,":",pGenotmp[gind])) #print tracks
       Gprob[ rbind(c(track,gind)) ] <<- pGenotmp[gind] #store calculations
      } 
    } #end for each genotype
  }
  #RUN RECURSION
  calcGprob(mkvecTMP=mkvec)
  #sum(Gprob) #PROB OF G OUTCOME MUST BE 1  
  
  return(list(G=G,Gprob=Gprob)) #return list
}


