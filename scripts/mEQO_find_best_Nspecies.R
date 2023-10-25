# Function to run mEQO for a range of Nmax with a range from one to maximum of maxNmax
# compile the AICs for regressions for the trait of interest as follows:
# linear model for trait as function of the assemblage of microbiome members that were selected by mEQO_ga at each tested value of Nmax 
mEQO_find_best_Nspecies <- function(microbiome = microbiome, trait = trait, maxNmax = maxNmax) {
  
  # define AIC function according to 
  aic<-sapply(1:maxNmax,function(N){
    print(paste0(paste0(paste0(paste0("Testing AIC for ", 
                                      N), " out of "), maxNmax), " Nmax values"))
    tryCatch(
      {  
        if(requireNamespace("vegan")){library(vegan)}
        if(requireNamespace("ade4")){library(ade4)}
        trait_data <- trait
        assemblage<-EQO_ga("c",microbiome,trait, maxIter = 100, Nmax=N, monitor = F)$abundance
        extractAIC(lm(trait_data~assemblage))[2]
        return(AIC(lm(trait_data~assemblage))+2*(N-1))
      },
      error = function(e) {
        cat("Error in iteration", N, ":", conditionMessage(e), "\n")
        return(NULL)})
  })
}