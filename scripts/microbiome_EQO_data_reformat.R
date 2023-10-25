# join trait data with microbiome to get order same
microbiome_EQO = function(microbiome = microbiome, 
                          trait = trait_data){
  
  # load helper functions  
  source("scripts/transpose_df.R")
  source("scripts/check_sample_lists.R")
  source("scripts/remove_trait_na.R")
  
  # stop function if arguments are missing from function call
  if(missing(microbiome)){
    stop("Microbiome data in wide format is required:\ndata.frame with taxa as row.names and samples as columns")}
  if(missing(trait)){
    stop("Trait data required:\nuse mEQO_trait_data() function to reshape the focal trait from psmelt data")}
  
  # 1. transpose the microbiome input data    
  microbiome = transpose_df(microbiome)
  
  # 2. custom function to check if samples match for microbiome and trait data and then split both components
  microbiome = check_sample_lists(microbiome = microbiome, trait = trait)
  
  # 3. remove samples with NA in trait, if present
  microbiome_noNA = remove_trait_na(microbiome = microbiome)
  
  # 4. reformat output as required format for mEQO_ga
  
  microbiome = as(microbiome_noNA[[1]][,-1], "matrix")
  colnames(microbiome) = colnames(microbiome_noNA[[1]])[c(2:(ncol(microbiome_noNA[[1]])))]
  rownames(microbiome) = microbiome_noNA[[1]][,1]
  
  trait = microbiome_noNA[[2]]
  
  # 5. save matrix for microbiome input for mEQO and named vector for trait data
  if(!dir.exists("output_data")){dir.create("output_data")}
  saveRDS(object = microbiome, file = "output_data/EQO_input_microbiome_matrix.rds")
  saveRDS(object = trait, file = "output_data/EQO_input_trait_named_vector.rds")
  
  # 6. return output
  return(list(microbiome = microbiome, trait = trait))
  
} # END