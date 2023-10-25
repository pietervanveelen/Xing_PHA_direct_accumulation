# check the sample lists of traits and microbiome data on lengths

# custom functions for making sample lists for comparing samples between microbiome data and trait data
samples1_char = function(microbiome){
  if(is.numeric(microbiome$rowname)){stop("rowname is numeric, but should contain sample IDs as character strings or factor levels")}
  if(is.factor(microbiome$rowname)){
    samples1 = sort(as.character(microbiome$rowname))}
  else {
    samples1 = sort(microbiome$rowname)}
  return(samples1)
}

samples2_char = function(trait){
  if(is.numeric(trait$Sample)){stop("Sample is a numeric vector, but should contain sample IDs as character strings or factor levels")}
  if(is.factor(trait$Sample)){
    samples2 = sort(as.character(trait$Sample))}
  else {
    samples2 = sort(trait$Sample)}
  return(samples2)
}

check_sample_lists = function(microbiome = microbiome, trait = trait)  {
  
  if(missing(microbiome)){
    stop("Microbiome data missing for sample list check vs. trait data sample list")}
  if(missing(trait)){
    stop("Trait data missing for sample list check vs. trait data sample list")}
  
  # remove samples for which trait data are missing
  if(!identical(samples1_char(microbiome = microbiome),
                samples2_char(trait = trait))){
    stop("Trait data and microbiome data are not of equal lengths:\n ensure identical samples ids")
  } else {
    microbiome = dplyr::left_join(microbiome, trait, by = c("rowname" = "Sample"))
  }
  return(microbiome = microbiome) 
}

