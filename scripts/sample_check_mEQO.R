# custom functions for making sample lists for comparing samples between microbiome data and trait data
samples1_char = function(microbiome){
  if(is.factor(microbiome$rowname)){
    samples1 = sort(as.character(microbiome$rowname))}
  else {
    samples1 = sort(microbiome$rowname)}
  return(samples1)
}

samples2_char = function(trait){
  if(is.factor(trait$Sample)){
    samples2 = sort(as.character(trait$Sample))}
  else {
    samples2 = sort(trait$Sample)}
  return(samples2)
}