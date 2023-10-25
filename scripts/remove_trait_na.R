# if present, remove NA data in trait data from microbiome data before proceeding to mEQO runs
remove_trait_na <- function(microbiome) {  
  
  microbiome_df = data.frame(microbiome, check.names = F)
  row.names(microbiome_df) <- microbiome$rowname
  trait_arg <- names(microbiome_df)[ncol(microbiome_df)]
  trait_var <- as.numeric(microbiome_df[[trait_arg]])
  names(trait_var) <- microbiome_df$rowname
  
  if (anyNA(trait_var)) {
    
    NA_count = sum(is.na(trait_var))
    NA_samples = names(which(is.na(trait_var)))
    cat(paste0("Trait data contained ", NA_count, " NA sample(s): ", NA_samples, " removed."))
    EQO_input <- microbiome_df[!is.na(trait_var), -ncol(microbiome_df)]
    trait_clean <- trait_var[!is.na(trait_var)]
  } else {
    EQO_input <- microbiome_df[, -ncol(microbiome_df)]
    trait_clean <- trait_var}
  
  return(list(microbiome = EQO_input, trait = trait_clean))
}