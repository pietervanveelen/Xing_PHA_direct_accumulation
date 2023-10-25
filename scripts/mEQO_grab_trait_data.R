# subset the trait(s) of interest
mEQO_trait_data <- function(ps_data = ps_data, trait = trait, check_names = F){
  
  trait_arg = substitute(trait)
  
  if(check_names == TRUE){
    trait_data = 
      ps_data %>% 
      as_tibble() %>% 
      mutate(Sample = if_else(grepl(".*[^0-9][0-9]{2}$", Sample),
                              Sample,
                              str_c(str_sub(Sample, end = -2),"0", str_sub(Sample, -1)))) %>% 
      mutate(across(.cols = everything(), .fns = ~str_replace(.x, ",", "."))) %>% 
      mutate(Sample = factor(Sample)) %>% 
      select(Sample, trait_arg) %>% 
      distinct() %>% 
      arrange(Sample)
  } else {
    trait_data = 
      ps_data %>% 
      as_tibble() %>% 
      mutate(across(.cols = everything(), .fns = ~str_replace(.x, ",", "."))) %>% 
      mutate(Sample = factor(Sample)) %>% 
      select(Sample, trait_arg) %>% 
      distinct() %>% 
      arrange(Sample)
  }
  
  return(trait_data)
}