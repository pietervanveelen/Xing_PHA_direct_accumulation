# find best Nmax using AIC, by constraining Nmax to a minimal number of species
best_Nmax = function(aic_input = aic_input, maxNmax = maxNmax, trait = trait, minimal_Nspecies = 3) {
  
  # define which trait is run
  if(missing(trait)){stop("Trait name is missing: provide a character string of the trait of interest")}
  trait_name = trait
  aic = as.list(aic_input)
  # define which N species ended with AIC = NULL because of error logical subscript too long    
  N_null = which(sapply(aic,is.null))
  
  # remove the NULL from the data by filtering
  aic_filt <- data.frame(AIC = do.call(rbind,
                                       Filter(Negate(is.null), aic)))
  
  # determine best Nmax, but should be at least 3
  aic_best_nmax = 
    tibble(aic_filt, N=seq(1:maxNmax)[!seq(1:maxNmax) %in% N_null]) %>% 
    filter(N>(minimal_Nspecies-1)) %>% 
    slice_min(AIC) %>% 
    pull(N)
  
  cat(paste("Best Nmax = ", aic_best_nmax," (with minimal number of species set at ", minimal_Nspecies,")"))
  
  # plot the AIC as function of Nmax
  plot_aic_nmax = 
    tibble(aic_filt, N=seq(1:maxNmax)[!seq(1:maxNmax) %in% N_null]) %>% 
    filter(N > (minimal_Nspecies-1)) %>% 
    ggplot(aes(x = N, y=AIC)) +
    geom_line(color="darkred", linewidth = 2) +
    geom_point(shape = 21, color = "white", fill = "darkred", size = 5) +
    scale_x_continuous(breaks = seq(1:10)) +
    labs(title = "Lowest AIC denotes optimal number of species",
         x="Maximal Number of Genera (Nmax)",
         y = glue::glue("AIC\n{trait_name}")) +
    theme_classic()
  
  # print(plot)
  print(plot_aic_nmax)
  
  # save plot
  if(!dir.exists("figures")){dir.create("figures")}
  ggsave(plot = plot_aic_nmax, filename = glue::glue("figures/{trait_name}_AIC_vs_Nmax.pdf"))
  
  return(aic_best_nmax)
} # END