## Clean taxonomy information from non-informative tags in phyloseq
## clean taxonomy tags with no information
library(phyloseq)
library(tidyverse)
library(glue)

tax_clean <- function(psdata){
  
  psdata_name <- deparse(substitute(psdata))
  # save uncleaned psdata
  if(!dir.exists("output_data")){dir.create("output_data")}
  saveRDS(psdata, file = paste0("output_data/",psdata_name,"_uncleaned.rds"))
  
  # specify NA taxon name tags to last known taxon names
  tax.clean <- data.frame(tax_table(psdata))
  
  tax.clean2 = 
    tax.clean %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(Kingdom = str_replace(Kingdom, pattern = "d__", replacement = "")) %>% 
    mutate(across(everything(),
                  ~ str_replace_all(.x,  "Ambiguous_taxa|metagenome|uncultured archeaon|uncultured bacterium|uncultured prokaryote|uncultured soil bacterium|uncultured rumen bacterium|uncultured compost bacterium|uncultured organism|uncultured$",
                                    replacement = NA_character_
                  ))) %>%
    replace(is.na(.), NA_character_) %>% 
    mutate(Phylum = if_else(is.na(Phylum), paste0("Phylum of ", Kingdom), Phylum),
           Class = if_else(is.na(Class), paste0("Class of ", Phylum), Class),
           Order = if_else(is.na(Order), paste0("Order of ", Class), Order),
           Family = if_else(is.na(Family), paste0("Family of ", Order), Family),
           Genus = if_else(is.na(Genus), paste0("Genus of ", Family), Genus),
           Species = if_else(is.na(Species), paste0("Species of ", Genus), Species)
    ) %>% 
    mutate(across(.cols = Kingdom:Species, .fns = ~if_else(str_detect(.,'\\bof\\b.*\\bof\\b'), paste0(word(., 1)," ", word(., 2)," ", word(., -1)), .)))
  
  # put cleaned tax_table into phyloseq object
  tax_table(psdata) <- tax_table(as.matrix(tax.clean2))
  
  #save psdata after cleaning as RDS object
  saveRDS(psdata, file = glue::glue("output_data/",{psdata_name},"_cleaned.rds", .sep = ""))
  print("taxonomy table cleaned and saved as .rds object in output_data")
  return(psdata)
}


