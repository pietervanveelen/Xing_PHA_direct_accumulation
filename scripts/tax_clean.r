## Clean taxonomy information from non-informative tags in phyloseq
## clean taxonomy tags with no information
library(phyloseq)
library(tidyverse)

tax_clean <- function(psdata){
  
  psdata_name <- deparse(substitute(psdata))
  # save uncleaned psdata
  if(!dir.exists("output_data")){dir.create("output_data")}
  saveRDS(psdata, file = paste0("output_data/",psdata_name, "_uncleaned.rds"))
  
  # specify NA taxon name tags to last known taxon names
  tax.clean <- data.frame(tax_table(psdata))
  
  tax.clean = tax.clean %>% mutate_all(funs(str_replace(., "Ambiguous_taxa", "")))
  tax.clean = tax.clean %>% mutate_all(funs(str_replace(., "metagenome", "")))
  tax.clean = tax.clean %>% mutate_all(funs(str_replace(., "uncultured archeaon", "")))
  tax.clean = tax.clean %>% mutate_all(funs(str_replace(., "uncultured bacterium", "")))
  tax.clean = tax.clean %>% mutate_all(funs(str_replace(., "uncultured", "")))
  tax.clean[is.na(tax.clean)] <- ""
  
  # function replace name tags from [https://github.com/joey711/phyloseq/issues/850]
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,2] == ""){
      Kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
      tax.clean[i, 2:7] <- Kingdom
    } else if (tax.clean[i,3] == ""){
      Phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
      tax.clean[i, 3:7] <- Phylum
    } else if (tax.clean[i,4] == ""){
      Class <- paste("Class_", tax.clean[i,3], sep = "")
      tax.clean[i, 4:7] <- Class
    } else if (tax.clean[i,5] == ""){
      Order <- paste("Order_", tax.clean[i,4], sep = "")
      tax.clean[i, 5:7] <- Order
    } else if (tax.clean[i,6] == ""){
      Family <- paste("Family_", tax.clean[i,5], sep = "")
      tax.clean[i, 6:7] <- Family
    } else if (tax.clean[i,7] == ""){
      tax.clean$Species[i] <- paste("Genus_",tax.clean$Genus[i], sep = "_")
    }
  }

  # put cleaned tax_table into phyloseq object
  tax_table(psdata) <- as.matrix(tax.clean)

  #save psdata after cleaning as RDS object
  saveRDS(psdata, file = paste("output_data/",psdata_name, "_cleaned.rds"))
  print("taxonomy table cleaned and saved as .rds object in output_data")
  print(psdata)
}
