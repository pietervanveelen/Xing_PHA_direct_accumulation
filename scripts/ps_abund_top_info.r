# Order taxa from phyloseq on abundance is phyloseq dataset

# psdata = phyloseq object
# top_nr = the number of taxa to be output in descending order of abundance

ps_abund_top_info <- function(psdata, top_nr)
{
  otu = otu_table(psdata) %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "taxon")
  abund <- otu %>% select(-taxon) %>% rowSums()
  otu <- otu %>% mutate(abund = abund) %>% 
    select(taxon, abund)
  tax = tax_table(psdata) %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "taxon")
  otu_tax = inner_join(tax, otu, by = "taxon")
  otu_tax %>% 
    arrange(desc(abund)) %>% 
    as_tibble(rownames = "rank") %>% 
    select(rank, taxon, abund, everything()) %>% 
    slice_head(n = top_nr)
}



