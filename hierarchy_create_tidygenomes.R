library(tidyverse)
library(tidygenomes)

din_genuslevel <- "input_v2_g2110/hierarchy_orthofinder_genuslevel"
din_specieslevel <- "input_v2_g2110/hierarchy_orthofinder_specieslevel"

genes_genuslevel <- read_pangenome(din_genuslevel)
lgc_genuslevel <- as_tidygenomes(genes_genuslevel)

dins_species <- list.dirs(din_specieslevel, recursive = F)

for (din_species in dins_species) {
  
  species <- str_extract(din_species, "[^/]+$")
  print(species)
  genes_specieslevel <- read_pangenome(din_species) %>%
    mutate(orthogroup = str_c(!! species, orthogroup, sep = "_"))
  lgc_specieslevel <- as_tidygenomes(genes_specieslevel)
  lgc_specieslevel$genomes$species <- species
  
  lgc_genuslevel <- 
    lgc_genuslevel %>%
    inflate_pangenome(lgc_specieslevel, species = species)
  
}

lgc_hier <- lgc_genuslevel %>% filter_genomes(! is.na(species))

save(lgc_hier, file = "parsed_v2_g2110/lgc_hierarchy_tidygenomes.rda")
