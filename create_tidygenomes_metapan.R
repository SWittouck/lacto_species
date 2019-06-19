library(tidyverse)
library(tidygenomes)
# devtools::load_all("../../tidygenomes/tidygenomes/")

# paths to input files
din_genuslevel <- "input_v2_g2110/hierarchy_orthofinder_genuslevel"
din_specieslevel <- "input_v2_g2110/hierarchy_orthofinder_specieslevel"
fin_lgc_repr <- "parsed_v3/lgc_tidygenomes_representatives.rda"

# create genus level pangenome
genes_genuslevel <- read_pangenome(din_genuslevel)
lgc_genuslevel <- as_tidygenomes(genes_genuslevel)

# inflate genus level pangenome with species level pangenomes
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

# filter out new species
lgc_genuslevel <- lgc_genuslevel %>% filter_genomes(! is.na(species))

# collapse species to gain metapangenome
lgc_metapan <- lgc_genuslevel %>% collapse_species()

# remove some objects to free memory 
rm(genes_genuslevel, genes_specieslevel, lgc_genuslevel, lgc_specieslevel)

# load tidygenomes object with representative genomes
load(fin_lgc_repr)

# change genome names in lgc_repr from accession numbers to species names 
lgc_repr <- update_genomes(lgc_repr, new_name = str_replace(species_full, " ", "_"))

# remove the pangenome from lgc_repr
lgc_repr$genes <- NULL 
lgc_repr$orthogroups <- NULL

# add the remaining lgc_repr data (pairs, tree, phylogroups) to lgc_metapan
lgc_metapan <- lgc_metapan %>% add_tidygenomes(lgc_repr)

# save lgc_metapan object
save(lgc_metapan, file = "parsed_v2_g2110/lgc_tidygenomes_metapan.rda")
