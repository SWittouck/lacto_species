library(tidyverse)
library(tidyorthogroups)
library(ggtree)

lgc <- create_tidyorthogroups("input_v3/representatives_orthofinder/", clade = "lgc")

genomes <- read_csv("input_v3/genomes_clusters.csv")
clusters <- read_csv("parsed_v3/clusters_all_named.csv")
phylogroups <- read_csv("parsed_v3/phylogroups.csv")
tree <- read.tree("parsed_v3/tree_rooted.tree")

clusters <- 
  clusters %>%
  mutate_at("species", str_replace, "Leuconostoc", "Leuc\\.") %>%
  mutate_at("species", str_replace, "(?<=[A-Z])[a-z]{5,}", "\\.")

genomes <- left_join(genomes, clusters)

tree_species <- tree

tree_species$tip.label <-
  tibble(genome = tree_species$tip.label) %>%
  left_join(genomes) %>%
  pull(species)

ape::write.tree(tree_species, "parsed_v3/tree_rooted_species.tree")

species_phylogroups <- 
  phylogroups %>%
  expand_clades(species_type, species_peripheral, tree_species) %>%
  rename(species = genome)

write_csv(species_phylogroups, "parsed_v3/species_phylogroups.csv")

lgc <-
  lgc %>%
  add_genome_tibble(genomes) %>%
  add_genome_tibble(species_phylogroups) %>%
  mutate_genomes(
    species_ordered = factor(species, levels = !! tree_species$tip.label)
  ) %>%
  mutate_genomes(
    is_type_of_phylogroup = species %in% !! phylogroups$species_type
  )

save(lgc, file = "parsed_v3/lgc_representatives_tidyorthogroups.rda")
