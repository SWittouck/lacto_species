library(tidyverse)
library(tidygenomes)

# read and preprocess data

genes <- read_pangenome("input_v3/representatives_orthofinder")
tree <- ape::read.tree("parsed_v3/tree_rooted.tree")
genomes_clusters <- read_csv("input_v3/genomes_clusters.csv")
clusters_species <- read_csv("parsed_v3/clusters_all_named.csv")
phylogroups <- read_csv("parsed_v3/phylogroups.csv")
pairs <- read_csv("input_v3/genome_pairs.csv.zip")

clusters_species <- 
  clusters_species %>%
  mutate_at("species", str_replace, "Leuconostoc", "Leuc\\.") %>%
  mutate_at("species", str_replace, "(?<=[A-Z])[a-z]{5,}", "\\.")

phylogroups <-
  phylogroups %>%
  rename(genome_type = species_type, genome_peripheral = species_peripheral)

# construct tidygenomes object

lgc <-
  as_tidygenomes(genes) %>%
  add_tidygenomes(tree) %>%
  add_genome_metadata(genomes_clusters) %>%
  add_genome_metadata(clusters_species, by = "cluster") %>%
  add_phylogroups(phylogroups, genome_identifier = species) %>%
  add_tidygenomes(pairs)

save(lgc, file = "parsed_v3/lgc_representatives_tidygenomes.rda")