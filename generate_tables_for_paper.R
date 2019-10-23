library(tidyverse)

din_input <- "input_v3"
din_parsed <- "parsed_v3"
dout <- "results_v3/paper_species_taxonomy"

# table 1 - species mergers and splits

clusters_exclusivity <-
  read_csv(paste0(din_parsed, "/clusters_exclusivity.csv"))

read_csv(paste0(din_input, "/taxonomy/splits_and_mergers.csv")) %>%
  left_join(clusters_exclusivity) %>%
  select(cluster, merger_split, names, min_cni_within) %>%
  arrange(merger_split, desc(min_cni_within)) %>% 
  write_csv(paste0(dout, "/table_1_splits_and_mergers2.csv"))

# table 2 - genome clusters without type genomes 

clusters_zerotypegenomes
read_csv(paste0(din_input, "/taxonomy/clusters_zerotypegenomes.csv"))%>% 
  write_csv(paste0(dout, "/table_2_clusters_zerotypegenomes.csv"))

# table S1 - genomes

genomes_parks <-
  read_csv(paste0(din_input, "/genomes_parks.csv")) %>%
  select(
    genome = accession, checkm_completeness, checkm_contamination, gtdb_species
  )

list(
  genomes_quality = paste0(din_input, "/genome_table.csv"),
  genomes_clusters = paste0(din_input, "/genomes_clusters.csv"),
  genomes_ncbi = paste0(din_input, "/genomes_assembly_reports.csv")
) %>%
  map(read_csv) %>%
  reduce(full_join) %>%
  select(
    genome, legen_completeness = completeness, legen_redundancy = redundancy, 
    legen_cluster = cluster, 
    ncbi_species = species, ncbi_isolation_source = isolation_source
  ) %>%
  left_join(
    paste0(din_parsed, "/clusters_all_named.csv") %>%
      read_csv() %>%
      select(legen_cluster = cluster, legen_species = species)
  ) %>%
  left_join(genomes_parks) %>%
  select(
    genome, legen_completeness, legen_redundancy, legen_cluster, legen_species,
    ncbi_species, ncbi_isolation_source, checkm_completeness, 
    checkm_contamination, gtdb_species
  ) %>%
  write_csv(paste0(dout, "/table_S1_genomes.csv"))

# table S2 - type genomes 

file.copy(paste0(din_input, "/lgc_genomes_type.csv"), paste0(dout, "/table_S2_type_genomes.csv"))

# table S3 - clusters

list(
  clusters_named = paste0(din_parsed, "/clusters_all_named.csv"),
  clusters_exclusivity = paste0(din_parsed, "/clusters_exclusivity.csv")
) %>%
  map(read_csv) %>%
  reduce(full_join) %>%
  select(
    cluster, species, n_genomes = n, min_cni_within, max_cni_between, 
    closest_neighbor = closest_cni, exclusive = exclusive_cni
  ) %>%
  write_csv(paste0(dout, "/table_S3_clusters.csv"))
