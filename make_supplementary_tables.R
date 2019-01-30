library(tidyverse)

din_input <- "input_v3"
din_parsed <- "parsed_v3"
dout <- "results_v3"

# genomes 

list(
  genomes_quality = paste0(din_input, "/genome_table.csv"),
  genomes_clusters = paste0(din_input, "/genomes_clusters.csv"),
  genomes_ncbi = paste0(din_input, "/genomes_assembly_reports.csv")
) %>%
  map(read_csv) %>%
  reduce(full_join) %>%
  select(genome, completeness, redundancy, ncbi_species = species, cluster, species, isolation_source) %>%
  left_join(
    paste0(din_parsed, "/clusters_all_named.csv") %>%
      read_csv() %>%
      select(cluster, species)
  ) %>%
  write_csv(paste0(dout, "/table_S1_genomes.csv"))

# type genomes 

file.copy(paste0(din_input, "/lgc_genomes_type.csv"), paste0(dout, "/table_S2_type_genomes.csv"))

# clusters 

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
