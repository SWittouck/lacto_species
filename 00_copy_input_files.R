library(tidyverse)

# remark: change din to the location of the data folder of the
# lgc_evolutionary_genomics pipeline
din <- "~/rhambo/media/harddrive/stijn/projects_phd/legen_pipeline/data_v3"
dout <- "input_v3"

# copy files from legen pipeline

c(
  "/quality_control/genome_table.csv",
  "/scgs/candidate_scg_table.csv",
  "/scgs/score_table.csv.zip",
  "/similarities/genome_pairs.csv.zip",
  "/genome_clusters/genomes_clusters.csv",
  "/genome_clusters/exclusivity.csv",
  "/taxonomy/clusters_all_named.csv",
  "/taxonomy/clusters_zerotypegenomes.csv",
  "/taxonomy/genomes_assembly_reports.csv",
  "/taxonomy/splits_and_mergers.csv",
  "/representatives_v3_2/phylogeny/RAxML_bipartitions.lgc"
) %>%
  tibble(path = .) %>%
  mutate(
    from = str_c(din, path), 
    to = str_c(dout, str_extract(path, "/[^/]+"))
  ) %>%
  {
    walk(.$to, dir.create, showWarnings = F)
    walk2(.$from, .$to, file.copy)
  }

# download GTDB (Parks et al.) taxonomy

url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv"
genomes_parks_full <- read_tsv(url)

genomes <- read_csv(paste0(dout, "/genome_table.csv"))

gtdb_ranks <- c(
  "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order", 
  "gtdb_family", "gtdb_genus", "gtdb_species"
)

genomes_parks <-
  genomes_parks_full %>%
  mutate(accession = str_remove(accession, "^.{3}")) %>%
  separate(gtdb_taxonomy, sep = ";", into = gtdb_ranks) %>%
  mutate_at(gtdb_ranks, str_remove, "^.{3}") %>%
  filter(accession %in% !! genomes$genome)

write_csv(genomes_parks, paste0(dout, "/genomes_parks.csv"))
