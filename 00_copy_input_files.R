library(tidyverse)

# remark: change din to the location of the data folder of the
# lgc_evolutionary_genomics pipeline
din <- "~/rhambo/media/harddrive/lgc_comparative_pangenomics/data_v3"
dout <- "input_v3"

if (! dir.exists(dout)) dir.create(dout)

# copy files from legen pipeline

c(
  "/quality_control/genome_table.csv",
  "/scgs/candidate_scg_table.csv",
  "/scgs/score_table.csv.zip",
  "/similarities/genome_pairs.csv.zip",
  "/genome_clusters/genomes_clusters.csv",
  "/genome_clusters/exclusivity.csv",
  "/genome_clusters/phylogeny/RAxML_bipartitions.lgc",
  "/taxonomy/blast/hits_to_type_16s_genes.tsv",
  "/taxonomy/lgc_genomes_type.csv",
  "/taxonomy/lgc_names.rda",
  "/taxonomy/rrnas_16S.txt",
  "/taxonomy/genomes_assembly_reports.csv"
) %>%
  map_chr(~ paste0(din, .)) %>%
  map(~ file.copy(., dout))

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
