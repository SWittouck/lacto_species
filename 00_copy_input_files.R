library(tidyverse)

# remark: change fin to the location of the data folder of the
# lgc_evolutionary_genomics pipeline
fin <- "~/rhambo/media/harddrive/lgc_comparative_pangenomics/data_v3"
fout <- "input_v3"

if (! dir.exists(fout)) dir.create(fout)

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
  map_chr(~ paste0(fin, .)) %>%
  map(~ file.copy(., fout))
