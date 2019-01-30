library(tidyverse)
library(ggtree) # for tree tidying, annotation and visualization
library(phytools) # for tree manipulation

# Input/output directories 
din <- "input_v3"
parsed <- "parsed_v3"
dout <- "results_v3"

# Create output directories
if (! dir.exists(parsed)) dir.create(parsed, recursive = T)
if (! dir.exists(dout)) dir.create(dout, recursive = T)

# Read necessary files
tree <-
  paste0(din, "/RAxML_bipartitions.lgc") %>%
  read.tree() 
genomes_clusters <- 
  paste0(din, "/genomes_clusters.csv") %>%
  read_csv()
clusters <-
  paste0(parsed, "/clusters_all_named.csv") %>%
  read_csv() 

# Add cluster information to the genomes table
genomes_clusters <-
  genomes_clusters %>%
  left_join(clusters) %>%
  mutate(species_type = if_else(str_detect(species, "New species|Unidentified species"), "new", "normal"))

# Find the MRCA node of the clade that contains L. mellifer and L. concavus;
# this is the outgroup clade
out_exs_names <- c("Lactobacillus mellifer", "Lactobacillus concavus")
out_exs_accessions <- genomes_clusters %>%
  filter(species %in% out_exs_names) %>%
  pull(genome) 
out_exs_tips <- which(tree$tip.label %in% out_exs_accessions)
out_mrca_node <- findMRCA(tree, tip = out_exs_tips)

# Add an artifical outgroup to the branch leading to the outgroup clade
tree <- bind.tip(
  tree, tip.label = "outgroup", edge.length = 2, 
  where = out_mrca_node, position = 0.05
)

# Reroot the tree on the branch leading to the outgroup
tree <- reroot(tree, node.number = which(tree$tip.label == "outgroup"), position = 1)

# Small helper function
ge <- function(x, y) x >= y 

# Add phylogroups to the genomes table
genomes_clusters_extended <- 
  genomes_clusters %>%
  left_join(tribble(
    ~ species, ~ phylogroup,
    "Lactobacillus floricola", "floricola group",
    "Lactobacillus amylophilus", "amylophilus group",
    "Lactobacillus delbrueckii", "delbrueckii group",
    "Lactobacillus mellifer", "mellifer group",
    "Lactobacillus alimentarius", "alimentarius group",
    "Lactobacillus dextrinicus", "dextrinicus group",
    "Lactobacillus composti", "composti group",
    "Lactobacillus perolens", "perolens group",
    "Lactobacillus casei",  "casei group",
    "Lactobacillus selangorensis", "selangorensis group",
    "Lactobacillus sakei", "sakei group",
    "Lactobacillus coryniformis", "coryniformis group",
    "Lactobacillus algidus", "algidus group",
    "Lactobacillus salivarius", "salivarius group",
    "Lactobacillus plantarum", "plantarum group",
    "Lactobacillus rossiae", "rossiae group",
    "Lactobacillus vaccinostercus", "vaccinostercus group",
    "Lactobacillus reuteri", "reuteri group",
    "Lactobacillus collinoides", "collinoides group",
    "Lactobacillus brevis", "brevis group",
    "Lactobacillus kunkeei", "kunkeei group",
    "Lactobacillus fructivorans", "fructivorans group",
    "Lactobacillus buchneri", "buchneri group",
    "Pediococcus acidilactici", "Pediococcus",
    "Leuconostoc mesenteroides", "Leuconostoc",
    "Fructobacillus fructosus", "Fructobacillus",
    "Weissella viridescens", "Weissella",
    "Oenococcus oeni", "Oenococcus"
  )) %>%
  mutate(
    species_type =
      case_when(
        str_detect(species, "^New") ~ "new",
        str_detect(species, "^Unidentified") ~ "unidentified",
        ! is.na(phylogroup) ~ "phylogroup",
        TRUE ~ "normal"
      )
  ) %>%
  mutate(species_type = if_else(is.na(phylogroup), species_type, "phylogroup"))
  
# Plot the tree
tree %>%
  ggtree(layout = "circular", col = "grey", size = 0.5, alpha = 1) %<+%
  genomes_clusters_extended +
  geom_tiplab2(
    aes(angle = angle, label = species, col = species_type),
    size = 1.5, align = T, linesize = 0.2, offset = 0.05
  ) +
  geom_point(aes(shape = label %>% as.integer() %>% ge(70) %>% as.character()), size = 0.8, color = "grey50") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_color_manual(values = c(
    "normal" = "#a6cee3", "phylogroup" = "#1f78b4", "new" = "#ff7f00", "unidentified" = "#33a02c"
  )) +
  xlim(c(0, 3.5)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA)
  )
ggsave(paste0(dout, "/figure_3.png"), bg = "white", units = "cm", width = 20, height = 20)
