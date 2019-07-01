library(tidyverse)
library(ggtree) # for tree tidying, annotation and visualization
library(phytools) # for tree manipulation

# Input/output directories 
din <- "input_v3"
parsed <- "parsed_v3"
dout_all <- "results_v3/04_tree"
dout_paper <- "results_v3/paper_species_taxonomy"

# Create output directories
if (! dir.exists(parsed)) dir.create(parsed, recursive = T)
if (! dir.exists(dout_all)) dir.create(dout_all, recursive = T)
if (! dir.exists(dout_paper)) dir.create(dout_paper, recursive = T)

# Read necessary files
tree <-
  paste0(din, "/RAxML_bipartitions.lgc") %>%
  ggtree::read.tree() 
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
  mutate(species_type = if_else(
    str_detect(species, "New species|Unidentified species"), "new", "normal"
  ))

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
tree <- root.phylo(tree, outgroup = "outgroup", resolve = T)
tree$edge.length[1] <- 1

# Save tree for later use
write.tree(tree, file = "parsed_v3/tree_rooted.tree")

# Small helper function
ge <- function(x, y) x >= y 

# Define phylogroups
phylogroups <-
  tribble(
    ~ species_type, ~ species_peripheral, ~ phylogroup,
    "L. floricola", "L. floricola", "floricola group",
    "L. amylophilus", "L. amylotrophicus", "amylophilus group",
    "L. delbrueckii", "L. iners", "delbrueckii group",
    "L. mellifer", "L. mellis", "mellifer group",
    "L. alimentarius", "L. terrae", "alimentarius group",
    "L. dextrinicus", "L. concavus", "dextrinicus group",
    "L. composti", "L. composti", "composti group",
    "L. perolens", "L. harbinensis", "perolens group",
    "L. casei", "L. sharpeae",  "casei group",
    "L. selangorensis", "L. selangorensis", "selangorensis group",
    "L. sakei", "L. fuchuensis", "sakei group",
    "L. coryniformis", "L. bifermentans", "coryniformis group",
    "L. algidus", "L. algidus", "algidus group",
    "L. salivarius", "L. vini", "salivarius group",
    "L. plantarum", "L. fabifermentans", "plantarum group",
    "L. rossiae", "L. siliginis", "rossiae group",
    "L. vaccinostercus", "L. oligofermentans", "vaccinostercus group",
    "L. reuteri", "L. mucosae", "reuteri group",
    "L. collinoides", "L. oryzae", "collinoides group",
    "L. brevis", "L. bambusae", "brevis group",
    "L. kunkeei", "L. ozensis", "kunkeei group",
    "L. fructivorans", "L. florum", "fructivorans group",
    "L. buchneri", "L. senioris", "buchneri group",
    "P. acidilactici", "P. inopinatus", "Pediococcus",
    "Leuc. mesenteroides", "Leuc. fallax", "Leuconostoc",
    "F. fructosus", "F. tropaeoli", "Fructobacillus",
    "W. viridescens", "W. soli", "Weissella",
    "O. oeni", "O. kitaharae", "Oenococcus",
    "C. intestini", "C. intestini", "Convivina"
  )

# Write phylogroups
write_csv(phylogroups, path = "parsed_v3/phylogroups.csv")

# Add phylogroups to the genomes table
genomes_clusters_extended <- 
  genomes_clusters %>%
  left_join(phylogroups %>% rename(species_short = species_type)) %>%
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
  ggtree(layout = "circular", col = "grey", size = 0.3) %<+%
  genomes_clusters_extended +
  geom_rect(xmin = 3.05, xmax = 3.85, ymin = 0, ymax = 300, fill = "grey95") +
  geom_tiplab2(
    aes(angle = angle, label = species_short, col = species_type),
    size = 1.5, align = T, linesize = 0.2, offset = 0.05
  ) +
  geom_point(
    aes(shape = label %>% as.integer() %>% ge(70) %>% as.character()), 
    size = 0.8, color = "grey50"
  ) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_color_manual(values = c(
    "normal" = "#a6cee3", "phylogroup" = "#1f78b4", 
    "new" = "#ff7f00", "unidentified" = "#33a02c"
  )) +
  xlim(c(0, 3.1)) 
ggsave(
  paste0(dout_all, "/tree.png"), 
  bg = "white", units = "cm", width = 17.4, height = 17.4, dpi = 300
)
ggsave(
  paste0(dout_paper, "/figure_3_tree.eps"), 
  bg = "white", units = "cm", width = 17.4, height = 17.4, fonts = c("sans")
)
