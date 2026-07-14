########################################
## Bacillota + Actinomycetota tree plot
########################################

library(ape)
library(data.table)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(phangorn)
library(stringr)
library(svglite)
library(ggstar)

setwd("/data/san/data1/users/david/ellen/e_faecium/tree/genomes")

tree_path <- "denovo_out_bacteria/infer/gtdbtk.bac120.decorated.tree"
tax_path <- "denovo_out_bacteria/infer/gtdbtk.bac120.decorated.tree-taxonomy"
custom_tax_path <- "custom_taxonomy.txt"

outer_fruit_tip_set <- c(
  "GCA_000742975.1_ASM74297v1_genomic",
  "GCA_008121495.1_ASM812149v1_genomic"
)
outer_fruit_color <- "#5BB99B"

endolysin_sensitive_tip_set <- c(
  "GCA_001594345.1_ASM159434v1_genomic",
  "GCA_008121495.1_ASM812149v1_genomic"
)

tree <- read.tree(tree_path)

tax <- fread(
  tax_path,
  header = FALSE,
  sep = "\t",
  col.names = c("user_genome", "classification")
) %>%
  mutate(
    classification = gsub("; ", ";", classification, fixed = TRUE),
    phylum_raw = str_match(classification, "(^|;)p__([^;]*)")[, 3],
    class_raw = str_match(classification, "(^|;)c__([^;]*)")[, 3],
    phylum_group = str_remove(phylum_raw, "_[A-Z].*$"),
    class_group = str_remove(class_raw, "_[A-Z].*$"),
    family = str_match(classification, "(^|;)f__([^;]*)")[, 3],
    genus = str_match(classification, "(^|;)g__([^;]*)")[, 3],
    species = str_match(classification, "(^|;)s__([^;]*)")[, 3],
    class_raw = if_else(is.na(class_raw) | class_raw == "", "unclassified", class_raw),
    class_group = if_else(is.na(class_group) | class_group == "", "unclassified", class_group)
  )

custom_tax <- fread(custom_tax_path) %>%
  mutate(
    is_custom = TRUE,
    custom_genus = na_if(str_match(classification, "(^|;)g__([^;]*)")[, 3], ""),
    custom_species = na_if(str_match(classification, "(^|;)s__([^;]*)")[, 3], ""),
    custom_tax_label = coalesce(custom_species, custom_genus, user_genome),
    custom_group = case_when(
      str_detect(classification, "Ruminococcus_B gnavus") ~ "Mediterraneibacter",
      str_detect(classification, "Mediterraneibacter") ~ "Mediterraneibacter",
      str_detect(classification, "Enterococcus") ~ "Enterococcus",
      str_detect(classification, "Streptococcus") ~ "Streptococcus",
      str_detect(classification, "Staphylococcus") ~ "Staphylococcus",
      str_detect(classification, "Blautia") ~ "Blautia",
      str_detect(classification, "Bacillus") ~ "Bacillus",
      str_detect(classification, "Clostridioides") ~ "Clostridioides",
      str_detect(classification, "Lactiplantibacillus") ~ "Lactiplantibacillus",
      str_detect(classification, "Listeria") ~ "Listeria",
      str_detect(classification, "Bifidobacterium") ~ "Bifidobacterium",
      str_detect(classification, "Faecalibacterium") ~ "Faecalibacterium",
      TRUE ~ custom_tax_label
    ),
    custom_color = case_when(
      str_detect(classification, "Ruminococcus_B gnavus") ~ "#C38ABE",
      str_detect(classification, "Mediterraneibacter") ~ "#C38ABE",
      str_detect(classification, "Enterococcus") ~ "#277382",
      str_detect(classification, "Streptococcus") ~ "#5BB99B",
      str_detect(classification, "Staphylococcus") ~ "#29932F",
      str_detect(classification, "Blautia") ~ "#E6A0C4",
      str_detect(classification, "Bacillus") ~ "#000000",
      str_detect(classification, "Clostridioides") ~ "#335080",
      str_detect(classification, "Lactiplantibacillus") ~ "#F4C7E4",
      str_detect(classification, "Listeria") ~ "#561671",
      str_detect(classification, "Bifidobacterium") ~ "#002e3a",
      str_detect(classification, "Faecalibacterium") ~ "#cb8686",    
      TRUE ~ "#111111"
    )
  ) %>%
  select(user_genome, is_custom, custom_group, custom_color)

plot_tax <- tax %>%
  filter(phylum_group %in% c("Bacillota", "Actinomycetota")) %>%
  add_count(class_raw, name = "class_count") %>%
  filter(class_count > 150) %>%
  left_join(custom_tax, by = "user_genome") %>%
  mutate(
    is_custom = coalesce(is_custom, FALSE),
    custom_group = if_else(is_custom, custom_group, NA_character_),
    custom_color = if_else(is_custom, custom_color, NA_character_)
  )

class_counts <- tax %>%
  filter(phylum_group %in% c("Bacillota", "Actinomycetota")) %>%
  count(phylum_group, class_raw, class_group, name = "class_count") %>%
  mutate(kept = class_count > 150) %>%
  arrange(desc(class_count), phylum_group, class_raw)

tree_plot <- keep.tip(tree, intersect(tree$tip.label, plot_tax$user_genome))
plot_tax <- plot_tax %>%
  filter(user_genome %in% tree_plot$tip.label)

missing_tax_tips <- setdiff(tree_plot$tip.label, plot_tax$user_genome)
if (length(missing_tax_tips) > 0L) {
  stop("Tree tips missing taxonomy rows: ", paste(missing_tax_tips, collapse = ", "))
}

actinomycetota_tips <- intersect(
  plot_tax$user_genome[plot_tax$phylum_group == "Actinomycetota"],
  tree_plot$tip.label
)
if (length(actinomycetota_tips) < 2L) {
  stop("Need at least two Actinomycetota tips to root the tree.")
}
tree_plot <- root(tree_plot, outgroup = actinomycetota_tips, resolve.root = TRUE)

fruit_tax <- plot_tax %>%
  mutate(fruit_class = class_raw) %>%
  select(
    user_genome,
    fruit_class,
    phylum_raw,
    phylum_group,
    class_raw,
    class_group,
    class_count,
    family,
    genus,
    species,
    classification,
    is_custom,
    custom_group,
    custom_color
  )

missing_fruit_tips <- setdiff(tree_plot$tip.label, fruit_tax$user_genome)
if (length(missing_fruit_tips) > 0L) {
  stop("Tree tips missing fruit rows: ", paste(missing_fruit_tips, collapse = ", "))
}

base_palette <- c(
  "Actinomycetes" = "#4C78A8",
  "Acidimicrobiia" = "#72B7B2",
  "Coriobacteriia" = "#54A24B",
  "Thermoleophilia" = "#B279A2",
  "Nitriliruptoria" = "#499894",
  "Rubrobacteria" = "#A0CBE8",
  "Clostridia" = "#335080",
  "Bacilli" = "#3D0149",
  "Bacilli_A" = "#9C4DCC",
  "Peptococcia" = "#7F5F86",
  "Negativicutes" = "#5BB99B",
  "Desulfitobacteriia" = "#9E7C9F",
  "Fermentithermobacillia" = "#C4A2C0",
  "Sulfobacillia" = "#BFA0C0"
)

observed_classes <- sort(unique(fruit_tax$fruit_class))
extra_classes <- setdiff(observed_classes, names(base_palette))
extra_palette <- setNames(character(0), character(0))
if (length(extra_classes) > 0L) {
  extra_palette <- setNames(
    grDevices::rainbow(length(extra_classes), s = 0.8, v = 0.9),
    extra_classes
  )
}
class_palette <- c(base_palette, extra_palette)
class_palette <- class_palette[observed_classes]

missing_palette_classes <- setdiff(unique(fruit_tax$fruit_class), names(class_palette))
if (length(missing_palette_classes) > 0L) {
  stop("Fruit classes missing palette colours: ", paste(missing_palette_classes, collapse = ", "))
}

bacilli_exact_tips <- plot_tax$user_genome[plot_tax$class_raw == "Bacilli"]
bacilli_exact_node <- getMRCA(tree_plot, bacilli_exact_tips)
bacilli_exact_descendant_tips <- tree_plot$tip.label[
  phangorn::Descendants(tree_plot, bacilli_exact_node, type = "tips")[[1]]
]

bacilli_group_tips <- plot_tax$user_genome[plot_tax$class_group == "Bacilli"]
bacilli_group_node <- getMRCA(tree_plot, bacilli_group_tips)
bacilli_group_descendant_tips <- tree_plot$tip.label[
  phangorn::Descendants(tree_plot, bacilli_group_node, type = "tips")[[1]]
]

fruit_tax <- fruit_tax %>%
  mutate(
    fruit_color = unname(class_palette[fruit_class]),
    in_exact_bacilli_mrca = user_genome %in% bacilli_exact_descendant_tips,
    in_group_bacilli_mrca = user_genome %in% bacilli_group_descendant_tips
  )

outer_fruit_tax <- fruit_tax %>%
  filter(user_genome %in% outer_fruit_tip_set) %>%
  transmute(
    user_genome,
    outer_marker = "Selected tips"
  )

tip_meta <- fruit_tax %>%
  transmute(
    label = user_genome,
    is_custom,
    custom_group,
    custom_color,
    endolysin_response = if_else(
      user_genome %in% endolysin_sensitive_tip_set,
      "Sensitive",
      "Insensitive"
    )
  )

tip_lab_meta <- tip_meta %>%
  left_join(
    fruit_tax %>% select(label = user_genome, fruit_class),
    by = "label"
  )

custom_tip_palette_table <- tip_meta %>%
  filter(is_custom) %>%
  distinct(custom_group, custom_color) %>%
  arrange(custom_group)
custom_tip_palette <- setNames(
  custom_tip_palette_table$custom_color,
  custom_tip_palette_table$custom_group
)

class_palette_table <- tibble(
  fruit_class = names(class_palette),
  fruit_color = unname(class_palette)
)

fwrite(fruit_tax, "gtdb_tree_plot_taxonomy.tsv", sep = "\t")
fwrite(class_counts, "gtdb_tree_class_counts.tsv", sep = "\t")
fwrite(class_palette_table, "gtdb_tree_class_palette.tsv", sep = "\t")
fwrite(
  fruit_tax %>% filter(phylum_raw != phylum_group | class_raw != class_group),
  "gtdb_tree_suffix_normalized_tips.tsv",
  sep = "\t"
)
fwrite(
  fruit_tax %>% filter(in_exact_bacilli_mrca, class_raw != "Bacilli"),
  "bacilli_exact_mrca_non_bacilli_tips.tsv",
  sep = "\t"
)
fwrite(
  fruit_tax %>% filter(in_group_bacilli_mrca, class_group != "Bacilli"),
  "bacilli_group_mrca_non_bacilli_tips.tsv",
  sep = "\t"
)

p <- ggtree(tree_plot, layout = "circular") %<+% tip_meta +
  geom_tree(linewidth = 0.1, color = "#b0b0b0") +
  geom_tippoint(
    aes(
      subset = is_custom,
      fill = custom_group,
      shape = endolysin_response,
      color = endolysin_response
    ),
    size = 8,
    stroke = 2
  ) +
  scale_fill_manual(
    name = "Genera of interest",
    values = custom_tip_palette,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(shape = 21, size = 5, color = "black")
    )
  ) +
  scale_shape_manual(
    name = "Endolysin",
    values = c("Insensitive" = 21, "Sensitive" = 23),
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(fill = "grey70", size = 5, color = c("black", "red"))
    )
  ) +
  scale_color_manual(
    values = c("Insensitive" = "black", "Sensitive" = "red"),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_fruit(
    data = fruit_tax %>% select(user_genome, fruit_class),
    geom = geom_tile,
    mapping = aes(y = user_genome, fill = fruit_class),
    width = 0.18,
    offset = 0.03
  ) +
  scale_fill_manual(
    name = "Class",
    values = class_palette,
    drop = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  # ggnewscale::new_scale_fill() +
  # geom_fruit(
  #   data = outer_fruit_tax,
  #   geom = geom_tile,
  #   mapping = aes(y = user_genome, fill = outer_marker),
  #   width = 0.08,
  #   offset = 0.04
  # ) +
  # scale_fill_manual(
  #   name = "Outer ring",
  #   values = c("Selected tips" = outer_fruit_color),
  #   guide = guide_legend(nrow = 1, byrow = TRUE)
  # ) +
  theme_tree() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = grid::unit(0.8, "cm")
  ) 

ggsave("gtdb_tree_custom.png", p, width = 45, height = 30, units = "cm")
ggsave("gtdb_tree_custom.svg", p, width = 45, height = 30, units = "cm")

p_right_legend <- p +
  guides(
    fill = guide_legend(ncol = 1, byrow = TRUE),
    shape = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(fill = "grey70", size = 5, color = c("black", "red"))
    )
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = grid::unit(1, "cm"),
    legend.spacing.y = grid::unit(0.12, "cm")
  )

ggsave("gtdb_tree_custom_legend_right.png", p_right_legend, width = 45, height = 35, units = "cm", dpi = 300)
ggsave("gtdb_tree_custom_legend_right.svg", p_right_legend, width = 45, height = 35, units = "cm", dpi = 300)



########################################
## plot with fill tiplab as the class colour
########################################

p2_base <- ggtree(tree_plot, layout = "circular") %<+% tip_lab_meta
p2_tip_segments <- p2_base$data %>%
  filter(isTip) %>%
  transmute(
    label,
    fruit_class,
    x,
    xend = max(p2_base$data$x, na.rm = TRUE) + 0.2,
    y,
    yend = y
  )

p2 <- p2_base +
  geom_segment(
    data = p2_tip_segments,
    aes(x = x, xend = xend, y = y, yend = yend, color = fruit_class),
    linewidth = 3.5,
    lineend = "butt",
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_color_manual(
    name = "Class",
    values = class_palette,
    drop = FALSE,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 3, size = 0)
    )
  ) +
  geom_tree(linewidth = 0.1, color = "#b0b0b0") +
  ggnewscale::new_scale_color() +
  geom_tippoint(
    aes(
      subset = is_custom,
      fill = custom_group,
      shape = endolysin_response,
      color = endolysin_response
    ),
    size = 8,
    stroke = 2
  ) +
  scale_fill_manual(
    name = "Genera of interest",
    values = custom_tip_palette,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(shape = 21, size = 5, color = "black", label = "")
    )
  ) +
  scale_shape_manual(
    name = "Endolysin",
    values = c("Insensitive" = 21, "Sensitive" = 23),
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(fill = "grey70", size = 5, color = c("black", "red"))
    )
  ) +
  scale_color_manual(
    values = c("Insensitive" = "black", "Sensitive" = "red"),
    guide = "none"
  ) +
  # ggnewscale::new_scale_fill() +
  theme_tree() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = grid::unit(0.8, "cm")
  ) 


p2_right_legend <- p2 +
  guides(
    shape = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(fill = "grey70", size = 5, color = c("black", "red"))
    )
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = grid::unit(1, "cm"),
    legend.spacing.y = grid::unit(0.12, "cm")
  )

ggsave("gtdb_tree_custom_legend_right_tiplab.png", p2_right_legend, width = 45, height = 35, units = "cm", dpi = 300)
ggsave("gtdb_tree_custom_legend_right_tiplab.svg", p2_right_legend, width = 45, height = 35, units = "cm", dpi = 300)



message("Tips plotted: ", length(tree_plot$tip.label))
message("Bacillota tips plotted: ", sum(plot_tax$phylum_group == "Bacillota"))
message("Actinomycetota outgroup tips plotted: ", length(actinomycetota_tips))
message("Rooted on Actinomycetota: ", is.rooted(tree_plot))
message("Classes coloured in geom_fruit: ", length(class_palette))
message("Classes kept with >150 genomes: ", sum(class_counts$kept))
message("Classes excluded with <=150 genomes: ", sum(!class_counts$kept))
message("Bacilli_A tips have their own colour: ", sum(fruit_tax$class_raw == "Bacilli_A"))
message("Custom tip points drawn: ", sum(tip_meta$is_custom))
message("Outer geom_fruit tips drawn: ", nrow(outer_fruit_tax))
message("Non-Bacilli tips inside exact Bacilli MRCA: ", sum(fruit_tax$in_exact_bacilli_mrca & fruit_tax$class_raw != "Bacilli"))
message("Non-Bacilli-group tips inside Bacilli-group MRCA: ", sum(fruit_tax$in_group_bacilli_mrca & fruit_tax$class_group != "Bacilli"))
