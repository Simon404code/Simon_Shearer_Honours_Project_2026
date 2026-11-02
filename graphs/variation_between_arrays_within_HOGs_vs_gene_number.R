#### Variation of tandem array size in HOGs 
### N0.HOG0000296 
### similarity_min	similarity_max	similarity_mean	similarity_range
### 97.5	100	99.6523	2.5

#HOG_ID(array)	Type	CNV(min-max)	avg_CNV	frequency	duplicate_rate	min_gene_length	max_gene_length	avg_gene_length
#N0.HOG0000296 OG0000000	Shell	2-58	7.4	15/50	15/50	364	668	441.5217391304348

#HOG_ID	list_of_gene_IDs	number_of_Genes	number_of_Genomes	unique_chromosomes	avg_gene_length	tandem_array_rate_percent	classifier	avg_CNV
#N0.HOG0000296 OG0000000	Eg_RK320_H1.07G017700;Eg_RK436_H1.07G020300;Eg_RK715_H2.07G019800;Eg_RK84_H1.07G023000;Eg_RK05_H1.07G019700;Eg_RK05_H1.07G019800;Eg_RK25_H1.07G019400;Eg_RK25_H1.07G019800;Eg_RK25_H1.07G020300;Eg_RK25_H1.07G020400;Eg_RK25_H1.07G020600;Eg_RK25_H1.07G020700;Eg_RK25_H1.07G021900;Eg_RK25_H1.07G022700;Eg_RK25_H2.07G019400;Eg_RK25_H2.07G019800;Eg_RK25_H2.07G020700;Eg_RK320_H2.07G023500;Eg_RK320_H2.07G024500;Eg_RK388_H1.07G022500;Eg_RK388_H1.07G022700;Eg_RK436_H2.07G028000;Eg_RK436_H2.07G028200;Eg_RK436_H2.07G029000;Eg_RK498_H1.07G024200;Eg_RK498_H1.07G024300;Eg_RK498_H2.07G018100;Eg_RK498_H2.07G018200;Eg_RK498_H2.07G018300;Eg_RK498_H2.07G018400;Eg_RK498_H2.07G018500;Eg_RK498_H2.07G018800;Eg_RK498_H2.07G018900;Eg_RK498_H2.07G019000;Eg_RK498_H2.07G019100;Eg_RK498_H2.07G019200;Eg_RK498_H2.07G020800;Eg_RK498_H2.07G020900;Eg_RK498_H2.07G021000;Eg_RK498_H2.07G021100;Eg_RK561_H1.07G018700;Eg_RK561_H1.07G019200;Eg_RK561_H1.07G019800;Eg_RK561_H1.07G020000;Eg_RK687_H1.07G020500;Eg_RK687_H1.07G020800;Eg_RK687_H1.07G021200;Eg_RK687_H1.07G021500;Eg_RK687_H2.07G018400;Eg_RK687_H2.07G018500;Eg_RK687_H2.07G018600;Eg_RK687_H2.07G018700;Eg_RK687_H2.07G018900;Eg_RK687_H2.07G019000;Eg_RK687_H2.07G019100;Eg_RK687_H2.07G019200;Eg_RK687_H2.07G019300;Eg_RK687_H2.07G019400;Eg_RK687_H2.07G019500;Eg_RK687_H2.07G019600;Eg_RK687_H2.07G019700;Eg_RK687_H2.07G019800;Eg_RK687_H2.07G019900;Eg_RK687_H2.07G020000;Eg_RK687_H2.07G020100;Eg_RK687_H2.07G020200;Eg_RK687_H2.07G020300;Eg_RK687_H2.07G020400;Eg_RK687_H2.07G020500;Eg_RK687_H2.07G020600;Eg_RK687_H2.07G020700;Eg_RK687_H2.07G020900;Eg_RK687_H2.07G021000;Eg_RK687_H2.07G021100;Eg_RK687_H2.07G021200;Eg_RK687_H2.07G021300;Eg_RK687_H2.07G021400;Eg_RK687_H2.07G021500;Eg_RK687_H2.07G021600;Eg_RK687_H2.07G021700;Eg_RK687_H2.07G021800;Eg_RK687_H2.07G021900;Eg_RK687_H2.07G022000;Eg_RK687_H2.07G022100;Eg_RK687_H2.07G022200;Eg_RK687_H2.07G022300;Eg_RK687_H2.07G022400;Eg_RK687_H2.07G022500;Eg_RK687_H2.07G022600;Eg_RK687_H2.07G022700;Eg_RK687_H2.07G022800;Eg_RK687_H2.07G022900;Eg_RK687_H2.07G023000;Eg_RK687_H2.07G023100;Eg_RK687_H2.07G023200;Eg_RK687_H2.07G023300;Eg_RK687_H2.07G023400;Eg_RK687_H2.07G023500;Eg_RK687_H2.07G023600;Eg_RK687_H2.07G023700;Eg_RK687_H2.07G023800;Eg_RK687_H2.07G023900;Eg_RK687_H2.07G024000;Eg_RK687_H2.07G024100;Eg_RK687_H2.07G024200;Eg_RK687_H2.07G024300;Eg_RK792_H1.07G021900;Eg_RK792_H1.07G022000;Eg_RK792_H2.07G021400;Eg_RK792_H2.07G021700;Eg_RK86_H2.07G020800;Eg_RK86_H2.07G022200;Eg_RK86_H2.07G022400;Eg_RM02_H2.07G018900;Eg_RM02_H2.07G019300	115	19	1	441.5217391304348	78.94736842105263	minority shell	7.4


library(tidyverse)

# --------------------------------------------------
# 1. Prepare tandem-array data
# --------------------------------------------------

arrays <- tandem_array_HOG_summary %>%
  filter(
    !is.na(HOG_ID),
    !is.na(genome),
    !is.na(arrayID),
    !is.na(genes_in_array),
    genes_in_array >= 2
  ) %>%
  mutate(
    HOG_ID = as.character(HOG_ID),
    genome = as.character(genome),
    genes_in_array = as.numeric(genes_in_array)
  ) %>%
  distinct(HOG_ID, genome, arrayID, .keep_all = TRUE)

# --------------------------------------------------
# 2. Find the largest individual array per HOG per genome
# --------------------------------------------------

hog_genome_sizes <- arrays %>%
  group_by(HOG_ID, genome) %>%
  summarise(
    largest_array_size = max(genes_in_array, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------------------------------
# 3. Calculate the number of contributing genomes
#    and the array-size range for each HOG
# --------------------------------------------------

hog_variation <- hog_genome_sizes %>%
  group_by(HOG_ID) %>%
  summarise(
    number_of_genomes = n_distinct(genome),
    smallest_array = min(largest_array_size, na.rm = TRUE),
    largest_array = max(largest_array_size, na.rm = TRUE),
    variance_range = largest_array - smallest_array,
    .groups = "drop"
  ) %>%
  filter(number_of_genomes > 1)

# --------------------------------------------------
# 4. Create the dot plot
# --------------------------------------------------

p <- ggplot(
  hog_variation,
  aes(
    x = number_of_genomes,
    y = variance_range
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.8,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks()
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks()
  ) +
  labs(
    title = "Variation in Tandem-Array Sizes Across the Pangenome",
    x = "Number of Genomes Contributing to a HOG",
    y = "The largest Array-Size difference within a HOG"
    ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0)
  )

print(p)

# --------------------------------------------------
# 5. Save the plot and the data used to make it
# --------------------------------------------------

ggsave(
  "tandem_array_variation_vs_genome_count.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "tandem_array_variation_vs_genome_count.pdf",
  plot = p,
  width = 9,
  height = 6
)

write_csv(
  hog_variation,
  "tandem_array_variation_by_HOG.csv"
)
