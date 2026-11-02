library(ggplot2)
library(tidyverse)
library(ggbeeswarm)

# --------------------------------------------------
# 1. Prepare the data
# --------------------------------------------------
array_totals <- tribble(
  ~genes_in_array, ~total_arrays,
  2, 22871,
  3, 2104,
  4, 420,
  5, 147,
  6, 60,
  7, 37,
  8, 20,
  9, 10,
  10, 10,
  11, 6,
  12, 5,
  13, 5,
  14, 3,
  15, 7
  #  16, 3,
  #  17, 2,
  #  18, 0,
  #  19, 1,
  #  20, 3,
  #  21, 3,
  #  22, 0,
  #  23, 2,
  #  24, 0,
  #  25, 1,
  #  26, 0,
  #  27, 0,
  #  28, 0,
  #  29, 0,
  #  30, 0,
  #  31, 0,
  #  32, 0,
  #  33, 0,
  #  34, 0,
  #  35, 0,
  #  36, 0,
  #  37, 1,
  #  38, 0,
  #  39, 0,
  #  40, 1,
  #  41, 0,
  #  42, 0,
  #  43, 0,
  #  44, 0,
  #  45, 0,
  #  46, 0,
  #  47, 0,
  #  48, 0,
  #  49, 0,
  #  50, 0,
  #  51, 0,
  #  52, 0,
  #  53, 0,
  #  54, 0,
  # 55, 0,
  #  56, 0,
  #  57, 0,
  #  58, 1
) %>%
  filter(total_arrays > 0)

array_totals <- array_totals %>%
  mutate(
    percentage = total_arrays / sum(total_arrays) * 100
  )



arrays <- tandem_array_HOG_summary %>%
  filter(
    !is.na(genome),
    !is.na(arrayID),
    !is.na(genes_in_array),
    !is.na(Subpopulation),
    genes_in_array >= 2
  ) %>%
  mutate(
    genome = as.character(genome),
    Subpopulation = as.factor(Subpopulation),
    genes_in_array = as.integer(genes_in_array)
  ) %>%
  distinct(genome, arrayID, .keep_all = TRUE)

# Check whether haplotypes map to multiple subpopulations
subpopulation_check <- arrays %>%
  distinct(genome, Subpopulation) %>%
  count(genome, name = "n_subpopulations") %>%
  filter(n_subpopulations > 1)

if (nrow(subpopulation_check) > 0) {
  warning(
    "Some haplotypes map to multiple subpopulations. ",
    "Inspect subpopulation_check."
  )
}

# --------------------------------------------------
# 2. Count arrays of each size in each haplotype
# --------------------------------------------------

haplotype_counts <- arrays %>%
  count(
    genome,
    Subpopulation,
    genes_in_array,
    name = "number_of_arrays"
  )

# --------------------------------------------------
# 3. Create plot_data and exclude zero counts
# --------------------------------------------------

plot_data <- haplotype_counts %>%
  filter(
    number_of_arrays > 0,
    genes_in_array %in% array_totals$genes_in_array
  ) %>%
  mutate(
    genes_in_array = factor(
      genes_in_array,
      levels = array_totals$genes_in_array
    )
  )

# --------------------------------------------------
# 4. Create boxplots with custom x-axis labels
# --------------------------------------------------

# Ensure array sizes are ordered numerically and match array_totals
size_levels <- as.character(array_totals$genes_in_array)

plot_data <- plot_data %>%
  filter(genes_in_array %in% array_totals$genes_in_array) %>%
  mutate(
    genes_in_array = factor(
      genes_in_array,
      levels = array_totals$genes_in_array
    )
  )

# Create x-axis labels: percentage, total count, and array size
x_labels <- setNames(
  paste0(
    round(array_totals$percentage), "%    ",
    "(n = ",
    format(array_totals$total_arrays, big.mark = ","),
    ")    ",
    array_totals$genes_in_array
  ),
  size_levels
)

p <- ggplot(
  plot_data,
  aes(
    x = genes_in_array,
    y = number_of_arrays
  )
) +
  geom_boxplot(
    aes(group = genes_in_array),
    outlier.shape = NA,
    fill = "grey92",
    colour = "grey35",
    width = 0.65
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(colour = Subpopulation),
    width = 0.18,
    size = 2.2,
    alpha = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "North" = "#1B9E77",
      "South" = "#D95F02",
      "Mackay" = "#FF69B4",
      "Mondi Tree Improvement Research" = "#C0C0C0",
      "Merensky Timber" = "#C0C0C0"
    ),
    labels = NULL
  ) +
  scale_x_discrete(
    labels = x_labels,
    drop = TRUE
  ) +
  labs(
    title = "Tandem-array counts across Eucalyptus haplotypes",
    subtitle = "Non-zero counts; logarithmic y-axis",
    x = "  Array size (genes per array), percentage of total, and total count",
    y = "Number of arrays per haplotype (log10 scale)",
    colour = "Subpopulation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.title.x = element_text(
      margin = margin(t = 12)
    ),
    legend.position = "right"
  )

# Display the plot
print(p)

# --------------------------------------------------
# 5. Save the graph and CSV table
# --------------------------------------------------

ggsave(
  filename = "tandem_array_boxplots_subpopulation_log10.png",
  plot = p,
  width = 12,
  height = 7,
  dpi = 300
)

write_csv(
  plot_data,
  "tandem_array_counts_by_haplotype_subpopulation_no_zeros.csv"
)

# Print output directory
getwd()
