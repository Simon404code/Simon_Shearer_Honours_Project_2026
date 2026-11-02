
# --------------------------------------------------
# 4. Add total array counts beneath x-axis labels
# --------------------------------------------------

library(tidyverse)
library(ggbeeswarm)

# Supplied total number of arrays for each array size
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
  15, 7,
  16, 3,
  17, 2,
  18, 0,
  19, 1,
  20, 3,
  21, 3,
  22, 0,
  23, 2,
  24, 0,
  25, 1,
  26, 0,
  27, 0,
  28, 0,
  29, 0,
  30, 0,
  31, 0,
  32, 0,
  33, 0,
  34, 0,
  35, 0,
  36, 0,
  37, 1,
  38, 0,
  39, 0,
  40, 1,
  41, 0,
  42, 0,
  43, 0,
  44, 0,
  45, 0,
  46, 0,
  47, 0,
  48, 0,
  49, 0,
  50, 0,
  51, 0,
  52, 0,
  53, 0,
  54, 0,
  55, 0,
  56, 0,
  57, 0,
  58, 1
) %>%
  filter(total_arrays > 0)

array_totals <- array_totals %>%
  mutate(
    percentage = total_arrays / sum(total_arrays) * 100
  )


# Ensure array sizes are ordered numerically
size_levels <- as.character(array_totals$genes_in_array)

# Keep only sizes with non-zero supplied totals
plot_data_section4 <- plot_data %>%
  mutate(genes_in_array = as.character(genes_in_array)) %>%
  filter(genes_in_array %in% size_levels) %>%
  mutate(
    genes_in_array = factor(
      genes_in_array,
      levels = size_levels
    )
  )

# Labels show size above and total below
x_labels <- setNames(
  paste0(
    round(array_totals$percentage), "%   ",
    "(n = ",
    format(array_totals$total_arrays, big.mark = ","),
    ")   ",
    array_totals$genes_in_array
  ),
  as.character(array_totals$genes_in_array)
)

# Build plot
p <- ggplot(
  plot_data_section4,
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
    )
  ) +
  scale_x_discrete(
    labels = x_labels,
    drop = TRUE
  ) +
  labs(
    title = "Tandem-array Counts Across Eucalyptus Haplotypes",
    subtitle = paste(
      "Points represent haplotypes.",
      "Totals below each size are supplied array counts."
    ),
    x = "Array size (genes per array) and total array count",
    y = "Number of arrays per haplotype",
    colour = "Subpopulation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      lineheight = 1.1
    ),
    axis.title.x = element_text(margin = margin(t = 12)),
    legend.position = "right"
  )
print(p)

# Save a wide, high-resolution graph
ggsave(
  filename = "tandem_array_boxplots_subpopulation_with_totals.png",
  plot = p,
  width = 18,
  height = 8,
  dpi = 300
)

# Also save as PDF for publication-quality scaling
ggsave(
  filename = "tandem_array_boxplots_subpopulation_with_totals.pdf",
  plot = p,
  width = 18,
  height = 8
)
