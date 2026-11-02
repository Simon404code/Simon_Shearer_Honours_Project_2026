
library(tidyverse)
library(ggplot2)
library(scales)

# ============================================================
# TANDEM-ARRAY SIZE VARIATION WITHIN HOGs
#
# X-axis:
#   Number of genomes containing the HOG.
#
# Y-axis:
#   Difference between the largest and smallest tandem-array
#   sizes among genomes containing a qualifying array.
#
# Points:
#   Blue = Variation > 0
#   Pink = Zero variation
#
# Black line:
#   Cumulative percentage of HOGs by HOG genome frequency.
#
# Dotted line:
#   50% cumulative frequency reference.
#
# Required objects:
#   combBed_v5_complete.2
#   tandem_array_HOG_summary
# ============================================================


# ------------------------------------------------------------
# 1. Check required objects and columns
# ------------------------------------------------------------

if (!exists("combBed_v5_complete.2")) {
  stop("Object 'combBed_v5_complete.2' was not found.")
}

if (!exists("tandem_array_HOG_summary")) {
  stop("Object 'tandem_array_HOG_summary' was not found.")
}

required_combBed <- c("globHOG", "genome")

required_tandem <- c(
  "HOG_ID",
  "genome",
  "genes_in_array"
)

missing_combBed <- setdiff(
  required_combBed,
  colnames(combBed_v5_complete.2)
)

missing_tandem <- setdiff(
  required_tandem,
  colnames(tandem_array_HOG_summary)
)

if (length(missing_combBed) > 0) {
  stop(
    "Missing columns in combBed: ",
    paste(missing_combBed, collapse = ", ")
  )
}

if (length(missing_tandem) > 0) {
  stop(
    "Missing columns in tandem_array_HOG_summary: ",
    paste(missing_tandem, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2. Calculate HOG presence across all genomes
# ------------------------------------------------------------
#
# Extract only the HOG ID from globHOG.
#
# Example:
# "N0.HOG0000782 OG0000003"
#
# Becomes:
# "N0.HOG0000782"
# ------------------------------------------------------------

hog_presence <- combBed_v5_complete.2 %>%
  mutate(
    HOG_ID = str_extract(
      str_trim(as.character(globHOG)),
      "^\\S+"
    ),
    genome = str_trim(as.character(genome))
  ) %>%
  filter(
    !is.na(HOG_ID),
    HOG_ID != "",
    !is.na(genome),
    genome != ""
  ) %>%
  distinct(HOG_ID, genome) %>%
  count(
    HOG_ID,
    name = "number_of_genomes"
  )


# ------------------------------------------------------------
# 3. Prepare tandem-array information
# ------------------------------------------------------------
#
# Only qualifying tandem arrays are included:
# genes_in_array >= 2
#
# If a HOG has multiple arrays in one genome, retain the
# largest array for that HOG in that genome.
# ------------------------------------------------------------

hog_arrays <- tandem_array_HOG_summary %>%
  mutate(
    HOG_ID = str_trim(as.character(HOG_ID)),
    genome = str_trim(as.character(genome)),
    genes_in_array = as.numeric(genes_in_array)
  ) %>%
  filter(
    !is.na(HOG_ID),
    HOG_ID != "",
    !is.na(genome),
    genome != "",
    !is.na(genes_in_array),
    genes_in_array >= 2
  ) %>%
  group_by(HOG_ID, genome) %>%
  summarise(
    largest_array = max(genes_in_array),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 4. Calculate array-size variation within each HOG
# ------------------------------------------------------------
#
# Variation = Largest array - Smallest array
#
# Only genomes containing a qualifying tandem array
# are considered for this calculation.
# ------------------------------------------------------------

hog_variation <- hog_arrays %>%
  group_by(HOG_ID) %>%
  summarise(
    genomes_with_arrays = n_distinct(genome),
    smallest_array = min(largest_array),
    largest_array = max(largest_array),
    
    variation = largest_array - smallest_array,
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 5. Combine HOG presence and array-size variation
# ------------------------------------------------------------
#
# X-axis uses the total number of genomes containing
# each HOG, not just genomes with tandem arrays.
# ------------------------------------------------------------

plot_data <- hog_presence %>%
  inner_join(
    hog_variation,
    by = "HOG_ID"
  ) %>%
  mutate(
    variation_category = if_else(
      variation == 0,
      "Zero variation",
      "Variation > 0"
    )
  ) %>%
  filter(
    number_of_genomes >= 2
  )


# ------------------------------------------------------------
# 6. Calculate cumulative HOG frequency
# ------------------------------------------------------------

cumulative_data <- plot_data %>%
  count(
    number_of_genomes,
    name = "HOG_count"
  ) %>%
  arrange(number_of_genomes) %>%
  mutate(
    cumulative_HOGs = cumsum(HOG_count),
    
    total_HOGs = sum(HOG_count),
    
    cumulative_frequency = (
      cumulative_HOGs / total_HOGs
    ) * 100
  )


# ------------------------------------------------------------
# 7. Scale cumulative frequency for secondary axis
# ------------------------------------------------------------

max_variation <- max(
  plot_data$variation,
  na.rm = TRUE
)

if (max_variation == 0) {
  max_variation <- 1
}

cumulative_data <- cumulative_data %>%
  mutate(
    cumulative_y = (
      cumulative_frequency / 100
    ) * max_variation
  )


# ------------------------------------------------------------
# 8. Calculate 50% reference-line position
# ------------------------------------------------------------

reference_50 <- max_variation * 0.50


# ------------------------------------------------------------
# 9. Create the graph
# ------------------------------------------------------------

p <- ggplot() +
  
  # Variation points
  geom_point(
    data = plot_data,
    aes(
      x = number_of_genomes,
      y = variation,
      colour = variation_category
    ),
    size = 2.2,
    alpha = 0.60
  ) +
  
  # 50% cumulative frequency dotted reference line
  geom_hline(
    yintercept = reference_50,
    colour = "black",
    linetype = "dotted",
    linewidth = 0.9
  ) +
  
  # Cumulative frequency line
  geom_line(
    data = cumulative_data,
    aes(
      x = number_of_genomes,
      y = cumulative_y,
      group = 1
    ),
    colour = "black",
    linewidth = 0.9
  ) +
  
  # Cumulative frequency points
  geom_point(
    data = cumulative_data,
    aes(
      x = number_of_genomes,
      y = cumulative_y
    ),
    colour = "black",
    size = 1.8
  ) +
  
  # Point colours
  scale_colour_manual(
    name = "Array-size variation",
    values = c(
      "Variation > 0" = "#377EB8",
      "Zero variation" = "#E78AC3"
    )
  ) +
  
  # Left and right y-axes
  scale_y_continuous(
    name = "Variation in number of genes in a \n Tandem gene array in a HOG",
    
    limits = c(0, max_variation),
    
    breaks = pretty(
      c(0, max_variation),
      n = 6
    ),
    
    sec.axis = sec_axis(
      ~ . / max_variation * 100,
      name = "Cumulative Frequency of HOGs (%)",
      breaks = seq(0, 100, by = 10)
    )
  ) +
  
  # X-axis
  scale_x_continuous(
    name = "Number of Haplotypes present in a HOG",
    breaks = seq(5, 50, by = 5),
    limits = c(1.5, 50.5)
  ) +
  
  # Theme
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.title.y.right = element_text(
      angle = 270,
      margin = margin(l = 12)
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    )
  )


# ------------------------------------------------------------
# 10. Display graph
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 11. Export high-resolution graph
# ------------------------------------------------------------

# High-resolution PNG
ggsave(
  "tandem_array_variation_vs_HOG_presence.png",
  plot = p,
  width = 14,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# High-quality PDF (vector format)
ggsave(
  "tandem_array_variation_vs_HOG_presence.pdf",
  plot = p,
  width = 14,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)



# ------------------------------------------------------------
# 12. Export underlying tables
# ------------------------------------------------------------

write_csv(
  plot_data,
  "tandem_array_variation_by_HOG_presence.csv"
)

write_csv(
  cumulative_data,
  "tandem_array_cumulative_frequency.csv"
)


# ------------------------------------------------------------
# 13. Print summary and diagnostics
# ------------------------------------------------------------

cat("\n============================================\n")
cat("TANDEM-ARRAY HOG VARIATION SUMMARY\n")
cat("============================================\n")

cat(
  "Total HOGs present in combBed:",
  nrow(hog_presence),
  "\n"
)

cat(
  "HOGs with qualifying tandem arrays:",
  nrow(plot_data),
  "\n"
)

cat(
  "HOGs with zero array-size variation:",
  sum(plot_data$variation == 0),
  "\n"
)

cat(
  "HOGs with non-zero array-size variation:",
  sum(plot_data$variation > 0),
  "\n"
)

cat(
  "Maximum array-size variation:",
  max(plot_data$variation),
  "\n"
)

cat(
  "Maximum number of genomes containing a HOG:",
  max(plot_data$number_of_genomes),
  "\n"
)

cat(
  "50% reference line (left y-axis position):",
  reference_50,
  "\n"
)

cat("\nOutput files created:\n")
cat("- tandem_array_variation_vs_HOG_presence.png\n")
cat("- tandem_array_variation_vs_HOG_presence.pdf\n")
cat("- tandem_array_variation_by_HOG_presence.csv\n")
cat("- tandem_array_cumulative_frequency.csv\n")

cat("\n============================================\n")
