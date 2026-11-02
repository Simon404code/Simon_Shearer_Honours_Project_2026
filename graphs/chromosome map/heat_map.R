
library(tidyverse)
library(viridis)

# ==================================================
# 1. CHROMOSOME LENGTHS (bp)
# ==================================================

chromosome_lengths <- tribble(
  ~chr,     ~length_bp,
  "Chr01",  44764794,
  "Chr02",  57174188,
  "Chr03",  64190686,
  "Chr04",  37932624,
  "Chr05",  67041112,
  "Chr06",  55026902,
  "Chr07",  58892594,
  "Chr08",  68844007,
  "Chr09",  36944981,
  "Chr10",  38023239,
  "Chr11",  41944707
)

# Bin size: 400 kb
bin_size <- 400000

output_dir <- "tandem_gene_chromosome_karyotype"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ==================================================
# 2. CHECK INPUT
# ==================================================

if (!exists("tandem_array_HOG_summary")) {
  stop(
    "The object 'tandem_array_HOG_summary' was not found. ",
    "Load the dataset before running this script."
  )
}

required_columns <- c(
  "chr", "genome", "arrayID",
  "start", "end", "genes_in_array"
)

missing_columns <- setdiff(
  required_columns,
  names(tandem_array_HOG_summary)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

cat("Chromosomes in the original dataset:\n")
print(sort(unique(as.character(tandem_array_HOG_summary$chr))))

# ==================================================
# 3. CLEAN AND FILTER THE DATA
# ==================================================

data <- tandem_array_HOG_summary %>%
  mutate(
    chr = str_to_title(
      str_trim(as.character(chr))
    ),
    chr = str_replace(
      chr,
      "^Chr(\\d)$",
      "Chr0\\1"
    ),
    start = suppressWarnings(as.numeric(start)),
    end = suppressWarnings(as.numeric(end)),
    genes_in_array = suppressWarnings(
      as.numeric(genes_in_array)
    )
  )

# Keep only the 11 chromosomes in the length table.
# Scaffolds and unplaced sequences are excluded.

excluded_chromosomes <- setdiff(
  unique(na.omit(data$chr)),
  chromosome_lengths$chr
)

cat("\nChromosomes/scaffolds excluded:\n")
print(excluded_chromosomes)

data <- data %>%
  filter(
    chr %in% chromosome_lengths$chr,
    !is.na(genome),
    !is.na(arrayID),
    !is.na(start),
    !is.na(end),
    !is.na(genes_in_array),
    genes_in_array > 0
  ) %>%
  mutate(
    midpoint = (start + end) / 2
  ) %>%
  distinct(
    genome, chr, arrayID,
    .keep_all = TRUE
  ) %>%
  left_join(
    chromosome_lengths,
    by = "chr"
  ) %>%
  filter(
    start >= 1,
    end >= start,
    midpoint <= length_bp
  )

if (nrow(data) == 0) {
  stop(
    "No valid tandem arrays remain after filtering. ",
    "Check chromosome names and coordinate columns."
  )
}

cat("\nNumber of tandem arrays retained:", nrow(data), "\n")
cat(
  "Chromosomes represented:",
  n_distinct(data$chr),
  "\n"
)

# ==================================================
# 4. COUNT TANDEM GENES IN EACH 400 KB BIN
# ==================================================

binned_counts <- data %>%
  mutate(
    bin_start = (
      floor((midpoint - 1) / bin_size) *
        bin_size
    ) + 1
  ) %>%
  group_by(chr, bin_start) %>%
  summarise(
    tandem_gene_count = sum(genes_in_array),
    number_of_arrays = n(),
    .groups = "drop"
  )

# ==================================================
# 5. GENERATE EVERY BIN, INCLUDING EMPTY BINS
# ==================================================

all_bins <- chromosome_lengths %>%
  rowwise() %>%
  mutate(
    bin_start = list(
      seq(1, length_bp, by = bin_size)
    )
  ) %>%
  ungroup() %>%
  unnest(bin_start) %>%
  mutate(
    bin_end = pmin(
      bin_start + bin_size - 1,
      length_bp
    ),
    bin_midpoint = (bin_start + bin_end) / 2,
    position_Mb = bin_midpoint / 1e6,
    bin_width_Mb = (bin_end - bin_start + 1) / 1e6
  ) %>%
  select(
    chr, length_bp,
    bin_start, bin_end,
    bin_midpoint, position_Mb,
    bin_width_Mb
  )

plot_data <- all_bins %>%
  left_join(
    binned_counts,
    by = c("chr", "bin_start")
  ) %>%
  mutate(
    tandem_gene_count = replace_na(
      tandem_gene_count, 0
    ),
    number_of_arrays = replace_na(
      number_of_arrays, 0L
    ),
    chr = factor(
      chr,
      levels = rev(chromosome_lengths$chr)
    )
  )

# Export the complete binned dataset.
write_csv(
  plot_data,
  file.path(
    output_dir,
    "tandem_gene_counts_per_400kb_bin.csv"
  )
)
# ==================================================
# 6. CREATE THE KARYOTYPE-STYLE HEATMAP
# ==================================================

max_position_Mb <- max(
  chromosome_lengths$length_bp
) / 1e6

karyotype_plot <- ggplot(
  plot_data,
  aes(
    x = position_Mb,
    y = chr,
    fill = tandem_gene_count
  )
) +
  geom_tile(
    aes(width = bin_width_Mb),
    height = 0.78
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    trans = "sqrt",
    name = "Tandem genes\nper 400 kb bin",
    na.value = "grey90"
  ) +
  scale_x_continuous(
    name = "Genomic position (Mb)",
    limits = c(0, max_position_Mb),
    breaks = seq(
      0,
      ceiling(max_position_Mb / 10) * 10,
      by = 10
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    name = "Chromosome",
    drop = FALSE,
    expand = expansion(add = 0.5)
  ) +
  labs(
    title = "Distribution of tandem genes across chromosomes"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold",
      colour = "black",
      size = 10
    ),
    axis.text.x = element_text(
      colour = "black"
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,
      size = 9
    ),
    plot.margin = margin(
      t = 12, r = 15, b = 12, l = 12
    )
  )

# Display the plot in RStudio.
print(karyotype_plot)
# ==================================================
# 7. SAVE THE KARYOTYPE
# ==================================================

ggsave(
  filename = file.path(
    output_dir,
    "tandem_gene_karyotype.png"
  ),
  plot = karyotype_plot,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "tandem_gene_karyotype.pdf"
  ),
  plot = karyotype_plot,
  width = 14,
  height = 6,
  units = "in",
  bg = "white"
)

# ==================================================
# 8. REPORT OUTPUTS
# ==================================================

cat(
  "\nKaryotype heatmap complete.\n",
  "Output directory: ",
  normalizePath(output_dir),
  "\nFiles generated:\n",
  "- tandem_gene_karyotype.png\n",
  "- tandem_gene_karyotype.pdf\n",
  "- tandem_gene_counts_per_400kb_bin.csv\n"
)
