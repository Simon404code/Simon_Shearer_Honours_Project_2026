# ============================================================
# SECOND GRAPH: STACKED BAR CHART
# Zero-variation vs non-zero-variation HOGs
# ============================================================

stacked_data <- plot_data %>%
  mutate(
    variation_category = if_else(
      variation == 0,
      "Zero variation",
      "Variation"
    )
  ) %>%
  count(number_of_genomes, variation_category) %>%
  complete(
    number_of_genomes = 2:50,
    variation_category = c("Zero variation", "Variation"),
    fill = list(n = 0)
  ) %>%
  group_by(number_of_genomes) %>%
  mutate(
    total_HOGs = sum(n),
    label_y = total_HOGs
  ) %>%
  ungroup()

p_stacked <- ggplot(
  stacked_data,
  aes(
    x = factor(number_of_genomes),
    y = n,
    fill = variation_category
  )
) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = if_else(n > 0, as.character(n), "")),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    colour = "black"
  ) +
  geom_text(
    data = stacked_data %>%
      distinct(number_of_genomes, total_HOGs),
    aes(
      x = factor(number_of_genomes),
      y = total_HOGs,
      label = total_HOGs
    ),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 2.8
  ) +
  scale_fill_manual(
    values = c(
      "Zero variation" = "#E78AC3",
      "Variation" = "#377EB8"
    )
  ) +
  scale_x_discrete(
    breaks = as.character(2:50),
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Number of genomes containing the HOG",
    y = "Number of HOGs",
    fill = "Variation category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 8
    ),
    legend.position = "top",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

# Display the second graph
print(p_stacked)

# Export the second graph
ggsave(
  "tandem_array_zero_vs_nonzero_variation_stacked.png",
  plot = p_stacked,
  width = 14,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "tandem_array_zero_vs_nonzero_variation_stacked.pdf",
  plot = p_stacked,
  width = 14,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)