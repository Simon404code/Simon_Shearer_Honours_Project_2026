library(dplyr)
library(ggplot2)
library(ggrepel)

make_pie <- function(df, title, start = 0) {
  total <- sum(df$n)
  
  df <- df %>%
    mutate(
      Percentage = n / total * 100,
      Label = paste0(round(Percentage, 1), "%"),
      legend = factor(
        paste0(category, " (n = ", format(n, big.mark = ","), ")"),
        levels = paste0(category, " (n = ", format(n, big.mark = ","), ")")
      )
    )
  
  ggplot(df, aes(x = 1, y = n, fill = legend)) +
    geom_col(width = 1, colour = "white") +
    coord_polar(theta = "y", start = start, clip = "off") +
    geom_text_repel(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      show.legend = FALSE
    ) +
    xlim(0.5, 2) +
    theme_void() +
    labs(
      title = title,
      subtitle = paste0("Total: ", format(total, big.mark = ",")),
      fill = "Category"
    )
}

# 1. Tandem HOG distribution
df1 <- data.frame(
  category = c("shell", "core", "soft-core", "singletons"),
  n = c(10984, 692, 198, 487)
)
make_pie(df1, "Distribution of tandem array carrying HOGs")

# 2. All HOGs vs tandem HOGs
df2 <- data.frame(
  category = c("HOGs without tandem arrays", "HOGs with tandem arrays"),
  n = c(111577, 12361)
)
make_pie(df2, "Distribution of HOGs with tandem arrays")

# 3. All genes vs tandem genes
df3 <- data.frame(
  category = c("Genes not in tandem arrays", "Genes in tandem arrays"),
  n = c(1979054, 56206)
)
make_pie(df3, "Distribution of genes across the pangenome")