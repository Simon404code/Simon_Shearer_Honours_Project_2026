
library(tidyverse)
library(writexl)

# ============================================================
# GENE AND TANDEM-ARRAY COUNTS PER HAPLOTYPE
#
# Total genes:
#   Source: combBed
#
# Genes in tandem arrays:
#   Source: all_genes_in_tandem_HOGs_2
#
# Distinct tandem arrays:
#   Source: all_genes_in_tandem_HOGs_2
#
# Output: one row per haplotype
# ============================================================


# ------------------------------------------------------------
# 1. Check required data frames and columns
# ------------------------------------------------------------

if (!exists("combBed_v5_complete.2")) {
  stop("The object 'combBed' is not loaded.")
}

if (!exists("all_genes_in_tandem_HOGs_2")) {
  stop("The object 'all_genes_in_tandem_HOGs_2' is not loaded.")
}

required_combBed <- c("genome", "id")

required_tandem <- c(
  "HOG_ID",
  "Gene_ID",
  "Genome",
  "array_ID"
)

missing_combBed <- setdiff(required_combBed, names(combBed_v5_complete.2))

missing_tandem <- setdiff(
  required_tandem,
  names(all_genes_in_tandem_HOGs_2)
)

if (length(missing_combBed) > 0) {
  stop(
    "combBed is missing columns: ",
    paste(missing_combBed, collapse = ", ")
  )
}

if (length(missing_tandem) > 0) {
  stop(
    "all_genes_in_tandem_HOGs_2 is missing columns: ",
    paste(missing_tandem, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2. Count all annotated genes per haplotype
# ------------------------------------------------------------
# Count unique gene IDs in the full gene table.

total_genes <- combBed_v5_complete.2 %>%
  transmute(
    haplotype = str_trim(as.character(genome)),
    gene_id = str_trim(as.character(id))
  ) %>%
  filter(
    !is.na(haplotype),
    haplotype != "",
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct(haplotype, gene_id) %>%
  count(haplotype, name = "total_genes")


# ------------------------------------------------------------
# 3. Prepare the tandem-HOG gene table
# ------------------------------------------------------------
# Use only genes explicitly listed in
# all_genes_in_tandem_HOGs_2.

tandem_hog_genes <- all_genes_in_tandem_HOGs_2 %>%
  transmute(
    HOG_ID = str_trim(as.character(HOG_ID)),
    gene_id = str_trim(as.character(Gene_ID)),
    haplotype = str_trim(as.character(Genome)),
    array_id = str_trim(as.character(array_ID))
  ) %>%
  filter(
    !is.na(gene_id),
    gene_id != "",
    !is.na(haplotype),
    haplotype != "",
    !is.na(array_id),
    array_id != "",
    array_id != "NA"
  ) %>%
  distinct(haplotype, gene_id, array_id, .keep_all = TRUE)


# ------------------------------------------------------------
# 4. Count distinct genes in tandem arrays per haplotype
# ------------------------------------------------------------
# A gene is counted once per haplotype, even if it occurs
# more than once in the input table.

tandem_genes <- tandem_hog_genes %>%
  distinct(haplotype, gene_id) %>%
  count(
    haplotype,
    name = "genes_in_tandem_arrays"
  )


# ------------------------------------------------------------
# 5. Count distinct tandem arrays per haplotype
# ------------------------------------------------------------
# An array is uniquely identified by its haplotype and
# array_ID. Each array is counted once.

distinct_arrays <- tandem_hog_genes %>%
  distinct(haplotype, array_id) %>%
  count(
    haplotype,
    name = "distinct_tandem_arrays"
  )


# ------------------------------------------------------------
# 6. Combine all counts
# ------------------------------------------------------------
# Retain haplotypes with no qualifying tandem arrays.
# Their tandem-array counts are set to zero.

haplotype_summary <- total_genes %>%
  left_join(tandem_genes, by = "haplotype") %>%
  left_join(distinct_arrays, by = "haplotype") %>%
  mutate(
    genes_in_tandem_arrays = replace_na(
      genes_in_tandem_arrays, 0L
    ),
    distinct_tandem_arrays = replace_na(
      distinct_tandem_arrays, 0L
    ),
    percentage_genes_in_tandem_arrays = if_else(
      total_genes > 0,
      round(
        100 * genes_in_tandem_arrays / total_genes,
        2
      ),
      0
    )
  ) %>%
  arrange(haplotype) %>%
  select(
    haplotype,
    total_genes,
    genes_in_tandem_arrays,
    distinct_tandem_arrays,
    percentage_genes_in_tandem_arrays
  )


# ------------------------------------------------------------
# 7. Check results
# ------------------------------------------------------------

cat("\n========== SUMMARY ==========\n")
cat("Number of haplotypes:", nrow(haplotype_summary), "\n")

cat(
  "Total annotated genes across haplotypes:",
  sum(haplotype_summary$total_genes),
  "\n"
)

cat(
  "Genes in tandem HOG arrays:",
  sum(haplotype_summary$genes_in_tandem_arrays),
  "\n"
)

cat(
  "Distinct tandem arrays:",
  sum(haplotype_summary$distinct_tandem_arrays),
  "\n"
)

cat("\n========== RESULTS BY HAPLOTYPE ==========\n")
#print(haplotype_summary, n = Inf)


# ------------------------------------------------------------
# 8. Save results
# ------------------------------------------------------------

write_csv(
  haplotype_summary,
  "gene_and_tandem_HOG_array_counts_per_haplotype.csv"
)
cat(
  "CSV file saved to:",
  file.path(
    getwd(),
    "gene_and_tandem_HOG_array_counts_per_haplotype.csv"
  ),
  "\n"
)

View(haplotype_summary)
