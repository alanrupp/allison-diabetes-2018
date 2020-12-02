# - Analysis of count data from Allison et al. Diabetes 2018 ------------------

library(tidyverse)

# - Read in count data and metadata -------------------------------------------
count_files <- list.files("counts", full.names = TRUE)
samples <- str_extract(count_files, "Sample_[0-9]+")
counts <- map_dfc(count_files, ~ read_tsv(.x, skip = 4, col_names = FALSE))
counts <- select(counts, `X1...1`, starts_with("X2"))
counts <- set_names(counts, "gene_id", samples)
counts <- column_to_rownames(as.data.frame(counts), "gene_id")

metadata <- read.csv("samples.csv", row.names = "Sample_ID")
metadata$Pair <- letters[rep(seq(nrow(metadata)/2), each = 2)]
metadata$Cells <- factor(metadata$Cells, levels = c("Sup", "Bead"))

# - Enrichment ----------------------------------------------------------------
# looking for enrichment using only the WT PBS samples
wt_pbs_samples <- rownames(metadata)[
  str_detect(metadata$Treatment, "PBS") & metadata$Mutant == "WT"
]

controls <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts[, wt_pbs_samples],
  colData = metadata[wt_pbs_samples, ],
  ~ Pair + Cells
) %>%
  DESeq2::DESeq()

enrichment <- DESeq2::results(controls) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

expression <- DESeq2::fpm(controls) %>% as.data.frame()
wt_expression <- rowMeans(expression)

# - Regulation ----------------------------------------------------------------
# set up a single model to compare expression by given groups
metadata$Group <- paste(metadata$Mutant, metadata$Treatment)
bead_samples <- rownames(metadata)[metadata$Cells == "Bead"]
regulation <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts[, bead_samples],
  colData = metadata[bead_samples, ],
  ~ Group
) %>%
  DESeq2::DESeq()

# get results for given contrasts
contrasts <- list(
  "3hLep" = c("Group", "WT 3hrLeptin", "WT 3hrPBS"),
  "10hLep" = c("Group", "WT 10hrLeptin", "WT 10hrPBS"),
  "ob" = c("Group", "ob/ob 10hrPBS", "WT 10hrPBS"),
  "obLep" = c("Group", "ob/ob 10hrLeptin", "ob/ob 10hrPBS")
)
results <- map(contrasts, ~ DESeq2::results(regulation, contrast = .x))

# - Recapitulate Figure 1B ----------------------------------------------------
# put genes regulated in any conditions into a single matrix with log2FoldChange
de_genes <- map(results, ~ rownames_to_column(as.data.frame(.x), "gene_id")) %>%
  bind_rows() %>%
  filter(padj < 0.05) %>%
  filter(!duplicated(gene_id)) %>%
  pull(gene_id)

regulated <- do.call(cbind, map(results, ~ .x[de_genes, "log2FoldChange"]))

pheatmap::pheatmap(regulated)

# - Recapitulate the tables of fold-change values -----------------------------
# combine all enrichment, expression, and fold-change values into a single df
combined <- data.frame(
  "Expression" = wt_expression,
  "Enrichment" = enrichment$log2FoldChange
) %>%
  bind_cols(map_dfc(results, ~ .x[, "log2FoldChange"]))

# recapitulate Table 1
table1_genes <- c("Serpina3i", "Socs3", "Serpina3h", "Traf3ip3", "Irgm2",
                  "Prokr2", "Serpina3c", "Serpina3n", "Irf9", "Bcl3", "Vwa5a",
                  "Serpina3m", "Stat3", "Asb4", "Vwce", "Rasl11a", "Junb",
                  "Atf3", "Ly6a", "Vwf", "Acer2", "Drd1a", "1700024P16Rik")

# convert gene_id to gene_name
genes <- read_csv("tx2gene.csv")
genes <- filter(genes, !duplicated(gene_id)) %>% select(starts_with("gene"))

combined %>% rownames_to_column("gene_id") %>%
  left_join(select(genes, gene_name, gene_id), ., by = "gene_id") %>%
  select(-gene_id) %>%
  filter(gene_name %in% table1_genes) %>%
  mutate(gene_name = factor(gene_name, levels = table1_genes)) %>%
  arrange(gene_name)

# - Recapitulate Table 2 ------------------------------------------------------
table2_genes <- c("Socs3", "Atf3", "Junb", "Arid5a", "Etv6")

combined %>% rownames_to_column("gene_id") %>%
  left_join(select(genes, gene_name, gene_id), ., by = "gene_id") %>%
  select(gene_name, `3hLep`) %>%
  filter(gene_name %in% table2_genes) %>%
  mutate(gene_name = factor(gene_name, levels = table2_genes)) %>%
  arrange(gene_name)

# - Recapitulate Table 3 -----------------------------------------------------
table3_genes <- c("Agrp", "Apoa1", "Ghrh", "Nmb", "Npy", "Rbp4", "Nts",
                  "Cartpt", "Edn3", "Gal", "Kiss1", "Ucn", "Pomc")

combined %>% rownames_to_column("gene_id") %>%
  left_join(select(genes, gene_name, gene_id), ., by = "gene_id") %>%
  select(-gene_id, -Enrichment, -Expression) %>%
  filter(gene_name %in% table3_genes) %>%
  mutate(gene_name = factor(gene_name, levels = table3_genes)) %>%
  arrange(gene_name)

# - Write count data to file for GEO -----------------------------------------
counts %>% rownames_to_column("gene_id") %>%
  left_join(select(genes, gene_name, gene_id), ., by = "gene_id") %>%
  write.csv(., "counts.csv", row.names = FALSE)
