# Computational Assignment: ENSEMBL Gene Expression Analysis with R

# Overview
In this assignment, you will apply the R programming skills you've learned through the Data Carpentry curriculum to analyze and visualize gene expression data from ENSEMBL. You will download a dataset, perform statistical analysis, and create meaningful visualizations to interpret the results.

# Learning Objectives
By completing this assignment, you will:

* Practice accessing and downloading genomic data from public repositories
* Apply data manipulation techniques using dplyr
* Perform statistical analysis on gene expression data
* Create informative visualizations using ggplot2
* Interpret biological significance from computational results

# Instructions

# Part 1: Data Acquisition
Use the biomaRt package to access the ENSEMBL database, retrieve gene expression data for a set of genes across different tissues in humans, and create a well-structured dataframe from the retrieved data.

```
# Installation of required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("biomaRt", "ggplot2", "dplyr", "tidyr"))

# Loading packages
library(biomaRt)
library(ggplot2)
library(dplyr)
library(tidyr)

# Connect to ENSEMBL database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Example code to get you started - you will need to modify this
# Get a list of available attributes to understand what data you can retrieve
attributes <- listAttributes(ensembl)
head(attributes)

# Retrieve gene expression data for 20 genes involved in the immune response
# You need to decide which genes to analyze and which tissues to include
immune_genes <- c("IL6", "TNF", "IL1B", "IL10", "IFNG", 
                 "TLR4", "CXCL8", "IL2", "IL4", "CD4", 
                 "CD8A", "FOXP3", "CTLA4", "PD1", "NFKB1", 
                 "STAT1", "STAT3", "IRF3", "IRF7", "MX1")

# Get ENSEMBL IDs for these genes
gene_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "hgnc_symbol",
                 values = immune_genes,
                 mart = ensembl)

# Now retrieve expression data
# Note: You'll need to retrieve appropriate data from ENSEMBL
```
# Part 2: Data Preprocessing and Exploration
Clean and restructure your data as needed, create summary statistics for gene expression across tissues, identify genes with high variability across tissues, and handle missing values appropriately.

```
# Example code for data preprocessing - you will need to extend this
# Assuming you now have a dataframe called 'expression_data'

# Check for missing values
missing_values <- expression_data %>%
  summarize_all(~sum(is.na(.)))

# Basic summary statistics
expression_summary <- expression_data %>%
  group_by(gene) %>%
  summarize(
    mean_expression = mean(expression, na.rm = TRUE),
    median_expression = median(expression, na.rm = TRUE),
    min_expression = min(expression, na.rm = TRUE),
    max_expression = max(expression, na.rm = TRUE),
    variance = var(expression, na.rm = TRUE),
    std_dev = sd(expression, na.rm = TRUE)
  ) %>%
  arrange(desc(variance))  # Sort by variance to identify highly variable genes

# Print the top 5 most variable genes
head(expression_summary, 5)
```

# Part 3: Statistical Analysis
Perform appropriate statistical tests to compare gene expression between tissues, identify significantly differentially expressed genes, calculate correlations between gene expression patterns, and create a heatmap of gene expression correlation.

```
# Example code for statistical analysis
# You will need to extend this based on your specific dataset

# ANOVA to test for differences across tissues for each gene
anova_results <- expression_data %>%
  group_by(gene) %>%
  do(anova_result = summary(aov(expression ~ tissue, data = .)))

# Extract p-values from ANOVA results and adjust for multiple testing
# (This is simplified - you'll need to expand this code)
p_values <- sapply(anova_results$anova_result, function(x) x[[1]]$`Pr(>F)`[1])
adjusted_p <- p.adjust(p_values, method = "BH")

# Create a dataframe with significant genes
significant_genes <- data.frame(
  gene = anova_results$gene,
  p_value = p_values,
  adjusted_p = adjusted_p
) %>%
  filter(adjusted_p < 0.05)

# Calculate correlation matrix for gene expression
gene_expression_wide <- expression_data %>%
  pivot_wider(names_from = gene, values_from = expression)

correlation_matrix <- cor(gene_expression_wide[,-1], use = "complete.obs")
```

# Part 4: Data Visualization
Create boxplots showing expression distribution across tissues for top genes, generate a heatmap of gene expression patterns, create a PCA plot to visualize relationships between tissues based on gene expression, and make your visualizations publication-ready with proper labels and formatting.

```
# Example code for visualization - you will need to extend this

# Boxplot of expression levels for top 5 most variable genes
top_genes <- expression_summary$gene[1:5]

expression_data %>%
  filter(gene %in% top_genes) %>%
  ggplot(aes(x = tissue, y = expression, fill = tissue)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of Top Variable Genes Across Tissues",
       x = "Tissue",
       y = "Expression Level")

# Heatmap of gene expression
library(pheatmap)

# Prepare data for heatmap (you'll need to adjust this for your data)
expression_matrix <- expression_data %>%
  pivot_wider(names_from = tissue, values_from = expression) %>%
  column_to_rownames("gene")

pheatmap(expression_matrix,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Gene Expression Heatmap")

# PCA plot
pca_data <- expression_data %>%
  pivot_wider(names_from = gene, values_from = expression) %>%
  column_to_rownames("tissue")

pca_result <- prcomp(pca_data, scale = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$tissue <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, label = tissue)) +
  geom_point(size = 3, aes(color = tissue)) +
  geom_text(hjust = 0.5, vjust = -0.5) +
  theme_bw() +
  labs(title = "PCA of Gene Expression Across Tissues",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)"))
```

# Deliverables
1. R script containing all your code, well-documented with comments
2. R Markdown file that includes:
- Introduction explaining your research question and approach
- Methods describing your data acquisition and analysis steps
- Results section with statistics and visualizations
- Discussion interpreting your findings
- Conclusion summarizing key insights

# Evaluation Criteria
* Correct implementation of data retrieval from ENSEMBL (15%)
* Appropriate data preprocessing and exploration (20%)
* Sound statistical analysis (25%)
* Quality of visualizations (25%)
* Code organization, documentation, and reproducibility (15%)

# Submission Instructions
* Submit your R script and R Markdown files through Canvas in a single double compressed (tar + gzip = tarball).
* Late submissions will be subject to the policy in the syllabus.
