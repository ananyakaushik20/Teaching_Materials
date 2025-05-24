# limma_dge_example.R

#Note: data used is simulated data 

# Load necessary package
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")
library(limma)

# -----------------------------
# 1. Simulate gene expression data
# -----------------------------
set.seed(42)

n_genes <- 1000
n_samples <- 6

# Simulated expression matrix: genes x samples
expression_data <- matrix(rnorm(n_genes * n_samples, mean = 5, sd = 1), 
                          nrow = n_genes, ncol = n_samples)
rownames(expression_data) <- paste0("Gene", 1:n_genes)
colnames(expression_data) <- paste0("Sample", 1:n_samples)

# -----------------------------
# 2. Assign sample groups
# -----------------------------
group <- factor(c("Control", "Control", "Control", "Treated", "Treated", "Treated"))

# -----------------------------
# 3. Add artificial differential expression
# -----------------------------
# Make first 50 genes upregulated in treated samples
expression_data[1:50, 4:6] <- expression_data[1:50, 4:6] + 2  # log2 fold change ~ 2

# -----------------------------
# 4. Create design matrix
# -----------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# -----------------------------
# 5. Fit linear model and apply contrasts
# -----------------------------
fit <- lmFit(expression_data, design)

contrast.matrix <- makeContrasts(Treated - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 6. Get top differentially expressed genes
# -----------------------------
top_genes <- topTable(fit2, number = 20, adjust.method = "BH")
print(top_genes)
