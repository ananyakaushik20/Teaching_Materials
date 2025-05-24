# egdeR_voom_dge_from_file.R

# Load packages

library(limma)
library(edgeR)

# -----------------------------
# 1. Load data
# -----------------------------
counts <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv")

# Ensure sample order matches
metadata <- metadata[match(colnames(counts), metadata$sample), ]
group <- factor(metadata$group)

# -----------------------------
# 2. Create DGEList object
# -----------------------------
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# -----------------------------
# 3. Design matrix
# -----------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# -----------------------------
# 4. voom transformation
# -----------------------------
v <- voom(dge, design, plot = TRUE)

# -----------------------------
# 5. Fit linear model and contrasts
# -----------------------------
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(Treated - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 6. Get DE results
# -----------------------------
top_genes <- topTable(fit2, number = 50, adjust.method = "BH")
write.csv(top_genes, "differential_expression_results.csv")
print(head(top_genes, 10))
