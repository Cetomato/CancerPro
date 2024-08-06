# Load libraries
library(data.table)
library(ComplexHeatmap)

# Read the pathology data from a TSV file
pgs <- fread("./pathology.tsv")

# Get the list of unique cancer types from the data
cancers <- names(table(pgs$Cancer))

# Extract prognostic favorable genes for each cancer type
prognostic_favorable <- lapply(cancers, function(cancer) 
  pgs[Cancer == cancer & !is.na(`prognostic - favorable`), .(Gene, `Gene name`, `Cancer`, `prognostic - favorable`)])
# Set prognosis to 1 for favorable genes
prognostic_favorable <- lapply(prognostic_favorable, function(x) x[, prognosis := 1])

# Select relevant columns and discard others for favorable genes
for (i in 1:length(prognostic_favorable)) {  
  prognostic_favorable[[i]] <- prognostic_favorable[[i]][, .(`Gene`, `Gene name`, `prognosis`)]
}

# Extract prognostic unfavorable genes for each cancer type
prognostic_unfavorable <- lapply(cancers, function(cancer) 
  pgs[Cancer == cancer & !is.na(`prognostic - unfavorable`), .(Gene, `Gene name`, `Cancer`, `prognostic - unfavorable`)])
# Set prognosis to -1 for unfavorable genes
prognostic_unfavorable <- lapply(prognostic_unfavorable, function(x) x[, prognosis := -1])

# Select relevant columns and discard others for unfavorable genes
for (i in 1:length(prognostic_unfavorable)) {  
  prognostic_unfavorable[[i]] <- prognostic_unfavorable[[i]][, .(`Gene`, `Gene name`, `prognosis`)]
}

# Combine favorable and unfavorable prognostic genes for each cancer type
prognostic_com <- list()
for (i in 1:length(cancers)) {
  prognostic_com[[i]] <- rbind(prognostic_favorable[[i]], prognostic_unfavorable[[i]])
}

# Rename the prognosis column to the respective cancer type with underscores
for (i in 1:length(cancers)) {
  setnames(prognostic_com[[i]], "prognosis", gsub(" ", "_", cancers[i]))
}

# Merge all the prognostic data into one data.table
prognostic_com_dt <- Reduce(function(x, y) merge(x, y, by = c("Gene", "Gene name"), all = TRUE), prognostic_com)

# Convert the data.table to a data.frame and set row names
prognostic_com_df <- as.data.frame(prognostic_com_dt)
rownames(prognostic_com_df) <- prognostic_com_df$`Gene name`
prognostic_com_df <- prognostic_com_df[, -c(1, 2)]

# Count the number of genes with non-missing values for each cancer type
cancer_gs <- apply(prognostic_com_df, 2, function(x) length(na.omit(x)))
cancer_gs

# Filter cancers with non-zero gene counts
cancer_sub <- names(cancer_gs[cancer_gs != 0])
prognostic_com_df <- prognostic_com_df[, cancer_sub]

# Sort genes by the number of cancers in which they are unfavorable
unfavorable_gene_sort <- sort(apply(prognostic_com_df, 1, function(x) sum(na.omit(unlist(x)) == -1)), decreasing = TRUE)
unfavorable_gene_sort[1:30]
sum(unfavorable_gene_sort > 3)  # Count of genes that are unfavorable in more than 3 cancers
sum(unfavorable_gene_sort >= 5)  # Count of genes that are unfavorable in 5 or more cancers
unfavorable_top <- unfavorable_gene_sort[1:sum(unfavorable_gene_sort >= 5)]
unfavorable_top

# Write the top unfavorable genes to a CSV file
write.table(data.frame(unfavorable_top), file = "./unfavorable_top.csv", row.names = TRUE, sep = ",")

# Sort genes by the number of cancers in which they are favorable
favorable_gene_sort <- sort(apply(prognostic_com_df, 1, function(x) sum(na.omit(x) == 1)), decreasing = TRUE)
favorable_gene_sort[1:30]
sum(favorable_gene_sort > 3)  # Count of genes that are favorable in more than 3 cancers
sum(favorable_gene_sort >= 5 )  # Count of genes that are favorable in 5 or more cancers
favorable_gene_sort[1:sum(favorable_gene_sort >= 5)]

# Write the favorable genes to a CSV file
write.table(data.frame(favorable_gene_sort[1:sum(favorable_gene_sort > 3)]), file = "favorable_gene_sort_3.csv", sep = ",")

# Count the number of favorable and unfavorable genes for each cancer type
favorable_n <- apply(prognostic_com_df, 2, function(x) sum(na.omit(x) == 1))
unfavorable_n <- apply(prognostic_com_df, 2, function(x) sum(na.omit(x) == -1))
biomarker_gene_count <- cbind(data.frame(favorable_n), data.frame(unfavorable_n))

# Write the biomarker gene counts to a CSV file
write.table(biomarker_gene_count, file = "biomarker_gene_count.csv", row.names = TRUE, sep = ",")

# Identify prostate cancer-specific poor prognosis genes
colnames(prognostic_com_df)
prostate_gs <- rownames(prognostic_com_df)[which(prognostic_com_df$prostate_cancer == -1)]

# Find genes that are poor prognosis in other cancers
other_minus1 <- apply(prognostic_com_df, 1, function(x) any(x[-14] == -1))
other_minus1[is.na(other_minus1)] <- FALSE

# Identify prostate-specific poor prognosis genes
prostate_specific <- setdiff(prostate_gs, rownames(prognostic_com_df)[other_minus1])
fwrite(data.frame(prostate_specific), file = "prostate_specific.csv")

# Identify endometrial cancer-specific poor prognosis genes
colnames(prognostic_com_df)
endometrial_gs <- rownames(prognostic_com_df)[which(prognostic_com_df$endometrial_cancer == -1)]

# Find genes that are poor prognosis in other cancers
other_minus1 <- apply(prognostic_com_df, 1, function(x) any(x[-5] == -1))
other_minus1[is.na(other_minus1)] <- FALSE

# Identify endometrial-specific poor prognosis genes
endometrial_specific <- setdiff(endometrial_gs, rownames(prognostic_com_df)[other_minus1])
fwrite(data.frame(endometrial_specific), file = "endometrial_specific.csv")

# Select genes that are either favorable or unfavorable in more than 3 cancers
genes_sub <- unique(c(names(favorable_gene_sort)[favorable_gene_sort > 3], names(unfavorable_gene_sort)[unfavorable_gene_sort > 3]))

# Create a subset of the prognostic data frame for these genes
prognostic_com_df_sub <- prognostic_com_df[genes_sub, ]
prognostic_com_df_sub[is.na(prognostic_com_df_sub)] <- 0

# Create a column annotation for the heatmap
colanno <- data.frame(Cancer = colnames(prognostic_com_df_sub))
rownames(colanno) <- colnames(prognostic_com_df_sub)

# Load the ComplexHeatmap library and plot the heatmap with advanced options
suppressMessages(library(ComplexHeatmap, verbose = FALSE))
library(pheatmap)
p1 <- ComplexHeatmap::pheatmap(
  as.matrix(prognostic_com_df_sub),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  cutree_cols = 5,
  cutree_rows = 5,
  color = hcl.colors(3, "Blue-Red", rev = FALSE),
  fontsize_number = 9,
  fontsize_row = 7,
  column_names_side = "top",
  angle_col = "90"
)

# Save the heatmap to a PDF file
dev.off()
pdf(file = "Fig_1_pan-cancer_prognostic_genes.pdf", width = 7, height = 26, onefile = TRUE)
plot(p1)
dev.off()
