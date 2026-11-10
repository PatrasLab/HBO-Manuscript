library(DESeq2)
library(sva)
library(readr)
library(dplyr)

####Pre-processing####

wd <- "/my/directory/file/path/"

#create counts matrix
wd <- "/mount/patras/clare/hbo-bulk-rna-seq/"

counts <- read_delim(file=paste0(
    wd, "mouse_counts.txt"),
             delim='\t')

counts <- counts |> select(!gene_chr:tf_family) |> select(!gene_id)

gene_names <- counts$gene_name

counts$gene_name <- NULL

counts_matrix <- as.matrix(counts[,1:ncol(counts)])
rownames(counts_matrix) <- gene_names

#generate metadata df
sample_names = colnames(counts_matrix[,1:ncol(counts_matrix)])

meta <- data.frame(
  row.names = sample_names,
  sex = case_when(
    grepl("MM_", sample_names) ~ "male",
    grepl("FM_", sample_names) ~ "female",
    TRUE ~ "unknown"
  ),
  treatment = case_when(
    grepl("_E_", sample_names) ~ "E",
    grepl("_H_", sample_names) ~ "H",
    grepl("_M_", sample_names) ~ "M",
    TRUE ~ "unknown"
  ),
  timepoint = case_when(
    grepl("_4h_", sample_names) ~ "4h",
    grepl("_24h_", sample_names) ~ "24h",
    TRUE ~ "unknown"
  )
)

table(meta$sex, meta$treatment) 


#create DESeq dataset

dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix), 
                                colData = meta, 
                                design = ~ 1) #don't care about design for this purpose

#filter out any sparse genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

####Normalization & PCA####

vst_data <- vst(dds, blind = TRUE)

head(vst_data)

vst_df <- as.data.frame(assay(vst_data))

vst_df$gene_name <- rownames(vst_df)

write.table(vst_df, paste0(
    wd, "mouse_uti_counts_norm_vst.tsv"),
             sep='\t', row.names = FALSE, col.names=TRUE, quote=FALSE)


#taking top 500 most variable genes (comparable to DESeq2 plotPCA)
mat <- assay(vst_data)
#mat <- assay(rlog_data)
ntop <- 500
top <- order(rowVars(mat), decreasing = TRUE)[1:ntop]

#run PCA with centering, not scaling (comparable to DESeq2 plotPCA)
pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)

#pull info
n_pc <- 6
scores <- as.data.frame(pca$x[, 1:n_pc])
scores$sample <- rownames(scores)
scores <- scores[, c("sample", paste0("PC", 1:n_pc))]

# variance explained
pct_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
print(round(pct_var[1:n_pc], 1))

# export
write.table(scores, paste0(
    wd, "mouse_uti_deseq2_PCA_coordinates.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

write.table(data.frame(PC = paste0("PC", 1:n_pc),
                       percent_variance = pct_var[1:n_pc]),
            paste0(wd, "mouse_uti_deseq2_PCA_variance_explained.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


#PCA coordinates with metadata
meta <- data.frame(
  sample = sample_names,
  sex = case_when(
    grepl("MM_", sample_names) ~ "male",
    grepl("FM_", sample_names) ~ "female",
    TRUE ~ "unknown"
  ),
  treatment = case_when(
    grepl("_E_", sample_names) ~ "E",
    grepl("_H_", sample_names) ~ "H",
    grepl("_M_", sample_names) ~ "M",
    TRUE ~ "unknown"
  ),
  timepoint = case_when(
    grepl("_4h_", sample_names) ~ "4h",
    grepl("_24h_", sample_names) ~ "24h",
    TRUE ~ "unknown"
  )
)

scores_meta <- inner_join(scores, meta, by="sample")

write.table(scores_meta, paste0(
    wd, "mouse_uti_deseq2_PCA_coordinates_with_metadata.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)


