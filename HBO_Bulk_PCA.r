library(DESeq2)
library(sva)
library(readr)
library(dplyr)

####Pre-processing####

wd <- "/my/directory/file/path/"

#create counts matrices
counts_df_51 <- read_delim(file=paste0(
    wd, "hbo051_gene_counts.csv"),
             delim=',')

counts_df_51 <- counts_df_51 |> select(!gene_chr:tf_family) |> select(!gene_id)

counts_df_24_41 <- read_delim(file=paste0(
    wd, "hbo024_041_gene_counts.csv"),
             delim=',')

#merge
counts_df <- full_join(counts_df_51, counts_df_24_41, by = "gene_name")

gene_names <- counts_df$gene_name
counts_df$gene_name <- NULL

counts_matrix <- as.matrix(counts_df[,2:ncol(counts_df)])
rownames(counts_matrix) <- gene_names

#create metadata table
sample_names = colnames(counts_matrix[,1:ncol(counts_matrix)])

meta <- data.frame(
  row.names = sample_names,
  line = case_when(
    grepl("H24", sample_names) ~ "H24",
    grepl("H41", sample_names) ~ "H41",
    grepl("H51", sample_names) ~ "H51",
    TRUE ~ "unknown"
  ),
  treatment = case_when(
    grepl("_EH", sample_names) ~ "EH",
    grepl("_E[0-9]", sample_names) ~ "EC",
    grepl("_U", sample_names) ~ "U",
    grepl("_H", sample_names) ~ "H",
    grepl("_M|_BL", sample_names) ~ "Media",
    TRUE ~ "unknown"
  )
)

#run DESeq / create DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix), 
                                colData = meta, 
                                design = ~ 1) #don't care about design for this purpose

#filter out any sparse genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

####Normalization & PCA per organoid line####

#HBO024

dds_sub <- dds[, dds$line %in% "H24"]
dds_sub <- dds_sub[, dds_sub$treatment %in% c("EC","EH","H","U")]

# re-filter genes for the subset, then transform
dds_sub <- dds_sub[rowSums(counts(dds_sub)) >= 10, ]

#run vst

vst_data <- vst(dds_sub, blind = TRUE)

head(vst_data)

vst_df <- as.data.frame(assay(vst_data))

vst_df$gene_name <- rownames(vst_df)

write.table(vst_df, paste0(
    wd, "hbo_counts_24_no_media_norm_vst.tsv"),
             sep='\t', row.names = FALSE, col.names=TRUE, quote=FALSE)
            

mat <- assay(vst_data)
ntop <- 500
top <- order(rowVars(mat), decreasing = TRUE)[1:ntop]

pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)

n_pc <- 6
scores <- as.data.frame(pca$x[, 1:n_pc])
scores$sample <- rownames(scores)
scores <- scores[, c("sample", paste0("PC", 1:n_pc))]

pct_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# export
write.table(scores, paste0(wd, "hbo_24_deseq2_no_media_PCA_coordinates.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

write.table(data.frame(PC = paste0("PC", 1:n_pc),
                       percent_variance = pct_var[1:n_pc]),
            paste0(wd, "hbo_24_deseq2_no_media_PCA_variance_explained.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#HBO041

dds_sub <- dds[, dds$line %in% "H41"]
dds_sub <- dds_sub[, dds_sub$treatment %in% c("EC","EH","H","U")]

# re-filter genes for the subset, then transform
dds_sub <- dds_sub[rowSums(counts(dds_sub)) >= 10, ]

#run vst

vst_data <- vst(dds_sub, blind = TRUE)

head(vst_data)

vst_df <- as.data.frame(assay(vst_data))

vst_df$gene_name <- rownames(vst_df)

write.table(vst_df, paste0(
    wd, "hbo_counts_41_no_media_norm_vst.tsv"),
             sep='\t', row.names = FALSE, col.names=TRUE, quote=FALSE)
            

mat <- assay(vst_data)
ntop <- 500
top <- order(rowVars(mat), decreasing = TRUE)[1:ntop]

pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)

n_pc <- 6
scores <- as.data.frame(pca$x[, 1:n_pc])
scores$sample <- rownames(scores)
scores <- scores[, c("sample", paste0("PC", 1:n_pc))]

pct_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# export
write.table(scores, paste0(wd, "hbo_41_deseq2_no_media_PCA_coordinates.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

write.table(data.frame(PC = paste0("PC", 1:n_pc),
                       percent_variance = pct_var[1:n_pc]),
            paste0(wd, "hbo_41_deseq2_no_media_PCA_variance_explained.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#HBO051

dds_sub <- dds[, dds$line %in% "H51"]
dds_sub <- dds_sub[, dds_sub$treatment %in% c("EC","EH","H","U")]

# re-filter genes for the subset, then transform
dds_sub <- dds_sub[rowSums(counts(dds_sub)) >= 10, ]

#run vst

vst_data <- vst(dds_sub, blind = TRUE)

head(vst_data)

vst_df <- as.data.frame(assay(vst_data))

vst_df$gene_name <- rownames(vst_df)

write.table(vst_df, paste0(
    wd, "hbo_counts_51_no_media_norm_vst.tsv"),
             sep='\t', row.names = FALSE, col.names=TRUE, quote=FALSE)
            

mat <- assay(vst_data)
ntop <- 500
top <- order(rowVars(mat), decreasing = TRUE)[1:ntop]

pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)

n_pc <- 6
scores <- as.data.frame(pca$x[, 1:n_pc])
scores$sample <- rownames(scores)
scores <- scores[, c("sample", paste0("PC", 1:n_pc))]

pct_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# export
write.table(scores, paste0(wd, "hbo_51_deseq2_no_media_PCA_coordinates.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

write.table(data.frame(PC = paste0("PC", 1:n_pc),
                       percent_variance = pct_var[1:n_pc]),
            paste0(wd, "hbo_51_deseq2_no_media_PCA_variance_explained.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
