library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)


# load gene matrix from cellranger output
data_dir <- c('/mount/patras/clare/sc-seq/hbo/hbo1/run_cellranger_count/run_count_hbo1_65195/outs/filtered_feature_bc_matrix',
              '/mount/patras/clare/sc-seq/hbo/hbo1/run_cellranger_count/run_count_hbo1_65194/outs/filtered_feature_bc_matrix')
data <- Read10X(data.dir = data_dir)

# create seurat object
hbo1 <- CreateSeuratObject(counts = data, project = "hbo1", min.cells = 5, min.features = 500)
hbo1

#assign a treatment variable
hbo1$treatment.numeric <- Idents(hbo1)
x <- Idents(hbo1)
x <- dplyr::recode(x,`1`="Mock",`2`="UPEC")
hbo1$treatment <- x
Idents(hbo1) <- "treatment"

##Quality control##

#create a column that denotes the percent of transcripts from mitochondrial genes
hbo1[["percent.mt"]] <- PercentageFeatureSet(hbo1, pattern = "^MT-")
head(hbo1$percent.mt)

#same thing for ribosomal genes
hbo1[["percent.ribosomal"]] <- PercentageFeatureSet(hbo1,pattern="^RP[LS]")
head(hbo1$percent.ribosomal)

##Normalize data##

hbo1 <- SCTransform(hbo1, vars.to.regress = "percent.mt", verbose = TRUE)

##Visualize QC metrics##

#individual metrics
VlnPlot(hbo1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal"), ncol = 4, log=FALSE, alpha=0.4)

#relationships between metrics
FeatureScatter(hbo1, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=0.5,jitter=TRUE)
FeatureScatter(hbo1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=0.5,jitter=TRUE)


num_cells <- ncol(hbo1);
num_features <- nrow(hbo1);

print(paste(num_cells, "cells", "x", num_features, "features observed"))

#28241 cells x 25065 features observed


hbo1[[]] |>
as_tibble(
  rownames="Cell.Barcode"
) -> qc.metrics

head(qc.metrics)

#scatter
qc.metrics |>
  arrange(percent.mt) |>
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Counts x features") +
  geom_hline(yintercept = 750) +
  geom_hline(yintercept = 2000) 

#histogram 1
qc.metrics |>
  ggplot(aes(log10(nCount_RNA))) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("UMI counts histogram") +
  geom_vline(xintercept = log10(1000))

#histogram 2
qc.metrics |>
  ggplot(aes(log10(nFeature_RNA))) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Feature counts histogram") +
  geom_vline(xintercept = log10(2000))

##Final filtering##

#filter out cells with <200 transcripts, >4000 transcripts, or >15% mitochondrial DNA
hbo1 <- subset(hbo1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)


#check # cells remaining
num_cells <- ncol(hbo1)
num_features <- nrow(hbo1)

print(paste(num_cells, "cells", "x", num_features, "features observed"))

#19652 cells x 25065 features observed

##Cell cycle scoring##

CellCycleScoring(
  hbo1, 
  s.features = cc.genes.updated.2019$s.genes, 
  g2m.features = cc.genes.updated.2019$g2m.genes, 
  set.ident = FALSE
) -> hbo1

head(hbo1[[]])


##Dimensionality reduction & clustering##

VizDimLoadings(hbo1, dims = 1:15, reduction = "pca")
DimHeatmap(hbo1, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(hbo1)
#proceed with 9 PCs


hbo1 <- RunPCA(hbo1, dims= 1:9, verbose = TRUE)
hbo1 <- RunUMAP(hbo1, dims = 1:9, verbose = TRUE)


hbo1 <- FindNeighbors(hbo1, dims = 1:9, verbose = TRUE)


#test and visualize clustering levels
hbo1 <- FindClusters(algorithm=4,resolution =0.1, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment') 

hbo1 <- FindClusters(algorithm=4,resolution =0.2, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment')

hbo1 <- FindClusters(algorithm=4,resolution =0.3, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment')

hbo1 <- FindClusters(algorithm=4,resolution =0.4, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment')

hbo1 <- FindClusters(algorithm=4,resolution =0.5, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment')

hbo1 <- FindClusters(algorithm=4,resolution =0.6, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment')

hbo1 <- FindClusters(algorithm=4,resolution =0.7, hbo1, verbose = TRUE)
DimPlot(hbo1, label = TRUE, split.by='treatment')

#proceed with resolution 0.3
hbo1 <- FindClusters(algorithm=4,resolution =0.3, hbo1, verbose = TRUE)


##Cell cycle phases by cluster and by infection status##

#visualize cell cycle scores on UMAP
DimPlot(hbo1,reduction="umap", group.by = "Phase", pt.size=0.6, alpha=1)
DimPlot(hbo1,reduction="umap", group.by = "Phase", split.by='treatment' ,pt.size=0.6, alpha=1)

#stacked barplot by cluster
hbo1@meta.data |>
  group_by(seurat_clusters,Phase) |>
  count() |>
  group_by(seurat_clusters) |>
  mutate(percent=100*n/sum(n)) |>
  ungroup() |>
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() + theme(text=element_text(size=20)) +
  ggtitle("Percentage of cell cycle phases per cluster")

#same plot faceted by infection status
hbo1@meta.data |>
  group_by(seurat_clusters, Phase, treatment) |>
  count() |>
  group_by(seurat_clusters) |>
  mutate(percent=100*n/sum(n)) |>
  ungroup() |>
  ggplot(aes(x=treatment,y=percent, fill=Phase)) +
  geom_col() + theme(text=element_text(size=30)) + facet_wrap(~ seurat_clusters, ncol=4) +
  ggtitle("Percentage of cell cycle phases per cluster by infection")


##Find markers##

FindAllMarkers(hbo1, only.pos = FALSE)
hbo1.markers <- FindAllMarkers(hbo1, only.pos = FALSE)

#group by cluster, order by log fold change, filter by p value
hbo1.markers.grouped.filt <- hbo1.markers |>
    group_by(cluster) |>
    dplyr::filter(p_val_adj > 0.05) |>
    dplyr::arrange(avg_log2FC)


##Expression of canonical bladder marker genes##

FeaturePlot(hbo1, features = c("KRT18", "KRT19", "KRT5", "KRT17", "KRT13", "KRT20", "UPK1A", "UPK1B",
    "UPK3B","UPK3A", "UPK2", "VIM", "S100A4", "COL3A1", "COL1A1", "COL1A2", "ACTA2",
    "TAGLN","DES", "CNN1", "ACTG2", "TPM2", "SELE", "PECAM1", "VCAM1", "CDH5",
    "LYZ","MS4A7", "CD14", "CD209", "CD3D", "MZB1", "CD79A", "GPM6A"),order=TRUE)
#The following requested variables were not found: VIM, COL3A1, COL1A2, DES, VCAM1, MS4A7, CD209 

#Markers with very little expression are KRT5, COL1A1, ACTA2, TAGLN, CNN1, ACTG2, TPM2, SELE, PECAM1, CDH5, LYZ, CD14, CD3D, MZB1, CD79A, GPM6A

#version that excludes the markers that are not present or are barely present
VlnPlot(hbo1, features = c("KRT18", "KRT19", "KRT5", "KRT17", "KRT13", "KRT20", "UPK1A", "UPK1B",
    "UPK3B","UPK3A", "UPK2", "S100A4"), ncol=6)


##Assign cell types##


new.cluster.ids <- c("Basal", "Non-proliferating", "Umbrella", 
                     "Non-proliferating", "Intermediate", "Intermediate",
                     "Umbrella", "Basal Progenitor")
names(new.cluster.ids) <- levels(hbo1)
hbo1 <- RenameIdents(hbo1, new.cluster.ids)

##Fig 1H## ##UMAP with cell types##

basalcolor <- "#FFCA99"
nonprocolor <- "#BA68C8"
umbrellacolor <- "#85B22C"
intermediatecolor <- "#51C3CC"
basalprocolor <- "#FF8E32"

p <- DimPlot(hbo1, reduction = "umap", label = FALSE, pt.size = 0.5, alpha=0.7, label.size=5, 
             cols=c(basalcolor, nonprocolor, umbrellacolor, intermediatecolor, basalprocolor)) + NoLegend()
p + annotate("label", x=2.9,y=7,label="Basal Progenitor",
                          fill="#FF8E32",size = 6,
                          label.size = 0,                 
                          label.r = unit(0.15, "lines")) +
    annotate("label", x=4.7,y=4.65,label="Non-proliferating",
                          fill="#BA68C8",size = 6,
                          label.size = 0,                 
                          label.r = unit(0.15, "lines")) +
    annotate("label", x=-5,y=5.1,label="Basal",
                          fill="#FFCA99",size = 6,
                          label.size = 0,                 
                          label.r = unit(0.15, "lines")) +
    annotate("label", x=4.6,y=-7,label="Umbrella",
                          fill="#85B22C",size = 6,
                          label.size = 0,                 
                          label.r = unit(0.15, "lines")) +
    annotate("label", x=-5,y=-5,label="Intermediate",
                          fill="#51C3CC",size = 6,
                          label.size = 0,                 
                          label.r = unit(0.15, "lines"))


##Fig 1I-L## ##UMAPs with expression##

#MKI67
FeaturePlot(hbo1, features = "MKI67", cols = c("#E6E6E3", "#004C54"), order=TRUE,
           pt.size=0.5, alpha=0.7)

#KRT20
FeaturePlot(hbo1, features = "KRT20",   cols = c("#E6E6E3", "#004C54"), order=TRUE,
           pt.size=0.5, alpha=0.7)

#UPK1A
FeaturePlot(hbo1, features = "UPK1A",   cols = c("#E6E6E3", "#004C54"), order=TRUE,
           pt.size=0.5, alpha=0.7)

#KRT19
FeaturePlot(hbo1, features = "KRT19",   cols = c("#E6E6E3", "#004C54"), order=TRUE,
           pt.size=0.5, alpha=0.7)


##Check IL17-C expression##

FetchData(hbo1, vars=c("IL17C", "seurat_clusters"), assay="SCT", layer="data") |> head(,50)
#The following requested variables were not found: IL17C 


##DEGs within clusters, infected over mock##

hbo1$celltype <- Idents(hbo1)
hbo1$celltype.inf <- paste(hbo1$celltype, hbo1$treatment, sep = " ")
Idents(hbo1) <- "celltype.inf"

#Basal cell infection response
basal.infect.response <- FindMarkers(hbo1, ident.1 = "Basal UPEC", ident.2 = "Basal Mock", verbose = TRUE)
basal.infect.response.sort.pos <- basal.infect.response |> arrange(desc(avg_log2FC)) |> filter(p_val_adj<=0.05)
basal.infect.response.sort.neg <- basal.infect.response |> arrange(avg_log2FC) |> filter(p_val_adj<=0.05)
head(basal.infect.response.sort.pos, n = 15)
head(basal.infect.response.sort.neg, n = 15)

basal.top.genes.pos <- rownames(basal.infect.response.sort.pos)[1:5]
basal.top.genes.neg <- rownames(basal.infect.response.sort.neg)[1:5]

basal.top10.pos.neg <- c(basal.top.genes.pos, basal.top.genes.neg)


#Basal Proliferating infection response
basalpro.infect.response <- FindMarkers(hbo1, ident.1 = "Basal Progenitor UPEC", ident.2 = "Basal Progenitor Mock", verbose = TRUE)
basalpro.infect.response.sort.pos <- basalpro.infect.response |> arrange(desc(avg_log2FC)) |> filter(p_val_adj<=0.05)
basalpro.infect.response.sort.neg <- basalpro.infect.response |> arrange(avg_log2FC) |> filter(p_val_adj<=0.05)
head(basalpro.infect.response.sort.pos, n = 15)
head(basalpro.infect.response.sort.neg, n = 15)

basalpro.top.genes.pos <- rownames(basalpro.infect.response.sort.pos)[1:5]
basalpro.top.genes.neg <- rownames(basalpro.infect.response.sort.neg)[1:5]

basalpro.top10.pos.neg <- c(basalpro.top.genes.pos, basalpro.top.genes.neg)

#Intermediate cells infection response
intermediate.infect.response <- FindMarkers(hbo1, ident.1 = "Intermediate UPEC", ident.2 = "Intermediate Mock", verbose = TRUE)
intermediate.infect.response.sort.pos <- intermediate.infect.response |> arrange(desc(avg_log2FC)) |> filter(p_val_adj<=0.05)
intermediate.infect.response.sort.neg <- intermediate.infect.response |> arrange(avg_log2FC) |> filter(p_val_adj<=0.05)
head(intermediate.infect.response.sort.pos, n = 15)
head(intermediate.infect.response.sort.neg, n = 15)

intermediate.top.genes.pos <- rownames(intermediate.infect.response.sort.pos)[1:5]
intermediate.top.genes.neg <- rownames(intermediate.infect.response.sort.neg)[1:5]

intermediate.top10.pos.neg <- c(intermediate.top.genes.pos, intermediate.top.genes.neg)

#Umbrella infection response
umbrella.infect.response <- FindMarkers(hbo1, ident.1 = "Umbrella UPEC", ident.2 = "Umbrella Mock", verbose = TRUE)
umbrella.infect.response.sort.pos <- umbrella.infect.response |> arrange(desc(avg_log2FC)) |> filter(p_val_adj<=0.05)
umbrella.infect.response.sort.neg <- umbrella.infect.response |> arrange(avg_log2FC) |> filter(p_val_adj<=0.05)
head(umbrella.infect.response.sort.pos, n = 15)
head(umbrella.infect.response.sort.neg, n = 15)

umbrella.top.genes.pos <- rownames(umbrella.infect.response.sort.pos)[1:5]
umbrella.top.genes.neg <- rownames(umbrella.infect.response.sort.neg)[1:5]

umbrella.top10.pos.neg <- c(umbrella.top.genes.pos, umbrella.top.genes.neg)

#Non-proliferating cells infection response
nonpro.infect.response <- FindMarkers(hbo1, ident.1 = "Non-proliferating UPEC", ident.2 = "Non-proliferating Mock", verbose = TRUE)
nonpro.infect.response.sort.pos <- nonpro.infect.response |> arrange(desc(avg_log2FC)) |> filter(p_val_adj<=0.05)
nonpro.infect.response.sort.neg <- nonpro.infect.response |> arrange(avg_log2FC) |> filter(p_val_adj<=0.05)
head(nonpro.infect.response.sort.pos, n = 15)
head(nonpro.infect.response.sort.neg, n = 15)

nonpro.top.genes.pos <- rownames(nonpro.infect.response.sort.pos)[1:5]
nonpro.top.genes.neg <- rownames(nonpro.infect.response.sort.neg)[1:5]

nonpro.top10.pos.neg <- c(nonpro.top.genes.pos, nonpro.top.genes.neg)

#combine all results
allcelltypes.top.genes.pos.neg <- c(basal.top10.pos.neg,basalpro.top10.pos.neg,
                                    intermediate.top10.pos.neg,umbrella.top10.pos.neg,
                                    nonpro.top10.pos.neg)
allcelltypes.top.genes.pos.neg


##Fig 2G data analysis & processing##

hbo1 <- ScaleData(hbo1, features = rownames(hbo1),assay="SCT")

avg.exp <- AverageExpression(hbo1, assays="SCT", layer="scale.data",
                             group.by=c("celltype","treatment"),
                             features= allcelltypes.top.genes.pos.neg,
                             verbose = TRUE)

avg.exp.mat <- as.matrix(avg.exp$SCT)

avg.exp.mat <- matrix(
  unlist(avg.exp.mat),
  nrow = nrow(avg.exp.mat),
  dimnames = dimnames(avg.exp.mat)
)

#reorder genes in the same order as the original gene list, retaining duplicate gene names
avg.exp.mat <- avg.exp.mat[allcelltypes.top.genes.pos.neg, , drop = FALSE]

#reorder matrix
ct_order <- c("Basal_Mock", "Basal_UPEC",
              "Basal Progenitor_Mock", "Basal Progenitor_UPEC",
               "Intermediate_Mock", "Intermediate_UPEC",
                      "Non-proliferating_Mock", "Non-proliferating_UPEC",
                      "Umbrella_Mock", "Umbrella_UPEC")
avg.exp.mat <- avg.exp.mat[, ct_order]


##Supplemental Table 2##

findmarkers_list <- list(
  "Basal" = basal.infect.response,
    "Basal Progenitor" = basalpro.infect.response,
  "Intermediate" = intermediate.infect.response,
    "Non-proliferating" = nonpro.infect.response,
    "Umbrella" = umbrella.infect.response
)

annot_mat <- matrix(
  "", 
  nrow = nrow(avg.exp.mat), 
  ncol = ncol(avg.exp.mat),
  dimnames = dimnames(avg.exp.mat)
)


for (i in seq_len(nrow(avg.exp.mat))) {
  gene <- allcelltypes.top.genes.pos.neg[i]
    
  for (j in seq_len(ncol(avg.exp.mat))) {
    colname <- colnames(avg.exp.mat)[j]
    
    if (grepl("_UPEC$", colname)) {
      celltype <- sub("_UPEC$", "", colname)
      
      if (celltype %in% names(findmarkers_list)) {
        df <- findmarkers_list[[celltype]]
        
        if (gene %in% rownames(df)) {
          pval <- df[gene, "p_val_adj"]
            annot_mat[i,j] <- pval

          }
        }
      }
    }
  }

mock_cols <- grep("_Mock$", colnames(annot_mat), value = TRUE)
annot_mat[, mock_cols] <- ""


##Supplemental Table 1##

#log1p corrected counts
avg.exp.ct.df <- as.data.frame(AverageExpression(hbo1, assays="SCT", layer="data",
                             group.by="celltype",
                             features= rownames(hbo1),
                             verbose = TRUE))


#scaled expression
avg.exp.ct.df.scale <- as.data.frame(AverageExpression(hbo1, assays="SCT", layer="scale.data",
                             group.by= "celltype",
                             features= rownames(hbo1),
                             verbose = TRUE))