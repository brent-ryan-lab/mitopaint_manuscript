# Title: dimensionality reduction mitopaint data (mean per well) v2
# Step: 6
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 21-06-2026

# load packages ####
library(data.table)
library(Seurat)
library(tidyverse)
# set variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- "unintegrated"
dims_use <- 1:50
k_param <- 15
res <- 1
perplexity <- 20
max_iter <- 4000
# load data ####
# load data as df
if (redu_state == "redu") {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,"_redu.csv", sep = ""), 
      header = TRUE)
  )
} else {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,".csv", sep = ""), 
      header = TRUE)
  )
}
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load metadata as meta
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_meta_", integrate_state, ".csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(meta) <- meta$V1
meta$V1 <- NULL
# remove any columns with NA/ non finite values
df <- df[, colSums(!is.finite(as.matrix(df))) == 0, drop = FALSE]
# remove any rows with NA/ non finite values
df <- df[apply(df, 1, function(x) all(is.finite(x))), , drop = FALSE]
# put data in Seurat object df.seurat ####
# assign data and meta to seurat object
df.seurat <- CreateSeuratObject(
  # transpose zscore matrix
  counts = t(as.matrix(df)),
  meta.data = meta,
  # assay is saved as MP for mitopaint
  assay = "MP"
)
# make mitopaint data default assay
DefaultAssay(df.seurat) <- "MP"
# manually copy counts to data layer to circumvent v5 seurat structure syntax
df.seurat <- SetAssayData(
  df.seurat,
  assay = "MP",
  layer = "data",
  new.data = GetAssayData(
    df.seurat,
    assay = "MP",
    layer = "counts"
  )
)
# scale data
df.seurat <- ScaleData(
  df.seurat,
  features = rownames(df.seurat)
)
# run pca ####
# make pca container to hold pca results
pca <- list()
# run pca on seurat
df.seurat <- RunPCA(
  df.seurat,
  # use all profiling features for the PCA
  features = rownames(df.seurat),
  # set seed so results are reproducible
  seed.use = 42
)
# find pca nearest neighbors for seurat
df.seurat <- FindNeighbors(
  df.seurat,
  # dims_use dictates how many PC dimensions to use
  # if dims_use = 50, then all PCs are used
  dims = dims_use,
  # k_param dictates how many nearest neighbors each point should have
  # this is dependent of number of technical well replicates in experiment
  k.param = k_param,
  reduction = "pca",
  graph.name = c("pca_nn", "pca_snn")
)
# save nn graphs in graphs
graphs <- list()
graphs$pca_nn <- df.seurat@graphs[["pca_nn"]]
graphs$pca_snn <- df.seurat@graphs[["pca_snn"]]
# find pca clusters for seurat, based on pca nearest neighbors
df.seurat <- FindClusters(
  df.seurat,
  # resolution dictates how granular the clustering is
  # higher res = more clusters, lower res = less clusters
  # this is tuneable based on how many subpopulations are expected
  resolution = res,
  graph.name = "pca_snn",
  cluster.name = "PCA_NN"
)
# store pca embeddings (row names = wells, col names = PCs) in a dataframe
pca$embeddings <-  as.data.frame(
  Embeddings(df.seurat, "pca")
)
# store pca loadings (row names = features, col names = PCs) in a matrix
# larger absolute value = larger contribution to PC for given feature
pca$loadings <- Loadings(df.seurat, "pca")
# store pca top features in a nested list
# finds top 10 features with highest neg/ pos loading for each PC
# same feature can be heavily weighted in multiple PCs
# this is used to understand which features are contributing to which PC
pca$top_features <- list()
# for loop goes through all calculated PCs (50)
for (i in 1:ncol(pca$loadings)) {
  pca$top_features[[i]] <- TopFeatures(
    object = df.seurat[["pca"]],
    dim = i,
    # asks for 20 top features in total 
    nfeatures = 20,
    # returns for both positive and negative (10x neg, 10x pos)
    balanced = TRUE
  )
}
# visually inspect top features for PC1 
# View(as.data.frame(pca$top_features[[1]]))
# convert pca top features into a single long df
# imap_dfr() function loops over the list and binds all rows at the end
# function(pc_...) gives the contents of the list, and list number (PC)
pca$top_features_df <- imap_dfr(pca$top_features, function(pc_features, pc_i) {
  pc_col <- as.integer(pc_i)
  # first makes df for positive and negative weighted features, then bind_rows()
  positive_df <- data.frame(
    # feature name
    feature = pc_features[["positive"]],
    # which PC it is in top 10 for
    PC = pc_col,
    # where in top 10 features it ranks
    ranking = seq_along(pc_features[["positive"]]),
    # positively weighted feature
    sign = "positive",
    # loading value 
    loading = pca$loadings[pc_features[["positive"]], pc_col],
    row.names = NULL
  )
  negative_df <- data.frame(
    # feature name
    feature = pc_features[["negative"]],
    # which PC it is in top 10 for
    PC = pc_col,
    # where in top 10 features it ranks
    ranking = seq_along(pc_features[["negative"]]),
    # negatively weighted feature
    sign = "negative",
    # loading value 
    loading = pca$loadings[pc_features[["negative"]], pc_col],
    row.names = NULL
  )
  bind_rows(positive_df, negative_df)
})
# pca variance and sd
# sd reflects the spread of observations along given PC
# pca$sd is a double object
pca$sd <- df.seurat[["pca"]]@stdev
# variance is sd squared, to give variance explained by each pc
# pca$var_pc is a double object
pca$var_pc <- pca$sd^2
# pca$var is a data table summaraising pca$sd and pca$var_pc
pca$var <- data.frame(
  # which PC
  PC = paste0("PC", seq_along(pca$sd)),
  # sd of given PC
  SD = pca$sd,
  # variance explained by given pc
  Variance = pca$var_pc,
  # variance explained by given pc as a percent of total PCA variance
  Percent_Variance = 100 * pca$var_pc / sum(pca$var_pc),
  # cumulative percent, useful for knowing how many PCs to retain
  # and visualising PC plateau
  Cumulative_Percent_Variance = cumsum(100 * pca$var_pc / sum(pca$var_pc))
)
# View rough pca dimplot just to visually inspect data
DimPlot(
  df.seurat,
  reduction = "pca",
  # coloured by compound, concentration, batch and PCA clusters
  group.by = c("Compound",
               "Concentration",
               "Batch",
               "PCA_NN")
)
# run tsne ####
# make tsne container to hold tsne results
tsne <- list()
# run tsne on seurat
df.seurat <- RunTSNE(
  df.seurat,
  # tsne is NOT run on existing dimensionality reduction (eg pca)
  dims = NULL,
  reduction = NULL,
  # tsne is calculated from entire original feature space (all feature rows)
  features = rownames(df.seurat),
  # perplexity = how many nearby points to a point to preserve a relationship
  # lower perplexity = preserve local, higher perplexity = preserve global
  perplexity = perplexity,
  # max_iter = how many iterations tsne spends to optimise points
  # higher interations = more optimisation, but takes longer computing
  max_iter = max_iter,
  reduction.name = "tsne",
  # set seed so results are reproducible
  seed.use = 42
)
# find tsne nearest neighbors for seurat
df.seurat <- FindNeighbors(
  df.seurat,
  # dims dictates to use the two tsne dimensions
  dims = 1:2,
  # k_param dictates how many nearest neighbors each point should have
  # this is dependent of number of technical well replicates in experiment
  k.param = k_param,
  reduction = "tsne",
  graph.name = c("tsne_nn", "tsne_snn")
)
# save nn graphs in graphs
graphs$tsne_nn <- df.seurat@graphs[["tsne_nn"]]
graphs$tsne_snn <- df.seurat@graphs[["tsne_snn"]]
# find tsne clusters for seurat, based on tsne nearest neighbors
df.seurat <- FindClusters(
  df.seurat,
  # resolution dictates how granular the clustering is
  # higher res = more clusters, lower res = less clusters
  # this is tuneable based on how many subpopulations are expected
  resolution = res,
  graph.name = "tsne_snn",
  cluster.name = "tSNE_NN"
)
# store tsne embeddings (row names = wells, col names = tsne) in a dataframe
tsne$embeddings <- as.data.frame(
  Embeddings(df.seurat, "tsne")
)
# View rough tsne dimplot just to visually inspect data
DimPlot(
  df.seurat,
  reduction = "tsne",
  # coloured by compound, concentration, batch and tsne clusters
  group.by = c("Compound",
               "Concentration",
               "Batch",
               "tSNE_NN")
)
# run umap ####
# make umap container to hold umap results
umap <- list()
# run umap on seurat
df.seurat <- RunUMAP(
  df.seurat,
  # umap is NOT run on existing dimensionality reduction (eg pca)
  dims = NULL,
  reduction = NULL,
  # umap is calculated from entire original feature space (all feature rows)
  features = rownames(df.seurat),
  reduction.name = "umap",
  # set seed so results are reproducible
  seed.use = 42
)
# find umap nearest neighbors for seurat
df.seurat <- FindNeighbors(
  df.seurat,
  # dims dictates to use the two umap dimensions
  dims = 1:2,
  # k_param dictates how many nearest neighbors each point should have
  # this is dependent of number of technical well replicates in experiment
  k.param = k_param,
  reduction = "umap",
  graph.name = c("umap_nn", "umap_snn")
)
# save nn graphs in graphs
graphs$umap_nn <- df.seurat@graphs[["umap_nn"]]
graphs$umap_snn <- df.seurat@graphs[["umap_snn"]]
# find umap clusters for seurat, based on umap nearest neighbors
df.seurat <- FindClusters(
  df.seurat,
  # resolution dictates how granular the clustering is
  # higher res = more clusters, lower res = less clusters
  # this is tuneable based on how many subpopulations are expected
  resolution = res,
  graph.name = "umap_snn",
  cluster.name = "UMAP_NN"
)
# store umap embeddings (row names = wells, col names = umap) in a dataframe
umap$embeddings <- as.data.frame(
  Embeddings(df.seurat, "umap")
)
# View rough umap dimplot just to visually inspect data
DimPlot(
  df.seurat,
  reduction = "umap",
  # coloured by compound, concentration, batch and umap clusters
  group.by = c("Compound",
               "Concentration",
               "Batch",
               "UMAP_NN")
)
# save meta ####
write.csv(df.seurat@meta.data,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv", sep = "")
)
# save data ####
write.csv(pca$embeddings,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_embeddings.csv", sep = "")
)
write.csv(pca$var,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_var.csv", sep = "")
)
write.csv(pca$loadings,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_loadings.csv", sep = "")
)
write.csv(pca$top_features_df,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_top_features.csv", sep = "")
)
write.csv(tsne$embeddings,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_tsne_embeddings.csv", sep = "")
)
write.csv(umap$embeddings,
          paste(
            "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_embeddings.csv", sep = "")
)
# save nn and snn edge list
dir.create("data/processed/graphs", recursive = TRUE, showWarnings = FALSE)
iwalk(graphs, function(g, graph_name) {
  edge_df <- Matrix::summary(g)
  edge_df$from <- rownames(g)[edge_df$i]
  edge_df$to <- colnames(g)[edge_df$j]
  write.csv(
    edge_df[, c("from", "to", "x")],
    paste0(
      "data/processed/",
      file_name, "_",
      integrate_state, "_",
      redu_state, "_",
      graph_name, "_edges.csv"
    ),
    row.names = FALSE
  )
})
rm(list = ls())
