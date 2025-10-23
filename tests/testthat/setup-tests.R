human_genes <- c("hsa:5594", "hsa:5894", "hsa:5604", "hsa:2002")
compounds_1 <- c("C00002", "C00008", "C00076")
compounds_2 <- c("C00022")
compounds_3 <- c("C00076")

res_genes <- data.frame(
  KEGGID = human_genes,
  log2FoldChange = rnorm(length(human_genes), mean = 0, sd = 1)
)
res_metabo_1 <- data.frame(
  KEGG = compounds_1,
  log2FC = rnorm(length(compounds_1), mean = 0, sd = 1)
)
res_metabo_2 <- data.frame(
  KEGG_ids = compounds_2,
  log2FoldChange = rnorm(length(compounds_2), mean = 0, sd = 1)
)
res_metabo_3 <- data.frame(
  KEGG_ids = compounds_3,
  log2FoldChange = rnorm(length(compounds_3), mean = 0, sd = 1)
)

de_results_list <- list(
  trans = list(
    de_table = res_genes,
    value_column = "log2FoldChange",
    feature_column = "KEGGID"
  ),
  metabo_1 = list(
    de_table = res_metabo_1,
    value_column = "log2FC",
    feature_column = "KEGG"
  ),
  metabo_2 = list(
    de_table = res_metabo_2,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  )
)

de_results_list_with_repetition <- list(
  genes = list(
    de_table = res_genes,
    value_column = "log2FoldChange",
    feature_column = "KEGGID"
  ),
  metabo_1 = list(
    de_table = res_metabo_1,
    value_column = "log2FC",
    feature_column = "KEGG"
  ),
  metabo_2 = list(
    de_table = res_metabo_2,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  ),
  metabo_3 = list(
    de_table = res_metabo_3,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  )
)

kgml_path <- system.file("extdata", "hsa04010.xml", package = "MOVIDA")

nodes_df_path <- system.file("extdata", "hsa04010_nodes.csv", package = "MOVIDA")
edges_df_path <- system.file("extdata", "hsa04010_edges.csv", package = "MOVIDA")

# kgml_path <- file.path("isnt", "extdata", "hsa04010.xml")

# nodes_df_path <- file.path("isnt", "extdata", "hsa04010_nodes.csv")
# edges_df_path <- file.path("isnt", "extdata", "hsa04010_edges.csv")
#  path: hsa04010
# write.table(nodes_df, "inst/extdata/hsa04010_nodes.csv",  sep = ";")
