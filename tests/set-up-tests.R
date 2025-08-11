# Helper function to create fake row names based on type
generate_row_names <- function(type, n, name) {
  switch(type,
    symbol = paste0(name, "_GENE", seq_len(n)),
    ensemblid = paste0("ENSG", sprintf("%011d", seq_len(n))),
    chebiid = paste0("CHEBI:", 10000 + seq_len(n)),
    uniprotid = paste0("P", sprintf("%05d", 10000 + seq_len(n))),
    inchikey = paste0("INCHIKEY_", substr(digest::md5sum(as.character(seq_len(n))), 1, 14)),
    paste0(name, "_gene", seq_len(n)) # default
  )
}

# Helper function to create a fake SummarizedExperiment object
create_fake_se <- function(nrow, ncol, name, rowname_type = "symbol") {
  assay <- stats::rnorm(nrow * ncol)
  assay <- matrix(assay, nrow = nrow, ncol = ncol)
  row_names <- generate_row_names(rowname_type, nrow, name)
  rowData <- S4Vectors::DataFrame(gene_id = row_names)
  rownames(rowData) <- row_names  # Set gene_id as rownames
  # Add a fake group column to colData
  sample_names <- paste0(name, "_sample", seq_len(ncol))
  group <- rep(c("A", "B"), length.out = ncol)
  colData <- S4Vectors::DataFrame(
    sample_id = sample_names,
    group = group
  )
  rownames(colData) <- sample_names  # Set sample_id as rownames
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay), rowData = rowData, colData = colData)
}

# Create three fake SE objects with different row name types
se_trans <- create_fake_se(10, 5, "SE1", rowname_type = "ensemblid")
se_prot <- create_fake_se(8, 4, "SE2", rowname_type = "uniprotid")
se_metabo <- create_fake_se(12, 6, "SE3", rowname_type = "chebiid")

# Create fake differential expression (DE) results
create_fake_de_results <- function(se, name, return_type = c("DESeqResults", "MArrayLM")) {

  return_type <- match.arg(return_type)
  n <- SummarizedExperiment::nrow(se)
  res_df <- data.frame(
    baseMean = stats::runif(n, min = 10, max = 1000),
    log2FoldChange = stats::rnorm(n, mean = 0, sd = 2),
    lfcSE = stats::runif(n, min = 0.1, max = 1),
    stat = stats::rnorm(n),
    pvalue = stats::runif(n, min = 0, max = 1),
    padj = stats::p.adjust(stats::runif(n), method = "BH")
  )
  rownames(res_df) <- SummarizedExperiment::rowData(se)$gene_id

  if (return_type == "DESeqResults") {
    return(DESeq2::DESeqResults(res_df))
  } else if (return_type == "MArrayLM") {
    # Create a minimal MArrayLM object (from limma)
    fit <- list()
    fit$coefficients <- matrix(res_df$log2FoldChange, ncol = 1)
    rownames(fit$coefficients) <- rownames(res_df)
    colnames(fit$coefficients) <- name
    fit$stdev.unscaled <- matrix(res_df$lfcSE, ncol = 1)
    fit$sigma <- rep(1, n)
    fit$df.residual <- rep(10, n)
    fit$p.value <- matrix(res_df$pvalue, ncol = 1)
    rownames(fit$p.value) <- rownames(res_df)
    colnames(fit$p.value) <- name
    class(fit) <- "MArrayLM"
    return(fit)
  }
}

de_results1 <- create_fake_de_results(se_trans, "A_vs_B")
de_results2 <- create_fake_de_results(se_prot, "A_vs_B")
de_results3 <- create_fake_de_results(se_metabo, "A_vs_B")

# Create fake enrichment results
create_fake_enrichment_results <- function(name, n_terms = 5) {
  data.frame(
    GO.ID = paste0("GO:", sprintf("%07d", 8150 + seq_len(n_terms))),
    Term = paste("Fake GO term", seq_len(n_terms)),
    Annotated = sample(20:100, n_terms, replace = TRUE),
    Significant = sample(1:20, n_terms, replace = TRUE),
    Expected = runif(n_terms, min = 1, max = 20),
    `Rank in p.value_classic` = seq_len(n_terms),
    p.value_elim = stats::runif(n_terms, min = 0, max = 1),
    p.value_classic = stats::runif(n_terms, min = 0, max = 1),
    genes = sapply(seq_len(n_terms), function(i) paste0("GENE", sample(1:100, sample(2:5, 1)), collapse = ",")),
    stringsAsFactors = FALSE
  )
}

enrichment_results1 <- create_fake_enrichment_results("A_vs_B")
enrichment_results2 <- create_fake_enrichment_results("A_vs_B")
enrichment_results3 <- create_fake_enrichment_results("A_vs_B")

dde_trans <- DeeDeeExperiment::DeeDeeExperiment(se_trans)
dde_prot <- DeeDeeExperiment::DeeDeeExperiment(se_prot)
dde_metabo <- DeeDeeExperiment::DeeDeeExperiment(se_metabo)

dde_trans <- DeeDeeExperiment::addDEA(dde_trans, de_results1)
dde_trans <- DeeDeeExperiment::renameDEA(dde_trans, "de_results1", "A_vs_B")

dde_prot <- DeeDeeExperiment::addDEA(dde_prot, de_results2)
dde_prot <- DeeDeeExperiment::renameDEA(dde_prot, "de_results2", "A_vs_B")

dde_metabo <- DeeDeeExperiment::addDEA(dde_metabo, de_results3)
dde_metabo <- DeeDeeExperiment::renameDEA(dde_metabo, "de_results3", "A_vs_B")


dde_trans <- DeeDeeExperiment::addFEA(dde_trans, enrichment_results1, de_name = "A_vs_B")
dde_trans <- DeeDeeExperiment::renameFEA(dde_trans, "enrichment_results1", "A_vs_BFEA")
dde_trans <- DeeDeeExperiment::linkDEAandFEA(dde_trans, "A_vs_B", "A_vs_BFEA", force = FALSE)

dde_prot <- DeeDeeExperiment::addFEA(dde_prot, enrichment_results2, de_name = "A_vs_B")
dde_prot <- DeeDeeExperiment::renameFEA(dde_prot, "enrichment_results2", "A_vs_BFEA")
dde_prot <- DeeDeeExperiment::linkDEAandFEA(dde_prot, "A_vs_B", "A_vs_BFEA", force = FALSE)

dde_metabo <- DeeDeeExperiment::addFEA(dde_metabo, enrichment_results3, de_name = "A_vs_B")
dde_metabo <- DeeDeeExperiment::renameFEA(dde_metabo, "enrichment_results3", "A_vs_BFEA")
dde_metabo <- DeeDeeExperiment::linkDEAandFEA(dde_metabo, "A_vs_B", "A_vs_BFEA", force = FALSE)

# Create shared metadata for all samples
shared_metadata <- data.frame(
  sample_id = paste0("Sample", seq_len(max(
    SummarizedExperiment::ncol(dde_trans),
    SummarizedExperiment::ncol(dde_prot),
    SummarizedExperiment::ncol(dde_metabo)
  ))),
  group = rep(c("A", "B"), length.out = max(
    SummarizedExperiment::ncol(dde_trans),
    SummarizedExperiment::ncol(dde_prot),
    SummarizedExperiment::ncol(dde_metabo)
  )),
  batch = rep(1, max(
    SummarizedExperiment::ncol(dde_trans),
    SummarizedExperiment::ncol(dde_prot),
    SummarizedExperiment::ncol(dde_metabo)
  )),
  stringsAsFactors = FALSE
)

movida_list <- list(
  se_prot = se_prot,
  se_metabo = se_metabo,
  se_trans = se_trans,
  organism = "Mm",
  metadata = NULL
)

movida_list_shared <- list(
  se_prot = se_prot,
  se_metabo = se_metabo,
  se_trans = se_trans,
  organism = "Mm",
  metadata = shared_metadata
)