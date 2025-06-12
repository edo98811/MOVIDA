# test get inchi key
source("utils/create_annotations.R")
query_result <- uniprot_inchi_query("WQZGKKKJIJFFOK-GASJEMHNSA-N", get_ensembl = TRUE, get_uniprot = TRUE)


se_trans <- readRDS("/Users/edoardofilippi/Development/Projects/MOVIDA/Example_data/Transcriptomics/se.rds")
se_prot <- readRDS("/Users/edoardofilippi/Development/Projects/MOVIDA/Example_data/Proteomics/se.rds")

query_result$protein_ids <- query_result$protein_ids[query_result$protein_ids %in% rownames(se_prot)]
query_result$gene_ids <- query_result$gene_ids[query_result$gene_ids %in% rownames(se_trans)]
