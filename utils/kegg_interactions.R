library(KEGGREST)

# Get reactions for glucose (C00031)
rxns <- keggLink("reaction", "cpd:C00031")

# Get enzymes linked to these reactions
enzymes <- unlist(lapply(names(rxns), function(r) {
  keggLink("enzyme", r)
}))

# Get human genes for those enzymes
human_genes <- unique(unlist(lapply(enzymes, function(e) {
  keggLink("hsa", e)
})))

# Translate to UniProt
UNIPROTs <- unlist(lapply(human_genes, function(g) {
  keggConv("uniprot", g)
}))