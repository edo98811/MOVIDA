# MOVIDA

A package to interactively explore multi-omics data.

## Installation

```R
install.packages("devtools")     
devtools::install_github("edoardofilippi/MOVIDA")
```
## Usage

To use the Kegg pathway visualization function:
```R
library(MOVIDA)
pathway <- "hsa04110"  # Example pathway ID
graph <- kegg_to_graph(pathway, organism = "hsa", return_type = "visNetwork")
graph
```



