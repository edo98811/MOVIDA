# MOVIDA

A package to interactively explore multi-omics data.

## Installation

```R
install.packages("devtools")     
devtools::install_github("edoardofilippi/MOVIDA")
```
## Usage

```R
library("MOVIDA")

movida_list <- list(
  dde_metabo = dde,
  organsim = "Mm"
)

MovidaApp(movida_list)
```



