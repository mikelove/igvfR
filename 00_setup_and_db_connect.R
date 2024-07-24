# load connection to IGVF arango db
library(reticulate)
use_python("/opt/homebrew/bin/python3.11")
a <- import("arango")
cl <- a$ArangoClient(hosts="https://db-dev.catalog.igvf.org")
db <- cl$db("igvf", username="guest", password="guestigvfcatalog")

# hack to save the symbol to gene map
if (!file.exists("sym2gene.rda")) {
  library(dplyr)
  library(tibble)
  library(EnsDb.Hsapiens.v86)
  sym2gene <- ensembldb::genes(EnsDb.Hsapiens.v86, 
                               return.type="DataFrame") |>
    as_tibble() |>
    dplyr::select(symbol, gene_id) |>
    dplyr::filter(!duplicated(symbol))
  save(sym2gene, file="sym2gene.rda")
} else {
  load("sym2gene.rda")
}
