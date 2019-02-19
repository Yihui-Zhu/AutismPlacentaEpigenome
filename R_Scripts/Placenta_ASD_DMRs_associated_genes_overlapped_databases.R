### load packages #################
library(pasilla)
library(Biobase)
library(DESeq)
library(GeneOverlap)

######################################
### Overlapping test #################
geneset1 <- Placenta_ASD_DMRs_associated_genes 
geneset2 <- other_database

## add two list together
geneset3 = c(geneset1, geneset2)

## universe <- length(unique(geneset3))
universe = #number of genes as background

common <- length(
  intersect(
    unique(geneset1),
    unique(geneset2)
  )
)

mat <- matrix(
  c(
    universe - length(union(geneset1, geneset2)),
    length(setdiff(geneset1, geneset2)),
    length(setdiff(geneset2, geneset1)),
    length(intersect(geneset1, geneset2))
  ),
  nrow=2
)

fr <- fisher.test(mat)
fr

overl <- newGeneOverlap(
  unique(geneset1),
  unique(geneset2),
  genome.size=universe)

overl <- testGeneOverlap(overl)

## show the overlap genes
getIntersection(overl)

