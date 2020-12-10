# Function definition for getting enrichment p-values
# Input: Selected, e.g., DE genes (genes.of.interest), all genes (gene.universe)
# Output: A vector of enrichment p-values for all KEGG pathways
# What it does: For each pathway, it tests (hypergeometric test) for
#	the association between two properties of a gene, 
#	being DE, being a member of the pathway.
# 	and calculates a p-value for association (enrichment).
getORApvals <- function(genes.of.interest, gene.universe) {

  # Gene universe
  all.geneIDs <- gene.universe

  # Define the function for hypergeometric test
  hyperg <- Category:::.doHyperGInternal
  hyperg.test <-
    function(pathway.genes, genes.of.interest, all.geneIDs, over=TRUE)
    {
      white.balls.drawn <- length(intersect(genes.of.interest, pathway.genes))
      white.balls.in.urn <- length(pathway.genes)
      total.balls.in.urn <- length(all.geneIDs)
      black.balls.in.urn <- total.balls.in.urn - white.balls.in.urn
      balls.pulled.from.urn <- length(genes.of.interest)
      hyperg(white.balls.in.urn, black.balls.in.urn,
           balls.pulled.from.urn, white.balls.drawn, over)
    }
  # get the p-values for all pathways
  pVals.by.pathway <-
    t(sapply(genes.by.pathway, hyperg.test, genes.of.interest, all.geneIDs))
    
  pVals.by.pathway
}

