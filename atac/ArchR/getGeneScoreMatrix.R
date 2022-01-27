
getGeneScoreMatrix <- function(
    ArchRProj,
    ClustersName = "ClustersGerm",
    file.path = NULL
) {
    GS_matrix <- getMatrixFromProject(ArchRProj, "GeneScoreMatrix")
    GS_mat <- GS_matrix@assays@data$GeneScoreMatrix %>% as.matrix()

    colnames(GS_mat) <- getCellColData(ArchRProj)[[ClustersName]]
    rownames(GS_mat) <- GS_matrix@elementMetadata$name

    GS_long <- reshape2:::melt.matrix(GS_mat)
    colnames(GS_long) <- c("symbol", "cell_type", "gene_score")

    GS_wide <- GS_long %>% 
        group_by(cell_type, symbol) %>% 
        summarize(gs = mean(gene_score)) %>% 
        spread(key = cell_type, value = gs)
    
    if (file.path) saveRDS(file.path)

    return (GE_wide)
}