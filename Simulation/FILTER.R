#' FILTER the HPO2Gene KnowledgeBase (H2GKB) to become disease-specific
#' 
#' This function filters the H2GKB to retain disease-specific phenotype information. This
#' is done by inputting a list of HPO IDs associated with a specific disease. The functio
#' 
#' @param disease_ids A list of Human Phenotype Ontology (HPO) IDs associated with a disease
#' @param heading Whether a heading is included with the list of HPO IDs (default: TRUE)
#' @param dir_H2GKB The folder path that leads to the H2GKB database if stored locally. 
#' Subdirectories should include "knowledgebase" and "skewness".
#' @param keep If dir_H2GKB is not stored locally, should a copy be kept?
#' @param top How many top scores to keep from each HPO ID (default: 1000)
#' @param filter After creating the gene-phenotype score matrix, what percent of scores should be kept?
#' 
#' @export
FILTER <- function(disease_ids, gene_list, heading = TRUE, dir_H2GKB, keep = TRUE, top = 1000, filter = 10) {
  # Download the H2GKB from online if not stored locally
  if(missing(dir_H2GKB)) {
    url <- "https://github.com/WGLab/Phen2Gene/releases/download/1.1.0/H2GKBs.zip"
    download.file(url, destfile = "Data/H2GKBs.zip")
    unzip("Data/H2GKBs.zip", exdir = "Data/H2GKBs/")
    dir <- "Data/H2GKBs/"
    unlink("Data/H2GKBs.zip")
  } else if (endsWith(dir_H2GKB, "/")) {
    dir <- dir_H2GKB
  } else {
    dir <- paste0(dir_H2GKB, "/")
  }
  
  # Extract certain HPO IDs from H2GKB
  files_skewness <- list.files(paste0(dir,"skewness"), full.names = TRUE)
  HPO_ids <- disease_ids
  HPO_ids <- gsub(":","_", HPO_ids)
  
  H2GKB_ids <- unlist(lapply(files_skewness, function(x) stringr::str_sub(x, -10, -1)))
  H2GKB_ids <- intersect(H2GKB_ids, HPO_ids)
  
  files_knowledgebase <- unlist(lapply(H2GKB_ids, function(x) {
    paste0(dir, "Knowledgebase/", x, ".candidate_gene_list")
  }))
  files_skewness <- unlist(lapply(H2GKB_ids, function(x) {
    paste0(dir, "skewness/", x)
  }))
  
  HPO_data <- lapply(files_knowledgebase, data.table::fread, header = TRUE)
  skew_data <- unlist(lapply(files_skewness, data.table::fread, header = FALSE))
  
  if(missing(dir_H2GKB) && keep == FALSE) {
    unlink("H2GKBs/", recursive = TRUE)
  }
  
  gene_H2GKB <- unique(unlist(lapply(HPO_data, function(x) x$Gene)))
  gene_H2GKB <- intersect(gene_list, gene_H2GKB)
  
  HPO_data <- lapply(HPO_data, function(x) {
    x[x$Gene%in%gene_H2GKB]
  })
  
  # Keep top scores for each HPO ID
  HPO_ranked <- lapply(HPO_data, function(x) {
    x[x$Rank <= top]
  })
  
  # Weight HPO scores using skewness
  HPO_weighted <- lapply(1:length(HPO_ranked), function(x) {
    HPO_ranked[[x]]$Score <- HPO_ranked[[x]]$Score * skew_data[x]
    return(HPO_ranked[[x]])
  })
  
  HPO_weighted_sub <- lapply(HPO_weighted, function(x) x[,.(Gene, Score)])
  
  # Make sure each HPO ID has all genes (fill in 0s)
  gene_H2GKB <- data.table::as.data.table(gene_H2GKB)
  
  HPO_allgenes <- lapply(HPO_weighted_sub, function(x) {
    HPO_merge <- data.table::merge.data.table(x, gene_H2GKB, by.x = "Gene", by.y = "gene_H2GKB", all.y = TRUE)
    HPO_merge$Score[is.na(HPO_merge$Score)] <- 0
    return(HPO_merge)
  })
  
  # Generate gene-phenotype score matrix
  HPO_scores <- lapply(HPO_allgenes, function(x) x[, .(Score)])
  gene_score_matrix <- matrix(unlist(HPO_scores), nrow = length(HPO_scores), byrow = TRUE)
  
  # Generate a gene similarity matrix based on gene-phenotype scores
  gene_similarity <- SMUT::eigenMapMatMult(t(gene_score_matrix), gene_score_matrix) / (sqrt(tcrossprod(colSums((gene_score_matrix)^2))))
  colnames(gene_similarity) <- gene_H2GKB$gene_H2GKB
  rownames(gene_similarity) <- gene_H2GKB$gene_H2GKB
  
  # Filter the gene similarity matrix
  # cutoff <- quantile(gene_similarity, (1-filter/100))
  # gene_similarity[gene_similarity < cutoff] <- 0
  
  return(gene_similarity)
}
