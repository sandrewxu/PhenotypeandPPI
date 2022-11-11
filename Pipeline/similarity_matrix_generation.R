similarity_matrix_generation <- function(score_mat, type = "cosine") {
  if (type == "cosine") {
    gene_sim_mat_cosine <- crossprod(score_mat) / (sqrt(tcrossprod(colSums((score_mat)^2))))
    diag(gene_sim_mat_cosine) <- 0
    return(gene_sim_mat_cosine)
  } else if (type == "spearman") {
    gene_sim_mat_spearman <- cor(score_mat, method = c("spearman"))
    diag(gene_sim_mat_spearman) <- 0
    return(gene_sim_mat_spearman)
  } else if (type == "kendall") {
    require(pcaPP)
    gene_sim_mat_kendall <- cor.fk(score_mat)
    diag(gene_sim_mat_kendall) <- 0
    return(gene_sim_mat_kendall)
  } else if (type == "pearson") {
    gene_sim_mat_pearson <- cor(score_mat, method = c("pearson"))
    diag(gene_sim_mat_pearson) <- 0
    return(gene_sim_mat_pearson)
  }
}