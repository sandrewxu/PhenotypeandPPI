RWR_vector <- function(gene_known, adj_mat, alpha = 0.25, iter = 25, difference = 1e-6) {
  require(SMUT)
  adj_mat <- as.matrix(adj_mat)
  
  gene_known_ids <- which(row.names(adj_mat)%in%gene_known)
  
  df_known <- data.frame(gene = row.names(adj_mat), score = rep(0, nrow(adj_mat)))
  
  for (i in 1:length(gene_known_ids)) {
    df_known$score[gene_known_ids[i]] <- 1/length(gene_known_ids)
  }
  
  # normalize adj_mat
  adj_mat <- t(t(adj_mat)/colSums(adj_mat + .Machine$double.eps))
  
  # initialize random walk
  # p_t+1 = (1-alpha)*(adj)*p^t + r*p0
  p_t <- as.matrix(df_known$score)
  p_t1 <- p_t
  
  for (i in 1:iter) {
    p_t1 <- (alpha)*eigenMapMatMult(adj_mat, p_t) + (1-alpha)*(df_known$score)
    
    cur_diff <- sum(abs(p_t1 - p_t)) #L1 norm / manhattan distance
    print(c("iteration:", i, "difference:", cur_diff))
    
    if (cur_diff < difference) {
      return(data.frame(gene = df_known$gene, score = p_t1))
    } else {
      p_t <- p_t1
    }
  }
  
  return(data.frame(gene = df_known$gene, score = p_t1))
}
