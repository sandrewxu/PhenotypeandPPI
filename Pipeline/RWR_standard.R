RWR_standard <- function(network, adj_mat, alpha = 0.75, iter = 25, difference = 1e-6) {
  require(SMUT)
  network <- as.matrix(network)
  adj_mat <- as.matrix(adj_mat)
  
  # make sure network and adj_mat have same genes
  gene_inter <- intersect(row.names(adj_mat), row.names(network))
  network <- network[gene_inter, gene_inter]
  adj_mat <- adj_mat[gene_inter, gene_inter]
  
  # reorder network and adj_mat
  network <- network[, order(colnames(network))]
  network <- network[order(row.names(network)),]
  adj_mat <- adj_mat[, order(colnames(adj_mat))]
  adj_mat <- adj_mat[order(row.names(adj_mat)),]
  
  # normalize adj_mat
  adj_mat <- t(t(adj_mat)/colSums(adj_mat))
  
  # initialize random walk
  # p_t+1 = (1-alpha)*(adj)*p^t + r*p0
  p_t <- network
  p_t1 <- p_t
  
  for (i in 1:iter) {
    p_t1 <- (alpha)*eigenMapMatMult(adj_mat, p_t) + (1-alpha)*(network)
    
    cur_diff <- sum(abs(p_t1 - p_t)) #L1 norm / manhattan distance
    print(c("iteration:", i, "difference:", cur_diff))
    
    if (cur_diff < difference) {
      return(p_t1)
    } else {
      p_t <- p_t1
    }
  }
  
  return(p_t1)
}
