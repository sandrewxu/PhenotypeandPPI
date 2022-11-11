#' INTEGRATE phenotype similarity into the general human network
#' 
#' This function uses a random walk process to diffuse a phenotype 
#' similarity network onto a general human network
#' 
#' @param network General human network to be inputted
#' @param adj_matrix Adjacency matrix / phenotype similarity network
#' @param r regularization parameter -- the degree of the network that should be retained
#' @param iter max number of iterations
#' @param difference difference between two consecutive networks to end diffusion
#' 
#' @export
INTEGRATE <- function(network, adj_matrix, r = 0.25, iter = 40, difference = 1e-6) {
  # Download STRING if no network?
  
  # Convert to dense matrices
  network <- as.matrix(network)
  adj_matrix <- as.matrix(adj_matrix)
  
  # Intersect to same set of genes and same order
  gene_inter <- intersect(row.names(adj_matrix), row.names(network))
  network <- network[gene_inter, gene_inter]
  adj_matrix <- adj_matrix[gene_inter, gene_inter]
  network <- network[, order(colnames(network))]
  network <- network[order(row.names(network)),]
  adj_matrix <- adj_matrix[, order(colnames(adj_matrix))]
  adj_matrix <- adj_matrix[order(row.names(adj_matrix)),]
  
  # Normalize adj_matrix by columns
  adj_matrix <- t(t(adj_matrix)/colSums(adj_matrix))
  
  # Initialize random walk
  ## P_t+1 = (1-r)WP_t + rP0
  ## P0 is network
  p_t <- network
  p_t1 <- p_t
  
  for (i in 1:iter) {
    # SMUT::eigenMapMatMult is an RcppEigen implementation of 
    # standard matrix multiplication
    p_t1 <- (1-r)*SMUT::eigenMapMatMult(adj_matrix, p_t) + r*network
    
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
