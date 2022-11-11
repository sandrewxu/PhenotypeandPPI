#' REDUCE noise of the integrated phenotype and general human network
#' 
#' This denoising function adapts an algorithm from Wang et al, 2018
#' https://www.nature.com/articles/s41467-018-05469-x
#' 
#' @param network integrated phenotype and general human network to denoise
#' @param K K nearest neighbors of each node
#' @param alpha regularization parameter
#' 
#' @export
REDUCE <- function(network, K = min(20, ceiling(nrow(W_in)/10)), alpha = 0.9) {
  # Pre-processing
  W_in <- (network + Matrix::t(network))/2
  Matrix::diag(W_in) <- 0
  
  ## remove rows and columns of all 0
  nullrows <- Matrix::rowSums(W_in) == 0
  W0 <- W_in[!nullrows, !nullrows]
  
  # Diagonal matrix for rescaling at end
  scale <- Matrix::colSums(abs(W0))
  
  # Generate P, the first transition matrix
  P <- rownorm(as.matrix(W0)) # normalize by rows
  P <- KNN(P, min(K, nrow(P)-1)) # Localize -- keep k nearest neighbors
  
  # Generate T , the DSM transition matrix
  P <- P + diag(rowSums(abs(P))+1) # rescale based on diagonal
  T0 <- computeDSM(P)
  
  # Use eigenvalues to compute W_inf
  eig <- eigen(T0)
  U <- eig$vectors
  D <- eig$values
  D0 <- Re(D-.Machine$double.eps)
  
  D_inf <- (1-alpha)*D0/(1-alpha*D0^2)
  D <- diag(Re(D_inf))
  
  W_inf <- SMUT::eigenMapMatMult(U, SMUT::eigenMapMatMult(D, t(U)))
  
  # Rescale values based on diagonal
  div <- matrix(rep(1-diag(W_inf), nrow(W_inf)), nrow = nrow(W_inf))
  diag(W_inf) <- 0
  W_inf <- W_inf/div
  
  # Rescale matrix to original size
  W_inf <- Re(W_inf)
  scale <- matrix(diag(scale), nrow = length(scale))
  W_out <- SMUT::eigenMapMatMult(scale, W_inf)
  
  # Post-processing
  W_out[W_out<0] <- 0
  W_out <- (W_out+t(W_out))/2
  
  res <- matrix(0, nrow = nrow(network), ncol = ncol(network))
  res[!nullrows, !nullrows] <- W_out
  rownames(res) <- rownames(network)
  colnames(res) <- colnames(network)
  
  return(res)
}

rownorm <- function (network) {
  network <- network*nrow(network) + .Machine$double.eps
  
  D <- Matrix::rowSums(network) + .Machine$double.eps
  D <- 1/D
  D <- diag(D)
  
  res <- SMUT::eigenMapMatMult(D, network)
  return(res)
}

KNN <- function(P, K) {
  # order P by rows and record original indices
  A <- t(apply(P, 1, function(x) sort(x, decreasing = TRUE)))
  B <- t(apply(P, 1 ,function(x) order(-x)))
  
  # filter for K nearest neighbors
  res <- A[, 1:K]
  
  # create an index to repopulate the original matrix
  inds <- matrix(rep(1:nrow(P), K), ncol = 1)
  loc <- matrix(B[, 1:K], ncol = 1)
  sub2ind <- cbind(inds, loc)
  
  # remap filtered values to original matrix
  PNN_matrix1 <- matrix(0, nrow = nrow(P), ncol = ncol(P))
  PNN_matrix1[sub2ind] <- res
  # PNN_matrix <- (PNN_matrix1+t(PNN_matrix1))/2
  
  return(PNN_matrix1)
}

computeDSM <- function(P) {
  zeroindex <- Matrix::rowSums(P) == 0
  P <- rownorm(P)
  P_rt <- sqrt(colSums(abs(P))+.Machine$double.eps)
  P <- P/matrix(rep(P_rt, nrow(P)), nrow = nrow(P), byrow = TRUE)
  P <- SMUT::eigenMapMatMult(P, t(P))
  P[zeroindex, ] <- 0
  P[, zeroindex] <- 0
  
  return(P)
}
