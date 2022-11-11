require(Matrix)

denoise_NE <- function(W_in, order = 2, K, alpha = 0.9) {
  # if (!is.null(adjacency)) {
  #   # make sure network and adjacency have same colnames / dims
  #   gene_inter <- intersect(row.names(adjacency), row.names(W_in))
  #   W_in <- W_in[gene_inter, gene_inter]
  #   adjacency <- adjacency[gene_inter, gene_inter]
  #   
  #   # reorder network and adjacency
  #   W_in <- W_in[, order(colnames(W_in))]
  #   W_in <- W_in[order(row.names(W_in)),]
  #   adjacency <- adjacency[, order(colnames(adjacency))]
  #   adjacency <- adjacency[order(row.names(adjacency)),] 
  # }
  
  
  # if no K was inputted
  if(missing(K)) {
    K <- min(20, ceiling(nrow(W_in)/10))
  }
  
  # make diagonal 0s
  W_in1 <- W_in*(1-diag(nrow(W_in))) # same as diag() <- 0
  
  # remove rows and columns of all 0
  zeroindex <- Matrix::rowSums(W_in)==0
  W0 <- W_in[!zeroindex, !zeroindex]
  W <- dn(W0, 'ave')
  W <- (W + t(W))/2
  
  DD <- colSums(abs(W0))
  
  if(length(unique(c(W))) == 2) {
    P <- W
  } else {
    P <- dominateset(abs(W), min(K, nrow(W)-1)) * sign(W)
  }
  
  P <- P + diag(colSums(abs(t(P)))+1)
  # P <- as.matrix(P)
  P <- TransitionFields(P)
  # can i take the eigenvectors of a sparse matrix?
  eig <- eigen(P)
  U <- eig$vectors
  D <- eig$values
  d <- Re(D-.Machine$double.eps)
  
  d <- (1-alpha)*d/(1-alpha*(d^order))
  D <- diag(Re(d))
  W <- U%*%tcrossprod(D,U) # U%*%D%*%t(U)
  W <- (W*(1-diag(nrow(W))))/matrix(rep(1-diag(W), nrow(W)), nrow = nrow(W))
  D <- matrix(diag(DD), nrow = length(DD))
  W <- Re(W)
  W <- D%*%W
  W[W<0] <- 0
  W <- (W + t(W))/2
  
  return(W[!zeroindex, !zeroindex])
}

dn <- function(net, type) {
  # multiply each value in the matrix by the dimension
  net <- net*nrow(net)
  
  # convert the matrix to doubles
  # no need in R?
  
  # D is the sum of the rows of net
  D <- rowSums(net) + .Machine$double.eps
  
  if(type == 'ave') {
    D <- 1/D
    D <- Matrix::Matrix(diag(D))
    return(D %*% net)
  } else if (type == 'gph') {
    D <- 1/sqrt(D)
    D <- Matrix::Matrix(diag(D))
    return(D %*% (net %*% D))
  }
}

dominateset <- function(aff_matrix, NR_OF_KNN) {
  A <- t(apply(aff_matrix, 1, function(x) sort(x, decreasing = TRUE)))
  B <- t(apply(aff_matrix, 1 ,function(x) order(-x)))
  res <- A[, 1:NR_OF_KNN]
  # inds <- matrix(rep(1:nrow(aff_matrix), NR_OF_KNN), ncol = NR_OF_KNN)
  inds <- matrix(rep(1:nrow(aff_matrix), NR_OF_KNN), ncol = 1)
  loc <- matrix(B[, 1:NR_OF_KNN], ncol = 1)
  PNN_matrix1 <- matrix(0, nrow = nrow(aff_matrix), ncol = ncol(aff_matrix))
  sub2ind <- cbind(inds, loc)
  PNN_matrix1[sub2ind] <- res
  PNN_matrix <- (PNN_matrix1+t(PNN_matrix1))/2
}

TransitionFields <- function(W) {
  W <- as.matrix(W)
  zeroindex <- rowSums(W) == 0
  W <- W * nrow(W)
  W <- dn(W, 'ave')
  W <- as.matrix(W)
  W_rt <- sqrt(colSums(abs(W))+.Machine$double.eps)
  div <- matrix(rep(W_rt, nrow(W)), nrow = nrow(W), byrow = TRUE)
  W <- W/div
  W <- tcrossprod(W) # W%*%t(W)
  Wnew <- W
  Wnew[zeroindex, ] <- 0
  Wnew[, zeroindex] <- 0
  
  return(Wnew)
}
