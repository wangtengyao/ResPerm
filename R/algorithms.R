#' Naive permutation test
#' @description
#' Implements the naive permutation test for testing whether the first
#' regression coefficient associated with Y~X is zero.
#' @param X design matrix
#' @param Y response vector
#' @param K number of permutations
#' @return p-value of the test
#' @export
NaivePermutationTest <- function(X, Y, K=999){
  n <- dim(X)[1]; p <- dim(X)[2]
  U <- svd(X[, -1], nu=n)$u[, -(1:(p-1))]
  e.hat <- as.vector(t(U) %*% X[, 1])
  eps.hat <- as.vector(t(U) %*% Y)
  eps.hat.perm = replicate(K, sample(eps.hat))
  test.stat <- sum(e.hat * eps.hat)
  test.stat.perm <- colSums(eps.hat.perm * e.hat)
  pval <- (sum(abs(test.stat.perm) >= abs(test.stat)) + 1) / (K + 1)
  return(pval)
}

#' Generate permutations
#' @description
#' generate a subgroup of permutations of S_n of cardinality K
#' @param n total number of elements to permute
#' @param K size of the subgroup
#' @return matrix of permutation indices, each row represents a permutation
#' @export
generatePermutations <- function(n, K){
  ind <- head(sample(n), n - n %% (K+1))
  mx <- matrix(ind, nrow=K+1)
  perm_ind <- matrix(0, K, n)
  for (r in 1:K){
    perm_ind[r, ] <- plyr::mapvalues(1:n, mx, mx[c((r+1):(K+1), 1:r), ])
  }
  return(perm_ind)
}

#' Residual permutation test
#' @description
#' Implements the residual permutation test for testing whether the first
#' regression coefficient associated with Y~X is zero.
#' @param X design matrix
#' @param Y response vector
#' @param K number of permutations
#' @return p-values for RPT and RPT_em
#' @example example/example.R
#' @export
ResidualPermutationTest <- function(X, Y, K=99){
  Z <- X[, 1]; X <- X[, -1]
  n <- dim(X)[1]; p <- dim(X)[2]

  stat1 <- stat2 <- rep(0, K)
  perm_ind <- generatePermutations(n, K)
  
  for (r in 1:K){
    idx <- perm_ind[r, ]  # idx represents the associated permutation index of P_r
    
    # compute projection of Z onto the complement of column space of (X, X[idx,])
    Vz <- tryCatch({
     Vz<- lm.fit(cbind(X,X[idx,]), Z)$residuals
    }, error=function(e){"error"}) 
    
    # use svd if (X, X[idx,]) is singular
    if(identical(Vz, "error")){
      tmp <- svd(cbind(X, X[idx,]), nu=2*p)$u
      Vz <- Z - tmp %*% (t(tmp) %*% Z)
    }

    # compute inner product of Z with Y and Y[idx] on the orthogonal space of (X, X[idx,])
    stat1[r] <- as.numeric(sum(Vz * Y))
    stat2[r] <- as.numeric(sum(Vz * Y[idx]))
  }
  
  # compute pval of RPT_empirical and RPT respectively
  pval_RPT_em <- (sum(abs(stat2) >= abs(stat1)) + 1) / (K+1)
  pval_RPT <- (sum(abs(stat2) >= min(abs(stat1))) + 1) / (K+1)
  return(list(pval=pval_RPT, pval_empirical=pval_RPT_em))
}


##### debugging #####
if (sys.nframe()==0L){
  vector.normalise <- function(v){v / sqrt(sum(v^2))}
  n <- 600; p <- 100; k <- 5; coeff_sig <- 1; sigma <- 1; AR <- 0; 
  coeff_shape='gaussian'
  X <- matrix(rnorm(n*p), n)
  alpha <- vector.normalise(c(rnorm(k), rep(0, p-k)))
  beta <- vector.normalise(c(rnorm(k), rep(0, p-k)))
  e <- rt(n, df=2)
  eps <- rt(n, df=2)
  
  b <- 0.5 # coefficient of interest
  Z <- X %*% alpha + e
  Y <- X %*% beta + b * Z + eps
  
  X_aug <- cbind(Z, X)
  NaivePermutationTest(X_aug, Y, K=99)
  ResidualPermutationTest(X_aug, Y, K=99)
}

