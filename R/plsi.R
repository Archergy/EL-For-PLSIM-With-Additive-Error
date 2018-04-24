#semi-parametric estimation of PLSIM procedure by Xia and Hardle
require(pracma)
gplsiabc <- function(x, y, z, h)
{
  p <- dim(x)[2]
  q <- dim(z)[2]
  n <- length(y)
  pq <-  p + q
  
  onep <- matrix(1, p, 1)
  onen <- matrix(1, n, 1)
  B <- diag(1, p)
  nd <- 1
  
  a <- matrix(0, n, 1)
  h2 <- 2*n^(-2 / (p+4))
  invzz <- solve(t(z) %*% z) %*% t(z)
  eyep1 <- diag(1, p+1) / (n^2)
  for(iter in 1:p)
  {
    ye <- y - a
    theta <- invzz %*% ye
    ye <- ye - z %*% theta
    
    ab <- matrix(1, p, n)
    for(i in 1:n)
    {
      xi <- x - repmat(x[i,], n, 1)
      kernel <- matrix(exp(-rowSums((xi %*% B)^2) / h2), ncol=1)
      onexi <- cbind(onen, xi)
      xk <- onexi * repmat(kernel, 1, p+1)
      abi <- solve(t(xk) %*% onexi + eyep1) %*% t(xk) %*% ye
      ab[,i] <- abi[2:(p+1)]
      a[i] <- abi[1]
    }
    
    eig <- eigen(ab %*% t(ab))
    B0 <- eig$vector
    D <- eig$values
    I <- sort(D, index.return = T)$ix
    D <- sort(D)
    
    B <- B0
    for(k in 1:p)
    {
      B0[,k] <- B[, I[k]]
    }
    for(k in 1:nd)
    {
      B[,k] <- B0[,p-k+1]
    }
    B <- B[,1:max(1, p - iter)]
  }
  
  beta <- B
  
  h2 <- 2 * h * h
  eyepq <- diag(1, pq) / (n^2)
  for(k in 1:10)
  {
    tmp <- repmat(x %*% beta, 1, n)
    d <- (tmp - t(tmp))^2
    ker <- exp(-d/h2)
    ker <- ker / repmat(matrix(colSums(ker), ncol = 1), 1, n)
    md <- matrix(0, n, p*p)
    mc <- matrix(0, n, p)
    me <- matrix(0, n, p)
    
    D22 <- 0
    C2 <- 0
    mcz <- matrix(0, n, q)
    mez <- matrix(0, n, q)
    mdd <- matrix(0, n, p*q)
    for(i in 1:n)
    {
      tmp <- x - repmat(x[i,], n, 1)
      tmp1 <- matrix(repmat(ker[,i], 1, p),n) * tmp
      md[i,] <- matrix(t(tmp1) %*% tmp, 1, p*p)
      mc[i,] <- colSums(tmp1)
      me[i,] <- t(y) %*% tmp1
      
      z1 <- matrix(repmat(ker[,i], 1, q),n) * z
      D22 <- D22 + t(z1) %*% z
      mcz[i,] <- colSums(z1)
      C2 <- C2 + t(y) %*% z1
      
      mdd[i,] <- matrix(t(tmp) %*% z1, 1, q*p)
    }
    
    for(j in 1:5)
    {
      ye <- y - z %*% theta
      tmp <- repmat(x %*% beta, 1, n)
      d1 <- tmp - t(tmp)
      ker1 <- d1 * ker
      s2 <- rowSums(d1 * ker1)
      s1 <- rowSums(ker1)
      s <- rowSums(ker)
      d <- s2 * s - s1 * s1 + 1/n^2
      a <- (ker %*% ye * s2 - ker1 %*% ye * s1) / d
      b <- (ker1 %*% ye * s - ker %*% ye * s1) / d
      D <- matrix(t(b*b) %*% md, p, p)
      C <- t(b) %*% me - t(b * a) %*% mc
      
      Cz <- C2 - t(a) %*% mcz
      D12 <- matrix(t(b) %*% mdd, p, q)
      D <- rbind(cbind(D, D12), cbind(t(D12), D22))
      C <- cbind(C, Cz)
      
      beta <- solve(D + eyepq) %*% t(C)
      theta <- beta[(p+1) : pq]; beta <- beta[1:p]
      beta <- sign(beta[1]) * beta / as.numeric(sqrt(t(beta) %*% beta))
      
    }
  }
  
  return(list(g = a, beta = matrix(beta, ncol=1), 
              theta = matrix(theta, ncol = 1)))
}


