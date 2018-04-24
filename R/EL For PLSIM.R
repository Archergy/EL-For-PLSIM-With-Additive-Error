require(locpol)
Wn <- function(t, x, beta, bw, tilde = FALSE, ker = EpaK)
{
  t <- drop(t)
  n <- nrow(x)
  Kh <- function(x)
  {
    ker(x/bw) / bw
  }
  
  S_n <- function(l)
  {
    t((x %*% beta - t)^l) %*% Kh(x %*% beta - t) / n
  }
  
  Un <-  numeric(n)
  for(i in 1:n)
  {
    Un[i] <- Kh(x[i,] %*% beta - t) * (S_n(2) - (x[i,] %*% beta - t)*S_n(1))
  }
  if(all(Un == 0))
    return(rep_len(0, n))
  
  if(tilde)
  {
    Un_tilde <- numeric(n)
    for(i in 1:n)
    {
      Un_tilde[i] <- Kh(x[i,] %*% beta - t) * ((x[i,] %*% beta - t)*S_n(0) - S_n(1))
    }
    
    return(Un_tilde / sum(Un))
  }
  else
  {
    return(Un / sum(Un))
  }
}
Wn <- compiler::cmpfun(Wn)


eta <- function(y, x, z, beta, theta, ker = EpaK)
{
  n <- length(y)
  p <- dim(x)[2]
  q <- dim(z)[2]
  h <- regCVBwSelC(x %*% beta, y - z %*% theta, deg = 1, kernel = ker)
  h <- h*n^(-2/15)
  
  loclin <- locPolSmootherC(x %*% beta, y - z %*% theta, xeval = x %*% beta,
                            bw = h, deg = 1, kernel = ker)
  
  gval <- loclin$beta0
  g_dval <- loclin$beta1
  
  # mu1 <- function(t) t(Wn(t, x, beta, h, ker = ker)) %*% x
  # mu2 <- function(t) t(Wn(t, x, beta, h, ker = ker)) %*% z
  Weight <- sapply(x %*% beta, Wn, x = x, beta = beta, bw = h, ker = ker)
  mu1 <- t(Weight) %*% x
  mu2 <- t(Weight) %*% z

  gamma <- beta[-1]
  Jacobi <- rbind(-gamma/sqrt(1-sum(gamma^2)), diag(1, p-1))
  #Jacobi <- cbind(c(1,0), Jacobi)
  
  etaval <- matrix(0, nrow = n, ncol = ncol(Jacobi) + q)
  
  epsilon <- as.numeric(y - gval - z %*% theta)
  etaval <- cbind(((epsilon * g_dval * (x - mu1))) %*% Jacobi, epsilon * (z - mu2))
  return(na.omit(etaval))
  
}

eta <- compiler::cmpfun(eta)

