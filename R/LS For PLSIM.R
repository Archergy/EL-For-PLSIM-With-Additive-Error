NR_iter <- function(x, y, z, beta_h, theta_h) #Newton-Raphson迭代
{
  #options(warn = -1)
  p = dim(x)[2]
  q = dim(z)[2]
  n <- length(y)
    
  h <- regCVBwSelC(x %*% beta_h, y - z %*% theta_h, deg = 1, kernel = EpaK)# * n^(-2/15)
  #print(h)
  Kh <- function(u)
  {
    ifelse(abs(u/h) <= 1, 3/4 * (1 - (u/h)^2) / h, 0)
  }
  
  
  Kh_d <- function(u)
  {
    ifelse(abs(u/h) <= 1, -3*u / (2 * h^3), 0)
  }
  
  Kh_dd <- function(u)
  {
    ifelse(abs(u/h) <= 1, -3 / (2 * h^3), 0)
  }
  
  Kst <- function(xi, s, t)
  {
    ui <- as.numeric(xi %*% beta_h)
    u <- as.numeric(x %*% beta_h)
    t1 <- as.numeric(Kh(u - ui))
    t2 <- as.numeric((u - ui)^s)
    t3 <- as.numeric((y - z %*% matrix(theta_h))^t)
    
    sum(t1 * t2 * t3)
  }
  K00 <- function(xi) Kst(xi, s=0, t=0)
  K01 <- function(xi) Kst(xi, s=0, t=1)
  K10 <- function(xi) Kst(xi, s=1, t=0)
  K11 <- function(xi) Kst(xi, s=1, t=1)
  K20 <- function(xi) Kst(xi, s=2, t=0)
  K21 <- function(xi) Kst(xi, s=2, t=1)
  
  A <- function(xi)
  {
    K20(xi) * K01(xi) - K10(xi) * K11(xi)
  }
  
  B <- function(xi)
  {
    K00(xi) * K20(xi) - K10(xi) * K10(xi)
  }
  
  
  g <- function(xi)
  {
    A(xi) / B(xi)
  }
  
  pKs1_ptheta <- function(xi, s)
  {
    u <- (x - rep(1,n) %*% t(xi)) %*% beta_h
    t1 <- Kh(u)
    t2 <- u^s
    t(t1 * t2) %*% (-z)
  }
  pK01_ptheta <- function(xi) pKs1_ptheta(xi, s=0)
  pK11_ptheta <- function(xi) pKs1_ptheta(xi, s=1)
  pK21_ptheta <- function(xi) pKs1_ptheta(xi, s=2)
  
  
  
  pKst_pgamma <- function(xi, s, t)
  {
    xmxi <- x - rep(1,n) %*% t(xi)
    u <- xmxi %*% beta_h
    
    t1 <- (y - z %*% theta_h)^t
    t2 <- Kh_d(u) * u^s
    t3 <- ifelse(rep(s == 0, n), 0, s * Kh(u) * u^(s-1))
    
    t(t1 * (t2 + t3)) %*% xmxi %*% Jacobi
  }
  pK00_pgamma <- function(xi) pKst_pgamma(xi, s=0, t=0)
  pK01_pgamma <- function(xi) pKst_pgamma(xi, s=0, t=1)
  pK10_pgamma <- function(xi) pKst_pgamma(xi, s=1, t=0)
  pK11_pgamma <- function(xi) pKst_pgamma(xi, s=1, t=1)
  pK20_pgamma <- function(xi) pKst_pgamma(xi, s=2, t=0)
  pK21_pgamma <- function(xi) pKst_pgamma(xi, s=2, t=1)
  
  p2Ks1_pgamma_pthetaT <- function(xi, s)
  {
    xmxi <- x - rep(1,n) %*% t(xi)
    u <- as.numeric(xmxi %*% beta_h)
    
    t1 <- Kh_d(u) * u^s
    t2 <- ifelse(rep(s == 0, n), 0, s * Kh(u) * u^(s-1))
    
    t((t1 + t2) * (-z)) %*% xmxi %*% Jacobi
  }
  p2K01_pgamma_pthetaT <- function(xi) p2Ks1_pgamma_pthetaT(xi, s=0) 
  p2K11_pgamma_pthetaT <- function(xi) p2Ks1_pgamma_pthetaT(xi, s=1) 
  p2K21_pgamma_pthetaT <- function(xi) p2Ks1_pgamma_pthetaT(xi, s=2) 
  
  
  p2Kst_pgamma2 <- function(xi, s, t) {
    xmxi <- x - rep(1,n) %*% t(xi)
    u <- xmxi %*% beta_h
    v <- xmxi %*% Jacobi
    
    t8 <-  (diag(1, p-1) * (1-sum(gamma^2)) + gamma %*% t(gamma)) / (1-sum(gamma^2))^(3/2)
    
    t1 <- (y - z %*% theta_h)^t
    t2 <- Kh_dd(u) * u^s
    t3 <- ifelse(rep(s == 0, n), 0, 2 * s * Kh_d(u) * u^(s-1))
    t4 <- ifelse(rep(s<2, n), 0, 2 * Kh(u))
    t5 <- xmxi %*% Jacobi
    
    t6 <- Kh_d(u) * u^s
    t7 <- ifelse(rep(s == 0, n), 0, s * Kh(u) * u^(s-1))
    
    t(t1 * (t2 + t3 + t4) * t5) %*% t5 + t(t1 * (t6 + t7)) %*% xmxi[,1] * t8
  }
  p2K00_pgamma2 <- function(xi) p2Kst_pgamma2(xi, s=0, t=0)
  p2K01_pgamma2 <- function(xi) p2Kst_pgamma2(xi, s=0, t=1)
  p2K10_pgamma2 <- function(xi) p2Kst_pgamma2(xi, s=1, t=0)
  p2K11_pgamma2 <- function(xi) p2Kst_pgamma2(xi, s=1, t=1)
  p2K20_pgamma2 <- function(xi) p2Kst_pgamma2(xi, s=2, t=0)
  p2K21_pgamma2 <- function(xi) p2Kst_pgamma2(xi, s=2, t=1)
  
  C <- function(xi)
  {
    K20(xi) *  pK01_ptheta(xi) - K10(xi) * pK11_ptheta(xi)
  }
  
  pA_pgamma <- function(xi)
  {
    t1 <- pK20_pgamma(xi) * K01(xi) + K20(xi) * pK01_pgamma(xi)
    t2 <- pK10_pgamma(xi) * K11(xi) + K10(xi) * pK11_pgamma(xi)
    
    t1 - t2
  }
  
  pB_pgamma <- function(xi)
  {
    t1 <- pK00_pgamma(xi) * K20(xi) + K00(xi) * pK20_pgamma(xi)
    t2 <- 2 * K10(xi) * pK10_pgamma(xi)
    
    t1 - t2
  }
  
  pCT_pgamma <- function(xi)
  {
    t1 <- t(pK01_ptheta(xi)) %*% pK20_pgamma(xi) + K20(xi) * p2K01_pgamma_pthetaT(xi)
    t2 <- t(pK11_ptheta(xi)) %*% pK10_pgamma(xi) + K10(xi) * p2K11_pgamma_pthetaT(xi)
    
    t1 - t2
  }
  
  D <- function(xi)
  {
    t1 <- pA_pgamma(xi) * B(xi)
    t2 <- A(xi) * pB_pgamma(xi)
    
    t1 - t2
  }
  
  p2A_pgamma2 <- function(xi)
  {
    t1 <- p2K20_pgamma2(xi) * K01(xi) + t(pK20_pgamma(xi)) %*% pK01_pgamma(xi)
    t2 <- t(pK01_pgamma(xi)) %*% pK20_pgamma(xi) + K20(xi) * p2K01_pgamma2(xi)
    t3 <- p2K10_pgamma2(xi) * K11(xi) + t(pK10_pgamma(xi)) %*% pK11_pgamma(xi)
    t4 <- t(pK11_pgamma(xi)) %*% pK10_pgamma(xi) + K10(xi) * p2K11_pgamma2(xi)
    
    (t1 + t2) - (t3 + t4)
  }
  
  p2B_pgamma2 <- function(xi)
  {
    t1 <- p2K00_pgamma2(xi) * K20(xi) + t(pK00_pgamma(xi)) %*% pK20_pgamma(xi)
    t2 <- t(pK20_pgamma(xi)) %*% pK00_pgamma(xi) + K00(xi) * p2K20_pgamma2(xi)
    t3 <- 2 * (t(pK10_pgamma(xi)) %*% pK10_pgamma(xi) + K10(xi) * p2K10_pgamma2(xi))
    
    t1 + t2 - t3
  }
  
  pDT_pgamma <- function(xi)
  { 
    t1 <- p2A_pgamma2(xi) * B(xi) + t(pA_pgamma(xi)) %*% pB_pgamma(xi)
    t2 <- t(pB_pgamma(xi) %*% pA_pgamma(xi)) + A(xi) * p2B_pgamma2(xi)
    
    t1 - t2
  }
  
  pg_ptheta <- function(xi)
  {
    C(xi) / B(xi)
  }
  
  pg_pgamma <- function(xi)
  {
    D(xi) / (B(xi)^2)
  }
  
  p2g_pgamma_pthetaT <- function(xi)
  {
    t1 <- pCT_pgamma(xi) * B(xi) - t(C(xi)) %*% pB_pgamma(xi)
    
    t1 / (B(xi)^2)
  }
  
  p2g_pgamma2 <- function(xi)
  {
    t1 <- pDT_pgamma(xi) * (B(xi)^2) - 2 * B(xi) * t(D(xi)) %*% pB_pgamma(xi)
    
    t1 / (B(xi)^4)
  }

  
  
  gamma <- beta_h[-1,]
  Jacobi <- rbind(-gamma/sqrt(1-sum(gamma^2)), diag(1, p-1))
  
  for(iter in 1:5)
  {
    pQ_ptheta <- matrix(0, nrow = 1, ncol = q)
    pQ_pgamma <- matrix(0, nrow = 1, ncol = p-1)
    p2Q_ptheta2 <- matrix(0, nrow = q, ncol = q)
    p2Q_pgamma_pthetaT <- matrix(0, nrow = q, ncol = p-1)
    p2Q_pgamma2 <- matrix(0, nrow = p-1, ncol = p-1)
    Hessian <- matrix(0, nrow = p-1+q, ncol = p-1+q)
    
    for(i in 1:n)
    {
      xi <- x[i,]
      t1 <- as.numeric(y[i] - z[i,] %*% theta_h - g(xi))
      t2 <- t(z[i,]) + pg_ptheta(xi)
      if(is.na(t1)) next
      
      pQ_ptheta <- pQ_ptheta + 2 * t1 * (-t2)
      pQ_pgamma<- pQ_pgamma + 2 * t1 * (-pg_pgamma(xi))
      
      p2Q_pgamma_pthetaT <- p2Q_pgamma_pthetaT + 2 * (t(t2) %*% pg_pgamma(xi) + t1 * (-p2g_pgamma_pthetaT(xi)))
      p2Q_ptheta2 <- p2Q_ptheta2 + 2 * t(t2) %*% t2
      p2Q_pgamma2 <- p2Q_pgamma2 + 2 * (t(pg_pgamma(xi)) %*% pg_pgamma(xi) + t1 * (-p2g_pgamma2(xi)))
    }
    
    Hessian[1:q, 1:q] <- p2Q_ptheta2
    Hessian[1:q, (q+1):(q+p-1)] <- p2Q_pgamma_pthetaT
    Hessian[(q+1):(q+p-1), 1:q] <- t(p2Q_pgamma_pthetaT)
    Hessian[(q+1):(q+p-1), (q+1):(q+p-1)] <- p2Q_pgamma2
    #Hessian <- Hessian + diag(1, p+q-1)
    #print(Hessian)
    #(eigen(Hessian))
    Gradient <- t(cbind(pQ_ptheta, pQ_pgamma))
    #print(Gradient)
    step <- -solve(Hessian) %*% Gradient

    if(max(abs(step)) < 10^(-5))
      return(list(theta_h = theta_h, beta_h = beta_h))
    #print(step)
    zeta <- c(theta_h, gamma) + step
    #print(zeta)
    theta_h <- matrix(zeta[1:q], ncol = 1)
    gamma <- zeta[(q+1):(q+p-1)]
    Jacobi <- rbind(-gamma/sqrt(1-sum(gamma^2)), diag(1, p-1))
    beta_h <- matrix(c(sqrt(1 - sum(gamma^2)), gamma), ncol = 1)
    print(c(theta_h, beta_h))
    #h <- regCVBwSelC(x %*% beta_h, y - z %*% theta_h, deg = 1, kernel = EpaK)# * n^(-2/15)
    #cat('h=', h, '\n')
  }
  return(list(theta_h = theta_h, beta_h = beta_h))
}
NR_iter_c <- compiler::cmpfun(NR_iter)




#estimation of unknown function
g <- function(t, x = X, y = e_YU, y_t = Y_tilde, z = e_ZU, 
              z_t = Z_tilde, theta_h = theta_hat, beta_h = beta_hat, bw = h)
{
  Kh <- function(x)
  {
    EpaK(x/bw) / bw
  }
  
  n <- length(y)
  Kst <- function(t, l1, l2)
  {
    u <- as.numeric(x %*% beta_h)
    t1 <- as.numeric(Kh(u - t))
    t2 <- as.numeric((u - t)^l1)
    t3 <- as.numeric((y - z %*% theta_h)^l2)
    
    sum(t1 * t2 * t3)
  }
  K00 <- function(t) Kst(t, l1=0, l2=0)
  K01 <- function(t) Kst(t, l1=0, l2=1)
  K10 <- function(t) Kst(t, l1=1, l2=0)
  K11 <- function(t) Kst(t, l1=1, l2=1)
  K20 <- function(t) Kst(t, l1=2, l2=0)
  K21 <- function(t) Kst(t, l1=2, l2=1)
  
  A <- K20(t) * K01(t) - K10(t) * K11(t)
  B <- K00(t) * K20(t) - K10(t) * K10(t)
  
  
  A / B + 1/n * sum(y_t - z_t %*% theta)
}












