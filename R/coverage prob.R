library(locpol)
library(emplik)
source('R\\EL For PLSIM.R')
source('R\\LS For PLSIM.R')
set.seed(123)
nit <- 500
n <- 100
N <- n * nit
X <- matrix(runif(2*N), nrow = N)
Z <- rnorm(N, 0, 0.4)
Eps <- rnorm(N, mean = 0, sd = 0.2)
U <- matrix(runif(N), nrow = N) #confounding variable
phi <- U^4 - 1/5
psi <- U^2 - 1/3

Y <- myfun(X, Z, beta, theta) + Eps
Z_tilde <- matrix(Z + psi, nrow = N)
Y_tilde <- matrix(Y + phi, nrow = N)

count <- 0
count0 <- 0

for(i in 1:nit)
{
  index <- (n * (i - 1) + 1) : (n * i)
  Y_tildei <- Y_tilde[index, , drop = F]
  Z_tildei <- Z_tilde[index, , drop = F]
  Ui <- U[index, , drop = F]
  Xi <- X[index, , drop = F]
  
  h1 <- sd(Ui)*(N^(-1/5)) # thumbBw(Ui, Y_tildei, deg = 0, kernel = gaussK)
  E_YU <- locPolSmootherC(Ui, Y_tildei, xeval = Ui, bw = h1, deg = 0, kernel = gaussK)$beta0
  E_ZU <- locPolSmootherC(Ui, Z_tildei, xeval = Ui, bw = h1, deg = 0, kernel = gaussK)$beta0
  
  e_YU <- Y_tildei - E_YU
  e_ZU <- Z_tildei - E_ZU
  
  etaval <- eta(e_YU, Xi, e_ZU, beta, theta, ker = EpaK)
  etaval0 <- eta(Y_tildei, Xi, Z_tildei, beta, theta, ker = EpaK)
  if(el.test(etaval, rep(0, ncol(etaval)))$Pval > 0.05) count <- count + 1
  if(el.test(etaval0, rep(0, ncol(etaval0)))$Pval > 0.05) count0 <- count0 + 1
  if(i %% 10 == 0) print(i)
}

count/nit
count0/nit
