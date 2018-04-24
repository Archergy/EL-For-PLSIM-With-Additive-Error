library(locpol)

a <- sqrt(3/2) - 1.645/sqrt(12)
b <- sqrt(3/2) + 1.645/sqrt(12)
beta <- matrix(c(1, 1) / sqrt(2), ncol = 1)
theta <- matrix(1, ncol = 1)
p <- length(beta)
q <- length(theta)

#g0 <- function(t) sin(t^2)
#g0 <- function(t) exp(t)
#g0 <- function(t) t^2
g0 <- function(t) sin(pi*(t-a) / (b-a))

myfun <- function(x, z, beta, theta)
{
  g0(x %*% beta) + z %*% theta
}


#generate simulation data
set.seed(3)
N <- 100
X <- matrix(runif(2*N), nrow = N)
Z <- matrix(rnorm(N, 0, 0.4), nrow = N)
#Z <- matrix(rep(c(0,1), times = N/2), nrow = N)
Eps <- rnorm(N, mean = 0, sd = 0.2)
U <- runif(N) #confounding variable
phi <- U^3 - 1/4
psi <- U^2 - 1/3

Y <- myfun(X, Z, beta, theta) + Eps
Z_tilde <- matrix(Z + psi)
Y_tilde <- matrix(Y + phi)


#eliminate the effects caused by additive error
h1 <- sd(U)*(N^(-1/5)) #thumbBw(U, Y_tilde, deg = 0, kernel = gaussK)
E_YU <- locPolSmootherC(U, Y_tilde, xeval = U, bw = h1, deg = 0, kernel = gaussK)$beta0
E_ZU <- locPolSmootherC(U, Z_tilde, xeval = U, bw = h1, deg = 0, kernel = gaussK)$beta0

e_YU <- Y_tilde - E_YU
e_ZU <- Z_tilde - E_ZU


#EL For PLSIM
library(emplik)
source('R\\EL For PLSIM.R')
etaval <- eta(e_YU, X, e_ZU, beta, theta, ker = EpaK)
el.test(etaval, rep(0, ncol(etaval)))[c('-2LLR','Pval')]

k <- 15
beta2 <- seq(0.50, 0.85, length.out = k)
beta1 <- sqrt(1-beta2^2)
theta1 <- seq(0.8, 1.2, length.out = k)
grids <- cbind(beta1, beta2, theta1) #用于绘制置信域的网格点

elgrid <- matrix(0, nrow = k, ncol = k)
elgrid0 <- elgrid

for(i in 1:k){
  for(j in 1:k){

    etaval <- eta(e_YU, X, e_ZU, grids[i,1:2], grids[j,3], ker = EpaK)
    elgrid[i, j] <- el.test(etaval, rep(0, ncol(etaval)))$`-2LLR`
    
    etaval0 <- eta(Y_tilde, X, Z_tilde, grids[i,1:2], grids[j,3], ker = EpaK)
    elgrid0[i, j] <- el.test(etaval0, rep(0, ncol(etaval0)))$`-2LLR`
  }
  print(i)
}

contour(grids[,2], grids[,3], elgrid, nlevels = 1, levels = qchisq(0.95, df = 2),
        drawlabels = F, lty = 1, ylab = expression(theta), xlab = expression(beta[2]))
contour(grids[,2], grids[,3], elgrid0, nlevels = 1, levels = qchisq(0.95, df = 2),
        drawlabels = F, lty = 4, add = T)
points(x = beta[2], y = theta, pch = 3)



#LS For PLSIM
source('R\\LS For PLSIM.R')

lm.fit <- lm(e_YU ~ e_ZU + X) #用线性回归得到参数的初始值
theta0 <-  matrix(lm.fit$coefficients[2], ncol = 1)
beta0 <- matrix(lm.fit$coefficients[3:4], ncol = 1)
beta0 <- beta0 / sqrt(sum(beta0^2))

theta0;beta0
theta;beta
h_star <- regCVBwSelC(X %*% beta0, e_YU - e_ZU %*% theta0, deg = 1, kernel = EpaK)


source("R\\plsi.R")
plsim <- gplsiabc(X, e_YU, e_ZU, h_star)
(beta0 = plsim$beta); (theta0 = plsim$theta)
result <- NR_iter_c(X, e_YU, e_ZU, beta0, theta0)
theta_hat <- result$theta_h
beta_hat <- result$beta_h


h <- regCVBwSelC(X %*% beta_hat, e_YU - e_ZU %*% theta_hat, deg = 1, kernel = EpaK) * N^(-2/15)
loclin <- locPolSmootherC(X %*% beta_hat, e_YU - e_ZU %*% theta_hat, 
                          xeval = X %*% beta_hat, bw = h, deg = 1, kernel = EpaK)
g_star <- loclin$beta0
g_star_d <- loclin$beta1


gamma <- beta_hat[-1, , drop = F]
Jacobi <- rbind(-gamma/sqrt(1-sum(gamma^2)), diag(1, p-1))
Jacobi <- cbind(c(1,rep_len(0,p-1)), Jacobi)


Weight <- sapply(X %*% beta_hat, Wn, x = X, beta = beta_hat, bw = h, tilde = FALSE, ker = EpaK)
mu1 <- t(Weight) %*% X
mu2 <- t(Weight) %*% e_ZU
mu <- cbind(g_star_d * (mu1 %*% Jacobi), mu2)
Lambda <- cbind(g_star_d * (X %*% Jacobi), e_ZU)
temp <- Lambda - mu
epsilon <- as.numeric(e_YU - g_star - e_ZU %*% theta_hat)

V <- t(temp) %*% diag(epsilon^2) %*% temp
V <- V / N


library(ellipse)
plot(ellipse(V[-1,-1], centre = c(gamma, theta_hat), level = 0.95), type = 'l', 
     ylab = expression(theta), xlab = expression(beta[2]), lty = 2)
points(beta[2], theta, pch = 3)
contour(grids[,2], grids[,3], elgrid, nlevels = 1, levels = qchisq(0.95, df = 2), labels = 0.95, 
        col = 1, tcl = 0.5, drawlabels = F, add = T)
legend('topright', legend = "True value", pch = 3)
legend('topleft', legend = c('EL', 'LS'), lty = c(1,2))


#save.image(file = 'R\\s4.RData')
#load(file = 'R\\s4.RData')
