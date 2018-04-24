#leave-one-out cross validation

hgrid <- seq(h_star/5, h_star*5, length.out = 10)
loo_cv <- function(h_grid){
  n <- length(y)
  k <- length(h_grid)
  loo_error <- rep_len(0, k)
  for(j in 1:k){
    for(i in 1:n){
      plsim <- gplsiabc(X[-i,], e_YU[-i], e_ZU[-i, , drop = F], h_grid[j])
      beta0 = plsim$beta
      theta0 = plsim$theta
      t <- X[-i, ] %*% beta0
      g_val <- plsim$g
      spl_fun <- splinefun(t, g_val)
      gi <- e_YU[i] - e_ZU[i, , drop = F] %*% theta0
      loo_error[j] <- loo_error[j] + (spl_fun(X[i,] %*% beta0) - gi)^2
    }
  }
  return(loo_error)
}
#undebug(loo_cv)
cv_error <- loo_cv(h_grid = hgrid)
plot(cv_error,type = 'b')
h_star <- hgrid[which.min(cv_error)]

