# Prob thet each diag is greater than 1
set_diag <- function(lower_bound = 1, prob = 0.95){
  target <- optimise(f = function(lambda){abs((1 - pexp(lower_bound, lambda/2)) - prob)}, lower = 0, upper = 2)
  return(list("lambda" = target$minimum, "Criteria" = target$objective))
}

# Prob that all diags are greater than 1
#set_diag_all <- function(p, lower_bound = 1, prob = 0.95){
#  target <- optimise(f = function(lambda){abs((1 - pexp(lower_bound, lambda/2))^p - prob)}, lower = 0, upper = 2)
#  return(list("lambda" = target$minimum, "Criteria" = target$objective))
#}

#cdf is $1-exp(-\lambda/2 x)$ so we want that $(1-1-exp(-\lambda/2))^p = exp(-lambda/2 p) = 0.95$ which means that $lambda = -2log(0.95)/p$ which is correct!
set_diag_one <- function(lower_bound = 1, prob = 0.95){
  lambda <- -2/(lower_bound)*log(prob)
  return(lambda)
}


set_diag_all <- function(p, lower_bound = 1, prob = 0.95){
  lambda <- -2/(lower_bound*p)*log(prob)
  return(lambda)
}


# Laplace spike-and-slab generation
rSSGaussian <- function(n, w_slab, s_slab){
  slab <- sample(c(0, 1), n, replace = TRUE, prob = c(1-w_slab, w_slab))
  return((1-slab)*0 + slab*rnorm(n, 0, s_slab))
}

rOmega_diagExp_SSGaussian <- function(N, p, w_slab, s_slab, lambda, truncDiag = FALSE){
  Omega <- array(0, dim = c(N, p, p))
  if(truncDiag){
    l <- pexp(1, lambda/2)
  }
  for(n in 1:N){
    if(truncDiag){
      u <- runif(p, min = l, max = 1)
      diag(Omega[n,,]) <- qexp(u, lambda/2) / 2 # /2 to make the symmerty easier
    } else{
      diag(Omega[n,,]) <- rexp(p, lambda/2) / 2 # /2 to make the symmerty easier
    }
    Omega[n,,][upper.tri(Omega[n,,], diag = FALSE)] <- rSSGaussian(p*(p-1)/2, w_slab, s_slab)
    Omega[n,,] <- Omega[n,,] + t(Omega[n,,])
  }
  return(Omega)
}

probPD_diagExp_SSGaussian <- function(N, p, w_slab, s_slab, lambda, truncDiag = FALSE, verbose = FALSE){
  pd_indicator <- rep(NA, N)
  for(n in 1:N){
    Omega <- rOmega_diagExp_SSGaussian(N = 1, p, w_slab, s_slab, lambda, truncDiag)
    pd_indicator[n] <- min(eigen(Omega[1,,])$values) > 0
    if(verbose & (n %% (N/10)) == 0){
      cat("\r", n/N*100, "% done")
      flush.console() 
      #https://stackoverflow.com/questions/5953718/overwrite-current-output-in-the-r-console
    }
  }
  return(mean(pd_indicator))
}

# maximum s_slab achieveing certain probPD
gridsearch_probPD_diagExp_SSGaussian <- function(N, p, w_slab, s_slab_grid, lambda, probPD, truncDiag = FALSE, verbose = FALSE){
  N_try <- length(s_slab_grid)
  probPD_eval <- rep(NA, N_try)
  for(i in 1:N_try){
    probPD_eval[i] <- probPD_diagExp_SSGaussian(N, p, w_slab, s_slab = s_slab_grid[i], lambda, truncDiag, verbose)
    if(verbose){
      cat("\r", "try", i, "done", "\n")
    }
    
  }
  return(list("s_slab" = s_slab_grid[which.min(abs(probPD_eval - probPD))], "criteria" = min(abs(probPD_eval - probPD))))
}


