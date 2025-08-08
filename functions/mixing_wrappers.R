
ssgraph_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda, s_spike = 2e-02){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  ssgraph_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  ssgraph_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  ssgraph_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  ssgraph_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  ssgraph_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  for(n in 1:length(iterations_vec)){
    
    ## ssgraph
    ### init 1
    timer_ssgraph1.start <- Sys.time()
    fit.ssgraph1 <- ssgraph(data_y, method="ggm", iter=iterations_vec_ssgraph[n], burnin = warmup, g.prior= w_slab
                            , var1 = s_spike^2, var2 = s_slab^2
                            , lambda = lambda, save=TRUE, verbose=TRUE, sig.start = Sigma_init1)
    timer_ssgraph1.end <- Sys.time()
    ### init 2
    timer_ssgraph2.start <- Sys.time()
    fit.ssgraph2 <- ssgraph(data_y, method="ggm", iter=iterations_vec_ssgraph[n], burnin = warmup, g.prior= w_slab
                            , var1 = s_spike^2, var2 = s_slab^2
                            , lambda = lambda, save=TRUE, verbose=TRUE, sig.start = Sigma_init2)
    timer_ssgraph2.end <- Sys.time()
    
    ssgraph_time[n, 1] <- difftime(timer_ssgraph1.end, timer_ssgraph1.start, units = "secs")
    ssgraph_time[n, 2] <- difftime(timer_ssgraph2.end, timer_ssgraph2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(fit.ssgraph1$K_hat - Theta_sim)
    ssgraph_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.ssgraph1$p_links)
    ssgraph_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.p)
    
    ## chain 2
    e2 <- abs(fit.ssgraph2$K_hat - Theta_sim)
    ssgraph_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.ssgraph2$p_links)
    ssgraph_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.p)
    
    ## between chains
    e12 <- abs(fit.ssgraph2$K_hat - fit.ssgraph1$K_hat)
    ssgraph_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.ssgraph1$p_links - fit.ssgraph2$p_links)
    ssgraph_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.p)
    
  }
  
  return(list(
    "Theta_MAE" = ssgraph_Theta_MAE,
    "p_MAE" = ssgraph_p_MAE,
    "chain_Theta" = ssgraph_chain_Theta_MAE,
    "chain_p" = ssgraph_chain_p_MAE,
    "Time" = ssgraph_time
  ))
  
}

bdgraph_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab){
  
  bdgraph_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  bdgraph_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  bdgraph_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  bdgraph_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  bdgraph_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  for(n in 1:length(iterations_vec)){
    
    ## bdgraph
    ### init 1
    timer_bdgraph1.start <- Sys.time()
    fit.bdgraph1 <- bdgraph(data_y, method="ggm", iter=iterations_vec[n]*ceiling(sqrt(p)), burnin = warmup, g.prior= w_slab, save=FALSE, verbose=TRUE, g.start = Omega_init1!=0, df.prior = 3)
    timer_bdgraph1.end <- Sys.time()
    ### init 2
    timer_bdgraph2.start <- Sys.time()
    fit.bdgraph2 <- bdgraph(data_y, method="ggm", iter=iterations_vec[n]*ceiling(sqrt(p)), burnin = warmup, g.prior= w_slab, save=FALSE, verbose=TRUE, g.start = Omega_init2!=0, df.prior = 3)
    timer_bdgraph2.end <- Sys.time()
    
    bdgraph_time[n, 1] <- difftime(timer_bdgraph1.end, timer_bdgraph1.start, units = "secs")
    bdgraph_time[n, 2] <- difftime(timer_bdgraph2.end, timer_bdgraph2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(fit.bdgraph1$K_hat - Theta_sim)
    bdgraph_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.bdgraph1$p_links)
    bdgraph_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.p)
    
    ## chain 2
    e2 <- abs(fit.bdgraph2$K_hat - Theta_sim)
    bdgraph_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.bdgraph2$p_links)
    bdgraph_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.p)
    
    ## between chains
    e12 <- abs(fit.bdgraph2$K_hat - fit.bdgraph1$K_hat)
    bdgraph_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.bdgraph1$p_links - fit.bdgraph2$p_links)
    bdgraph_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.p)
    
  }
  
  return(list(
    "Theta_MAE" = bdgraph_Theta_MAE,
    "p_MAE" = bdgraph_p_MAE,
    "chain_Theta" = bdgraph_chain_Theta_MAE,
    "chain_p" = bdgraph_chain_p_MAE,
    "Time" = bdgraph_time
  ))
  
}

bdgraph.mpl_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab){
  
  
  bdgraph.mpl_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  bdgraph.mpl_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  bdgraph.mpl_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  bdgraph.mpl_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  bdgraph.mpl_time <- array(NA, dim = c(length(iterations_vec), 2))
  bdgraph.mpl.omega_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  S <- t(data_y)%*%data_y
  n_train <- nrow(data_y)
  
  for(n in 1:length(iterations_vec)){
    
    ## bdgraph.mpl
    ### init 1
    timer_bdgraph.mpl1.start <- Sys.time()
    fit.bdgraph.mpl1 <- bdgraph.mpl(data_y, method="ggm", iter=iterations_vec[n]*p, burnin = warmup, g.prior= w_slab, save=TRUE, verbose=TRUE, g.start = Omega_init1!=0)
    timer_bdgraph.mpl1.end <- Sys.time()
    ### init 2
    timer_bdgraph.mpl2.start <- Sys.time()
    fit.bdgraph.mpl2 <- bdgraph.mpl(data_y, method="ggm", iter=iterations_vec[n]*p, burnin = warmup, g.prior= w_slab, save=TRUE, verbose=TRUE, g.start = Omega_init2!=0)
    timer_bdgraph.mpl2.end <- Sys.time()
    
    bdgraph.mpl_time[n, 1] <- difftime(timer_bdgraph.mpl1.end, timer_bdgraph.mpl1.start, units = "secs")
    bdgraph.mpl_time[n, 2] <- difftime(timer_bdgraph.mpl2.end, timer_bdgraph.mpl2.start, units = "secs")
    
    fit.bdgraph.mpl1$p_links <- plinks(fit.bdgraph.mpl1, round = 3)
    fit.bdgraph.mpl2$p_links <- plinks(fit.bdgraph.mpl2, round = 3)
    
    # Sampling Omegas - NEED TO SET save = TRUE, had save = FALSE previously
    ## have reduced iterations by p here
    timer_bdgraph.mpl.omega1.start <- Sys.time()
    fit.bdgraph.mpl1$K_hat <- bdgraph.mpl_sampleOmega(fit.bdgraph.mpl1, iter=iterations_vec[n], p, S, n_train, store_samples = FALSE)$Omega_hat
    timer_bdgraph.mpl.omega1.end<- Sys.time()
    
    timer_bdgraph.mpl.omega2.start <- Sys.time()
    fit.bdgraph.mpl2$K_hat <- bdgraph.mpl_sampleOmega(fit.bdgraph.mpl2, iter=iterations_vec[n], p, S, n_train, store_samples = FALSE)$Omega_hat
    timer_bdgraph.mpl.omega2.end<- Sys.time()
    
    bdgraph.mpl.omega_time[n, 1] <- difftime(timer_bdgraph.mpl.omega1.end, timer_bdgraph.mpl.omega1.start, units = "secs")
    bdgraph.mpl.omega_time[n, 2] <- difftime(timer_bdgraph.mpl.omega2.end, timer_bdgraph.mpl.omega2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(fit.bdgraph.mpl1$K_hat - Theta_sim)
    bdgraph.mpl_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.bdgraph.mpl1$p_links)
    bdgraph.mpl_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.p)
    
    ## chain 2
    e2 <- abs(fit.bdgraph.mpl2$K_hat - Theta_sim)
    bdgraph.mpl_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.bdgraph.mpl2$p_links)
    bdgraph.mpl_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.p)
    
    ## between chains
    e12 <- abs(fit.bdgraph.mpl2$K_hat - fit.bdgraph.mpl1$K_hat)
    bdgraph.mpl_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.bdgraph.mpl1$p_links - fit.bdgraph.mpl2$p_links)
    bdgraph.mpl_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.p)
    
  }
  
  return(list(
    "Theta_MAE" = bdgraph.mpl_Theta_MAE,
    "p_MAE" = bdgraph.mpl_p_MAE,
    "chain_Theta" = bdgraph.mpl_chain_Theta_MAE,
    "chain_p" = bdgraph.mpl_chain_p_MAE,
    "Time" = bdgraph.mpl_time,
    "Time.omega" = bdgraph.mpl.omega_time
  ))
  
}

Gibbs_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  Gibbs_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  Gibbs_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  Gibbs_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  Gibbs_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  Gibbs_noRB_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  Gibbs_noRB_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  Gibbs_noRB_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  Gibbs_noRB_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  Gibbs_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  updates_per_iter <- p
  updates_per_column <- p
  
  for(n in 1:length(iterations_vec)){
    
    ## Gibbs
    ### init 1
    timer_Gibbs1.start <- Sys.time()
    fit.Gibbs1 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init1, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda))
    timer_Gibbs1.end <- Sys.time()
    
    b <- coef(fit.Gibbs1)
    fit.Gibbs1$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.Gibbs1$p_links[upper.tri(fit.Gibbs1$p_links, diag = TRUE)] <- b[,6]
    fit.Gibbs1$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.Gibbs1$p_links_noRB[upper.tri(fit.Gibbs1$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.Gibbs1$postSample) != 0)
    fit.Gibbs1$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.Gibbs1$postmean_noRB[upper.tri(fit.Gibbs1$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.Gibbs1$postSample))
    
    rm(b)
    
    ### init 2
    timer_Gibbs2.start <- Sys.time()
    fit.Gibbs2 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init2, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda))
    timer_Gibbs2.end <- Sys.time()
    
    b <- coef(fit.Gibbs2)
    fit.Gibbs2$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.Gibbs2$p_links[upper.tri(fit.Gibbs2$p_links, diag = TRUE)] <- b[,6]
    fit.Gibbs2$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.Gibbs2$p_links_noRB[upper.tri(fit.Gibbs2$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.Gibbs2$postSample) != 0)
    fit.Gibbs2$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.Gibbs2$postmean_noRB[upper.tri(fit.Gibbs2$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.Gibbs2$postSample))
    
    rm(b)
    
    Gibbs_time[n, 1] <- difftime(timer_Gibbs1.end, timer_Gibbs1.start, units = "secs")
    Gibbs_time[n, 2] <- difftime(timer_Gibbs2.end, timer_Gibbs2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(icov(fit.Gibbs1) - Theta_sim)
    Gibbs_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.noRB <- abs(fit.Gibbs1$postmean_noRB - Theta_sim)
    Gibbs_noRB_Theta_MAE[n, 1] <- mean(e1.noRB[upper.tri(e1.noRB, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.Gibbs1$p_links)
    Gibbs_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    e1.p.noRB <- abs((Theta_sim!=0) - fit.Gibbs1$p_links_noRB)
    Gibbs_noRB_p_MAE[n, 1] <- mean(e1.p.noRB[upper.tri(e1.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.noRB)
    rm(e1.p)
    rm(e1.p.noRB)
    
    ## chain 2
    e2 <- abs(icov(fit.Gibbs2) - Theta_sim)
    Gibbs_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.noRB <- abs(fit.Gibbs2$postmean_noRB - Theta_sim)
    Gibbs_noRB_Theta_MAE[n, 2] <- mean(e2.noRB[upper.tri(e2.noRB, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.Gibbs2$p_links)
    Gibbs_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    e2.p.noRB <- abs((Theta_sim!=0) - fit.Gibbs2$p_links_noRB)
    Gibbs_noRB_p_MAE[n, 2] <- mean(e2.p.noRB[upper.tri(e2.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.noRB)
    rm(e2.p)
    rm(e2.p.noRB)
    
    ## between chains
    e12 <- abs(icov(fit.Gibbs2) - icov(fit.Gibbs1))
    Gibbs_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.noRB <- abs(fit.Gibbs1$postmean_noRB - fit.Gibbs2$postmean_noRB )
    Gibbs_noRB_chain_Theta_MAE[n] <- mean(e12.noRB[upper.tri(e12.noRB, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.Gibbs1$p_links - fit.Gibbs2$p_links)
    Gibbs_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    e12.p.noRB <- abs(fit.Gibbs1$p_links_noRB - fit.Gibbs2$p_links_noRB)
    Gibbs_noRB_chain_p_MAE[n] <- mean(e12.p.noRB[upper.tri(e12.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.noRB)
    rm(e12.p)
    rm(e12.p.noRB)
    
  }
  
  return(list(
    "Theta_MAE" = Gibbs_Theta_MAE,
    "p_MAE" = Gibbs_p_MAE,
    "chain_Theta" = Gibbs_chain_Theta_MAE,
    "chain_p" = Gibbs_chain_p_MAE,
    "Theta_MAE_noRB" = Gibbs_noRB_Theta_MAE,
    "p_MAE_noRB" = Gibbs_noRB_p_MAE,
    "chain_Theta_noRB" = Gibbs_noRB_chain_Theta_MAE,
    "chain_p_noRB" = Gibbs_noRB_chain_p_MAE,
    "Time" = Gibbs_time
  ))
  
}

BDMH_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  BDMH_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  BDMH_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  BDMH_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  BDMH_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  BDMH_noRB_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  BDMH_noRB_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  BDMH_noRB_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  BDMH_noRB_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  BDMH_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  updates_per_iter <- p
  updates_per_column <- ceiling(sqrt(p))
  
  pbirth <- 0.75
  pdeath <- 0.5*(1-pbirth)
  
  for(n in 1:length(iterations_vec)){
    
    ## BDMH
    ### init 1
    timer_BDMH1.start <- Sys.time()
    fit.BDMH1 <- modelSelectionGGM(data_y, sampler='birthdeath', Omegaini=Omega_init1, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), pbirth=pbirth, pdeath=pdeath)
    timer_BDMH1.end <- Sys.time()
    
    b <- coef(fit.BDMH1)
    fit.BDMH1$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.BDMH1$p_links[upper.tri(fit.BDMH1$p_links, diag = TRUE)] <- b[,6]
    fit.BDMH1$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.BDMH1$p_links_noRB[upper.tri(fit.BDMH1$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.BDMH1$postSample) != 0)
    fit.BDMH1$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.BDMH1$postmean_noRB[upper.tri(fit.BDMH1$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.BDMH1$postSample))
    
    rm(b)
    
    ### init 2
    timer_BDMH2.start <- Sys.time()
    fit.BDMH2 <- modelSelectionGGM(data_y, sampler='birthdeath', Omegaini=Omega_init2, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), pbirth=pbirth, pdeath=pdeath)
    timer_BDMH2.end <- Sys.time()
    
    b <- coef(fit.BDMH2)
    fit.BDMH2$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.BDMH2$p_links[upper.tri(fit.BDMH2$p_links, diag = TRUE)] <- b[,6]
    fit.BDMH2$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.BDMH2$p_links_noRB[upper.tri(fit.BDMH2$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.BDMH2$postSample) != 0)
    fit.BDMH2$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.BDMH2$postmean_noRB[upper.tri(fit.BDMH2$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.BDMH2$postSample))
    
    rm(b)
    
    BDMH_time[n, 1] <- difftime(timer_BDMH1.end, timer_BDMH1.start, units = "secs")
    BDMH_time[n, 2] <- difftime(timer_BDMH2.end, timer_BDMH2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(icov(fit.BDMH1) - Theta_sim)
    BDMH_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.noRB <- abs(fit.BDMH1$postmean_noRB - Theta_sim)
    BDMH_noRB_Theta_MAE[n, 1] <- mean(e1.noRB[upper.tri(e1.noRB, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.BDMH1$p_links)
    BDMH_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    e1.p.noRB <- abs((Theta_sim!=0) - fit.BDMH1$p_links_noRB)
    BDMH_noRB_p_MAE[n, 1] <- mean(e1.p.noRB[upper.tri(e1.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.noRB)
    rm(e1.p)
    rm(e1.p.noRB)
    
    ## chain 2
    e2 <- abs(icov(fit.BDMH2) - Theta_sim)
    BDMH_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.noRB <- abs(fit.BDMH2$postmean_noRB - Theta_sim)
    BDMH_noRB_Theta_MAE[n, 2] <- mean(e2.noRB[upper.tri(e2.noRB, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.BDMH2$p_links)
    BDMH_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    e2.p.noRB <- abs((Theta_sim!=0) - fit.BDMH2$p_links_noRB)
    BDMH_noRB_p_MAE[n, 2] <- mean(e2.p.noRB[upper.tri(e2.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.noRB)
    rm(e2.p)
    rm(e2.p.noRB)
    
    ## between chains
    e12 <- abs(icov(fit.BDMH2) - icov(fit.BDMH1))
    BDMH_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.noRB <- abs(fit.BDMH1$postmean_noRB - fit.BDMH2$postmean_noRB )
    BDMH_noRB_chain_Theta_MAE[n] <- mean(e12.noRB[upper.tri(e12.noRB, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.BDMH1$p_links - fit.BDMH2$p_links)
    BDMH_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    e12.p.noRB <- abs(fit.BDMH1$p_links_noRB - fit.BDMH2$p_links_noRB)
    BDMH_noRB_chain_p_MAE[n] <- mean(e12.p.noRB[upper.tri(e12.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.noRB)
    rm(e12.p)
    rm(e12.p.noRB)
    
  }
  
  return(list(
    "Theta_MAE" = BDMH_Theta_MAE,
    "p_MAE" = BDMH_p_MAE,
    "chain_Theta" = BDMH_chain_Theta_MAE,
    "chain_p" = BDMH_chain_p_MAE,
    "Theta_MAE_noRB" = BDMH_noRB_Theta_MAE,
    "p_MAE_noRB" = BDMH_noRB_p_MAE,
    "chain_Theta_noRB" = BDMH_noRB_chain_Theta_MAE,
    "chain_p_noRB" = BDMH_noRB_chain_p_MAE,
    "Time" = BDMH_time
  ))
  
}

LIT_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  LIT_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  LIT_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  LIT_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  LIT_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  LIT_noRB_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  LIT_noRB_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  LIT_noRB_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  LIT_noRB_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  LIT_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  updates_per_iter <- p
  updates_per_column <- 1
  
  for(n in 1:length(iterations_vec)){
    
    ## LIT
    ### init 1
    timer_LIT1.start <- Sys.time()
    fit.LIT1 <- modelSelectionGGM(data_y, sampler='LIT', Omegaini=Omega_init1, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda))
    timer_LIT1.end <- Sys.time()
    
    b <- coef(fit.LIT1)
    fit.LIT1$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.LIT1$p_links[upper.tri(fit.LIT1$p_links, diag = TRUE)] <- b[,6]
    fit.LIT1$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.LIT1$p_links_noRB[upper.tri(fit.LIT1$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.LIT1$postSample) != 0)
    fit.LIT1$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.LIT1$postmean_noRB[upper.tri(fit.LIT1$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.LIT1$postSample))
    
    rm(b)
    
    ### init 2
    timer_LIT2.start <- Sys.time()
    fit.LIT2 <- modelSelectionGGM(data_y, sampler='LIT', Omegaini=Omega_init2, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda))
    timer_LIT2.end <- Sys.time()
    
    b <- coef(fit.LIT2)
    fit.LIT2$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.LIT2$p_links[upper.tri(fit.LIT2$p_links, diag = TRUE)] <- b[,6]
    fit.LIT2$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.LIT2$p_links_noRB[upper.tri(fit.LIT2$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.LIT2$postSample) != 0)
    fit.LIT2$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.LIT2$postmean_noRB[upper.tri(fit.LIT2$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.LIT2$postSample))
    
    rm(b)
    
    LIT_time[n, 1] <- difftime(timer_LIT1.end, timer_LIT1.start, units = "secs")
    LIT_time[n, 2] <- difftime(timer_LIT2.end, timer_LIT2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(icov(fit.LIT1) - Theta_sim)
    LIT_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.noRB <- abs(fit.LIT1$postmean_noRB - Theta_sim)
    LIT_noRB_Theta_MAE[n, 1] <- mean(e1.noRB[upper.tri(e1.noRB, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.LIT1$p_links)
    LIT_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    e1.p.noRB <- abs((Theta_sim!=0) - fit.LIT1$p_links_noRB)
    LIT_noRB_p_MAE[n, 1] <- mean(e1.p.noRB[upper.tri(e1.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.noRB)
    rm(e1.p)
    rm(e1.p.noRB)
    
    ## chain 2
    e2 <- abs(icov(fit.LIT2) - Theta_sim)
    LIT_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.noRB <- abs(fit.LIT2$postmean_noRB - Theta_sim)
    LIT_noRB_Theta_MAE[n, 2] <- mean(e2.noRB[upper.tri(e2.noRB, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.LIT2$p_links)
    LIT_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    e2.p.noRB <- abs((Theta_sim!=0) - fit.LIT2$p_links_noRB)
    LIT_noRB_p_MAE[n, 2] <- mean(e2.p.noRB[upper.tri(e2.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.noRB)
    rm(e2.p)
    rm(e2.p.noRB)
    
    ## between chains
    e12 <- abs(icov(fit.LIT2) - icov(fit.LIT1))
    LIT_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.noRB <- abs(fit.LIT1$postmean_noRB - fit.LIT2$postmean_noRB )
    LIT_noRB_chain_Theta_MAE[n] <- mean(e12.noRB[upper.tri(e12.noRB, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.LIT1$p_links - fit.LIT2$p_links)
    LIT_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    e12.p.noRB <- abs(fit.LIT1$p_links_noRB - fit.LIT2$p_links_noRB)
    LIT_noRB_chain_p_MAE[n] <- mean(e12.p.noRB[upper.tri(e12.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.noRB)
    rm(e12.p)
    rm(e12.p.noRB)
    
  }
  
  return(list(
    "Theta_MAE" = LIT_Theta_MAE,
    "p_MAE" = LIT_p_MAE,
    "chain_Theta" = LIT_chain_Theta_MAE,
    "chain_p" = LIT_chain_p_MAE,
    "Theta_MAE_noRB" = LIT_noRB_Theta_MAE,
    "p_MAE_noRB" = LIT_noRB_p_MAE,
    "chain_Theta_noRB" = LIT_noRB_chain_Theta_MAE,
    "chain_p_noRB" = LIT_noRB_chain_p_MAE,
    "Time" = LIT_time
  ))
  
}

GIMH_Gibbs_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda, tempering = 0.75){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  GIMH_Gibbs_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_Gibbs_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_Gibbs_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  GIMH_Gibbs_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  GIMH_Gibbs_noRB_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_Gibbs_noRB_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_Gibbs_noRB_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  GIMH_Gibbs_noRB_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  GIMH_Gibbs_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  updates_per_iter <- p#
  updates_per_column <- ceiling(sqrt(p))
  prob_global <- 1 # so eevry column udpate is a full row form the parallel
  truncratio <- -1 # presumbaly this means no truncation?
  
  for(n in 1:length(iterations_vec)){
    
    ## GIMH_Gibbs
    ### init 1
    timer_GIMH_Gibbs1.start <- Sys.time()
    fit.GIMH_Gibbs1 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init1, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='regression', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), prob_global=prob_global, tempering=tempering, truncratio=truncratio)
    timer_GIMH_Gibbs1.end <- Sys.time()
    
    b <- coef(fit.GIMH_Gibbs1)
    fit.GIMH_Gibbs1$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_Gibbs1$p_links[upper.tri(fit.GIMH_Gibbs1$p_links, diag = TRUE)] <- b[,6]
    fit.GIMH_Gibbs1$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_Gibbs1$p_links_noRB[upper.tri(fit.GIMH_Gibbs1$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_Gibbs1$postSample) != 0)
    fit.GIMH_Gibbs1$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_Gibbs1$postmean_noRB[upper.tri(fit.GIMH_Gibbs1$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_Gibbs1$postSample))
    
    rm(b)
    
    ### init 2
    timer_GIMH_Gibbs2.start <- Sys.time()
    fit.GIMH_Gibbs2 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init2, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='regression', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), prob_global=prob_global, tempering=tempering, truncratio=truncratio)
    timer_GIMH_Gibbs2.end <- Sys.time()
    
    b <- coef(fit.GIMH_Gibbs2)
    fit.GIMH_Gibbs2$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_Gibbs2$p_links[upper.tri(fit.GIMH_Gibbs2$p_links, diag = TRUE)] <- b[,6]
    fit.GIMH_Gibbs2$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_Gibbs2$p_links_noRB[upper.tri(fit.GIMH_Gibbs2$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_Gibbs2$postSample) != 0)
    fit.GIMH_Gibbs2$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_Gibbs2$postmean_noRB[upper.tri(fit.GIMH_Gibbs2$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_Gibbs2$postSample))
    
    rm(b)
    
    GIMH_Gibbs_time[n, 1] <- difftime(timer_GIMH_Gibbs1.end, timer_GIMH_Gibbs1.start, units = "secs")
    GIMH_Gibbs_time[n, 2] <- difftime(timer_GIMH_Gibbs2.end, timer_GIMH_Gibbs2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(icov(fit.GIMH_Gibbs1) - Theta_sim)
    GIMH_Gibbs_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.noRB <- abs(fit.GIMH_Gibbs1$postmean_noRB - Theta_sim)
    GIMH_Gibbs_noRB_Theta_MAE[n, 1] <- mean(e1.noRB[upper.tri(e1.noRB, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.GIMH_Gibbs1$p_links)
    GIMH_Gibbs_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    e1.p.noRB <- abs((Theta_sim!=0) - fit.GIMH_Gibbs1$p_links_noRB)
    GIMH_Gibbs_noRB_p_MAE[n, 1] <- mean(e1.p.noRB[upper.tri(e1.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.noRB)
    rm(e1.p)
    rm(e1.p.noRB)
    
    ## chain 2
    e2 <- abs(icov(fit.GIMH_Gibbs2) - Theta_sim)
    GIMH_Gibbs_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.noRB <- abs(fit.GIMH_Gibbs2$postmean_noRB - Theta_sim)
    GIMH_Gibbs_noRB_Theta_MAE[n, 2] <- mean(e2.noRB[upper.tri(e2.noRB, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.GIMH_Gibbs2$p_links)
    GIMH_Gibbs_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    e2.p.noRB <- abs((Theta_sim!=0) - fit.GIMH_Gibbs2$p_links_noRB)
    GIMH_Gibbs_noRB_p_MAE[n, 2] <- mean(e2.p.noRB[upper.tri(e2.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.noRB)
    rm(e2.p)
    rm(e2.p.noRB)
    
    ## between chains
    e12 <- abs(icov(fit.GIMH_Gibbs2) - icov(fit.GIMH_Gibbs1))
    GIMH_Gibbs_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.noRB <- abs(fit.GIMH_Gibbs1$postmean_noRB - fit.GIMH_Gibbs2$postmean_noRB )
    GIMH_Gibbs_noRB_chain_Theta_MAE[n] <- mean(e12.noRB[upper.tri(e12.noRB, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.GIMH_Gibbs1$p_links - fit.GIMH_Gibbs2$p_links)
    GIMH_Gibbs_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    e12.p.noRB <- abs(fit.GIMH_Gibbs1$p_links_noRB - fit.GIMH_Gibbs2$p_links_noRB)
    GIMH_Gibbs_noRB_chain_p_MAE[n] <- mean(e12.p.noRB[upper.tri(e12.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.noRB)
    rm(e12.p)
    rm(e12.p.noRB)
    
  }
  
  return(list(
    "Theta_MAE" = GIMH_Gibbs_Theta_MAE,
    "p_MAE" = GIMH_Gibbs_p_MAE,
    "chain_Theta" = GIMH_Gibbs_chain_Theta_MAE,
    "chain_p" = GIMH_Gibbs_chain_p_MAE,
    "Theta_MAE_noRB" = GIMH_Gibbs_noRB_Theta_MAE,
    "p_MAE_noRB" = GIMH_Gibbs_noRB_p_MAE,
    "chain_Theta_noRB" = GIMH_Gibbs_noRB_chain_Theta_MAE,
    "chain_p_noRB" = GIMH_Gibbs_noRB_chain_p_MAE,
    "Time" = GIMH_Gibbs_time
  ))
  
}

GIMH_BDMH_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda, tempering = 0.75){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  GIMH_BDMH_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_BDMH_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_BDMH_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  GIMH_BDMH_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  GIMH_BDMH_noRB_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_BDMH_noRB_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  GIMH_BDMH_noRB_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  GIMH_BDMH_noRB_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  GIMH_BDMH_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  updates_per_iter <- p#
  updates_per_column <- ceiling(sqrt(p))
  prob_global <- 1 # so eevry column udpate is a full row form the parallel
  truncratio <- -1 # presumbaly this means no truncation?
  
  pbirth <- 0.75
  pdeath <- 0.5*(1-pbirth)
  
  for(n in 1:length(iterations_vec)){
    
    ## GIMH_BDMH
    ### init 1
    timer_GIMH_BDMH1.start <- Sys.time()
    fit.GIMH_BDMH1 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init1, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='regression', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), prob_global=prob_global, tempering=tempering, truncratio=truncratio, pbirth=pbirth, pdeath=pdeath)
    timer_GIMH_BDMH1.end <- Sys.time()
    
    b <- coef(fit.GIMH_BDMH1)
    fit.GIMH_BDMH1$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_BDMH1$p_links[upper.tri(fit.GIMH_BDMH1$p_links, diag = TRUE)] <- b[,6]
    fit.GIMH_BDMH1$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_BDMH1$p_links_noRB[upper.tri(fit.GIMH_BDMH1$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_BDMH1$postSample) != 0)
    fit.GIMH_BDMH1$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_BDMH1$postmean_noRB[upper.tri(fit.GIMH_BDMH1$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_BDMH1$postSample))
    
    rm(b)
    
    ### init 2
    timer_GIMH_BDMH2.start <- Sys.time()
    fit.GIMH_BDMH2 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init2, niter=iterations_vec[n], burnin=warmup, scale=FALSE, global_proposal='regression', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), prob_global=prob_global, tempering=tempering, truncratio=truncratio, pbirth=pbirth, pdeath=pdeath)
    timer_GIMH_BDMH2.end <- Sys.time()
    
    b <- coef(fit.GIMH_BDMH2)
    fit.GIMH_BDMH2$p_links <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_BDMH2$p_links[upper.tri(fit.GIMH_BDMH2$p_links, diag = TRUE)] <- b[,6]
    fit.GIMH_BDMH2$p_links_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_BDMH2$p_links_noRB[upper.tri(fit.GIMH_BDMH2$p_links_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_BDMH2$postSample) != 0)
    fit.GIMH_BDMH2$postmean_noRB <- matrix(NA, nrow = p, ncol = p)
    fit.GIMH_BDMH2$postmean_noRB[upper.tri(fit.GIMH_BDMH2$postmean_noRB, diag = TRUE)] <- colMeans(as.matrix(fit.GIMH_BDMH2$postSample))
    
    rm(b)
    
    GIMH_BDMH_time[n, 1] <- difftime(timer_GIMH_BDMH1.end, timer_GIMH_BDMH1.start, units = "secs")
    GIMH_BDMH_time[n, 2] <- difftime(timer_GIMH_BDMH2.end, timer_GIMH_BDMH2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(icov(fit.GIMH_BDMH1) - Theta_sim)
    GIMH_BDMH_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.noRB <- abs(fit.GIMH_BDMH1$postmean_noRB - Theta_sim)
    GIMH_BDMH_noRB_Theta_MAE[n, 1] <- mean(e1.noRB[upper.tri(e1.noRB, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.GIMH_BDMH1$p_links)
    GIMH_BDMH_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    e1.p.noRB <- abs((Theta_sim!=0) - fit.GIMH_BDMH1$p_links_noRB)
    GIMH_BDMH_noRB_p_MAE[n, 1] <- mean(e1.p.noRB[upper.tri(e1.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.noRB)
    rm(e1.p)
    rm(e1.p.noRB)
    
    ## chain 2
    e2 <- abs(icov(fit.GIMH_BDMH2) - Theta_sim)
    GIMH_BDMH_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.noRB <- abs(fit.GIMH_BDMH2$postmean_noRB - Theta_sim)
    GIMH_BDMH_noRB_Theta_MAE[n, 2] <- mean(e2.noRB[upper.tri(e2.noRB, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.GIMH_BDMH2$p_links)
    GIMH_BDMH_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    e2.p.noRB <- abs((Theta_sim!=0) - fit.GIMH_BDMH2$p_links_noRB)
    GIMH_BDMH_noRB_p_MAE[n, 2] <- mean(e2.p.noRB[upper.tri(e2.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.noRB)
    rm(e2.p)
    rm(e2.p.noRB)
    
    ## between chains
    e12 <- abs(icov(fit.GIMH_BDMH2) - icov(fit.GIMH_BDMH1))
    GIMH_BDMH_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.noRB <- abs(fit.GIMH_BDMH1$postmean_noRB - fit.GIMH_BDMH2$postmean_noRB )
    GIMH_BDMH_noRB_chain_Theta_MAE[n] <- mean(e12.noRB[upper.tri(e12.noRB, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.GIMH_BDMH1$p_links - fit.GIMH_BDMH2$p_links)
    GIMH_BDMH_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    e12.p.noRB <- abs(fit.GIMH_BDMH1$p_links_noRB - fit.GIMH_BDMH2$p_links_noRB)
    GIMH_BDMH_noRB_chain_p_MAE[n] <- mean(e12.p.noRB[upper.tri(e12.p.noRB, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.noRB)
    rm(e12.p)
    rm(e12.p.noRB)
    
  }
  
  return(list(
    "Theta_MAE" = GIMH_BDMH_Theta_MAE,
    "p_MAE" = GIMH_BDMH_p_MAE,
    "chain_Theta" = GIMH_BDMH_chain_Theta_MAE,
    "chain_p" = GIMH_BDMH_chain_p_MAE,
    "Theta_MAE_noRB" = GIMH_BDMH_noRB_Theta_MAE,
    "p_MAE_noRB" = GIMH_BDMH_noRB_p_MAE,
    "chain_Theta_noRB" = GIMH_BDMH_noRB_chain_Theta_MAE,
    "chain_p_noRB" = GIMH_BDMH_noRB_chain_p_MAE,
    "Time" = GIMH_BDMH_time
  ))
  
}

regression.pl_mixing <- function(data_y, Theta_sim, iterations_vec, Omega_init1, Omega_init2, w_slab, s_slab, lambda){
  
  Sigma_init1 <- pd.solve(Omega_init1)
  Sigma_init2 <- pd.solve(Omega_init2)
  
  regression.pl_Theta_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  regression.pl_p_MAE <- array(NA, dim = c(length(iterations_vec), 2))
  regression.pl_chain_Theta_MAE <- array(NA, dim = c(length(iterations_vec)))
  regression.pl_chain_p_MAE <- array(NA, dim = c(length(iterations_vec)))
  
  regression.pl_time <- array(NA, dim = c(length(iterations_vec), 2))
  
  
  for(n in 1:length(iterations_vec)){
    
    ## regression.pl
    ### init 1
    timer_regression.pl1.start <- Sys.time()
    fit.regression.pl1 <- regression.pl_GGM(Y=data_y, iter=iterations_vec[n], burnin = warmup, w_slab=w_slab, s_slab=s_slab, lambda=lambda, Omega_init=Omega_init1, store_samples = FALSE)
    timer_regression.pl1.end <- Sys.time()
    ### init 2
    timer_regression.pl2.start <- Sys.time()
    fit.regression.pl2 <- regression.pl_GGM(Y=data_y, iter=iterations_vec[n], burnin = warmup, w_slab=w_slab, s_slab=s_slab, lambda=lambda, Omega_init=Omega_init2)
    timer_regression.pl2.end <- Sys.time()
    
    regression.pl_time[n, 1] <- difftime(timer_regression.pl1.end, timer_regression.pl1.start, units = "secs")
    regression.pl_time[n, 2] <- difftime(timer_regression.pl2.end, timer_regression.pl2.start, units = "secs")
    
    ## chain 1
    e1 <- abs(fit.regression.pl1$Omega_hat - Theta_sim)
    regression.pl_Theta_MAE[n, 1] <- mean(e1[upper.tri(e1, diag = TRUE)], na.rm = TRUE)
    e1.p <- abs((Theta_sim!=0) - fit.regression.pl1$p_links)
    regression.pl_p_MAE[n, 1] <- mean(e1.p[upper.tri(e1.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e1)
    rm(e1.p)
    
    ## chain 2
    e2 <- abs(fit.regression.pl2$Omega_hat - Theta_sim)
    regression.pl_Theta_MAE[n, 2] <- mean(e2[upper.tri(e2, diag = TRUE)], na.rm = TRUE)
    e2.p <- abs((Theta_sim!=0) - fit.regression.pl2$p_links)
    regression.pl_p_MAE[n, 2] <- mean(e2.p[upper.tri(e2.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e2)
    rm(e2.p)
    
    ## between chains
    e12 <- abs(fit.regression.pl2$Omega_hat - fit.regression.pl1$Omega_hat)
    regression.pl_chain_Theta_MAE[n] <- mean(e12[upper.tri(e12, diag = TRUE)], na.rm = TRUE)
    e12.p <- abs(fit.regression.pl1$p_links - fit.regression.pl2$p_links)
    regression.pl_chain_p_MAE[n] <- mean(e12.p[upper.tri(e12.p, diag = FALSE)], na.rm = TRUE)
    
    rm(e12)
    rm(e12.p)
    
  }
  
  return(list(
    "Theta_MAE" = regression.pl_Theta_MAE,
    "p_MAE" = regression.pl_p_MAE,
    "chain_Theta" = regression.pl_chain_Theta_MAE,
    "chain_p" = regression.pl_chain_p_MAE,
    "Time" = regression.pl_time
  ))
  
}
