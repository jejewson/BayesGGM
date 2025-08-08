GLASSO_inference <- function(data_y, Theta_sim, criteria){
  
  timer_GLASSO.start <- Sys.time()
  Omega_hat <- as.matrix(OmegaHat_GLASSO(data_y, criteria=paste0('glasso-', criteria), maxiter=5))
  timer_GLASSO.end <- Sys.time()
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  
  return(list(
    "Time" = difftime(timer_GLASSO.end, timer_GLASSO.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "TP" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(Omega_hat[upper.tri(Omega_hat, diag = FALSE)] != 0)),
    "TN" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(Omega_hat[upper.tri(Omega_hat, diag = FALSE)] == 0)),
    "FP" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(Omega_hat[upper.tri(Omega_hat, diag = FALSE)] != 0)),
    "FN" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(Omega_hat[upper.tri(Omega_hat, diag = FALSE)] == 0))
  ))
  
}

ssgraph_inference <- function(data_y, Theta_sim, iterations, warmup, w_slab, s_slab, lambda, s_spike = 2e-02){
  p <- ncol(data_y)
  
  ## initalisation
  Omegaini1 <- mombf:::initialEstimateGGM(data_y, Omegaini='glasso-ebic')
  Omega_init1 <- (as.matrix(Omegaini1) + t(as.matrix(Omegaini1)))/2
  Sigma_init1 <- pd.solve(Omega_init1)
  
  ## ssgraph
  ### init 1
  timer_ssgraph1.start <- Sys.time()
  fit.ssgraph1 <- ssgraph(data_y, method="ggm", iter=iterations+warmup, burnin = warmup, g.prior= w_slab
                          , var1 = s_spike^2, var2 = s_slab^2
                          , lambda = lambda, save=TRUE, verbose=TRUE, sig.start = Sigma_init1)
  timer_ssgraph1.end <- Sys.time()
  
  Omega_hat <- fit.ssgraph1$K_hat
  p_hat <- fit.ssgraph1$p_links
  G_hat <- bayes_graph_select(p_hat, thresh = 0.95, FDR = TRUE)
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  e1.p <- abs((Theta_sim!=0) - fit.ssgraph1$p_links)
  
  actual <- as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)
  predicted_probs <- p_hat[upper.tri(p_hat, diag = FALSE)]
  ROC <- pROC::roc(actual, predicted_probs)
  
  #coverage_temp <- matrix(NA, nrow = p, ncol = p)
  #for(j in 1:p){
  #  for(k in j:p){
  #    coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(fit.ssgraph1$Omega_samples[,j,k], prob = 0.025)) & (Theta_sim[j,k] <= quantile(fit.ssgraph1$Omega_samples[,j,k], prob = 0.975)))
  #  }
  #}
  
  return(list(
    "Time" = difftime(timer_ssgraph1.end, timer_ssgraph1.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "p_MAE" = mean(e1.p[upper.tri(e1.p, diag = TRUE)]),
    "p_zeros_MAE" = mean(e1.p[upper.tri(e1.p, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "p_nonzeros_MAE" = mean(e1.p[which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    #"TP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"TN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    #"FP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"FN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    "TP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "TN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "FP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "FN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "AUC" = ROC$auc#,
    #"coverage_nonzeros" = mean(ssgraph_coverage_temp[(Theta_sim - diag(diag(Theta_sim))) != 0], na.rm = TRUE)
  ))
  
}


bdgraph_inference <- function(data_y, Theta_sim, iterations, warmup, w_slab){
  p <- ncol(data_y)
  
  ## initalisation
  Omegaini1 <- mombf:::initialEstimateGGM(data_y, Omegaini='glasso-ebic')
  Omega_init1 <- (as.matrix(Omegaini1) + t(as.matrix(Omegaini1)))/2
  Sigma_init1 <- pd.solve(Omega_init1)
  
  S <- t(data_y)%*%data_y
  n_train <- nrow(data_y)
  
  ## bdgraph
  ### init 1
  timer_bdgraph1.start <- Sys.time()
  fit.bdgraph1 <- bdgraph(data_y, method="ggm", iter=iterations*p+warmup*p, burnin = warmup*p, g.prior= w_slab, save=TRUE, verbose=TRUE, g.start = Omega_init1!=0, df.prior = 3)
  bdgraph_samples_temp <- bdgraph.mpl_sampleOmega(fit.bdgraph1, iter=iterations, p, S, n_train)
  #fit.bdgraph1$K_hat <- bdgraph_samples_temp$Omega_hat
  fit.bdgraph1$Omega_samples <- bdgraph_samples_temp$Omega_sample
  timer_bdgraph1.end <- Sys.time()
  
  fit.bdgraph1$p_links <- plinks(fit.bdgraph1, round = 3)
  
  Omega_hat <- fit.bdgraph1$K_hat
  p_hat <- fit.bdgraph1$p_links
  G_hat <- bayes_graph_select(p_hat, thresh = 0.95, FDR = TRUE)
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  e1.p <- abs((Theta_sim!=0) - fit.bdgraph1$p_links)
  
  actual <- as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)
  predicted_probs <- p_hat[upper.tri(p_hat, diag = FALSE)]
  ROC <- pROC::roc(actual, predicted_probs)
  
  coverage_temp <- matrix(NA, nrow = p, ncol = p)
  for(j in 1:p){
    for(k in j:p){
      coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(fit.bdgraph1$Omega_samples[,j,k], prob = 0.025)) & (Theta_sim[j,k] <= quantile(fit.bdgraph1$Omega_samples[,j,k], prob = 0.975)))
    }
  }
  
  return(list(
    "Time" = difftime(timer_bdgraph1.end, timer_bdgraph1.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "p_MAE" = mean(e1.p[upper.tri(e1.p, diag = TRUE)]),
    "p_zeros_MAE" = mean(e1.p[upper.tri(e1.p, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "p_nonzeros_MAE" = mean(e1.p[which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    #"TP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"TN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    #"FP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"FN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    "TP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "TN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "FP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "FN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "AUC" = ROC$auc,
    "coverage_nonzeros" = mean(coverage_temp[(Theta_sim - diag(diag(Theta_sim))) != 0], na.rm = TRUE)
  ))
  
}


bdgraph.mpl_inference <- function(data_y, Theta_sim, iterations, warmup, w_slab){
  p <- ncol(data_y)
  
  ## initalisation
  Omegaini1 <- mombf:::initialEstimateGGM(data_y, Omegaini='glasso-ebic')
  Omega_init1 <- (as.matrix(Omegaini1) + t(as.matrix(Omegaini1)))/2
  Sigma_init1 <- pd.solve(Omega_init1)
  
  S <- t(data_y)%*%data_y
  n_train <- nrow(data_y)
  
  ## bdgraph.mpl
  ### init 1
  timer_bdgraph.mpl1.start <- Sys.time()
  fit.bdgraph.mpl1 <- bdgraph.mpl(data_y, method="ggm", iter=iterations*p+warmup*p, burnin = warmup*p, g.prior= w_slab, save=TRUE, verbose=TRUE, g.start = Omega_init1!=0)
  bdgraph.mpl_samples_temp <- bdgraph.mpl_sampleOmega(fit.bdgraph.mpl1, iter=iterations, p, S, n_train, store_samples = TRUE)
  fit.bdgraph.mpl1$K_hat <- bdgraph.mpl_samples_temp$Omega_hat
  fit.bdgraph.mpl1$Omega_samples <- bdgraph.mpl_samples_temp$Omega_sample
  timer_bdgraph.mpl1.end <- Sys.time()
  
  fit.bdgraph.mpl1$p_links <- plinks(fit.bdgraph.mpl1, round = 3)
  
  Omega_hat <- fit.bdgraph.mpl1$K_hat
  p_hat <- fit.bdgraph.mpl1$p_links
  G_hat <- bayes_graph_select(p_hat, thresh = 0.95, FDR = TRUE)
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  e1.p <- abs((Theta_sim!=0) - fit.bdgraph.mpl1$p_links)
  
  actual <- as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)
  predicted_probs <- p_hat[upper.tri(p_hat, diag = FALSE)]
  ROC <- pROC::roc(actual, predicted_probs)
  
  coverage_temp <- matrix(NA, nrow = p, ncol = p)
  for(j in 1:p){
    for(k in j:p){
      coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(fit.bdgraph.mpl1$Omega_samples[,j,k], prob = 0.025)) & (Theta_sim[j,k] <= quantile(fit.bdgraph.mpl1$Omega_samples[,j,k], prob = 0.975)))
    }
  }
  
  return(list(
    "Time" = difftime(timer_bdgraph.mpl1.end, timer_bdgraph.mpl1.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "p_MAE" = mean(e1.p[upper.tri(e1.p, diag = TRUE)]),
    "p_zeros_MAE" = mean(e1.p[upper.tri(e1.p, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "p_nonzeros_MAE" = mean(e1.p[which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    #"TP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"TN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    #"FP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"FN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    "TP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "TN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "FP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "FN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "AUC" = ROC$auc,
    "coverage_nonzeros" = mean(coverage_temp[(Theta_sim - diag(diag(Theta_sim))) != 0], na.rm = TRUE)
  ))
  
}

Gibbs_inference <- function(data_y, Theta_sim, iterations, warmup, w_slab, s_slab, lambda){
  p <- ncol(data_y)
  
  ## initalisation
  Omegaini1 <- mombf:::initialEstimateGGM(data_y, Omegaini='glasso-ebic')
  Omega_init1 <- (as.matrix(Omegaini1) + t(as.matrix(Omegaini1)))/2
  Sigma_init1 <- pd.solve(Omega_init1)
  
  updates_per_iter <- p
  updates_per_column <- p
  
  ## Gibbs
  ### init 1
  timer_Gibbs1.start <- Sys.time()
  fit.Gibbs1 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init1, niter=iterations+warmup, burnin=warmup, scale=FALSE, global_proposal='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda))
  timer_Gibbs1.end <- Sys.time()
  
  b <- coef(fit.Gibbs1)
  fit.Gibbs1$p_links <- matrix(NA, nrow = p, ncol = p)
  fit.Gibbs1$p_links[upper.tri(fit.Gibbs1$p_links, diag = TRUE)] <- b[,6]
  fit.Gibbs1$Omega_samples <- array(NA, dim = c(iterations, p, p))
  for(j in 1:iterations){
    fit.Gibbs1$Omega_samples[j,,][upper.tri(fit.Gibbs1$Omega_samples[j,,], diag = TRUE)] <- fit.Gibbs1$postSample[j,]
  }
  
  Omega_hat <- icov(fit.Gibbs1)
  p_hat <- fit.Gibbs1$p_links
  G_hat <- bayes_graph_select(p_hat, thresh = 0.95, FDR = TRUE)
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  e1.p <- abs((Theta_sim!=0) - fit.Gibbs1$p_links)
  
  actual <- as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)
  predicted_probs <- p_hat[upper.tri(p_hat, diag = FALSE)]
  ROC <- pROC::roc(actual, predicted_probs)
  
  coverage_temp <- matrix(NA, nrow = p, ncol = p)
  for(j in 1:p){
    for(k in j:p){
      coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(fit.Gibbs1$Omega_samples[,j,k], prob = 0.025)) & (Theta_sim[j,k] <= quantile(fit.Gibbs1$Omega_samples[,j,k], prob = 0.975)))
    }
  }
  
  return(list(
    "Time" = difftime(timer_Gibbs1.end, timer_Gibbs1.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "p_MAE" = mean(e1.p[upper.tri(e1.p, diag = TRUE)]),
    "p_zeros_MAE" = mean(e1.p[upper.tri(e1.p, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "p_nonzeros_MAE" = mean(e1.p[which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    #"TP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"TN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    #"FP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"FN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    "TP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "TN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "FP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "FN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "AUC" = ROC$auc,
    "coverage_nonzeros" = mean(coverage_temp[(Theta_sim - diag(diag(Theta_sim))) != 0], na.rm = TRUE)
  ))
  
}

GIMH_Gibbs_inference <- function(data_y, Theta_sim, iterations, warmup, w_slab, s_slab, lambda, tempering){
  p <- ncol(data_y)
  
  ## initalisation
  Omegaini1 <- mombf:::initialEstimateGGM(data_y, Omegaini='glasso-ebic')
  Omega_init1 <- (as.matrix(Omegaini1) + t(as.matrix(Omegaini1)))/2
  Sigma_init1 <- pd.solve(Omega_init1)
  
  updates_per_iter <- p#
  updates_per_column <- ceiling(sqrt(p))
  prob_global <- 1 # so eevry column udpate is a full row form the parallel
  truncratio <- -1 # presumbaly this means no truncation?
  
  ## GIMH_Gibbs
  ### init 1
  timer_GIMH_Gibbs1.start <- Sys.time()
  fit.GIMH_Gibbs1 <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Omega_init1, niter=iterations+warmup, burnin=warmup, scale=FALSE, global_proposal='regression', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(s_slab^2), priorModel=modelbinomprior(w_slab), priorDiag=exponentialprior(lambda=lambda), prob_global=prob_global, tempering=tempering, truncratio=truncratio)
  timer_GIMH_Gibbs1.end <- Sys.time()
  
  b <- coef(fit.GIMH_Gibbs1)
  fit.GIMH_Gibbs1$p_links <- matrix(NA, nrow = p, ncol = p)
  fit.GIMH_Gibbs1$p_links[upper.tri(fit.GIMH_Gibbs1$p_links, diag = TRUE)] <- b[,6]
  fit.GIMH_Gibbs1$Omega_samples <- array(NA, dim = c(iterations, p, p))
  for(j in 1:iterations){
    fit.GIMH_Gibbs1$Omega_samples[j,,][upper.tri(fit.GIMH_Gibbs1$Omega_samples[j,,], diag = TRUE)] <- fit.GIMH_Gibbs1$postSample[j,]
  }
  
  Omega_hat <- icov(fit.GIMH_Gibbs1)
  p_hat <- fit.GIMH_Gibbs1$p_links
  G_hat <- bayes_graph_select(p_hat, thresh = 0.95, FDR = TRUE)
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  e1.p <- abs((Theta_sim!=0) - fit.GIMH_Gibbs1$p_links)
  
  actual <- as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)
  predicted_probs <- p_hat[upper.tri(p_hat, diag = FALSE)]
  ROC <- pROC::roc(actual, predicted_probs)
  
  coverage_temp <- matrix(NA, nrow = p, ncol = p)
  for(j in 1:p){
    for(k in j:p){
      coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(fit.GIMH_Gibbs1$Omega_samples[,j,k], prob = 0.025)) & (Theta_sim[j,k] <= quantile(fit.GIMH_Gibbs1$Omega_samples[,j,k], prob = 0.975)))
    }
  }
  
  return(list(
    "Time" = difftime(timer_GIMH_Gibbs1.end, timer_GIMH_Gibbs1.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "p_MAE" = mean(e1.p[upper.tri(e1.p, diag = TRUE)]),
    "p_zeros_MAE" = mean(e1.p[upper.tri(e1.p, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "p_nonzeros_MAE" = mean(e1.p[which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    #"TP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"TN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    #"FP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"FN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    "TP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "TN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "FP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "FN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "AUC" = ROC$auc,
    "coverage_nonzeros" = mean(coverage_temp[(Theta_sim - diag(diag(Theta_sim))) != 0], na.rm = TRUE)
  ))
  
}

regression.pl_inference <- function(data_y, Theta_sim, iterations, warmup, w_slab, s_slab, lambda){
  p <- ncol(data_y)
  
  ## initalisation
  Omegaini1 <- mombf:::initialEstimateGGM(data_y, Omegaini='glasso-ebic')
  Omega_init1 <- (as.matrix(Omegaini1) + t(as.matrix(Omegaini1)))/2
  Sigma_init1 <- pd.solve(Omega_init1)
  
  ## regression.pl
  ### init 1
  timer_regression.pl1.start <- Sys.time()
  fit.regression.pl1 <- regression.pl_GGM(Y=data_y, iter=iterations+warmup, burnin = warmup, w_slab=w_slab, s_slab=s_slab, lambda=lambda, Omega_init=Omega_init1, store_samples = TRUE)
  timer_regression.pl1.end <- Sys.time()
  
  Omega_hat <- fit.regression.pl1$Omega_hat
  p_hat <- fit.regression.pl1$p_links
  G_hat <- bayes_graph_select(p_hat, thresh = 0.95, FDR = TRUE)
  
  e1 <- abs(as.matrix(Omega_hat) - Theta_sim)
  e1.p <- abs((Theta_sim!=0) - fit.regression.pl1$p_links)
  
  actual <- as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)
  predicted_probs <- p_hat[upper.tri(p_hat, diag = FALSE)]
  ROC <- pROC::roc(actual, predicted_probs)
  
  coverage_temp <- matrix(NA, nrow = p, ncol = p)
  for(j in 1:p){
    for(k in j:p){
      if(j == k){
        coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(fit.regression.pl1$Omega_samples[,j,j], prob = 0.025)) & (Theta_sim[j,k] <= quantile(fit.regression.pl1$Omega_samples[,j,j], prob = 0.975)))
      } else{
        coverage_temp[j,k] <- as.numeric((Theta_sim[j,k] >= quantile(c(fit.regression.pl1$Omega_samples[,j,k], fit.regression.pl1$Omega_samples[,k,j]), prob = 0.025)) & (Theta_sim[j,k] <= quantile(c(fit.regression.pl1$Omega_samples[,j,k], fit.regression.pl1$Omega_samples[,k,j]), prob = 0.975)))
      }
    }
  }
  
  return(list(
    "Time" = difftime(timer_regression.pl1.end, timer_regression.pl1.start, units = "secs"),
    "Theta_MAE" = mean(e1[upper.tri(e1, diag = TRUE)]),
    "Theta_diag_MAE" = mean(diag(e1)),
    "Theta_zeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "Theta_nonzeros_MAE" = mean(e1[upper.tri(e1, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    "p_MAE" = mean(e1.p[upper.tri(e1.p, diag = TRUE)]),
    "p_zeros_MAE" = mean(e1.p[upper.tri(e1.p, diag = FALSE)][which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0)]),
    "p_nonzeros_MAE" = mean(e1.p[which(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0)]),
    #"TP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"TN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    #"FP0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] >= 0.5)),
    #"FN0.5" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(p_hat[upper.tri(p_hat, diag = FALSE)] < 0.5)),
    "TP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "TN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "FP0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] == 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 1)),
    "FN0.95" = sum(as.numeric(Theta_sim[upper.tri(Theta_sim, diag = FALSE)] != 0) & as.numeric(G_hat[upper.tri(G_hat, diag = FALSE)] == 0)),
    "AUC" = ROC$auc,
    "coverage_nonzeros" = mean(coverage_temp[(Theta_sim - diag(diag(Theta_sim))) != 0], na.rm = TRUE)
  ))
  
}
