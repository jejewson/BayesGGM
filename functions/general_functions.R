## regression.pl Sampling
regression.pl_GGM <- function(Y, iter, burnin, w_slab, s_slab, lambda, Omega_init, store_samples = FALSE, store_model_samples = FALSE){
  
  p <- ncol(Y)
  
  p_links <- matrix(1, nrow = p, ncol = p)
  Omega_hat <- matrix(0, nrow = p, ncol = p)
  Omega_95_U <- matrix(NA, nrow = p, ncol = p)
  Omega_95_L <- matrix(NA, nrow = p, ncol = p)
  
  if(store_samples){
    Omega_samples <- array(NA, dim = c((iter - burnin), p, p))
  }
  if(store_model_samples){
    Z_samples <- array(NA, dim = c((iter - burnin), p, p))
  }
  
  for (id_var in 1:p){
    
    X = Y[,-id_var]
    y = Y[,id_var]
    #phi = 1./ Theta_sim[[1]][id_var, id_var] - Atchade did some empirical Bayes for the diagonals
    
    fit.modelSelection_lr <- modelSelection(y=y, x=X, center=FALSE, scale=FALSE,
                                   priorVar = igprior(1, lambda/2),
                                   priorCoef=normalidprior(s_slab^2),
                                   priorModel=modelbinomprior(w_slab),
                                   niter=iter, burnin=burnin, 
                                   verbose=FALSE, enumerate = FALSE,
                                   deltaini = (Omega_init[id_var,-id_var] != 0),
                                   initSearch=='none')# modelSelection GGM also only really cares about what is 0
    
    p_links[id_var,id_var] <- 1
    p_links[id_var, -id_var] <- fit.modelSelection_lr$margpp
    
    atch_samples <- rnlp(msfit = fit.modelSelection_lr, niter=(iter - burnin), burnin=0)
    Omega_hat[id_var, id_var] <- mean(1/atch_samples[,p+1])
    Omega_hat[id_var, -id_var] <- - colMeans(atch_samples[,2:p] / matrix(atch_samples[,p+1], nrow = (iter - burnin),ncol = p-1, byrow = FALSE))
    
    
    Omega_95_U[id_var, id_var] <- quantile(1/atch_samples[,p+1], probs = 0.975)
    Omega_95_U[id_var, -id_var] <- apply(-atch_samples[,2:p] / matrix(atch_samples[,p+1], nrow = (iter - burnin),ncol = p-1, byrow = FALSE), MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
    Omega_95_L[id_var, id_var] <- quantile(1/atch_samples[,p+1], probs = 0.025)
    Omega_95_L[id_var, -id_var] <- apply(-atch_samples[,2:p] / matrix(atch_samples[,p+1], nrow = (iter - burnin),ncol = p-1, byrow = FALSE), MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
    
    if(store_model_samples){
      Z_samples[, id_var, id_var] <- 1
      for(n in 1:(iter - burnin)){
        Z_samples[n, id_var, -id_var] <- fit.modelSelection_lr$postSample[n,]
      }
    } else Z_samples <- NULL
    
    if(store_samples){
      Omega_samples[, id_var, id_var] <- 1/atch_samples[,p+1]
      for(n in 1:(iter - burnin)){
        Omega_samples[n, id_var, -id_var] <- -atch_samples[n,2:p]/atch_samples[n,p+1]
      }## difficult to do joint, but if we only care about marginals this is fine
    } else Omega_samples <- NULL
    
  }
  
  Omega_hat <- (Omega_hat + t(Omega_hat))/2
  p_links <- (p_links + t(p_links))/2
  
  return(list("p_links" = p_links, "Omega_hat" = Omega_hat, "Omega_95_L" = Omega_95_L, "Omega_95_U" = Omega_95_U, "Omega_samples" = Omega_samples, "Z_samples" = Z_samples))
}

# Posteriro Sampling for bdgraph or bdgraph.mpl
bdgraph.mpl_sampleOmega <- function(fit.bdgraph.mpl, iter, p, S, n_train, store_samples = FALSE){
  
  if(store_samples){
    Omega_sample <- array(NA, dim = c(iter, p, p))
  } else Omega_sample <- NULL
  Omega_hat <- matrix(0, nrow = p, ncol = p)
  
  for(n in 1:iter){
    #sample iteration from graph_weights
    ind <- sample(1:length(fit.bdgraph.mpl$sample_graphs), 1, prob = fit.bdgraph.mpl$graph_weights/sum(fit.bdgraph.mpl$graph_weights))
    G_sample <- matrix(0, nrow= p, ncol = p)
    G_sample[upper.tri(G_sample, diag = FALSE)] <- as.numeric(strsplit(fit.bdgraph.mpl$sample_graphs[ind], split = "")[[1]])
    G_sample <- G_sample + t(G_sample) + diag(1, p)
    if(store_samples){
      Omega_sample[n,,] <-  BDgraph::rgwish(n = 1, adj = G_sample, b = 3 + n_train, D = diag(1, p) + S)
      Omega_hat <- Omega_hat + Omega_sample[n,,]/iter
    } else{
      Omega_hat <- Omega_hat + BDgraph::rgwish(n = 1, adj = G_sample, b = 3 + n_train, D = diag(1, p) + S )/iter
    }

  }
  
  return(list("Omega_hat" = Omega_hat, "Omega_sample" = Omega_sample))
}

FDR <- function(FP, TP){
  return(ifelse((FP + TP) > 0, FP/(FP + TP), 0))
}

bayes_graph_select <- function(probs, thresh = 0.95, FDR = FALSE){
  p <- ncol(probs)
  
  if(!FDR){
    out <- matrix(as.numeric(probs >= thresh), nrow = p, ncol = p)
  } else if(FDR){
    probs <- probs[upper.tri(probs, diag = FALSE)]
    probs_decr <- sort(probs, decreasing = TRUE)
    av_prob_inclusion <- cumsum(probs_decr)/(1:(p*(p-1)/2))
    cut_off_ind <- max(which(av_prob_inclusion >= thresh))
    out <- diag(0.5, p)
    if(cut_off_ind >= 1){
      out[upper.tri(out, diag = FALSE)][order(probs, decreasing = TRUE)[1:cut_off_ind]] <- 1
    }
    out <- out + t(out)
  }
  
  return(out)
  
}


dmvnorm_prec <- function(x, sigmainv, logdet.sigmainv, mu = rep(0, ncol(sigmainv)), log = FALSE) {
  if (missing(logdet.sigmainv)) logdet.sigmainv= determinant(sigmainv,logarithm=TRUE)$modulus
  d= mahalanobis(x, center=mu, cov=sigmainv, inverted=TRUE)
  ans= -0.5 * (d + ncol(sigmainv) * log(2*pi) - logdet.sigmainv)
  if(!log) ans= exp(ans)
  return(ans)
}

glasso_getBIC = function(y, sfit, method='EBIC') {
  #Extraction of BIC/EBIC, copied from huge.select in package huge
  n= nrow(y); d= ncol(y)
  if (method == 'BIC') {
    gamma = 0
  } else if (method == 'EBIC') {
    gamma = 0.5
  }
  #ans = -n * sfit$loglik + log(n) * sfit$df + 4 * gamma * log(d) * sfit$df
  #Manual implementation (deprecated)    
  logl= npar= double(length(sfit$icov))
  for (i in 1:length(logl)) {
    logl[i]= sum(dmvnorm_prec(y, sigmainv=sfit$icov[[i]], log=TRUE))
    npar[i]= (sum(sfit$icov[[i]] != 0) - d)/2
  }
  if (method == 'BIC') {
    ans= -2*logl + npar * log(n)
  } else if (method == 'EBIC') {
    ans= -2*logl + npar * (log(n) + 4 * gamma * log(d))
  }
  return(ans)
}


OmegaHat_GLASSO = function(y, criteria, maxiter = 5) {
  
    if (criteria == 'glasso-bic') {
      method <- "BIC"
      #gamma = 0
    } else if (criteria == 'glasso-ebic') {
      #gamma = 0.5
      method <- "EBIC"
    }
    
    sfit= huge::huge(y, method="glasso", scr=TRUE, verbose=FALSE, nlambda = 10)
    #sfit2= huge::huge.select(sfit, criterion="ebic", ebic.gamma = gamma, verbose=FALSE)
    ebic= glasso_getBIC(y=y, sfit=sfit, method=method)
    
    #topbic= which.min(sfit2$ebic.score)
    topbic= which.min(ebic)
    found= (topbic != 1) & (topbic != 10)
    
    niter= 1
    while ((!found) & (niter<=maxiter)) {
      ## huge wants a decreasing sequence
      if((niter == 1) & (topbic == 1)){## No need for lambda bigger, 
        #break## but can zoom in
        rholist= seq(sfit$lambda[1], sfit$lambda[2], length=10)
      }
      if (topbic==10) {
        rholist= seq(sfit$lambda[10], sfit$lambda[10]/10, length=10)
      } else {
        rholist= seq(10*sfit$lambda[1], sfit$lambda[1], length=10)
      }
      
      sfit= huge::huge(y, lambda=rholist, method="glasso", scr=TRUE, verbose=FALSE)
      #sfit2= huge::huge.select(sfit, criterion="ebic", ebic.gamma = gamma, verbose=FALSE)
      ebic= glasso_getBIC(y=y, sfit=sfit, method=method)
      
      #topbic= which.min(sfit2$ebic.score)
      topbic= which.min(ebic)
      found= (topbic != 1) & (topbic != 10)
      niter= niter+1
      
    }
    
    ans= Matrix::Matrix(sfit$icov[[topbic]], sparse=TRUE)
    

  
  return(ans)
  
}

generate_G <- function(p,type, size = NULL, graph_prob = NULL){
  
  if (type == "random"){
    
    G <- diag(0.5, p)
    G[lower.tri(G, diag = FALSE)] <- sample(c(0, 1), p*(p-1)/2, replace = TRUE, prob = c(1-graph_prob, graph_prob))
    G <- G + t(G)
    
    
  } else if (type == "scale-free"){
    
    data <- bdgraph.sim(n = 1, p = p, vis = FALSE, graph = "scale-free",  b = 3, D = diag(1, p))
    
    G <- data$G
    
  } else if (type == "tridiagonal"){
    
    if (is.null(size)){
      size = p 
    }
    
    G <- diag(1, p)
    diag(G[,-1][1:size-1, 1:size-1]) <- rep(1,size-1)
    diag(G[-1,][1:size-1, 1:size-1]) <- rep(1,size-1)
    
  } else if (type == "multidiagonal"){
    
    G <- diag(1, p)
    for (i in 1:2){
      diag(G[,-c(1:i)]) <- rep(1,p-i)
      diag(G[-c(1:i),]) <- rep(1,p-i)
    }
    
  } else if (type == "blockdiagonal"){
    
    G <- kronecker(diag(p), matrix(1, nrow = size, ncol = size))[1:p, 1:p]
    
  } else {
    
    print("not implemented")
    break
  }
  
  return(G)
}

generate_Theta_y <- function(n, G, min_rho = 0.1, max_rho = 0.5){
  # generate partial correlations
  # generate diagonals - do a Gamma (reflect most greater than 1)
  
  p <- nrow(G)
  
  pd <- FALSE
  #while(!pd){
  Rho <- diag(0.5, p)
  Rho[lower.tri(Rho, diag = FALSE)] <- sample(c(-1, 1), p*(p-1)/2, replace = TRUE)*
    sample(seq(min_rho, max_rho, by = 0.1), p*(p-1)/2, replace = TRUE)*
    G[lower.tri(G, diag = FALSE)]## I had a mistaken / 2 here!
  Rho <- Rho + t(Rho)
  
  pd <- (min(eigen(Rho)$values) > 0)
  #}
  
  if(pd){
    sqrt_omega_jj <- sqrt(rgamma(p, shape = 3, rate = 1))
    Theta <- outer(sqrt_omega_jj, sqrt_omega_jj)*Rho
    
    y = rmvnorm(n = n, sigma = pd.solve(Theta))
    
    ggm <- list(G = G, Theta = Theta, y = y)
  } else{
    ggm <- NA
  }
  
  return(ggm)
}



traceplot.modelSelection <- function(modelSelection.object, p){
  N_MCMC <- modelSelection.object@.Data[[1]]@Dim[1]
  model_size <- rep(NA, N_MCMC)
  for(n in 1:N_MCMC){
    model_size[n] <- sum(modelSelection.object$postSample[n,] != 0) - p
  }
  return(model_size)
}

plotcoda.modelSelection <- function(modelSelection.object, p, seed = 1, N_plot=NULL){
  N_MCMC <- modelSelection.object@.Data[[1]]@Dim[1]
  if(is.null(N_plot)){
    N_plot <- min(100, (p*(p+1)/2))
  }
  inclusion_probs <- matrix(NA, nrow = N_MCMC, ncol = N_plot)
  set.seed(seed)
  plot_index <- sample(1:(p*(p+1)/2), N_plot, replace = FALSE)
  for(n in 1:N_MCMC){
    inclusion_probs[n,] <- colMeans(as.matrix(modelSelection.object$postSample[1:n,plot_index] != 0))
  }
  return(inclusion_probs)
}

traceplot_llp.modelSelection <- function(y, modelSelection.object, p, thinning = 1){
  N_MCMC <- modelSelection.object@.Data[[1]]@Dim[1]
  N_MCMC_thinned <- ceiling(N_MCMC/thinning)
  log_target <- rep(NA, N_MCMC_thinned)
  Omega <- matrix(0, nrow = p, ncol = p)
  for(n in 1:N_MCMC_thinned){
    # fill Omega from correct index
    Omega[upper.tri(Omega, diag = TRUE)] <- modelSelection.object$postSample[(n-1)*thinning + thinning,]
    Omega[lower.tri(Omega, diag = FALSE)] <- 0
    Omega <- Omega + t(Omega) - diag(diag(Omega))
    llp <- sum(LaplacesDemon::dmvnp(y, mu = rep(0, p), Omega, log=TRUE))
    prior_diag <- sum(dgamma(diag(Omega), shape = modelSelection.object@.Data[[7]]$priorDiag@priorPars[1],
                             rate = modelSelection.object@.Data[[7]]$priorDiag@priorPars[2]/2, log = TRUE))
    prior_zeros <- sum(Omega[upper.tri(Omega, diag = FALSE)] == 0)*log(1-modelSelection.object@.Data[[7]]$priorModel@priorPars$p) + 
      sum(Omega[upper.tri(Omega, diag = FALSE)] != 0)*log(modelSelection.object@.Data[[7]]$priorModel@priorPars$p)
    prior_nonzeros <- sum(dnorm(Omega[upper.tri(Omega, diag = FALSE)][which(Omega[upper.tri(Omega, diag = FALSE)] != 0)], 0, sqrt(modelSelection.object@.Data[[7]]$priorCoef@priorPars[1]), log = TRUE))
    log_target[n] <- llp + prior_diag + prior_zeros + prior_nonzeros
  }
  return(log_target)
}









regression_MSE_eval <- function(y_test, Omega_hat){
  n_test <- nrow(y_test)
  p <- ncol(y_test)
  
  out <- rep(0, n_test)
  for(i in 1:n_test){
    for(j in 1:p){
      out[i] <- out[i] + 1/p*(y_test[i,j] - y_test[i,-j]%*%(-Omega_hat[j,-j])/Omega_hat[j,j])^2
    }
  }
  return(mean(out))
}

regression_R2_eval <- function(y_test, Omega_hat){
  n_test <- nrow(y_test)
  p <- ncol(y_test)
  
  y_hat <- matrix(NA, nrow = n_test, ncol = p)
  
  out <- rep(0, n_test)
  for(i in 1:n_test){
    for(j in 1:p){
      y_hat[i,j] <- y_test[i,-j]%*%(-Omega_hat[j,-j])/Omega_hat[j,j]
    }
  }
  return(cor(as.vector(y_test), as.vector(y_hat))^2)
}

traceplot.modelSelection.regression.pl <- function(regression.pl.object, p){
  N_MCMC <- dim(regression.pl.object$Z_samples)[1]
  model_size <- rep(NA, N_MCMC)
  for(n in 1:N_MCMC){
    #model_size[n] <- (sum(regression.pl.object$Omega_samples[n,,] != 0) - p)/2
    model_size[n] <- (sum(regression.pl.object$Z_samples[n,,]) - p)/2
  }
  return(model_size)
}  

plotcoda.modelSelection.regression.pl <- function(regression.pl.object, p, seed = 1, N_plot=NULL){
  N_MCMC <- dim(regression.pl.object$Z_samples)[1]
  if(is.null(N_plot)){
    N_plot <- min(100, (p*(p+1)/2))
  }
  inclusion_probs <- matrix(NA, nrow = N_MCMC, ncol = N_plot)
  set.seed(seed)
  plot_index1 <- sample(1:p, N_plot, replace = TRUE)
  plot_index2 <- sample(1:p, N_plot, replace = TRUE)
  for(i in 1:N_plot){
    #inclusion_probs[,i] <- cumsum(regression.pl.object$Omega_samples[,plot_index1[i], plot_index2[i]] != 0)/(2*(1:N_MCMC)) + cumsum(regression.pl.object$Omega_samples[,plot_index2[i], plot_index1[i]] != 0)/(2*(1:N_MCMC)) 
    inclusion_probs[,i] <- cumsum(regression.pl.object$Z_samples[,plot_index1[i], plot_index2[i]])/(2*(1:N_MCMC)) + cumsum(regression.pl.object$Z_samples[,plot_index2[i], plot_index1[i]])/(2*(1:N_MCMC)) 
  }
  return(inclusion_probs)
}

BDgraph_traceplot <- function (bdgraph.obj, thin = 1, acf = FALSE, pacf = FALSE, main = NULL, 
                               ...) 
{
  if ((inherits(bdgraph.obj, "bdgraph")) | (inherits(bdgraph.obj, 
                                                     "ssgraph"))) {
    if (is.null(bdgraph.obj$all_graphs)) 
      stop("'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option 'save = TRUE'")
    sample_graphs = bdgraph.obj$sample_graphs
    all_graphs = bdgraph.obj$all_graphs
    graph_weights = bdgraph.obj$graph_weights
    sizesample_graphs = sapply(sample_graphs, function(x) length(which(unlist(strsplit(as.character(x), 
                                                                                       "")) == 1)))
    sizeall_graphs = sizesample_graphs[all_graphs]
    which_G_max = which(max(graph_weights) == graph_weights)
    size_selected_g = sizeall_graphs[which_G_max]
    sample_mcmc = sizeall_graphs
    if (is.null(main)) 
      main = "Trace of graph size"
    ylab = "Graph size"
  }
  else {
    if (!is.vector(bdgraph.obj)) 
      stop("'bdgraph.obj' must be an object of functions 'bdgraph()', 'bdgraph.mpl()', or 'ssgraph()' or a vector")
    sample_mcmc = bdgraph.obj
    if (is.null(main)) 
      main = "Trace of MCMC sample"
    ylab = ""
  }
  #if (acf == FALSE & pacf == FALSE) 
  #  op = graphics::par(mfrow = c(1, 1), pty = "s")
  #if (acf == TRUE & pacf == TRUE) 
  #  op = graphics::par(mfrow = c(2, 2), pty = "s")
  #if (acf == TRUE & pacf == FALSE) 
  #  op = graphics::par(mfrow = c(1, 2), pty = "s")
  #if (acf == FALSE & pacf == TRUE) 
  #  op = graphics::par(mfrow = c(1, 2), pty = "s")
  x_vec = (1:length(sample_mcmc))[seq(1, length(sample_mcmc), by = thin)]
  graphics::plot(x = x_vec, y = sample_mcmc[seq(1, length(sample_mcmc), by = thin)], type = "l", main = main, 
                 col = "black", ylab = ylab, xlab = "Iteration", ...)
  #if (!is.vector(bdgraph.obj)) 
  #  graphics::lines(x = x_vec, y = rep(size_selected_g, length(sample_mcmc[seq(1, length(sample_mcmc), by = thin)])), 
  #                  col = "blue")
  if (acf == TRUE) 
    acf(sample_mcmc, main = "ACF for graph size")
  if (pacf == TRUE) 
    pacf(sample_mcmc, main = "PACF for graph size")
  #graphics::par(op)
}

BDgraph_plotcoda <- function (bdgraph.obj, thin = NULL, thin_plot = 1, control = TRUE, main = NULL, 
                              verbose = TRUE, ...) 
{
  if ((inherits(bdgraph.obj, "bdgraph")) | (inherits(bdgraph.obj, 
                                                     "ssgraph"))) {
    if (is.null(bdgraph.obj$all_graphs)) 
      stop("'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option 'save = TRUE'")
    if (is.null(bdgraph.obj$all_graphs)) 
      stop("'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option 'save = TRUE'")
  }
  else {
    stop("'bdgraph.obj' must be an object of functions 'bdgraph()', 'bdgraph.mpl()', or 'ssgraph()'")
  }
  if (is.null(thin)) 
    thin = ceiling(length(bdgraph.obj$all_graphs)/1000)
  if (!is.numeric(thin)) 
    stop("'thin' must be a number")
  if (is.matrix(thin)) 
    stop("'thin' must be a number")
  sample_graphs = bdgraph.obj$sample_graphs
  p = nrow(bdgraph.obj$last_graph)
  qp = p * (p - 1)/2
  all_weights = bdgraph.obj$all_weights
  all_graphs = bdgraph.obj$all_graphs
  allG_new = all_graphs[c(thin * (1:floor(length(all_graphs)/thin)))]
  all_weights_new = all_weights[c(thin * (1:floor(length(all_weights)/thin)))]
  length_allG_new = length(allG_new)
  result = matrix(0, qp, length_allG_new)
  vec_result = 0 * result[, 1]
  for (g in 1:length_allG_new) {
    if (verbose == TRUE) {
      mes = paste(c("Calculation ... in progress : ", floor(100 * 
                                                              g/length_allG_new), "%"), collapse = "")
      cat(mes, "\r")
      utils::flush.console()
    }
    which_edge = which(unlist(strsplit(as.character(sample_graphs[allG_new[g]]), 
                                       "")) == 1)
    vec_result[which_edge] = vec_result[which_edge] + all_weights_new[g]
    result[, g] = vec_result/sum(all_weights_new[c(1:g)])
  }
  if (control) 
    if (p > 15) {
      #randomLinks = sample(x = 1:qp, size = (qp - 100), 
      #                     replace = FALSE)
      #result[randomLinks, ] = 0
      randomLinks = sample(x = 1:qp, size = 100, 
                           replace = FALSE)
      result <- result[randomLinks, ]
    }
  if (verbose == TRUE) {
    mes = paste(c("Calculation ... done.                        "), 
                collapse = "")
    cat(mes, "\r")
    cat("\n")
    utils::flush.console()
  }
  graphics::matplot(x = (thin * (1:length_allG_new)[seq(1, length_allG_new, by = thin_plot)]), y = t(result[,seq(1, length_allG_new, by = thin_plot)]), 
                    type = "l", lty = 1, col = "black", xlab = "Iteration", 
                    ylab = "Posterior edge probability")#, cex.lab = 1.3, cex.axis = 1.2)
  if (is.null(main)) 
    main = "Trace of the Posterior Probabilities of the Links."
  graphics::title(main = main)#, cex.main = 1.2)
}
