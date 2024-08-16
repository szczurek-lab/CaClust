#' Run or continue one chain of CaClust model sampling
#'
#' @param A Named matrix of alternate variant read counts in cells, of dimensions variants x cells
#' @param D Named matrix of reference variant read counts in cells, of dimensions variants x cells
#' @param Omega Matrix of estimated clonal profiles, of dimensions variants x clones
#' @param BCR Array of one-hot encoded nucleotides in cells BCR sequences, of dimensions cells x positions x nucleotides
#' @param Psi Prior for clone assignment of cells; defaults to uniform
#' @param alpha_0 Initial value for concentration parameter; NULL will set initial value to mean of the prior
#' @param alpha_prior Prior for the concentration parameter
#' @param clusters DEPRECIATED List of initial hyperclusters; NULL will result in initial hyperclusters of cells with identical BCR sequences.
#' @param t DEPRECIATED Initial assignment of cells to hyperclusters; must agree with clusters if provided
#' @param g_prior_strength Strength hyperparameter with which the overall BCR sequences distribution is included in the prior for profiles of BCR hypercluster, must be non-negative
#' @param g_prior_base Strength hyperparameter with which the symmetrical Dirichlet distribution is included in the prior for profiles of BCR hypercluster, must be non-negative
#' @param relax_rate_prior Prior for the relax rate of the input clonal profiles
#' @param prior0 Prior for the theta0 parameter
#' @param prior1 Prior for the theta1 parameter
#' @param keep_base_clone Does the Omega matrix contain first column with an unmutated reference clone, which should be kept through the inference; defaults to TRUE
#' @param min_iter Minimum number of sampling iterations performed by the model
#' @param max_iter Maximum number of sampling iterations performed by the model
#' @param buin_frac Fraction of sampling iterations to be considered as burn-in in calculating final estimates
#' @param write_skip Save and use only 1 out of every write_skip iterations, useful for memory considerations and autocorrelation
#' @param checkpoint If a run should be continued from a checkpoint, a named list with the last sampled values of: clusters, t, I, B, alpha0, theta0, theta1, C, relax_rate. Otherwise default NULL will result in a new model run.
#'
#' @return
#' @export
#'
#' @examples
#'
caclust_sampling <- function(A, D, Omega, BCR, Psi=NULL, 
                             alpha_0=10, alpha_prior=c(10, 1),
                             clusters=NULL, t=NULL,
                             g_prior_strength=0,
                             g_prior_base=0.01,
                             relax_rate_prior=c(1, 19),
                             prior0=c(0.002, 0.998), prior1=c(0.45, 0.55),
                             keep_base_clone=TRUE,
                             min_iter=1, max_iter=500, 
                             buin_frac=0.75, write_skip=1,
                             checkpoint=NULL) {
  library(magrittr)
  
  if(is.null(Psi)){
    Psi <- rep(1/ncol(Omega), ncol(Omega))
  }
  if(dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
     dim(A)[1] != dim(Omega)[1] || dim(Omega)[2] != length(Psi)){
    stop(paste0("A and D must have the same size;\n ",
                "A and Omega must have the same number of variants;\n",
                "Omega and Psi must have the same number of clones"))
  }
  
  
  ## preprocessing
  N <- dim(A)[1]             # number of variants
  M <- dim(A)[2]             # number of cells
  K <- dim(Omega)[2]         # number of clones
  L <- dim(BCR)[2]           # length of BCR seq
  #Omega is the matrix of inferred mutation profiles
  #Psi is prior of cluster->clone assignment
  
  # make flat BCR
  flat_BCR <- matrix(0, nrow=M, ncol=4*L)
  for(cell in 1:M){
    flat_BCR[cell,] <- BCR[cell,,]
  }
  
  # >>> initialise or get clustering >>>
  g <- apply(BCR, c(2,3), sum) / apply(BCR, 2, sum) * g_prior_strength + g_prior_base
  
  if( is.null(checkpoint) ){
    identical_ordering <- data.frame(flat_BCR) %>% as.list() %>% do.call( grouping, . )
    
    ends <- attr(identical_ordering, 'ends')
    starts <- c(0, ends) + 1
    
    t <- rep.int(0, M)
    clusters <- list()
    for( i in seq_along(ends) ){
      clusters[[i]] <- identical_ordering[ starts[i]:ends[i] ]
      t[ clusters[[i]] ] <- i
    }
    
    #t is cell to cluster assign
    #clusters is the list of clusters containing their members
    
    I <- extraDistr::rcat(M, Psi)
    
    flat_B <- matrix(0, nrow=M, ncol=4*L)
    
    # a cluster by cell matrix of assignment
    clust_cell <- matrix(0, nrow=M, ncol=M)
    
    clust_cell[ cbind(t, 1:M) ] <- 1
    
    # faster B sampling
    non_empty_clusters <- which( sapply(clusters, length)>0 )
    k <- length(non_empty_clusters)
    
    up_priors <- clust_cell[non_empty_clusters,] %*% flat_BCR %>% matrix( ncol=4 )
    up_priors <- up_priors + matrix( t(g), nrow=k*L, ncol=4, byrow=TRUE)
    
    flat_B[non_empty_clusters,] <- rdirichlet(k*L, up_priors) %>% matrix( nrow=k )
    
    if( is.null(alpha_0) ){
      alpha_0 <- alpha_prior[1] / alpha_prior[2]
    }
    
  }else{
    n_checkpoint_clusts <- length(checkpoint$clusters)
    clusters <- checkpoint$clusters
    clusters[[M+1]] <- NULL #expand cluster list to make room for up to M clusters
    t <- checkpoint$t
    
    I <- 1:M
    I[1:n_checkpoint_clusts] <- checkpoint$I
    
    flat_B <- matrix(0, nrow=M, ncol=4*L)
    flat_B[1:n_checkpoint_clusts,1:(4*L)] <- checkpoint$B
    
    alpha_0 <- checkpoint$alpha0
  }
  # <<< initialise or get clustering <<<
  
  
  # >>> initialise variables >>>
  C1 <- Omega
  C0 <- 1 - Omega
  A1 <- A                  #number of alteration reads
  B1 <- D - A              #number of reference reads
  
  A1[is.na(A1)] <- 0
  B1[is.na(B1)] <- 0
  
  #reads number list for each clone
  #these are lists with numbers of agreements/disagreements between A, D and C0, C1
  S1_list <- list()
  S2_list <- list()
  S3_list <- list()
  S4_list <- list()
  for (k in seq_len(K)) {
    S1_list[[k]] <- A1 * C0[, k]
    S2_list[[k]] <- B1 * C0[, k]
    S3_list[[k]] <- A1 * C1[, k]
    S4_list[[k]] <- B1 * C1[, k]
  }
  # <<< initialise variables <<<
  
  
  # >>> Prepare for sampling >>>
  idx_vec <- seq_len(N)
  idx_mat <- seq_len(N*M)
  
  n_element <- length(idx_vec)
  
  if (is.null(dim(prior1)) && length(prior1) == 2) {
    #two variable to a matrix
    prior1 <- t(matrix(rep(prior1, n_element), nrow = 2))
  }
  if (!is.matrix(prior1)) {
    stop("prior1 need to be a matrix of n_element x 2")
  }
  
  #scale the theta_0 prior, to be comparable to the number of total reference reads
  prior0 <- prior0 * sum(B1) / 2
  #scale the theta_i prior, to match the number of alterated reads at pos i 
  prior1 <- prior1 * (apply(A, 1, sum) + 1)
  
  save_it <- floor( max_iter/write_skip )
  
  logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
  assign_all_j <- matrix(0, nrow = save_it, ncol = M)
  
  theta0_all <- matrix(0, nrow = save_it, ncol = 1)
  theta1_all <- matrix(0, nrow = save_it, ncol = n_element)
  
  C_all <- matrix(0, nrow = save_it, ncol = N*K)
  relax_rate_all <- matrix(0, nrow = save_it, ncol = 1)
  
  alpha_0_all <- matrix(0, nrow = save_it, ncol = 1)
  
  adjacency_matrix <- matrix(0, nrow=M, ncol=M)
  
  if( is.null(checkpoint) ){
    relax_rate <- relax_rate_prior[1] / sum(relax_rate_prior)
    
    C <- Omega
  }else{
    relax_rate <- checkpoint$relax_rate
    
    C <- checkpoint$C
  }
  
  C_prior <- Omega
  C_prior[Omega == 1] <- 1 - relax_rate
  C_prior[Omega == 0] <- relax_rate
  if(keep_base_clone){
    C_prior[,1] <- 0
  }
  C_prior_oddlog <- log(C_prior) - log(1 - C_prior)
  Iden_mat <- matrix(0, nrow = M, ncol = K)
  # <<< Prepare for sampling <<<
  
  ## Initialise or get theta
  if( is.null(checkpoint) ){
    theta0 <- stats::rbeta(1, prior0[1], prior0[2])
    theta1 <- stats::rbeta(rep(1,n_element), prior1[,1], prior1[,2])
  }else{
    theta0 <- checkpoint$theta0
    theta1 <- checkpoint$theta1
  }
  
  
  ## Set parent env of all called functions to this env
  environment(sample_I) <- environment()
  environment(sample_c_and_rr) <- environment()
  environment(sample_theta) <- environment()
  environment(sample_alpha_0) <- environment()
  environment(sample_T) <- environment()
  environment(sample_B) <- environment()
  
  
  ## Gibbs sampling
  start_time <- proc.time()
  
  # try first adjusting profiles to the scRNA seq data
  for (it in 2:(max_iter*10)) {
    # Sample clonal profiles and update cell-clone assignment
    sample_I()
    assign_j <- I[t]
    
    # Update C and relax rate
    sample_c_and_rr()
    
    # Sample theta with assigned clones
    sample_theta()
  }
  
  for (it in 2:max_iter) {
    #do Gibbs sampling clustering update
    
    # Sample clonal profiles and update cell-clone assignment
    sample_I()
    assign_j <- I[t]
    
    # Update C and relax rate
    sample_c_and_rr()
    
    # Sample theta with assigned clones
    sample_theta()
    
    # Update clustering and update adjacency matrix
    sample_T()
    
    for( cluster in clusters[ which(sapply(clusters, length)>0) ] ){
      adjacency_matrix[ cluster, cluster ] <- adjacency_matrix[ cluster, cluster ] + 1
    }
    
    # Sample cluster BCR profiles
    sample_B()
    
    # Sample alpha_0
    sample_alpha_0()
    
    # Calculate logLikelihood
    logLik_all[it] <- get_logLik(A1, B1, C, assign_j, theta0, theta1) +
                      get_logLik_BCR(flat_BCR, flat_B, t)
    
    # Save only the 1/write_skip sample
    if( it %% write_skip == 0 ){
      curr <- it / write_skip
      
      assign_all_j[curr, ] <- assign_j
      
      theta0_all[curr, 1] <- theta0
      theta1_all[curr, ] <- theta1
      
      C_all[curr, ] <- C
      relax_rate_all[curr, 1] <- relax_rate
      
      alpha_0_all[curr, 1] <- alpha_0
    }
  }
  end_time <- proc.time()
  
  ## Return values
  n_buin = ceiling(it * buin_frac / write_skip)
  last_it = floor(it / write_skip)
  
  C_prob <- matrix(colMeans(C_all[n_buin:last_it, ]), nrow=N, ncol=K)
  
  theta0 <- mean( theta0_all[n_buin:last_it, ] )
  theta1 <- colMeans( as.matrix(theta1_all[n_buin:last_it, ]) )
  
  alpha_0 <- mean( alpha_0_all[n_buin:last_it,] )
  
  prob_mat_j <- matrix(0, nrow = M, ncol = K)
  for(cell in 1:M){
    prob_mat_j[cell,] <- sapply(1:K, function(k){
      sum(assign_all_j[n_buin:last_it,cell]==k)
    }) / (last_it-n_buin+1)
  }
  
  logLik_post <- get_logLik(A1, B1, C_prob, prob_mat_j, theta0, theta1)
  DIC <- devianceIC(logLik_all[n_buin:last_it], logLik_post)
  
  #collapse output to non-empty clusters
  non_empty_clusters <- which( sapply(clusters, length)>0 )
  
  clusters <- clusters[non_empty_clusters]
  flat_B <- flat_B[non_empty_clusters,]
  I <- I[non_empty_clusters]
  
  for( i in seq_along(clusters) ){
    t[ clusters[[i]] ] <- i
  }
  
  names(t) <- colnames(A)
  
  return_list <- list("theta0" = theta0, "theta1" = theta1,
                      "alpha0" = alpha_0,
                      "alpha0_all" = alpha_0_all[1:last_it],
                      "theta0_all" = as.matrix(theta0_all[1:last_it, ]),
                      "theta1_all" = as.matrix(theta1_all[1:last_it, ]),
                      "element" = idx_mat, "logLik" = logLik_all[2:it],
                      "assign_all_j" = assign_all_j[1:last_it,],
                      "prob_mat_j"= prob_mat_j,
                      "relax_rate" = mean(relax_rate_all[n_buin:last_it]),
                      "C_prob" = C_prob,
                      "C_all" = C_all[1:last_it, ],
                      "relax_rate_all" = relax_rate_all[1:last_it], 
                      "DIC"=DIC,
                      "clusters" = clusters,
                      "t"= t,
                      "I"= I,
                      "B"=flat_B,
                      'Adj'=adjacency_matrix,
                      'time'=end_time-start_time)
  
  return_list
}
