#Variables:
#C: The observed Omega matrix of clonal profiles
#C: The true C matrix of clonal profiles
#prob_mat: the likelihood of the A and D reads
#S1_list: for each clone: a number of reads w/ mutations in cells from clones w/o the mutation
#S2_list: w/o -||- w/o
#S3_list: w/ -||- w/
#S4_list: w/o -||- w/
#Psi: prior probability of cluster -> clone assignment
#A1: numbers of observed mutations
#B1: numbers of observed reference
#K: number of clones
#Q: number of clusters
#t: cluster -> clone assignment
#it: iteration number
rcatlp <- extraDistr::rcatlp
rdirichlet <- extraDistr::rdirichlet

# add small epsilon, to avoid -infty logs
epsilon <- 10 ** (-50)


loglike_RNA <- function(cells, C, theta1, theta0, A, D){
  alternate <- t( A[,cells] )
  reference <- t( D[,cells]-A[,cells] )
  
  alternate_probabilities <- theta0*(1-C) + theta1*C
  reference_probabilities <- 1-alternate_probabilities
  
  colSums( alternate %*% log(alternate_probabilities+epsilon) + 
           reference %*% log(reference_probabilities+epsilon) )
}

Geweke_Z <- function(X, first=0.1, last=0.5) {
  
  if (is.null(dim(X))){
    X <- as.matrix(X, ncol = 1)}
  N <- nrow(X)
  A <- X[1:floor(first*N), , drop = FALSE]
  B <- X[ceiling(last*N):N, , drop = FALSE]
  
  A_col_mean <- colMeans(A)
  B_col_mean <- colMeans(B)
  A_col_var <- rowSums((t(A) - A_col_mean)^2) / (nrow(A) - 1)#^2
  B_col_var <- rowSums((t(B) - B_col_mean)^2) / (nrow(A) - 1)#^2
  
  min_var <- 10^(-50)
  Z <- (A_col_mean - B_col_mean) / sqrt(A_col_var + B_col_var + min_var)
  
  Z
}


devianceIC <- function(logLik_all, logLik_post) {
  
  logLik_mean = mean(logLik_all)
  logLik_var = var(logLik_all)
  
  p_D_Spiegelhalter = -2 * logLik_mean -  (-2 * logLik_post)
  DIC_Spiegelhalter = -2 * logLik_post + 2 * p_D_Spiegelhalter
  
  p_D_Gelman = 2 * logLik_var
  DIC_Gelman = -2 * logLik_post + 2 * p_D_Gelman
  
  DIC = DIC_Gelman
  
  #cat(paste("DIC:", round(DIC, 2), 
  #          "D_mean:", round(-2 * logLik_mean, 2), 
  #          "D_post:", round(-2 * logLik_post, 2), 
  #          "logLik_var:", round(logLik_var, 2), "\n"))
  
  list("DIC" = DIC, 
       "logLik_var" = logLik_var, 
       "D_mean" = -2 * logLik_mean, 
       "D_post" = -2 * logLik_post, 
       "DIC_Gelman" = DIC_Gelman, 
       "DIC_Spiegelhalter" = DIC_Spiegelhalter)
}


get_logLik <- function(A1, B1, C, Assign, theta0, theta1) {
  if (is.null(dim(Assign)) || length(dim(Assign)) == 1) {
    Assign_prob <- matrix(0, length(Assign), ncol(C))
    for (i in seq_len(length(Assign))) {
      Assign_prob[i, Assign[i]] = 1
    }
  } else {
    Assign_prob <- Assign
  }
  
  prob_mat <- C %*% t(Assign_prob)
  
  Lik_mat <- (exp(log(theta1) * A1 + log(1 - theta1) * B1) * prob_mat + 
              exp(log(theta0) * A1 + log(1 - theta0) * B1) * (1 - prob_mat))
  
  logLik <- (sum(log(Lik_mat), na.rm = TRUE) + 
             sum(lchoose(A1 + B1, A1), na.rm = TRUE))
  
  if(is.infinite(logLik)){
    print(Lik_mat)
    stop()
  }
  
  logLik
}


get_logLik_BCR <- function(flat_BCR, flat_B, t){
  logLik_all <- sum(log( flat_B[t,][ as.logical(flat_BCR) ] ), na.rm=TRUE)
  
  if(is.infinite(logLik_all)){
    print('bcr err')
    
    stop()
  }
  logLik_all
}


sample_T <- function(){
  # Likelihoods of all cells' RNA in different clones
  alternate_probabilities <- theta0*(1-C) + theta1*C
  reference_probabilities <- 1-alternate_probabilities
  
  scRNA_logLik_mat <- t(A1) %*% log(alternate_probabilities+epsilon) + 
                      t(B1) %*% log(reference_probabilities+epsilon)
  scRNA_logLik_mat <- scRNA_logLik_mat - apply(scRNA_logLik_mat, 1, max)
  
  # Likelihoods of all cells' BCR in different clusters (including potential clusters)
  non_empty_clusters <- which( sapply( clusters, length ) > 0 )
  queued_clusters <- setdiff( 1:M, non_empty_clusters )
  
  BCR_logLik_all <- flat_BCR %*% log(t(flat_B + epsilon))
  
  for( cell in seq_len(M) ){
    # get cluster populations without current cell
    old_t <- t[cell]
    clusters[[ old_t ]] <<- setdiff( clusters[[old_t]], cell )
    
    cluster_sizes <- sapply( clusters, length )
    
    non_empty_clusters <- which( cluster_sizes > 0 )
    
    # if no potential clusters remain, queue as many potential clusters 
    # as possible, but no more than remaining cells to be sampled
    if( length(queued_clusters)==0 ){
      queued_clusters <- setdiff( 1:M, non_empty_clusters )
      q <- length(queued_clusters)
      
      queued_clusters <- queued_clusters[ 1:min(q, M+1-cell) ]
      q <- length(queued_clusters)
      
      repeated_prior <- matrix( t(g), nrow=q*L, ncol=4, byrow=TRUE)
      
      flat_B[queued_clusters,] <<- rdirichlet(q*L, repeated_prior) %>% matrix( nrow=q )
      I[queued_clusters] <<- sample( 1:K, q, prob=Psi, replace=TRUE )
      
      # we calculate the BCR likelihoods of the remaining cells to the queued clusters.
      # we don't transpose if we only queued one cluster
      if( q==1 ){
        BCR_logLik_all[cell:M, queued_clusters] <- flat_BCR[cell:M,] %*% 
                                                   log(flat_B[queued_clusters,]+epsilon)
      }else{
        BCR_logLik_all[cell:M, queued_clusters] <- flat_BCR[cell:M,] %*% 
                                                   log(t(flat_B[queued_clusters,]+epsilon))
      }
    }
    
    # CRP restaurant prior probabilities on clustering
    CRP_logPrior <- log( cluster_sizes[non_empty_clusters] / (alpha_0+M-1) )
    
    # Adding the possibility of joining a new cluster, 
    # extracting the first from the queued clusters
    new_cluster <- queued_clusters[ 1 ]
    queued_clusters <- queued_clusters[ -1 ]
    
    non_empty_clusters <- append( non_empty_clusters, new_cluster )
    
    CRP_logPrior <- append( CRP_logPrior, log(alpha_0/(alpha_0+M-1)) )
    
    # Likelihoods of cell's BCR in different clusters
    BCR_logLikes <- BCR_logLik_all[cell, non_empty_clusters]  
    BCR_logLikes <- BCR_logLikes - max(BCR_logLikes)
    
    # Likelihoods of cell's RNA in different clusters
    scRNA_logLikes <- scRNA_logLik_mat[ cell, I[non_empty_clusters] ]
    
    # Normalise
    log_p <- CRP_logPrior + BCR_logLikes + scRNA_logLikes
    p <- exp(log_p - max(log_p))
    
    # Sample the new cluster for the cell
    new_t <- sample( non_empty_clusters, 1, prob=p )
    
    if( new_t==new_cluster ){
      clusters[[ new_cluster ]] <<- cell
      t[cell] <<- new_t
    }else{
      clusters[[ new_t ]] <<- append( clusters[[new_t]], cell )
      t[cell] <<- new_t
    }
  }
}

# sampling BCR profiles
sample_B <- function(){
  # matrix of cluster-cell assignment
  clust_cell <- matrix(0, nrow=M, ncol=M)
  
  clust_cell[ cbind(t, 1:M) ] <- 1
  
  # faster B sampling
  non_empty_clusters <- which( sapply(clusters, length)>0 )
  k <- length(non_empty_clusters)
  
  up_priors <- clust_cell[non_empty_clusters,] %*% flat_BCR %>% matrix( ncol=4 )
  up_priors <- up_priors + matrix( t(g), nrow=k*L, ncol=4, byrow=TRUE)
  
  flat_B[non_empty_clusters,] <<- rdirichlet(k*L, up_priors) %>% matrix( nrow=k )
  
  # sample the profiles of potential clusters, for clustering sampling later
  empty_clusters <- setdiff( 1:M, non_empty_clusters )
  repeated_prior <- matrix( t(g), nrow=(M-k)*L, ncol=4, byrow=TRUE)
  
  flat_B[empty_clusters,] <<- rdirichlet((M-k)*L, repeated_prior) %>% matrix( nrow=(M-k) )
}

sample_I <- function(){
  non_empty_clusters <- which( sapply(clusters, length) > 0 )
  l_clusters <- length( non_empty_clusters )
  
  cell_cluster <- matrix(0, nrow=M, ncol=M)
  cell_cluster[ cbind(1:M, t[1:M]) ] <- 1 
  cell_cluster <- cell_cluster[,non_empty_clusters]
  
  clusters_alt_reads <- A1 %*% cell_cluster
  clusters_ref_reads <- B1 %*% cell_cluster
  
  alternate_probabilities <- theta0*(1-C) + theta1*C
  reference_probabilities <- 1-alternate_probabilities
  
  scRNA_clusters_logLik <- t(clusters_alt_reads) %*% log(alternate_probabilities+epsilon) + 
                           t(clusters_ref_reads) %*% log(reference_probabilities+epsilon)
  
  # rcatlp gives values from 0 to K-1
  I[ non_empty_clusters ] <<- rcatlp( l_clusters, scRNA_clusters_logLik ) + 1

  # sample the assignment of potential clusters for clustering sampling later
  empty_clusters <- setdiff( 1:M, non_empty_clusters )
  I[ empty_clusters ] <<- sample( 1:K, size=length(empty_clusters), 
                                  replace=TRUE, prob=Psi )
}


sample_c_and_rr <- function(){
  if (it > 0.1 * min_iter){
    diff0 <- sum(Omega == C)
    diff1 <- sum(Omega != C)
    relax_rate <<- r_rate <- stats::rbeta(1, 
                                          relax_rate_prior[1] + diff1,
                                          relax_rate_prior[2] + diff0)
    
    #C_prior is the prior probability of C=1 given C:=Omega
    C_prior <- Omega
    C_prior[Omega == 1] <- 1 - r_rate
    C_prior[Omega == 0] <- r_rate
    if(keep_base_clone){
      C_prior[,1] <- 0
    }
    C_prior_oddlog <- log(C_prior) - log(1 - C_prior)
  }
  
  #Iden_mat is the assignment of cells to clones (through clusters)
  Iden_mat[,] <- 0
  Iden_mat[cbind(1:M, I[t])] <- 1 
  
  # calculate log_probability matrix with genotype 0 and 1
  P0_mat <- A1 * log(theta0+epsilon) + B1 * log(1 - theta0+epsilon)
  P1_mat <- A1 * log(theta1+epsilon) + B1 * log(1 - theta1+epsilon)
  
  oddR_log <- P1_mat %*% Iden_mat  - P0_mat %*% Iden_mat 
  oddR_log <- oddR_log + C_prior_oddlog
  oddR_log[which(oddR_log > 50)] <- 50
  oddR_log[which(oddR_log < -50)] <- -50
  C_prob_tmp <- exp(oddR_log) / (exp(oddR_log) + 1)
  
  C <<- C_n <- matrix(stats::rbinom(N*K, size = 1, C_prob_tmp), nrow=N)
  
  for (k in seq_len(K)) {
    S1_list[[k]] <<- A1 * (1 - C_n[,k])
    S2_list[[k]] <<- B1 * (1 - C_n[,k])
    S3_list[[k]] <<- A1 * C_n[, k]
    S4_list[[k]] <<- B1 * C_n[, k]
  }
}

sample_theta <- function(){
  S1_wgt <- S2_wgt <- 0 # weighted S1
  S3_wgt <- S4_wgt <- matrix(0, nrow = N, ncol = M)
  for (k in seq_len(K)) {
    idx <- which(I[t] == k)
    S1_wgt <- S1_wgt + sum(S1_list[[k]][,idx], na.rm = TRUE)
    S2_wgt <- S2_wgt + sum(S2_list[[k]][,idx], na.rm = TRUE)
    S3_wgt[,idx] <- S3_wgt[,idx] + S3_list[[k]][,idx]
    S4_wgt[,idx] <- S4_wgt[,idx] + S4_list[[k]][,idx]
  }
  
  S3_wgt[,] <- rowSums(S3_wgt, na.rm = TRUE)
  S4_wgt[,] <- rowSums(S4_wgt, na.rm = TRUE)
  
  theta0 <<- stats::rbeta(1, 
                          prior0[1] + S1_wgt, 
                          prior0[2] + S2_wgt)
  theta1 <<- stats::rbeta(rep(1, n_element),
                          prior1[,1] + S3_wgt[idx_vec],
                          prior1[,2] + S4_wgt[idx_vec])
}

sample_alpha_0 <- function(){
  # k is the number of non-empty clusters
  k <- sum(sapply(clusters, length)>0)
  
  #sample new value of sigma
  xSigma <- rbeta(1, alpha_0+1, M)
  
  y <- alpha_prior[1]
  x <- alpha_prior[2]
  
  from_first <- (runif(1, max=y+k+M*(x-log(xSigma))) < y+k)
  if(from_first){
    alpha_0 <<- rgamma(1, shape=y+k, rate=x-log(xSigma))
  }else{
    alpha_0 <<- rgamma(1, shape=y+k-1, rate=x-log(xSigma))
  }
}


