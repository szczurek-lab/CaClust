source('caclust_sampling_gibbs_phases.R')
source('caclust_helper_funcs_fast.R')


# >>> save paths >>>
LIBRARY_PATH <- '../project_libraries'
# <<< save paths <<<


# >>> get dependencies from default library location or project library >>>
.libPaths(LIBRARY_PATH)
for( dep in c('igraph', 'extraDistr', 'matrixStats', 'doParallel', 'dplyr', 'magrittr') ){
  if( !require(dep, character.only = TRUE) ){
    library(dep, 
            lib.loc=LIBRARY_PATH,
            character.only = TRUE)
  }
}
# <<< get dependencies <<<


#read toy data
data <- readRDS('../toy_data.rds')


#these are required inputs for the CaClust model
A <- data$A #matrix of alternate read counts in cells
D <- data$D #matrix of reference read counts in cells
Omega <- data$Omega #estimates of clonal genotypes present in the sample
BCR <- data$BCR #encodings of cells' BCR sequences


#with these we can run a model chain 
chain <- caclust_sampling(A, D, Omega, BCR,
                          keep_base_clone=F)
#we turn the keep_base_clone flag off, as in the toy data all three clones exhibit variants


#the two most important results from the chain:
#probabilities of cell to clone assignment
assignment_probs <- chain$prob_mat_j
#probabilities of variant presence in clones
profile_probs <- chain$C_prob


#we can assign cells to their most probable clone of origin
assignment <- apply(assignment_probs, 1, which.max)
#and check how many cells were correctly assigned
sum( assignment==data$I[data$t] )


#we can take the most probable genotypes
genotypes_estimate <- profile_probs >= .5
#and check how many positions agree with the true hidden genotypes
sum( genotypes_estimate==data$C )


#if we believe the chain has not converged yet, we can resume its run from a checkpoint of the last sampled values
checkpoint <- chain[c('theta0', 'theta1',
                      't','B','clusters','I',
                      'alpha0', 'relax_rate')]
checkpoint$C <- matrix(chain$C_all[500,1:30], ncol=3)

chain_extension <- caclust_sampling(A, D, Omega, BCR,
                                    keep_base_clone=F,
                                    checkpoint=checkpoint)

#and then check again
assignment2 <- apply(chain_extension$prob_mat_j, 1, which.max)
sum( assignment2==data$I[data$t] )
