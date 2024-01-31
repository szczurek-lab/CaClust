
# >>> save paths >>>

LIBRARY_PATH <- '../project_libraries'

input_save <- "../input/"
output_save <- "../results/"

parameters_save <- "params_caclust.csv"

# <<< save paths <<<




# >>> parameters >>>

params <- read.csv( parameters_save, comment.char = '#' )

n_chain <- params$n_chain
n_cores <- params$n_cores

subject <- params$subject
n_clone <- params$n_clone

Psi <- NULL

alpha_0 <- 1
alpha_prior <- c(1, 1)

relax_rate_prior <- c(1, 4)

max_iter <- 1000
min_iter <- 1
write_skip <- 5
# <<< parameters >>>




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


# >>> read input data >>>

A <- read.csv( paste0(input_save, subject, '/A.csv'), 
               check.names = FALSE,
               row.names = 1 ) %>% as.matrix()
D <- read.csv( paste0(input_save, subject, '/D.csv'),
               check.names = FALSE,
               row.names = 1 ) %>% as.matrix()

BCR <- readRDS( paste0(input_save, subject, '/BCR.rds') )

Omega <- read.csv( paste0(input_save, subject, '/Omega_', n_clone, '.csv'),
                   row.names = 1 ) %>% as.matrix()

# <<< read input data <<<


## check input data
if ( !(all(rownames(A) == rownames(D))) ){
  stop("Variants for A and D are not identical.")
}
if ( !(all(colnames(A) == colnames(D))) ){
  stop("Cells for A and D are not identical.")
}

## Match exome-seq and scRNA-seq data
if ( !any(rownames(D) %in% rownames(Omega)) ){
  stop("No matches in variant names between C and D arguments.")
}

## match variants
common_vars <- intersect( rownames(Omega), rownames(D) )
A <- A[common_vars,, drop = FALSE]
D <- D[common_vars,, drop = FALSE]
Omega <- Omega[common_vars,, drop = FALSE]

## match BCR
common_cells <- intersect( colnames(D), dimnames(BCR)[[1]] )
A <- A[,common_cells, drop = FALSE]
D <- D[,common_cells, drop = FALSE]
BCR <- BCR[common_cells,,, drop = FALSE]

## pass data to specific functions
registerDoParallel(cores = n_cores)

## Run n_chain times independently
ids_list <- foreach::foreach(i=1:n_chain, .export = ls(environment())) %dopar% {
  source('caclust_sampling_gibbs_phases.R')
  source('caclust_helper_funcs_fast.R')
  
  caclust_sampling(A=A, D=D, Omega=Omega, BCR=BCR, 
		   alpha_0=alpha_0, alpha_prior=alpha_prior,
  		   max_iter=max_iter, write_skip=write_skip, min_iter=min_iter,
  		   relax_rate_prior=relax_rate_prior)
} 

output_dir <- paste0(output_save, subject)

dir.create(output_dir, showWarnings=FALSE)

saveRDS(ids_list, paste0(output_dir, '/', subject, '_gibbs_results.rds'))

