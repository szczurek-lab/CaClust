LIBRARY_PATH <- '../project_libraries'

# >>> install dependencies >>>
.libPaths(LIBRARY_PATH)
for( dep in c('igraph', 'extraDistr', 'matrixStats', 'doParallel', 'dplyr', 'magrittr') ){
  if( !require(dep, character.only = TRUE) ){
    install.packages(dep, lib=LIBRARY_PATH)
    library(dep, 
            lib.loc=LIBRARY_PATH,
            character.only = TRUE)
  }
}
# <<< install dependencies <<<

