# CaClust 

Thank you for using CaClust, a probabilistic graphical model of Follicular Lymphoma clonal structure. Here you will find all functions needed for model inference as well as examples on real and toy data.

# Structure

- scripts
- input.zip
- toy_data.rds

Directory scripts contains all functions needed for model inference as well as two examples of running CaClust on real and toy data. Input.zip contains data for the real world example and toy_data.rds contains data for the toy example.

## scripts

- scripts
  - caclust_helper_funcs_fast.R
    Contains sampling functions for all model variables needed in the main inference routine in caclust_sampling_gibbs_phases.R.
  - caclust_sampling_gibbs_phases.R
    Contains the function for running one chain of CaClust model inference, new or from checkpoint, with descriptions of all parameters and input.
  - install_dependencies.R
    Script for installing all CaClust dependencies in a directory within this project.
  - params_caclust.csv
    Subject parameter file for the real data example.
  - run_caclust_gibbs_phases.R
    Example of running CaClust chains on real data from input.zip, which must be unpacked to a directory specified by input_path.
  - run_toy_example.R
    Example of running CaClust on toy_data.rds, extracting cell assignment and clonal profiles, and restarting a model run from checkpoint.
