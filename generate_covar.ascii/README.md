# generate_covar

This program is used to generate a file holding covariance matrix information
which can be used to generate ensemble perturbations by second-order exact
sampling. The file reads files from the model (usually snapshots from a
model trajectory) and computes EOFs (empirical orthogonal functions). PDAF
provides the routine PDAF_eofcovar for this. The procedure is described
in the documentation at https://pdaf.awi.de/trac/wiki/EnsembleGeneration.

The actual program is in the file
generate_covar.F90

The other source code files provide the infrastructure as in the codes for
the actual data assimilation in the directories online/ and offline/.
The source code files are copied from the directory offline/. 
(Actually, only the files init_parallel_pdaf_offline.F90, initialize_grid_mod.F90,
model_pdaf_mod_offline.F90 are specific for the offline-coupled code. The files
assimilation_pdaf_mod.F90, parallel_pdaf_mod.F90, parser_mpi.F90, and
statevector_pdaf_mod.F90 are identical for the online ond offline codes)

