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
These source code files are copied from the directory offline/src. 
(Actually, only the files init_parallel_pdaf_offline.F90, initialize_grid_mod.F90,
and model_pdaf_mod_offline.F90 are specific for the offline-coupled code.
The other files are identical for the online ond offline codes)

## Compilation

To compile use
```
make generate_covar PDAF_ARCH=<arch> PDAFDIR=<mydir>
```
analogous to compiling the online example program.

## Running
Execute
```
./generate_covar
```
This will read the `ens*nc` files from ../inputs_2fields and generate the output file `covar.nc`. Running
```
./generate_covar -infile 'true'
```
will read the `true*nc` files and generate the covariance matrix file from this data.


