# 2dmodel_pdaf - advanced implementation of tutorial model with PDAF

This directory contains a variant of the tutorial model with 2 model fields coupled to PDAF. 
The model uses the domain-decomposition parallelization as in the parallel tutorial cases.
The code supports all ensemble-based methods (ensemble Kalman filters, particle filters, hybrid Kalman-nonlinear filter).

In this variant we follow the same structure as we use in real applications with complex models.
In particular, this variant uses an advanced scheme to handle the fields in the state vector.

Notes on the model: 
* The tutorial model is very simple and without real dynamics. Even though the model uses a 2-dimensional domain, it is practically a 1-dimensional case (which is just rotated by 45 degrees). This model allows to perform DA with very small ensembles with a minimum of 3 states. This is far lower than we would use in real cases, but it is useful to perform some experiments on a notebook computer. For real cases one would rather use an ensemble of at least 20 states.
* Because the model is practically 1-dimensional, using the localized ensemble filters does not improve the state estimate. However, using a localize filter shows the effects on the space of the resulting state estimate.
* The model does not read restart files. Thus, one cannot perform a realistic sequence of offline analysis steps. Only a single offline analysis step is done for demonstration.
  

## Structure of the examples

The directories containing the model and coupled model-PDAF codes are the following
* 2dmodel/ - contains the pure tutorial model without coupling to PDAF
* offline/ - contains the analysis program for offline-coupled assimilation with PDAF
* online/ - contains the model with online-coupling to PDAF

Additional directories are
* generate_covar/ - contains the code for a program generating a covariance matrix file from model outputs. This can be used to initialize an ensemble using second-order exact sampling (activated using '-type_ens_init 2' when running the DA cases.
* inputs_online_2fields/ - contains the input files, e.g. observation files and files to provide ensemble states. It also contains the script 'pdaf_tutorial_online.py' used to generate the input files.
* plotting/ - contains plot scripts
* test_online/ - contains the script 'run.sh' which shows an example of how to run the online-coupled case in separate directories



## Examples

### Compiling

For compiling use
```
make <target> PDAF_ARCH=<arch> PDAFDIR=<mydir>
```
where
* `<target>`: the name of the program you like to compile
* `<arch>`: The name of the make.arch file without .h as in other compilations with PDAF (see files in make.arch/ of PDAF)
* `<mydir>`: The path to the PDAF directory
Possible names for <target> are shown when just running the command line without specifying <target>.

### Overview of output files

Running the online or offline DA produces a number of output files. There are separate files for the first model field 'A' and the second model field 'B'. Below, we list the file just showing field 'A'
* `stateA_step*_for.txt / stateA_step*_ana.txt` - files holding the ensemble mean of field A of the model for the forecast ('for') or analysis ('ana').
* `ensA_XX_step*_for.txt / ensA_XX_step*_ana.txt` - files holding each single ensemble member state (with index 'XX') for the forecast ('for') or analysis ('ana').
* `varianceA_step*_for.txt / varianceA_step*_ana.txt` - files holding the ensemble variance of field A of the model for the forecast ('for') or analysis ('ana').

### Plotting

The directory plotting/ contains Python plot scripts. They can be used like
```
plot_diff.py FILE1 FILE2      # Plot the difference of two files and display their RMS difference
plot_file.py FILE             # Plot the field in one file
plot_obs.py OBSFILE           # Plot one observation file
```

### Examples for online coupled DA

The online coupled DA is implemented with PDAF's fully parallel mode. Thus, the number of processes has to be a multiple of the ensemble size. 

To run the global ESTKF filter with ensemble size 4 over 18 time steps with an analysis step at each second step use
```
mpirun -np 4 ./model_pdaf -dim_ens 4 -filtertype 6
```
To plot the difference of the final state estimate of this experiments with the true model state use
```
python ../plotting/plot_diff.py stateA_step18_ana.txt ../inputs_online_2fields/trueA_step18.txt
```

To run the same experiment as before but using domain-decomposition with 2 processes per model use
```
mpirun -np 8 ./model_pdaf -dim_ens 4 -filtertype 6
```
The outputs of this run should be the same as before. However, in the upper part of screen output you will see information about the domain decomposition and the dimension of the fields in the state vector will reflect the decomposition.

The domain-localized filter LESTKF can be run e.g. using
```
mpirun -np 4 ./model_pdaf -dim_ens 4 -filtertype 7 -cradius 5.0
```
Here, 'cradius' is used to set the localization radius to 5 grid points. Plotting the state estimate as before shows that the localization changes the shape the state estimate and leads to a larger error (this is a particularity of this model and the fact that the ensemble states have the same shape as the true state).

### Examples for offline coupled DA

The offline coupled DA performs one single analysis step. It is also implemented with parallelization using domain decomposition. However, one can run it with a single process.

To run the global ESTKF filter with ensemble size 4 at step 1
```
./offline_pdaf -dim_ens 4 -filtertype 6 -step 1
```
(Depending the MPI-library you use you might need to write 'mpirun -np 1' at the beginning of the command)
Here setting step=1 means that the observations at time step 1 are assimilated. The program is implemented to always read the same input files. Thus, in a more realistic case, one would change this in the source code or add a script that copies restart files from the model to the default names. 

## Twin experiments

### Generating synthetic observations for twin experiments

The example allows to generate files holding synthetic observations from a model forecast and to later use these files as inputs for twin experiments. As there are no real observations for this model, this is for demonstrating how to use this functionality.

The example is implemented such that when it is run using parallelization, separate files for each sub-domain will be written. Accordingly, when performed the actual twin experiment (see below) one has to use the same number of processes per model. The example further assumes that observations have been generated for each time step.

The example is further implemented so that synthetic observations can only be generated for one observation type at a time. Thus, to generate files for each of the two observation types in the example, one has to do two runs. In each case one performs the run with ensemble size 1:
```
mpirun -np 2 ./model_pdaf -dim_ens 1 -filtertype 100 -delt_obs 1 -write_state F -write_ens F
mpirun -np 2 ./model_pdaf -dim_ens 1 -filtertype 100 -delt_obs 1 -write_state F -write_ens F -assim_A F -assim_B T
```
In the second run, the observation type A is deactivated, while observation type B is switched on. These runs will generate 4 files, i.e. one file for each observation type and sub-domain:
```
obsA_syn_0000.nc, obsA_syn_0001.nc, obsB_syn_0000.nc, obsB_syn_0000.nc
```
Here the 4-digit number is the process rank which has written the file. The observations are written as a list of values. The information about observation locations or errors are initialized by the respective observation module when performing a twin experiment.

Note: The functionality to generate synthetic observation is only available for the online coupled case.

### Performing a twin experiment

When the files holding the synthetic observations have been generated one can use them in twin experiments. In this case, each observation module reads the regular observation files to initialize all values. Then the file holding the synthetic observations is read, which overwrites the values of the actual observations.

The twin experiments can be run assimilation any combination of the observation types as long as the necessary files holding the synthetic observations were generated before.

To run the ESTKF with ensemble size 4, assimilating observations of type A, and using 2 processes per model to be consistent with the generated synthetic observation files one can execute:
```
mpirun -np 8 ./model_pdaf -dim_ens 4 -filtertype 6 -delt_obs 2 -twin_experiment T
```
To run the analogous experiment assimilating synthetic observations of both types, one can execute:
```
mpirun -np 8 ./model_pdaf -dim_ens 4 -filtertype 6 -delt_obs 2 -twin_experiment T -assim_B T
```

Note: One can also perform twin experiments with the offline coupled case, as long as one has generated synthetic observation files before with the online-coupled program.