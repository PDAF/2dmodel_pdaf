# Example for running online example in separate directories

This directory presents an application example of the online-coupled
program `model_pdaf` running the ensemble in separate directories
as is recommended for complex models.

The run script run.sh is configured and tested with OpenMPI. Other variants
of MPI might have a different syntax.

To do a test run:
1. compile `model_pdaf` in the directory `online/`
2. execute `./run.sh`

The script `run.sh` will 
* generate the experiment directory `exp1/`
* generate numbered ensemble run directories inside `exp1/`
* copy the program model_pdaf into each directory
* link the directory ../inputs_online_2fields into the test directory
* create the multi-program mpirun command line and display it 
* start the parallel run with ensemble size 4

Outputs will be written in the directory `exp1/001/`.

Some options, e.g. ensemble size, type of DA method and localization radius,
can be set at the beginning of run.sh. One can adapt these to perform
different experiments. One can further adapt EXPDIR to store each experiment
in a different directory.

Each experiment is independent. Thus one can safely remove a directory
with e.g. `rm -rf exp1`.

>[!Important]
> If the script does not create any `.nc`-files and ends with `All nodes which are allocated for this job are already filled.`, this is likely because your machine does not have enough physical cores (8 are required)
> Add the `--oversubscribe` option to mpirun (Open MPI) to allow MPI to run more ranks than the number of available physical cores. This can severely degrade performance.
