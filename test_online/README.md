This directory presents an application example of the online-coupled
program model_pdaf running separate directories as is recommended
for complex models.

The run script run.sh is configures and tested with OpenMPI. Other variants
of MPI might have a different syntax.

To do a test run:
1. compile model_pdaf in the directory online/
2. execute ./run.sh

The script run.sh will 
- generate numbered ensemble run directories
- copy the program model_pdaf into each directory
- link the directory ../inputs_online_2fields into the test directory
- create the multi-program mpirun command line 
- start the paralle run

Outputs will be written in the directory 001/

The script cleanup.sh can be used to remove the ensemble run 
directories and link inputs_online_2fields.
