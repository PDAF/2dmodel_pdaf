#!/bin/sh

# PDAF configuration
NENS=4          # Ensemble size
FILTERTYPE=6    # e.g. 6=ESTKF, 7=LESTKF
CRADIUS=5.0     # Cut-off radius for localization

# Number of processes per model task
nproc=2

# Name of model-PDAF program
EXE='model_pdaf'

# Path where the program is stored
PATH_EXE='../online/'

#################################################################

# Create configuration string for the command line parser
OPTIONS="-dim_ens $NENS -filtertype $FILTERTYPE -cradius $CRADIUS"


#################################################################
# Link input directory

ln -s ../inputs_online_2fields .

#################################################################
# Create and prepare run directories

echo "Prepare ensemble directories"
for((ENS=1;ENS<=$NENS;ENS++))
do
  ENSDIR=`printf %03d $ENS`

  # Create directory
  mkdir $ENSDIR

  # Copy program
  cp ${PATH_EXE}/${EXE} $ENSDIR
done

#################################################################
# Create string for multi-program execution

var1="-np $nproc -wdir 001 $EXE $OPTIONS"
for((ENS=2;ENS<=$NENS;ENS++))
do
  ENSDIR=`printf %03d $ENS`

  var0="-np $nproc -wdir ${ENSDIR} $EXE $OPTIONS"
  var1="$var1 : $var0"
done

#################################################################
# Start the parallel program

mpirun $var1
