#!/bin/sh

# PDAF configuration
NENS=4          # Ensemble size
FILTERTYPE=6    # e.g. 6=ESTKF, 7=LESTKF
CRADIUS=5.0     # Cut-off radius for localization
FORGET=1.0      # Forgetting factor for inflation
TYPE_ENS_INIT=1 # Type of ensemle initialization

# Number of processes per model task
nproc=2

# Name of model-PDAF program
EXE='model_pdaf'

# Directory for the experiment
EXPDIR='exp1'

# Directory of inputs
# (relative to EXPDIR or absolute path)
INDIR='../../inputs_2fields'

# Path where the program is stored
# (relative to EXPDIR or absolute path)
PATH_EXE='../../online/'

#################################################################

# Create configuration string for the command line parser
OPTIONS="-dim_ens $NENS -filtertype $FILTERTYPE -cradius $CRADIUS -forget $FORGET -type_ens_init $TYPE_ENS_INIT"


#################################################################
# Create and prepare run directories

echo "-- Prepare ensemble directories"

mkdir $EXPDIR
cd $EXPDIR

for((ENS=1;ENS<=$NENS;ENS++))
do
  ENSDIR=`printf %03d $ENS`

  # Create directory
  mkdir ${ENSDIR}

  # Link program
  cp ${PATH_EXE}/${EXE} $ENSDIR
done

#################################################################
# Link input directory

ln -s $INDIR .

#################################################################
# Create string for multi-program execution

echo "-- Create command line"

var1="-np $nproc -wdir 001 $EXE $OPTIONS"
for((ENS=2;ENS<=$NENS;ENS++))
do
  ENSDIR=`printf %03d $ENS`

  var0="-np $nproc -wdir ${ENSDIR} $EXE $OPTIONS"
  var1="$var1 : $var0"
done
echo 'mpirun ' $var1

#################################################################
# Start the parallel program

echo "-- Run program"

mpirun $var1
