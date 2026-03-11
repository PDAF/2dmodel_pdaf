#!/bin/sh

# Script to perform series of data assimilation runs using the global
# filter ESTKF  with the online-coupled model 'model_pdaf' in the
# directory 'online/'. Varied are the type of ensemble initialization,
# the florgetting factor, and the ensemble size.
#
# For each case, a directory for the outputs is created and the output
# files are moved into this directory.
#
# Lars Nerger, AWI, 2026-03


FILTERTYPE=6     # 6=ESTKF
DIM_ENS=8

for TYPE_ENS_INIT in 1 2 3
do
  for FORGET in 1.0 0.9 0.8 0.7 0.6 0.5 0.4
  do
    for DIM_ENS in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
    do
      expdir=ESTKF_N${DIM_ENS}_etype${TYPE_ENS_INIT}_f${FORGET}
      mkdir $expdir
      echo Run  -dim_ens $DIM_ENS -filtertype $FILTERTYPE \
        -type_ens_init $TYPE_ENS_INIT -forget $FORGET 
      mpirun --oversubscribe -np $DIM_ENS ./model_pdaf -dim_ens $DIM_ENS -filtertype $FILTERTYPE \
        -type_ens_init $TYPE_ENS_INIT -forget $FORGET > out.txt
      mv *nc $expdir
      mv out.txt $expdir
    done
  done
done
