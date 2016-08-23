#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o log/
#PBS -e log/
#PBS -N sim1
#PBS -m ea
#PBS -t 1-50

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR


Rscript ${SRC}/simulation.R ${index} $PBS_ARRAYID
