#!/bin/bash
#PBS -l nodes=5:ppn=20
#PBS -l walltime=7:00:00
#PBS -o log/
#PBS -e log/
#PBS -N sim1
#PBS -m ea
#PBS -t 1-100

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR


Rscript ${SRC}/simulation.R ${index} $PBS_ARRAYID