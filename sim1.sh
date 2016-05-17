#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -q qwork
#PBS -o log/
#PBS -e log/
#PBS -N sim1
#PBS -m ea
#PBS -t 1-25
#PBS -V

module load bioinformatics/R/3.2.5

SRC=$HOME/coexpression/may2016simulation/sim1-protocol-mammouth

Rscript $SRC/simulation.R ${index} $PBS_ARRAYID