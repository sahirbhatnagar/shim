#~/bin/bash
#PBS -m ea
#PBS -l walltime=15:00:00
#PBS -l nodes=1:ppn=1
#PBS -q qwork
#PBS -r n
#PBS -N sim4
#PBS -o log/
#PBS -e log/
#PBS -V

module load bioinformatics/R/3.2.5

SRC=$HOME/coexpression/august2016simulation/linear

Rscript $SRC/figures-for-manuscript.R ${index} $PBS_ARRAYID
