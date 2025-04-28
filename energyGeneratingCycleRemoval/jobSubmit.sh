#!/bin/bash

#PBS -l walltime=23:59:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l place=free:group=board
#PBS -J 1-102
#PBS -r y
#PBS -N egcRemovalLP
#PBS -A GBAnonlinearSJ
#PBS -m abe
#PBS -M youremail@gmail.com

cd $PBS_O_WORKDIR
echo “JobID: “$PBS_JOBID

#module load intel/xe2013
#module load R/3.0.2

module load R/3.6.1

cat $PBS_NODEFILE


Rscript energyGeneratingCycleRemovalLP.R $PBS_ARRAY_INDEX &>logfiles/energyGeneratingCycleRemovalLP.$PBS_ARRAY_INDEX.$PBS_JOBID.out





