#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=rickerTestRun
#SBATCH --array=1-2


module load arcc/1.0 gcc/12.2.0 r/4.2.2

cd /project/modelscape/analyses/MissingTS/missing-data

config=Model_Runs/RickerConfig.txt

datFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
parFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
clsize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
saveFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
index1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
index2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)


echo "This is array task ${SLURM_ARRAY_TASK_ID} loading data file ${datFile} and parameter file ${parFile} using cluster size ${clsize} and saving file at ${saveFile}. We are doing data set numbers ${index1} through ${index2}" >> output.txt

#Rscript Model_Runs/modelruns_ricker.R ${datFile} ${parFile} ${clsize} ${saveFile} ${index1} ${index2} > Model_Runs/outputRickerA_${SLURM_ARRAY_TASK_ID}.txt