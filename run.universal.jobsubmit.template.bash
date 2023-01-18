#!/bin/bash
#FOR RUNNING ON SLURM (BIOWULF):
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 10g
#SBATCH --job-name=steveTest
#SBATCH -o logs/%x.%N.%j.out # STDOUT
#SBATCH -e logs/%x.%N.%j.out # STDERR
#
#FOR RUNNING ON SGE (CCAD):
###############################
# SGE settings Here
# Basically, if a line starts with "#$", then you can enter any
# qsub command line flags .  See the qsub(1) man page.
# Redirect the STDOUT and STDERR files:
#$ -o $JOB_ID_$JOB_NAME.stdout
#$ -j yes
#$ -S /bin/bash
#$ -q all.q
#$ -l h_vmem=10G
# Do some validation checking, and bail on errors
#$ -w e
# Operate in the current directory
#$ -cwd
# End SGE Settings
###############################
# Usable variables: JOB_TITLE, HOSTNAME, JOB_NAME, JOB_ID

if [ "${JOB_ID}" != "" ]; then
   echo "JOB_TITLE=$JOB_TITLE"
   echo "HOSTNAME=$HOSTNAME"
   echo "JOB_NAME=$JOB_NAME"
   echo "JOB_ID=$JOB_ID"
else
  if [ "${SLURM_JOB_ID}" != "" ]; then
   echo "SLURM_JOB_ID=$SLURM_JOB_ID"
   echo "SLURM_JOB_NAME=$SLURM_JOB_NAME"
   echo "SLURM_JOB_NODELIST=$SLURM_JOB_NODELIST"
   echo "SLURM_LOCALID=$SLURM_LOCALID"
   echo "SLURM_NODEID=$SLURM_NODEID"
   echo "SLURM_PROCID=$SLURM_PROCID"
   echo "SLURM_SUBMIT_DIR=$SLURM_SUBMIT_DIR"
   echo "SLURM_SUBMIT_HOST=$SLURM_SUBMIT_HOST"
  else
     echo "WARNING: JOB QUEUE ENGINE NOT DETECTED!"
     echo "ENV vars JOB_ID and SLURM_JOB_ID are not set!"
     echo "Are you running this just from bash? If so, that's fine."
  fi
fi

echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "$( date "+%Y.%m.%d.%H%M%S" )"
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"



echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "$( date "+%Y.%m.%d.%H%M%S" )"
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"

echo "SCRIPT_EXECUTION_SUCCESSFULLY_COMPLETED"



