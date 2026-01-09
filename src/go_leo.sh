#!/bin/bash
#SBATCH --account="IscrB_EXCEED"
#SBATCH --job-name="test"
#SBATCH --time=1:00:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --qos=boost_qos_dbg
#SBATCH --error=test.err


module load profile/candidate
module load nvhpc/25.3

./mhit36