#!/bin/bash
# Account, partition, QOS
#SBATCH --account=kscale
#SBATCH --partition=standard
#SBATCH --qos=standard

# Main job spec
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32G

# Output & monitoring
#SBATCH --job-name="rtw_GAL"
#SBATCH --output=/home/users/emg97/log/log_rtw_GAL.out
#SBATCH --error=/home/users/emg97/log/log_rtw_GAL.err

SCRIPTPATH="/home/users/emg97/emgScripts/"

echo "Start Job $SLURM_JOB_ID, array $SLURM_ARRAY_TASK_ID on $HOSTNAME"

cd $SCRIPTPATH
module load jaspy
module load jasmin-sci

# Run Python with unbuffered output so prints appear immediately in .out
python3 -u realtime_waves_xr.py