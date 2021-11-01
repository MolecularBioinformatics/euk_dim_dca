#!/bin/bash
#SBATCH --job-name=DCA
#SBATCH --account=nn9795k
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2GB

## Set up job environment:
set -o errexit  # Exit the script on any error
# set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2019.03
module load HMMER/3.2.1-intel-2018b
module load MUSCLE/3.8.31-intel-2018b
module load Julia/1.6.1-linux-x86_64

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9795k/conda/biopython 

# path to config files
CONF=/cluster/projects/nn9795k/selina/configs
# CONF=../configs

python run_workflow.py $CONF/config_3ep3.txt $CONF/paths_saga.txt all -r
