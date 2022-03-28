#!/bin/bash
#SBATCH --job-name=analyze_tet_A
#SBATCH -J download
#SBATCH --time=72:00:00
#SBATCH --account=fc_pickeringlab
#SBATCH -p savio2_htc
#SBATCH --ntasks=4
#SBATCH --mem=84g
#SBATCH --mail-type=END
#SBATCH --mail-user=shruteek@berkeley.edu
#SBATCH --output=%j_output.txt
#SBATCH --error=%j_error.txt
## Command(s) to run:
module load python
## export PATH=$PATH:/global/scratch/projects/fc_pickeringlab/tools/miniconda3/bin/
dir='/global/scratch/projects/fc_pickeringlab/Optimized-sgRNA-Design/'
module load python
source activate /global/scratch/projects/fc_pickeringlab/miniconda3/envs/bowtie
cd /global/scratch/projects/fc_pickeringlab/Optimized-sgRNA-Design/datafiles
startTime=$(date +%s)
bowtie -S -a -v 2 --verbose --large-index 101_R1_illumina_index --suppress 1,5,6,7 bowtie_test_guides.fastq bowtie_illumina_test.sam >& job.out
timePassed=$(expr $(date +%s) - $startTime)
echo $timePassed >> job.out