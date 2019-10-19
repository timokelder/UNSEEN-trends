#$ -S /bin/bash
#$ -l h_rt=2:00:00
#$ -q ded-parallelx.q
#$ -l h_vmem=3G
#$ -t 1
#$ -o //home/timok/timok/ensex/out/out_$JOB_NAME.$JOB_ID.$HOSTNAME.$TASK_ID
#$ -e //home/timok/timok/ensex/err/statistics/err_$JOB_NAME.$JOB_ID.$HOSTNAME.$TASK_ID
#$ -cwd

echo "Got $NSLOTS slots for job $SGE_TASK_ID."

module load R/R-3.4.3-met
#month=$SGE_TASK_ID
Rscript "predictability.R"

