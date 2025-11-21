#!/bin/bash
#SBATCH --job-name=aracne_bootstrap
#SBATCH --output=job_logs/4a_SHERLOCK3_TFanalysis_makenetwork_sherlock1_par_%j.out    
#SBATCH --error=job_logs/4a_SHERLOCK3_TFanalysis_makenetwork_sherlock1_par_j.err        
#SBATCH --array=1-100%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=100:00:00

ml Java/8-LTS  

cd output/tf_analysis
# Set variables
SEED=$SLURM_ARRAY_TASK_ID
EXPR="sherlock1/counts_cpm.tsv"
TF="humanTF.tsv"
OUTDIR="sherlock1/network_res"

# Create output directory if it doesn't exist
mkdir -p ${OUTDIR}

#java -Xmx5G -jar Aracne.jar -e sherlock1/counts_cpm.tsv -o sherlock1/network_res --tfs humanTF.tsv --pvalue 1E-8 --seed 1 --calculateThreshold

# Run ARACNe
java -Xmx10G -jar Aracne.jar \
  -e ${EXPR} \
  -o ${OUTDIR} \
  --tfs ${TF} \
  --pvalue 1E-8 \
  --seed ${SEED}

  
#Consolidate results
#java -Xmx5G -jar Aracne.jar -o brush/network_res --consolidate
