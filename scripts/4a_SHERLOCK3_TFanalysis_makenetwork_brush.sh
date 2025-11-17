#!/bin/bash   
#SBATCH --job-name=4a_SHERLOCK3_TFanalysis_makenetwork_brush       # Name for your job
#SBATCH --time 100:00:00                        # Runtime in minutes.
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb                      # Reserve 2 GB RAM for the job
#SBATCH --output=job_logs/4a_SHERLOCK3_TFanalysis_makenetwork_brush%j.out       # Standard out goes to this file
#SBATCH --error=job_logs/4a_SHERLOCK3_TFanalysis_makenetwork_brush%j.err        # Standard err goes to this file
#SBATCH --mail-user=k.phung@umcg.nl        # email you wish to be notified at
#SBATCH --mail-type=ALL                         # ALL will alert you of job beginning, completion, failure etc


#module purge
ml Java/8-LTS  

#java -version

cd output/tf_analysis

#brush ------------------------------------------------------------------------------------------------------------------------------
mkdir -p brush/network_res
java -Xmx5G -jar Aracne.jar -e brush/counts_cpm.tsv -o brush/network_res --tfs humanTF.tsv --pvalue 1E-8 --seed 1 --calculateThreshold
    for i in {1..100}
    do
          java -Xmx5G -jar Aracne.jar -e brush/counts_cpm.tsv -o brush/network_res --tfs humanTF.tsv --pvalue 1E-8 --seed $i
    done

#Consolidate results
java -Xmx5G -jar Aracne.jar -o brush/network_res --consolidate


