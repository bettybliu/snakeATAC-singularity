#!/bin/bash

##### 0. HOW TO RUN#	
#	verify the CONTAINER in section 1 below is valid
#	check requested #nodes, memory, time etc. in section 3 and 4
#	save and exit
#	run this shell script from command line with the following command
# bash run_snakemake.sh

##### 1. USER INPUTS
CONTAINER=/oak/stanford/groups/wjg/bliu/containers/atac.sif

##### 2. PREP
#	unlock working directory, dry run to check code, should see a lot of green text displayed
${CONTAINER} "snakemake --unlock -s Snakefile.py; snakemake -ns Snakefile.py"

##### 3. SINGLES
#	for computation-intensive tasks that require no info from other samples in the group 
#       (e.g. alignment, peak calling), split meta file into single-sample files and 
#       submit individual snakemake jobs to slurm. store the job IDs so the group analysis
#	only starts after all jobs are completed.
#	NOTE: the single/group split analysis was implemented because we had difficulties
#	calling sbatch from within the container. This also enables greater portability to
#	non-slurm computing clusters in the future, e.g. google cloud, aws) 
META=$(grep "METADATA_FILE = " snakeATAC_config.py |tr "'" "\n"| sed -n "2p")
rm -rf .tmp; mkdir .tmp

for ((NUM=2; NUM<=$(wc -l < $META); NUM++))
do 
    METAPATH=.tmp/tmp_meta_$(echo $((NUM-1))).txt
    sed -n "1p;${NUM}p" $META > ${METAPATH} 
    # wrap sbatch in another bash because sbatch exits shell immediately after job submission
    SNAKE_CMD="snakemake --nolock -T -p -j 10 -s Snakefile.py --config RULE_GROUP='single' META=${METAPATH}"
    bash -c "sbatch --parsable -p sfgf,wjg,biochem -n 8 -t 24:00:00 --mem-per-cpu 64g --wrap \"${CONTAINER} ${SNAKE_CMD}\" >> .tmp/tmp_joblist.txt"
done
echo "Submitted single analysis jobs:"
cat .tmp/tmp_joblist.txt

##### 4. GROUP
#	after all single-sample jobs are completed, run group analysis tasks
SNAKE_CMD_GROUP="snakemake -T -p -j 10 -s Snakefile.py --config RULE_GROUP='group' META=${META}"
sbatch --dependency=afterok:$(cat .tmp/tmp_joblist.txt|tr '\n' ',' | sed 's/,$/\n/') -p sfgf,wjg,biochem -n 8 -t 24:00:00 --mem-per-cpu 64g --wrap "${CONTAINER} '${SNAKE_CMD_GROUP}'"
