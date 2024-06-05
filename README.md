# Bulk ATAC Analysis in a Singularity Container
Author: Betty Liu <liubetty@stanford.edu>

> [!NOTE] 
> This repo is now maintained on the [public Greenleaf Lab repo](https://github.com/GreenleafLab/snakeATAC_singularity)

### Introduction
Welcome to the new era of bulk ATAC analysis in the Greenleaf Lab!!! This document describes how to run snakeATAC analysis using a singularity container. The Stanford Sherlock website has a [description of what a container means](https://www.sherlock.stanford.edu/docs/software/using/singularity/) and why Singulairty was chosen over Docker. The templates and slurm commands have only been tested with Sherlock, but with the containerization of all packages it is easy to port this to another HPC platform in the future. Sherlock has ```singularity``` enabled for all users by default so you can check out the commands by typing ```singularity help``` .

### Steps
1. Clone this repository into your analysis folder and navigate into the downloaded folder:
```
git clone https://github.com/bettybliu/snakeATAC-singularity.git snakeATAC
cd snakeATAC
```
2. Go to file ```snakeATAC_config.py```, change variables in the ```User Inputs``` section based on your experiment. 
3. Go to file ```meta.txt``` and change the name of experiments and path to your reads. (Optional) You could run the following command to automatically generate ```meta.txt``` based on fastq file names. If your fastq filenames don't contain ```_R1``` and ```_R2```, use the shell command ```rename``` to change the names. Check ```meta.txt``` after running the command to make sure the sample labels are correct.	
```
python snakeATAC_config.py
```
4. Go to file ```fastq_screen.conf```, check the paths of the genomes to screen against.
5. Go to file ```Snakefile.py``` and change the analysis you wish to perform by commenting/uncommenting the output file names in ```rule_group_dict```.
6. Run snakeATAC with the following command.
```
bash run_snakemake.sh
```  

### Using Singularity for Other Stuff
The following container was built to run snakeATAC, but can also be used to run direct commands. It has three conda environments: base, py35, py27. SnakeATAC uses both py35 and py27 environments. I have a downloaded copy in my oak folder that's used by this pipeline by default. You can download it to a different location using ```singularity pull --arch amd64 library://liubetty/default/atac:latest```. 

```CONTAINER=/oak/stanford/groups/wjg/bliu/containers/atac.sif```

```singularity shell ${CONTAINER}``` opens an intearctive shell within the container

```singularity exec ${CONTAINER} COMMAND OPTIONS``` to use any COMMAND within the container (conda is not initiated by default if you run the container this way)

```singularity exec ${CONTAINER} bash -c "source activate py35; COMMAND OPTIONS"``` to activate the conda environment py35 inside the container and run COMMAND

```${CONTAINER} COMMAND OPTIONS``` to use any COMMAND within the container with the conda environment py35 already activated -- this was achieved by adding the activation commands into the ```%runscript``` section of the singularity definition file during build, and when you run a container directly on the command line like this, it automatically sources the code in ```%runscript``` first. 
- ```${CONTAINER} "conda activate py27; COMMAND OPTIONS"``` to use any COMMAND within the py27 environment inside the container 
- Consider adding the following to your ```~/.bashrc```:```alias sing='/oak/stanford/groups/wjg/bliu/containers/atac.sif'```, then just do ```sing COMMAND OPTIONS```


### Snippet from ```run_snakemake.sh``` for quick reference
```
#!/bin/bash

##### 0. HOW TO RUN#
#       verify the CONTAINER in section 1 below is valid
#       check requested #nodes, memory, time etc. in section 3 and 4
#       save and exit
#       run this shell script from command line with the following command
# bash run_snakemake.sh

##### 1. USER INPUTS
CONTAINER=/oak/stanford/groups/wjg/bliu/containers/atac.sif

##### 2. PREP
#       unlock working directory, dry run to check code, should see a lot of green text
${CONTAINER} "snakemake --unlock -s Snakefile.py; snakemake -ns Snakefile.py"

##### 3. SINGLES
#       for computation-intensive tasks that require no info from other samples in the group
#       (e.g. alignment, peak calling), split meta file into single-sample files and
#       submit individual snakemake jobs to slurm. store the job IDs so the group analysis
#       only starts after all jobs are completed.
#       NOTE: the single/group split analysis was implemented because we had difficulties
#       calling sbatch from within the container. This also enables greater portability to
#       non-slurm computing clusters in the future, e.g. google cloud, aws)
META=$(grep "METADATA_FILE = " snakeATAC_config.py |tr "'\"" "\n"| sed -n "2p")
rm -rf .tmp; mkdir .tmp

for ((NUM=2; NUM<=$(wc -l < $META); NUM++))
do
    METAPATH=.tmp/tmp_meta_$(echo $((NUM-1))).txt
    sed -n "1p;${NUM}p" meta.txt > ${METAPATH}
    # wrap sbatch in another bash because sbatch exits shell immediately after job submission
    SNAKE_CMD="snakemake --nolock -T -p -j 10 -s Snakefile.py \
                --config RULE_GROUP='single' META=${METAPATH}"
    bash -c "sbatch --parsable -p sfgf,wjg,biochem -n 8 -t 24:00:00 --mem-per-cpu 64g \
        --wrap \"${CONTAINER} ${SNAKE_CMD}\" >> .tmp/tmp_joblist.txt"
done
echo "Submitted single analysis jobs:"
cat .tmp/tmp_joblist.txt

##### 4. GROUP
#       after all single-sample jobs are completed, run group analysis tasks
#       need to remove the snakeATAC.txt output generated from single samples analysis first
SNAKE_CMD_GROUP="rm -rf snakeATAC.txt; snakemake -T -p -j 10 -s Snakefile.py \
                    --config RULE_GROUP='group' META=${META}"
sbatch --dependency=afterok:$(cat .tmp/tmp_joblist.txt|tr '\n' ',' | sed 's/,$/\n/') \
    -p sfgf,wjg,biochem -n 8 -t 24:00:00 --mem-per-cpu 64g \
    --wrap "${CONTAINER} '${SNAKE_CMD_GROUP}'"
```

### Singularity definition file
```
Bootstrap: library
From: ubuntu:18.04

%setup


%files


%environment
    # clear any user defined R libraries
    export R_LIBS_USER=''

%post
    # install essential ubuntu packages
    apt-get update && apt-get -y upgrade
    apt-get -y install \
        build-essential \
        make \
        wget \
        git \
        zip \
        unzip \
        vim \
        locales \
        libglu1-mesa-dev \
        libglib2.0-0 \
        libxext6 \
        libsm6 \
        libxrender1 \
        libreadline6-dev \
        libz-dev \
        gawk

    locale-gen en_US.UTF-8
    ln -s /lib/x86_64-linux-gnu/libreadline.so.7.0 /lib/x86_64-linux-gnu/libreadline.so.6
    ln -s /usr/lib/x86_64-linux-gnu/libicui18n.so.60 /usr/lib/x86_64-linux-gnu/libicui18n.so.58
    rm -rf /var/lib/apt/lists/*
    apt-get clean

    # install Anaconda3 
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O anaconda.sh 
    /bin/bash anaconda.sh -bfp /usr/local/anaconda
    rm anaconda.sh

    # set conda path (temporary)
    export PATH="/usr/local/anaconda/bin:$PATH"

    # conda configuration of channels from .condarc file
    conda config --file /.condarc --add channels defaults
    conda config --file /.condarc --add channels conda-forge
    conda config --file /.condarc --add channels bioconda


    ########################################################
    ###################### Python 3.5 ######################
    ########################################################

    conda create --name py35 python=3.5

    # install basic python packages and command line packages
    conda install -n py35 \
        numpy=1.15.2 \
        scipy=1.1.0 \
        pandas=0.23.4 \
        matplotlib=3.0.0 \
        cython=0.28.5 \
        hdf5=1.12.0 \
        gsl=2.2 \
        openjdk=8.0.152 

    # install basic bioinformatics command line packages
    conda install -n py35 -c bioconda\
        deeptools=3.2.1 \
        bowtie2=2.3.4.3 \
        samtools=1.7 \
        pysam=0.14.1 \
        bedtools=2.30.0 \
        cutadapt=1.18 \
        fastqc=0.11.9 \
        preseq=2.0.3 \
        subread=2.0.1 \
        ngs-bits=2018_04 \
        fastq-screen=0.14.0 \
        ucsc-bedgraphtobigwig=357 \
        ucsc-bedclip=332 \
        snakemake=3.6.1 \
        parso=0.7.0 \
        ipython=6.5.0

    conda install -n py35 -c dranew bcl2fastq=2.19.0

    # picard jar
    wget "https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar" \
        -O /usr/local/bin/picard.jar
    
    # snakeATAC tools and resources
    mkdir /usr/local/snakeATAC
    wget "https://bettyliu.s3.us-west-1.amazonaws.com/snakeATAC/atac_tools.zip" \
        -O /usr/local/snakeATAC/atac_tools.zip
    unzip /usr/local/snakeATAC/atac_tools.zip -d /usr/local/snakeATAC/
    rm /usr/local/snakeATAC/atac_tools.zip

    ########################################################
    ###################### Python 2.7 ######################
    ########################################################

    # create a python2.7 conda environment
    conda create --name py27 python=2.7
    conda install -n py27 -c bioconda \
        macs2=2.1.4 \
        pysam=0.15.3 \
        cython=0.29.14

    conda install -n py27 \
        pandas=0.24.2 \
        r-base=3.6.1 \
        r-essentials=3.6.0

    # install custom packages only compatible with python2.7
    git clone https://github.com/GreenleafLab/NucleoATAC.git
    /usr/local/anaconda/envs/py27/bin/pip install NucleoATAC/
    rm -rf NucleoATAC/

    # clean up
    conda clean --tarballs

    # set conda path (permanent)
    echo ". /usr/local/anaconda/etc/profile.d/conda.sh" >> ${SINGULARITY_ENVIRONMENT}
    echo "conda activate" >> ${SINGULARITY_ENVIRONMENT}

    
%runscript
    echo "Running snakeATAC singularity container >>>"
    exec bash -c "source activate py35; $*"

%startscript
    

%test
    

%labels
    Author liubetty@stanford.edu
    Group greenleaf.stanford.edu
    Version v0.0.9.1

%help
    This is a singularity container for running bulk ATACseq data analysis 
    using snakemake (snakeATAC).
```
