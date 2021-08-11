# Configuration file 

import os
import sys
import glob
from collections import defaultdict

################## User Inputs ###################
SPECIES_GENOME = "hg38"
FASTQ_DIR = '/oak/stanford/groups/wjg/bliu/data/test/fastq'
METADATA_FILE = '/oak/stanford/groups/wjg/bliu/data/test/snakeATAC_sing/meta.txt'
FASTQ_SCREEN_CONF = "/oak/stanford/groups/wjg/bliu/data/test/snakeATAC_sing/fastq_screen.conf"

# unfortunately not everything in the shared lab resources folders have consistent naming
# 	so always double check if these files exist if you are working outside of hg19/hg38/mm9
BEDS_DICT = {"mm9": "/home/groups/wjg/lab/snakeATAC/resources/mm9/mm9_tss.ensembl.bed",
             "hg19": "/home/groups/wjg/lab/snakeATAC/resources/hg19/hg19_tss.ensembl.bed",
             "hg38": "/oak/stanford/groups/wjg/share/resources/hg38/hg38.tss.bed"}
BEDS = {"TSS": BEDS_DICT[SPECIES_GENOME]}
TSS = '/oak/stanford/groups/wjg/share/resources/hg38/hg38.tss.bed' if SPECIES_GENOME=="hg38" \
        else '/home/groups/wjg/lab/snakeATAC/resources/%s/%s_tss.ensembl.bed' % (SPECIES_GENOME, SPECIES_GENOME)

REFERENCE_FILE = '/home/groups/wjg/lab/genomes/%s/%s' % (SPECIES_GENOME, SPECIES_GENOME)
CHROM_SIZES = '/home/groups/wjg/lab/genomes/gSizes/%s.all.genomsize' % (SPECIES_GENOME)
GENOME_SIZE_DICT = {"mm9": 1.87e9, "sacCer3": 1.2e7, "hg19": 2.7e9, "hg38": 3.0e9}
EFFECTIVE_GENOME_SIZE = GENOME_SIZE_DICT[SPECIES_GENOME]
BLACKLIST = None # can use blacklist in chraccr analysis instead
################## End User Inputs ###############

# Adding paths to useful tools for running snakeATAC 
# 	DON'T change this section, the directories are fixed inside the singularity container file system
PICARD_JAR = '/usr/local/bin/picard.jar'
SNAKE_DIR = '/usr/local/snakeATAC/'
ATAC_TOOLS = os.path.join(SNAKE_DIR, 'atac_tools')

# metadata file
def make_meta(filename):
    r1_files = list(map(os.path.abspath, glob.glob(os.path.join(FASTQ_DIR, "*_R1*.f*"))))
    if (len(r1_files) < 1):
        sys.exit("No fastqs with _R1 found.")
    r2_files = [os.path.join(os.path.dirname(r1_file), os.path.basename(r1_file).replace('R1', 'R2')) for r1_file in
                r1_files]
    if all([os.path.isfile(r2_file) for r2_file in r2_files]) is False:
        sys.exit("Not all matching _R2 files found.")
    sample_labels = [os.path.basename(r1_file).split("_R1")[0] for r1_file in r1_files]
    with open(filename, 'w') as outfile:
        outfile.write("\t".join(["Name", "Read1", "Read2"]) + "\n")
        for sample_label, r1_file, r2_file in zip(sample_labels, r1_files, r2_files):
            outfile.write("\t".join([sample_label, r1_file, r2_file]) + "\n")

# make the meta data file
if __name__ == "__main__":
    make_meta(METADATA_FILE)
