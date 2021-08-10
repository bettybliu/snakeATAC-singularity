import pandas as pd
import os
import sys

include: "snakeATAC_config.py"

# changing metadata file variable to allow slurm submission of single samples
RULE_GROUP = config['RULE_GROUP'] if ('RULE_GROUP' in config.keys()) else 'all'
METADATA_FILE = config['META'] if ('META' in config.keys()) else METADATA_FILE

metadata = pd.read_table(METADATA_FILE, index_col = False)
sample_labels = metadata.Name.tolist()
output_suffix = sample_labels[0] if RULE_GROUP=="single" else "group"

rule_group_dict = {
    "single": [
        expand("output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam", sample_label = sample_labels),
        expand("output/plots/qc/insert_size/{sample_label}_insert_size_histogram.pdf",sample_label = sample_labels),
        expand("output/peaks/{sample_label}_summits.bed", sample_label = sample_labels),
        expand("output/plots/TSS/{sample_label}.TSS.insertion_profile.png", sample_label = sample_labels),
        expand("output/counts/{sample_label}.saf",sample_label = sample_labels),
        expand("output/plots/{bed_data}/{sample_label}.{bed_data}.Vplot.eps", bed_data = BEDS.keys(), sample_label = sample_labels),
        expand("output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html", sample_label = sample_labels),
        expand("output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html", sample_label = sample_labels),
        expand("output/bams/qc/complexity/{sample_label}.extrapolated_yield.txt", sample_label = sample_labels),
        expand("output/fastqs/qc/{sample_label}_R1_trimmed_screen.html", sample_label = sample_labels),
        expand("output/coverage_data/{sample_label}.insertion_track.bw",sample_label = sample_labels),
	],
    "group": [
        expand("output/beds/qc/in_peaks/{sample_label}.reads_in_peaks.txt",sample_label = sample_labels),
        "output/plots/qc/sample_correlation_plot.pdf",
        "output/bams/qc/compiled_flagstats.txt",
        "output/bams/qc/compiled_counts.txt",
        "output/plots/qc/compiled_idxstats.mito.pdf",
        "output/bams/qc/compiled_picard.dedup_metrics.txt",
        "output/peaks/all_samples/all_sample_peaks.narrowPeak.bed",
        "output/beds/qc/all_samples.nuc_matrix.txt",
        "output/coverage_data/qc/multiBamSummary_out_readCounts.tab",
	]
}
rule_group_dict['all'] = rule_group_dict['single'] + rule_group_dict['group']

rule all:
    input:
        rule_group_dict[RULE_GROUP]
    output:
        "success/snakeATAC_%s.txt" % output_suffix
    shell: 
        "echo $(date)  > {output};"
        "echo created by Evan Boyle and the Greenleaf lab >> {output};"
        "echo singularity container integration by Betty Liu >> {output}"

rule trim_adapters_seqpurge:
    input:
        left = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample_label]["Read1"],
        right = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample_label]["Read2"],
    output:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right ="output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"
    params:
        error_out_file = "error_files/{sample_label}_trim",
        run_time="2:30:00",
        cores="8",
        memory="6000",
        job_name="trimming"
    benchmark: 
        "benchmarks/trimming/{sample_label}.txt"
    threads: 8
    shell: 
        "SeqPurge -a1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a2 CTGTCTCTTATACACATCTGACGCTGCCGACGA -threads {threads} -out1 {output.left} -out2 {output.right} -in1 {input.left} -in2 {input.right}"

rule trim_adapters_cutadapt:
    input:
        left = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample_label]["Read1"],
        right = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample_label]["Read2"],
    output:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right ="output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz",
    params:
        error_out_file = "error_files/{sample_label}_trim",
        run_time="2:30:00",
        cores="1",
        memory="6000",
        job_name="trimming"
    benchmark: 
        "benchmarks/trimming/{sample_label}.txt"
    shell: 
        "cutadapt -a Trans2_rc=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A Trans1_rc=CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length 20 --overlap=5 -o {output.left} --paired-output {output.right} {input.left} {input.right}"

ruleorder: trim_adapters_seqpurge > trim_adapters_cutadapt

rule fastqc_unmapped_trimmed:
    input:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right = "output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"
    output:
        "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html",
        "output/fastqs/qc/{sample_label}_R2_trimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.zip",
        "output/fastqs/qc/{sample_label}_R2_trimmed_fastqc.zip"
    params:
        error_out_file = "error_files/{sample_label}_trim_fastqc",
        run_time="00:15:00",
        cores="1",
        memory="6000",
        job_name="fastqc"
    benchmark: 
        "benchmarks/fastqc/{sample_label}_trim.txt"
    shell: 
        "fastqc {input.left} {input.right} --outdir=" + "output/fastqs/qc/"

rule fastqc_unmapped_untrimmed:
    input:
        left = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample_label]["Read1"],
        right = lambda wildcards: metadata.loc[metadata.Name == wildcards.sample_label]["Read2"],
    output:
        lh = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html",
        rh = "output/fastqs/qc/{sample_label}_R2_untrimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        lz = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.zip",
        rz = "output/fastqs/qc/{sample_label}_R2_untrimmed_fastqc.zip",
        # if run is interrupted
        #templh = temp("output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.html"),
        #temprh = temp("output/fastqs/qc/" + input[1].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.html"),
        #templz = temp("output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.zip"),
        #temprz = temp("output/fastqs/qc/" + input[1].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.zip"),
    params:
        error_out_file = "error_files/{sample_label}_untrim_fastqc",
        run_time="00:15:00",
        cores="1",
        memory="6000",
        job_name="fastqc"
    benchmark: 
        "benchmarks/fastqc/{sample_label}_untrim.txt"
    run:
        # files must be deleted if a previous run failed -- they do not have sample_label in them and cannot be listed as temp files
        shell("rm -f output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.html"),
        shell("rm -f output/fastqs/qc/" + input[1].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.html"),
        shell("rm -f output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.zip"),
        shell("rm -f output/fastqs/qc/" + input[1].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.zip"),
        shell("fastqc {input.left} {input.right} --outdir=output/fastqs/qc/;"),
        shell("mv output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")  + "_fastqc.html {output.lh};"),
        shell("mv output/fastqs/qc/" + input[1].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")  + "_fastqc.html {output.rh};"),
        shell("mv output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")  + "_fastqc.zip {output.lz};"),
        shell("mv output/fastqs/qc/" + input[1].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")  + "_fastqc.zip {output.rz};"),

rule fastqscreen_trimmed:
    input:
        left = rules.trim_adapters_seqpurge.output.left,
        right = rules.trim_adapters_seqpurge.output.right,
        conf = FASTQ_SCREEN_CONF,
    output:
        "output/fastqs/qc/{sample_label}_R1_trimmed_screen.txt",
        "output/fastqs/qc/{sample_label}_R2_trimmed_screen.txt",
        "output/fastqs/qc/{sample_label}_R1_trimmed_screen.html",
        "output/fastqs/qc/{sample_label}_R2_trimmed_screen.html",
    threads: 4
    params:
        error_out_file = "error_files/{sample_label}_fastq_screen",
        run_time="00:15:00",
        cores="4",
        memory="6000",
        job_name="fastq_screen"
    benchmark: 
        "benchmarks/fastq_screen/{sample_label}_fastq_screen.txt"
    shell: 
        "fastq_screen --subset 500000 --outdir output/fastqs/qc/ --threads {threads} --conf {input.conf} --bowtie2 '-X 2000 --no-mixed --no-discordant' --aligner bowtie2 {input.left} {input.right}"

rule run_bowtie:
    input:
        idx = REFERENCE_FILE + ".1.bt2",
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right ="output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz",
    output:
        bam = "output/bams/unprocessed/{sample_label}.bam",
        idx = "output/bams/unprocessed/{sample_label}.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_bowtie",
        run_time = "4:59:00",
        cores = "8",
        memory = "8000",
        job_name = "bwt2"
    benchmark: "benchmarks/bowtie/{sample_label}.txt"
    threads: 8
    shell:  # -X 2000 # prevents mates separated by a lot
        "bowtie2 -X 2000 --threads {threads} --rg-id {wildcards.sample_label} --rg 'SM:{wildcards.sample_label}' -x " + REFERENCE_FILE + " -1 {input.left} -2 {input.right} | samtools view -b -S - | samtools sort -o output/bams/unprocessed/{wildcards.sample_label}.bam -; "
        "samtools index output/bams/unprocessed/{wildcards.sample_label}.bam; "

rule estimate_library_complexity:
    input:
        bam = rules.run_bowtie.output.bam
    output:
        lc = "output/bams/qc/complexity/{sample_label}.extrapolated_yield.txt",
        c = "output/bams/qc/complexity/{sample_label}.downsampled_yield.txt"
    params:
        error_out_file = "error_files/{sample_label}_estimate_lc",
        run_time = "1:00:00",
        cores = "1",
        memory = "8000",
        job_name = "lc_extrap"
    benchmark: "benchmarks/preseq/{sample_label}.txt"
    threads: 1
    shell: 
        "preseq lc_extrap -P -o {output.lc} -B {input.bam}; " +
        "preseq c_curve -P -s 100000 -o {output.c} -B {input.bam}"

rule calc_flagstats:
    input:
        "output/bams/unprocessed/{sample_label}.bam"
    output:
        "output/bams/qc/flagstats/{sample_label}.flagstat.txt" 
    params:
        error_out_file="error_files/flagstats",
        run_time="00:05:00",
        cores="1",
        memory="3000",
        job_name="flagstat"
    shell: 
        "samtools flagstat {input} | awk '{{print \"{wildcards.sample_label}\\t\" $0}}' > {output};"

rule calc_idxstats:
    input:
        "output/bams/unprocessed/{sample_label}.bam"
    output:
        "output/bams/qc/idxstats/{sample_label}.idxstats.txt" 
    params:
        error_out_file="error_files/idxstats",
        run_time="00:05:00",
        cores="1",
        memory="1000",
        job_name="idxstats"
    shell: 
        "samtools idxstats {input} | awk '{{print \"{wildcards.sample_label}\\t\" $0}}' > {output};"


rule plot_flagstats:
    input:
        expand("output/bams/qc/flagstats/{sample_label}.flagstat.txt", sample_label=sample_labels)
    output:
        table = "output/bams/qc/compiled_flagstats.txt",
        pdf = "output/plots/qc/compiled_flagstats.pdf"
    params:
        error_out_file="error_files/flagstat_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_flagstat"
    shell: 
        "awk 'BEGIN {{OFS = \"\\t\"; print \"sample_label\",\"total\",\"secondary\",\"supplementary\",\"duplicates\",\"mapped\",\"paired\",\"read1\",\"read2\",\"proper_pair\",\"both_mapped\",\"singletons\",\"separate_chr\",\"separate_chr_mapq_above5\"}} FNR == 1 && NR != 1 {{print \"\"}} FNR == 1 {{printf $1}} {{printf \"\\t\" $2 }} END {{print \"\"}} ' {input} > {output.table};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_boxplot.R {output.table} read_count {output.pdf}"

rule plot_idxstats:
    input:
        expand("output/bams/qc/idxstats/{sample_label}.idxstats.txt", sample_label=sample_labels)
    output:
        qc_table = "output/bams/qc/compiled_idxstats.txt",
        mito_table = "output/bams/qc/compiled_idxstats.mito_fraction.txt",
        qc_pdf = "output/plots/qc/counts/compiled_idxstats.counts.pdf",
        mito_pdf = "output/plots/qc/compiled_idxstats.mito.pdf",
    params:
        error_out_file="error_files/idxstats_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_idxstats"
    shell: 
        "awk 'BEGIN {{OFS = \"\\t\"; print \"sample_label\",\"chr\",\"ref_length\",\"mapped\",\"unmapped\"}} {{totals[$1] += $4}} $2 == \"chrM\" {{mito[$1] += $4}} {{print}} \
            END {{print \"sample_label\",\"total_reads\",\"mito_reads\",\"mito_percent\" > \"{output.mito_table}\"; for(s in totals) {{print s,totals[s],mito[s],mito[s]/totals[s] * 100 > \"{output.mito_table}\"}} }}' {input} > {output.qc_table};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.mito_table} sample_label mito_percent {output.mito_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.mito_table} sample_label total_reads {output.qc_pdf};"

# From JB nature paper: Reads mapping to the mitochondria, unmapped contigs and chromosome Y were removed and not considered
# so remove the extra stuff...
# consider another call to filter out mitochondrial reads?
# see this http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3646.html#methods
rule rm_mito: # chokes on sam file with weird headers
    input:
        bam = rules.run_bowtie.output.bam,
        idx = rules.run_bowtie.output.idx
    output:
        bam = "output/bams/noMT/{sample_label}.noMT.bam",
        idx = "output/bams/noMT/{sample_label}.noMT.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_remove_mitochondrial_reads",
        run_time = "00:30:00",
        cores = "1",
        memory = "4000",
        job_name = "rm_mt_reads"
    threads: 1
    shell: 
        "samtools idxstats {input.bam} | cut -f 1 | grep -v chrM | xargs samtools view -b {input.bam} > {output.bam}; " # something like this
        "samtools index {output.bam}"

rule filter_bams:
    input: 
        bam = rules.rm_mito.output.bam,
        idx = rules.rm_mito.output.idx
    output: 
        bam = "output/bams/filtered/{sample_label}.noMT.filtered.bam"
    params:
        error_out_file="error_files/{sample_label}_filtered_bams",
        mapq_threshold="20",
        run_time="00:30:00",
        cores="1",
        memory="8000",
        job_name="filter_bams"
    threads: 1
    run:
        # -F 1804: exclude flag, exludes unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
        # -f 2: flags to require, properly aligned
        # -q 30: exlude low MAPQ, set as parameter to adjust
        if BLACKLIST is None:
            shell("samtools view -F 1804 -f 2 -q {params.mapq_threshold} -b {input.bam} > {output}")
        else:
            shell("samtools view -F 1804 -f 2 -q {params.mapq_threshold} -b {input.bam} | bedtools intersect -v -abam - -b " + BLACKLIST + " -wa > {output}")

# perhaps just mark duplicates and remove downstream. For now we remove
# PICARD does not play nice with the read names we get from the SRA...
# PICARD is not good at handling memory, if crashes rerun with more memory...
# include below if not standard illumnia naming convention. Tries to extract read location information.
# READ_NAME_REGEX=null
rule rm_duplicates_picard: 
    input: # low mapping quality reads can also be removed
        bam = rules.filter_bams.output.bam,
    output:
        bam = "output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam",
        idx = "output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam.bai",
        raw_metrics = "output/picard/duplicates/raw/picard_dedup_metrics_{sample_label}.txt",
        parsed_metrics = "output/picard/duplicates/parsed/picard_dedup_metrics_{sample_label}.parsed.txt", 
    params:
        error_out_file =  "error_files/{sample_label}_picard_rmdup",
        run_time="01:00:00",
        cores="1",
        memory="20000",
        job_name="picard_rm_duplicate_reads"
    benchmark: "benchmarks/picard_MarkDuplicates/{sample_label}.txt"
    threads: 1
    shell:  # -Xms4g # this seems to get the process killed... # WE CAN INCLUDE READ_NAME INFO if we have illumina reads...
        "java -XX:ParallelGCThreads=3 -jar " + PICARD_JAR + " MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.raw_metrics} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT READ_NAME_REGEX=null; "
        "samtools index {output.bam}; " # and index
        "grep -A 1 ESTIMATED_LIBRARY_SIZE {output.raw_metrics} | tail -1 | cut -f 9-10 | xargs echo {wildcards.sample_label} | tr ' ' $'\t' > {output.parsed_metrics};"

rule count_bam_reads:
    input:
        "output/bams/unprocessed/{sample_label}.bam",
        "output/bams/noMT/{sample_label}.noMT.bam",
        "output/bams/filtered/{sample_label}.noMT.filtered.bam",
        "output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam",
    output:
        "output/bams/qc/counts/{sample_label}.counts.txt"
    params:
        error_out_file="error_files/bam_counts",
        run_time="00:10:00",
        cores="1",
        memory="3000",
        job_name="count_bam"
    shell: 
        "u=$(samtools flagstat output/bams/unprocessed/{wildcards.sample_label}.bam | head -1 | cut -f 1 -d ' ');"
        "M=$(samtools flagstat output/bams/noMT/{wildcards.sample_label}.noMT.bam | head -1 | cut -f 1 -d ' ');"
        "f=$(samtools flagstat output/bams/filtered/{wildcards.sample_label}.noMT.filtered.bam | head -1 | cut -f 1 -d ' ');"
        "d=$(samtools flagstat output/bams/deduped/{wildcards.sample_label}.noMT.filtered.deduped.bam | head -1 | cut -f 1 -d ' ');"
        "echo {wildcards.sample_label} $'\t' $u $'\t' $M $'\t' $f $'\t' $d > output/bams/qc/counts/{wildcards.sample_label}.counts.txt"

rule merge_bam_files:
    input:
        expand("bams/deduped/{sample_label}.noMT.filtered.deduped.bam", sample_label=sample_labels)
    output:
        bam = "bams/merged/merged.noMT.filtered.deduped.bam",
        idx = "bams/merged/merged.noMT.filtered.deduped.bam.bai"
    params:
        error_out_file="error_files/merge_bams",
        run_time="00:30:00",
        cores="8",
        memory="6000",
        job_name="merge_bams"
    benchmark: "benchmarks/merge_bams.txt"
    threads: 8
    shell: 
        "samtools merge -r -l 2 -@ {threads} {output.bam} {input};"
        "samtools index {output.bam}"

#awk 'BEGIN {command = "paste "} FILENAME != previous {command = command  "<(cut -f " NF " FILENAME ")"; previous=FILENAME }' output/beds/{sample_label}.insertions.bed.gz
#ls output/beds/WT-3h_S14_L001.insertions.bed.gz output/beds/Mz-3h_S10_L001.insertions.bed.gz output/beds/Kz-3h_S8_L001.insertions.bed.gz  | while read line; do echo -n "<(cut -f 11 " $line ") "; done
rule plot_duplicate_stats:
    input:
        expand("output/picard/duplicates/parsed/picard_dedup_metrics_{sample_label}.parsed.txt", sample_label=sample_labels)
    output:
        duplicate_table = "output/bams/qc/compiled_picard.dedup_metrics.txt",
        est_libsize_pdf = "output/plots/qc/compiled_picard_rmdup.est_libsize.pdf",
        duplicate_percent_pdf = "output/plots/qc/compiled_picard_rm_dup.duplicate_percent.pdf"
    params:
        error_out_file="error_files/bam_count_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_bam"
    shell: 
        "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"percent_duplication\",\"estimated_library_size\"}} {{print}}' {input} > {output.duplicate_table};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.duplicate_table} sample_label percent_duplication {output.duplicate_percent_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.duplicate_table} sample_label estimated_library_size {output.est_libsize_pdf};"
        

rule plot_bam_reads:
    input:
        expand("output/bams/qc/counts/{sample_label}.counts.txt", sample_label=sample_labels)
    output:
        fraction_table = "output/bams/qc/compiled_counts.fraction.txt",
        count_table = "output/bams/qc/compiled_counts.txt",
        disjoint_table = "output/bams/qc/compiled_counts.disjoint.txt",
        qc_fraction_pdf = "output/plots/qc/counts/compiled_counts.qc_fraction.pdf",
        total_count_pdf = "output/plots/qc/counts/compiled_counts.total.pdf",
        post_filter_count_pdf = "output/plots/qc/counts/compiled_counts.post_filter.pdf",
        post_dedup_count_pdf = "output/plots/qc/counts/compiled_counts.post_dedup.pdf",
        post_mito_count_pdf = "output/plots/qc/counts/compiled_counts.post_mitochondria.pdf",
        disjoint_count_pdf = "output/plots/qc/counts/compiled_counts.disjoint.pdf",
        mito_fraction_pdf = "output/plots/qc/counts/compiled_counts.mito_fraction.pdf",
        filter_fraction_pdf = "output/plots/qc/counts/compiled_counts.filter_fraction.pdf",
        duplicate_fraction_pdf = "output/plots/qc/counts/compiled_counts.duplicate_fraction.pdf"
    params:
        error_out_file="error_files/bam_count_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_bam"
    shell: 
        "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"total\",\"post_mitochondria\",\"post_filter\",\"post_dedup\"}} {{print $1,$2,$3,$4,$5}}' {input} > {output.count_table};"
        "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"mitochondria\",\"filtered\",\"duplicate\",\"informative\"}} {{print $1,$2 - $3,$3 - $4,$4 - $5,$5}}' {input} > {output.disjoint_table};"
        "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"mitochondria\",\"filtered\",\"duplicate\"}} {{print $1,($2-$3)/$2,($3-$4)/$3,($4-$5)/$4}}' {input} > {output.fraction_table};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_boxplot.R {output.fraction_table} fraction_reads_removed {output.qc_fraction_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.count_table} sample_label total {output.total_count_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.count_table} sample_label post_filter {output.post_filter_count_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.count_table} sample_label post_dedup {output.post_dedup_count_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.count_table} sample_label post_mitochondria {output.post_mito_count_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.fraction_table} sample_label filtered {output.filter_fraction_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.fraction_table} sample_label duplicate {output.duplicate_fraction_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_bargraph.R {output.fraction_table} sample_label mitochondria {output.mito_fraction_pdf};"
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_stackbargraph.R {output.disjoint_table} sample_label {output.disjoint_count_pdf};"

rule plot_insert_size_hist_snake:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx
    output:
        histogram_plot = "output/plots/qc/insert_size/{sample_label}_insert_size_histogram.pdf",
        histogram_data = "output/bams/qc/insert_size/{sample_label}_insert_size_histogram.data.txt"
    params:
        error_out_file="error_files/{sample_label}_insert_size_hist",
        run_time="20:00:00",
        cores="1",
        memory="6000",
        job_name="plot_insert_size"
    threads: 1
    shell: 
        "samtools view {input.bam} | cut -f 9 | awk '{{print ($1**2)**0.5}}' | sort -n | uniq -c | awk 'BEGIN {{print \"sample_label\\tinsert_size\\tcount\"}} {{print \"{wildcards.sample_label}\\t\" $2 \"\\t\" $1}}' > {output.histogram_data};" 
        "/usr/local/anaconda/envs/py27/bin/" + "Rscript --vanilla " + ATAC_TOOLS + "/qc_histogram.R {output.histogram_data} insert_size count 1 5 'Insert size' {output.histogram_plot};"

rule plot_insert_size_hist_picard:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx
    output:
        histogram_plot = "output/plots/qc/insert_size/{sample_label}_insert_size_histogram.pdf",
        histogram_data = "output/picard/insert_size/{sample_label}_insert_size_histogram.data.txt"
    params:
        error_out_file="error_files/{sample_label}_insert_size_hist",
        run_time="20:00:00",
        cores="1",
        memory="6000",
        job_name="plot_insert_size"
    threads: 1
    shell: 
        "java -jar " + PICARD_JAR + " CollectInsertSizeMetrics I={input.bam} O={output.histogram_data} H={output.histogram_plot} W=1000 STOP_AFTER=50000000" 

ruleorder: plot_insert_size_hist_snake > plot_insert_size_hist_picard

# make big wigs to then use deeptools to make plots centered on features
rule make_coverage_bigwig:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx
    output:
        "output/coverage_data/qc/{sample_label}.100bp_coverage.bw"
    params:
        error_out_file="error_files/{sample_label}_deeptools_make_big_wig",
        run_time="20:00:00",
        cores="8",
        memory="6000",
        job_name="bam2bigwig"
    threads: 8
    benchmark: "benchmarks/make_coverage_bigwig/{sample_label}_deeptools_make_bigwig.txt"
    shell: 
        # --ignoreForNormalization chrX chrY taken off because we should remove this before calculating peaks anyway
        "bamCoverage -b {input.bam} -o {output} --binSize 1 --normalizeUsing CPM --extendReads --exactScaling --numberOfProcessors {threads}; "

rule make_sample_correlation_matrix:
    input: 
         expand("output/coverage_data/qc/{sample_label}.100bp_coverage.bw", sample_label=sample_labels)
    output:
        matrix = "output/coverage_data/qc/multiBamSummary_out.npz",
        read_counts = "output/coverage_data/qc/multiBamSummary_out_readCounts.tab"
    params:
        error_out_file="error_files/sample_correlation_matrix",
        run_time="23:59:59",
        cores="8",
        memory="16000",
        job_name="dt_correlation_mat"
    benchmark: "benchmarks/sample_correlation_matrix/v1.txt"
    threads: 8
    shell: 
        "multiBigwigSummary bins --bwfiles {input} --binSize 10000  --labels " + " ".join(["'" + label + "'" for label in sample_labels]) + " --outFileName {output.matrix} --outRawCounts {output.read_counts} --numberOfProcessors {threads}; "
        
rule plot_sample_correlation:
    input:
        rules.make_sample_correlation_matrix.output.matrix
    output:
        "output/plots/qc/sample_correlation_plot.pdf"
    params:
        error_out_file="error_files/sample_correlation_plot",
        run_time="00:30:00",
        cores="1",
        memory="6000",
        job_name="dt_correlation_plot"
    benchmark: "benchmarks/sample_correlation_plot/v1.txt"
    shell: 
        "plotCorrelation --corData {input} --colorMap RdYlBu --plotNumbers --plotFile {output} --corMethod spearman --whatToPlot heatmap"

rule make_insertion_bed:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx
    output:
        bed = "output/beds/{sample_label}.insertions.bed.gz"
    params:
        error_out_file="error_files/{sample_label}_bam2bed",
        run_time="20:00:00",
        cores="1",
        memory="6000",
        job_name="bam2bed"
    threads: 1
    benchmark: "benchmarks/make_bed/{sample_label}_bam2bed.txt"
    shell: 
        "bedtools bamtobed -i {input.bam} | awk -F $\'\\t\' 'BEGIN {{OFS = FS}} $6 == \"+\" {{$2 = $2 + 4; $3 = $2 + 1; print}} $6 == \"-\" {{$3 = $3 - 5; $2 = $3 - 1; print $0}}' | sort -k1,1 -k2,2n | gzip -c > {output.bed};"

# BEFORE MACS we should
# 1. filter high-quality aligning reads: rm_low_quality_reads
# 2. remove duplicate reads using Picard: rm_duplicates
# 3. remove mitochondrial reads: rm_mito (might be necessary for MACS2?)
# 4. remove read pairs that aren't connected, mate unmapped, or different chromosomes: done in step 1
# 5. remove things that aren't mapped uniquely: done in step 1

rule run_MACS2_bam:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx # do we need the bam index? better safe
         # rules.rm_blacklist # if remove from blacklist first
    output:
        narrowPeak = "output/peaks/{sample_label}_peaks.narrowPeak",
        broadPeak = "output/peaks/{sample_label}_peaks.broadPeak",
        gappedPeak = "output/peaks/{sample_label}_peaks.gappedPeak",
        peak_xls = "output/peaks/{sample_label}_peaks.xls",
        peak_bed = "output/peaks/{sample_label}_summits.bed",
        peak_treat = "output/peaks/{sample_label}_treat_pileup.bdg",
        peak_control = "output/peaks/{sample_label}_control_lambda.bdg",        
    params:
        error_out_file = "error_files/{sample_label}_MACS2_bam",
        run_time = "00:59:59",
        cores = "1",
        memory = "8000",
        job_name = "macs2"
    benchmark: "benchmarks/macs2/{sample_label}.bam.txt"
    shell:
        "/usr/local/anaconda/envs/py27/bin/" + "macs2 callpeak -g " + str(EFFECTIVE_GENOME_SIZE) + " --name {wildcards.sample_label} --treatment {input.bam} --outdir output/peaks --format BAMPE --nomodel --broad --nolambda --keep-dup all -p 0.01 -B --SPMR;"
        "/usr/local/anaconda/envs/py27/bin/" + "macs2 callpeak -g " + str(EFFECTIVE_GENOME_SIZE) + " --name {wildcards.sample_label} --treatment {input.bam} --outdir output/peaks --format BAMPE --nomodel --call-summits --nolambda --keep-dup all -p 0.01 -B --SPMR;"         

rule run_MACS2_bed:
    input:
        bed = rules.make_insertion_bed.output.bed,
    output:
        narrowPeak = "output/peaks/{sample_label}_peaks.narrowPeak",
        broadPeak = "output/peaks/{sample_label}_peaks.broadPeak",
        gappedPeak = "output/peaks/{sample_label}_peaks.gappedPeak",
        peak_xls = "output/peaks/{sample_label}_peaks.xls",
        peak_bed = "output/peaks/{sample_label}_summits.bed",
        peak_treat = "output/peaks/{sample_label}_treat_pileup.bdg",
        peak_control = "output/peaks/{sample_label}_control_lambda.bdg",
    params:
        error_out_file = "error_files/{sample_label}_MACS2_bed",
        run_time = "00:59:59",
        cores = "1",
        memory = "8000",
        job_name = "macs2"
    benchmark: "benchmarks/macs2/{sample_label}.bed.txt"
    shell:
        "/usr/local/anaconda/envs/py27/bin/" + "macs2 callpeak -g " + str(EFFECTIVE_GENOME_SIZE) + " --name {wildcards.sample_label} --treatment {input.bed} --outdir output/peaks --format BED --shift -75 --extsize 150 --nomodel --broad --nolambda --keep-dup all -p 0.01 -B --SPMR;"
        "/usr/local/anaconda/envs/py27/bin/" + "macs2 callpeak -g " + str(EFFECTIVE_GENOME_SIZE) + " --name {wildcards.sample_label} --treatment {input.bed} --outdir output/peaks --format BED --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01 -B --SPMR;"

ruleorder: run_MACS2_bed > run_MACS2_bam

rule run_MACS2_all_sample_bed:
     input:
         beds = expand("output/beds/{sample_label}.insertions.bed.gz", sample_label=sample_labels)
     output:
         peak_file = "output/peaks/all_samples/all_sample_peaks.narrowPeak",
         broad_peak_file = "output/peaks/all_samples/all_sample_peaks.broadPeak",
         gapped_peak_file = "output/peaks/all_samples/all_sample_peaks.gappedPeak",
         peak_xls = "output/peaks/all_samples/all_sample_peaks.xls",
         peak_bed = "output/peaks/all_samples/all_sample_summits.bed"
     params:
         error_out_file = "error_files/all_sample_MACS2_bed",
         run_time = "00:59:59",
         cores = "1",
         memory = "8000",
         job_name = "macs2"
     benchmark: "benchmarks/macs2/all_sample.bed.txt"
     shell:
         "/usr/local/anaconda/envs/py27/bin/" + "macs2 callpeak --name all_sample --treatment <(zcat {input.beds}) --outdir output/peaks/all_samples --format BED --shift -75 --extsize 150 --nomodel --broad --nolambda --keep-dup all -p 0.01 -B --SPMR;"
         "/usr/local/anaconda/envs/py27/bin/" + "macs2 callpeak --name all_sample --treatment <(zcat {input.beds}) --outdir output/peaks/all_samples --format BED --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01 -B --SPMR;"
 
rule run_MACS2_bdgcmp:
    input:
        peak_treatment = "output/peaks/{sample_label}_treat_pileup.bdg",
        peak_control = "output/peaks/{sample_label}_control_lambda.bdg"
    output:
        peak_fe = "output/tracks/{sample_label}_FE.bdg",
        peak_loglr = "output/tracks/{sample_label}_logLR.bdg",
    benchmark: "benchmarks/output/peaks/{sample_label}.bdgcmp.txt"
    shell:
        "/usr/local/anaconda/envs/py27/bin/" + "macs2 bdgcmp -t {input.peak_treatment} -c {input.peak_control} -m FE --outdir output/tracks --o-prefix {wildcards.sample_label} -p 0.0001 > /dev/null;"
        "/usr/local/anaconda/envs/py27/bin/" + "macs2 bdgcmp -t {input.peak_treatment} -c {input.peak_control} -m logLR --outdir output/tracks --o-prefix {wildcards.sample_label} -p 0.0001 > /dev/null" 

rule convert_MACS2_signal_to_bigwig:
    input:
        peak_fe = "output/tracks/{sample_label}_FE.bdg",
        peak_loglr = "output/tracks/{sample_label}_logLR.bdg"
    output:
        peak_fe = "output/tracks/{sample_label}_FE.bw",
        peak_loglr = "output/tracks/{sample_label}_logLR.bw",
        peak_fe_filter = "output/tracks/{sample_label}_FE.filter.bdg",
        peak_loglr_filter = "output/tracks/{sample_label}_logLR.filter.bdg"
    shell: 
        "grep -Ff <(awk '{{print $1 \"\t\"}}' " + CHROM_SIZES + ") {input.peak_fe} | sort -k1,1 -k2,2n | bedtools slop -i stdin -g " + CHROM_SIZES + " -b 0 | bedClip stdin " + CHROM_SIZES + " {output.peak_fe_filter};"
        "bedGraphToBigWig {output.peak_fe_filter} " + CHROM_SIZES + " {output.peak_fe};"
        "grep -Ff (awk '{{print $1 \"\t\"}}' " + CHROM_SIZES + ") {input.peak_loglr} | sort -k1,1 -k2,2n | bedtools slop -i stdin -g " + CHROM_SIZES + " -b 0 | bedClip stdin " + CHROM_SIZES + " {output.peak_loglr_filter};"
        "bedGraphToBigWig {output.peak_loglr_filter} " + CHROM_SIZES + " {output.peak_loglr};"

rule make_vplots:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        bed = lambda wildcards: BEDS[wildcards.bed_data]
    output:
        "output/plots/{bed_data}/{sample_label}.{bed_data}.Vplot.eps",
        "output/plots/{bed_data}/{sample_label}.{bed_data}.InsertionProfile.eps",
        "output/plots/{bed_data}/{sample_label}.{bed_data}.InsertSizes.eps",
        "output/plots/{bed_data}/{sample_label}.{bed_data}.VMat"
    params:
        error_out_file="error_files/{sample_label}_vplot",
        run_time="00:30:00",
        cores="8",
        memory="6000",
        job_name="v_plot"
    benchmark: "benchmarks/vplots/{sample_label}.{bed_data}.txt"
    threads: 8
    shell: 
        "/usr/local/anaconda/envs/py27/bin/pyatac vplot --bed {input.bed} --bam {input.bam} --out output/plots/{wildcards.bed_data}/{wildcards.sample_label}.{wildcards.bed_data} --plot_extra --cores {threads}"

rule make_insertion_bw:
    input:
        bed = rules.make_insertion_bed.output.bed,
    output:
        bg = temp("output/beds/{sample_label}.insertions.bedGraph"),
        bw = "output/coverage_data/{sample_label}.insertion_track.bw"
    params:
        error_out_file="error_files/{sample_label}_insertion_bw",
        run_time="00:30:00",
        cores="1",
        memory="2000",
        job_name="insertion_bw"
    benchmark: "benchmarks/insertion_bw/{sample_label}.txt"
    threads: 1
    shell: 
        "zcat {input.bed} | awk 'BEGIN {{chr = 0; pos = 0 ; print \"track type=bedGraph\"}} NR == 1 {{ chr = $1; pos = $3; coverage = 1}} NR > 1 && chr == $1 && pos == $3 {{coverage += 1}} NR > 1 && (chr != $1 || pos != $3) {{print chr \"\\t\" pos - 1 \"\\t\" pos \"\\t\" coverage; chr = $1; pos = $3; coverage = 1}} END {{print chr \"\\t\" pos - 1 \"\\t\" pos \"\\t\" coverage}}' > output/beds/{wildcards.sample_label}.insertions.bedGraph; bedGraphToBigWig {output.bg} " + CHROM_SIZES + " {output.bw}" 
        
rule calc_insertion_matrix:
    input:
        bw = "output/coverage_data/{sample_label}.insertion_track.bw",
        bed = lambda wildcards: BEDS[wildcards.bed_data]
    output:
        "output/coverage_data/{bed_data}/{sample_label}.{bed_data}.insertion_matrix.txt.gz",
    threads: 8
    params:
        error_out_file="error_files/{sample_label}_vplot",
        run_time="00:30:00",
        cores="8",
        memory="6000",
        job_name="bw_insertion_matrix"
    benchmark: "benchmarks/vplots/{sample_label}.{bed_data}.txt"
    shell: 
        "computeMatrix reference-point -S {input.bw} -R {input.bed} -p {threads} -a 1000 -b 1000 -bs 10 --outFileName {output} --missingDataAsZero"

rule plot_insertion_py2:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        TSS = TSS
    output:
        "output/coverage_data/TSS/{sample_label}.TSS.insertion_profile.txt",
        "output/coverage_data/TSS/{sample_label}.TSS.insertion_matrix.txt.gz",
        "output/plots/TSS/{sample_label}.TSS.insertion_profile.png",
    threads: 4
    run:
        shell("/usr/local/anaconda/envs/py27/bin/" + "python " + ATAC_TOOLS + "/make_insertion_enrichment_plot.py --bam {input.bam} --bed  {input.TSS} --out_basename {wildcards.sample_label}.TSS --plt_dir output/plots/TSS --txt_dir output/coverage_data/TSS --sample_name {wildcards.sample_label} --cores {threads} --feature_name TSS")

rule calculate_all_sample_nuc_content:
     input:
         "output/peaks/all_samples/all_sample_peaks.narrowPeak"
     output:
         bed = "output/peaks/all_samples/all_sample_peaks.narrowPeak.bed",
         nuc = "output/peaks/all_samples/all_sample_peaks.narrowPeak.nuc.txt"
     benchmark: "benchmarks/macs2/nuc.all_samples.bed.txt"
     shell: 
         "awk '{{OFS=\"\t\"; print $1, $2, $3, $4, int($5)}}' {input} | bedtools slop -i stdin -g " + CHROM_SIZES + " -b 0 | bedClip stdin " + CHROM_SIZES + "  {output.bed}; "
         "bedtools nuc -fi " + REFERENCE_FILE + ".fa -bed {output.bed} | awk '(NR > 1)' > {output.nuc}"
 
rule calculate_reads_in_peaks:
    input:
        sample_bed = "output/beds/{sample_label}.insertions.bed.gz",
        all_sample_peaks = "output/peaks/all_samples/all_sample_peaks.narrowPeak"
    output:
        "output/beds/qc/in_peaks/{sample_label}.reads_in_peaks.txt"
    benchmark: "benchmarks/macs2/{sample_label}.reads_in_peaks.txt"
    shell: 
        "bedtools intersect -a {input.all_sample_peaks} -b {input.sample_bed} -c > {output}"
 
rule compile_all_sample_nuc_matrix:
    input:
        narrowpeak = "output/peaks/all_samples/all_sample_peaks.narrowPeak",
        count_files = expand("output/beds/qc/in_peaks/{sample_label}.reads_in_peaks.txt", sample_label = sample_labels),
        nuc_content = rules.calculate_all_sample_nuc_content.output.nuc
    output:
        "output/beds/qc/all_samples.nuc_matrix.txt"
    benchmark: "benchmarks/beds/all_samples.nuc_matrix.txt"
    shell: 
        "coverage_columns=$(ls {input.count_files}  | while read count_file; do echo -n '<(cut -f 11 ' $count_file ') '; done); "
        "echo chr start end name score strand signal p q peak at gc A C G T N O length {sample_labels} | tr ' ' '\t' > {output}; "
        "eval \"paste {input.narrowpeak} <(cut -f 6- {input.nuc_content}) $coverage_columns >> {output}\""

rule FRiP_on_narrow:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        narrowPeak = "output/peaks/{sample_label}_peaks.narrowPeak"
    output:
        saf = "output/counts/{sample_label}.saf",
        counts = "output/counts/{sample_label}_narrowPeaks.counts"
    shell: 
        "awk 'OFS=\"\t\" {{print $1\"-\"$2+1\"-\"$3\"\t\"$1\"\t\"$2+1\"\t\"$3, \"+\"}}' {input.narrowPeak} > {output.saf};"
        "featureCounts -p -a {output.saf} -F SAF -o {output.counts} {input.bam}"
