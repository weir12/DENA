import json
import yaml
from os.path import join, basename, dirname


configfile: 'config.yaml'
thread=config['threads']
ref=config['ref']
out_dir=config['output_basedir']

SAMPLE=["WT","KO"]
ID=[config["WT"]['id'],config["KO"]['id']]
#Rules ------------------------------
rule all:
    input:
        expand("logs/{sample}/tombo_annotate_raw_with_fastqs_{id}.log",zip,sample=SAMPLE,id=ID),
        expand("logs/{sample}/tombo_squiggle_{id}.log",zip,sample=SAMPLE,id=ID),
        'logs/Selection_of_candidate_sites.log',
        expand(join(out_dir,"{sample}_align_result","{id}.sort.bam"),zip,sample=SAMPLE,id=ID)



rule Selection_of_candidate_sites:
     input:
         differ_res=config['differ_res']
     output:
         join(out_dir,"candidate_sites_result","candidate_sites.txt")
     log:
         'logs/Selection_of_candidate_sites.log'
     shell:
         """
         bash ./step1_select_m6a_sites/find_motif.sh -i {input.differ_res} -o {output} -a {ref}
         """

rule tombo_annotate_raw_with_fastqs:
    input:
        lambda wildcards:config[wildcards.sample]['fastq']
    params:
        se_su=lambda wildcards:config[wildcards.sample]['sequenceing_summary'],
        fast5=lambda wildcards:config[wildcards.sample]['fast5']
    output:
        "logs/{sample}/tombo_annotate_raw_with_fastqs_{id}.log"
    shell:
        """
        tombo preprocess annotate_raw_with_fastqs --fast5-basedir {params.fast5} --fastq-filenames {input} \
        --sequencing-summary-filenames {params.se_su} --processes {thread} --overwrite  --basecall-group "Basecall_1D_001" \
        1>{output} 2>&1
        """

rule Sequence_alignment_of_FASTQ_files:
     input:
         fastq=lambda wildcards:config[wildcards.sample]['fastq']
     output:
         join(out_dir,"{sample}_align_result","{id}.sort.bam")
     shell:
         """
         cat {input.fastq} |minimap2 -t {thread}  -ax map-ont  --secondary=no {ref} -  |samtools sort -@ {thread} -T {wildcards.id} -O bam -o  {output} -  && samtools index -@ {thread} {output}
         """


rule tombo_squiggle:
    input:
        "logs/{sample}/tombo_annotate_raw_with_fastqs_{id}.log"
    params:
        ref={ref},
        fast5=lambda wildcards:config[wildcards.sample]['fast5']
    log:
        "logs/{sample}/tombo_squiggle_{id}.log"
    shell:
        """
        tombo resquiggle --processes {thread} --ignore-read-locks --max-scaling-iterations 5 \
        --rna --basecall-group "Basecall_1D_001" --num-most-common-errors 5 --include-event-stdev \
        --overwrite --signal-length-range 0 500000 {params.fast5} {params.ref} 1> 
        """

rule extract_alignment_features:
   input:
        s=join(out_dir,"candidate_sites_result","candidate_sites.txt"),
        b=join(out_dir,"{sample}_align_result","{id}.sort.bam")
    output:
    	directory(join(out_dir,"features"))
    shell:
    	"""
		python ./step2_feature_extraction/label_reads.py --sites_file {input.s}  --bam {input.b}--label {wildcards.sample}
		
    	"""
rule extract_signal_features:
	input:
		


