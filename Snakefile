"""
Author: Y. Ahmed-Braimah
--- Snakemake workflow to process single-end RNA-Seq.
"""

import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Global config files: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']

# Full path to an uncompressed GTF file with known gene annotations.
GTF = config['GTF']

# Full path to a folder where output files will be created.
OUT_DIR = config['OUT_DIR']

# Samples and their corresponding filenames.
FILES = json.load(open(config['SAMPLES_JSON'])) # The "samples.json" file can be created with
SAMPLES = sorted(FILES.keys())                  # the "make_json_SE/PE_samples.py" script

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all: 
    input: 
        join(OUT_DIR, 'Cuffmerge', 'merged.gtf'),
        expand(join(OUT_DIR, 'fastQC_initial', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC_final', '{sample}' + '_trimmed_R1_fastqc.html'), sample = SAMPLES),
        join(OUT_DIR, 'MultiQC_fastQC', 'multiqc_report.html'),
        expand(join(OUT_DIR, '{sample}', 'Cuffquant', 'abundances.cxb'), sample = SAMPLES),
        join(OUT_DIR, 'eXpress', 'genes.TMM.EXPR.matrix'),
        join(OUT_DIR, 'Cuffnorm', 'expression_data', 'run.info')


## Rule to check quality of raw reads
rule fastqc_init:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1']
    output:
        join(OUT_DIR, 'fastQC_initial', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, '{sample}', 'fastQC_initial.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'fastQC_initial.benchmark.tsv')
    message: 
        """--- Checking read quality of sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC_initial')):
            os.makedirs(join(OUT_DIR, 'fastQC_initial'))
        shell('fastqc'
                ' -o ' + join(OUT_DIR, 'fastQC_initial') + 
                ' {input.r1} > {log} 2>&1')

## Rule to trim raw fastq files
rule trimming:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1']
    output:
        r1 = join(OUT_DIR, '{sample}', '{sample}_trimmed_R1.fastq')
    log:
        join(OUT_DIR, '{sample}', 'sickle.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'sickle.benchmark.tsv')
    message: 
        """--- Trimming raw reads of sample "{wildcards.sample}" """
    run:
        shell('sickle se' 
                ' -f {input.r1}'
                ' -t sanger'
                ' -o {output.r1}'
                ' -q 20'
                ' -l 30 > {log} 2>&1')

## Rule to check quality of trimmed reads
rule fastqc_final:
    input:
        r1 = rules.trimming.output.r1
    output:
        join(OUT_DIR, 'fastQC_final', '{sample}' + '_trimmed_R1_fastqc.html')
    log:
        join(OUT_DIR, '{sample}', 'fastQC_final.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'fastQC_final.benchmark.tsv')
    message: 
        """--- Checking trimmed read quality of sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC_final')):
            os.makedirs(join(OUT_DIR, 'fastQC_final'))
        shell('fastqc'
                ' -o ' + join(OUT_DIR, 'fastQC_final') + 
                ' {input.r1} > {log} 2>&1')

## Rule to generate bowtie2 genome index file
rule index:
    input:
        dna = DNA
    output:
        index = join(dirname(DNA), rstrip(DNA, '.fa') + '.rev.1.bt2'),
        bt2i = join(dirname(DNA), rstrip(DNA, '.fa') + '.ok')
    log:
        join(dirname(DNA), 'bt2.index.log')
    benchmark:
        join(dirname(DNA), 'bt2.index.benchmark.tsv')
    message: 
        """--- Building bowtie2 genome index """
    run:
        shell('samtools faidx {input.dna}')
        shell('bowtie2-build {input.dna} ' + join(dirname(DNA), rstrip(DNA, '.fa')) + ' > {log} 2>&1')
        shell('touch ' + join(dirname(DNA), rstrip(DNA, '.fa') + '.ok'))

## Rule for mapping reads to the genome with Tophat
rule tophat:
    input:
        reads = rules.trimming.output.r1,
        idx = rules.index.output.bt2i
    output: 
        bam = join(OUT_DIR, '{sample}', 'tophat_out', 'accepted_hits.bam')
    params: 
        gtf = GTF
    log:
        join(OUT_DIR, '{sample}', 'tophat_out', 'tophat.map.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'tophat_out', 'tophat.map.benchmark.tsv')
    message: 
        """--- Mapping sample "{wildcards.sample}" with Tophat."""
    run: 
          shell('tophat'                                     
                ' -o ' + join(OUT_DIR, '{wildcards.sample}', 'tophat_out') + '/'    
                ' -G {params.gtf}'                                 
                ' -p 2 ' + join(dirname(DNA), rstrip(DNA, '.fa')) +                           
                ' {input.reads} > {log} 2>&1')

## Rule to collate fastQC outputs with multiQC
rule multiQC_fastQC:
    input:
        expand(join(OUT_DIR, 'fastQC_initial', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC_final', '{sample}' + '_trimmed_R1_fastqc.html'), sample = SAMPLES)
    output:
        file = join(OUT_DIR, 'MultiQC_fastQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC_fastQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC_fastQC', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC_fastQC """
    run:
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC_fastQC') + ' ' +
                join(OUT_DIR, 'fastQC_initial') + ' ' +
                join(OUT_DIR, 'fastQC_final') + 
                ' > {log} 2>&1')

## Rule for assembling rtansfrags with Cufflinks
rule cufflinks:
    input: 
        bam = rules.tophat.output.bam
    output: 
        gtf = join(OUT_DIR, '{sample}', 'cufflinks_out', 'transcripts.gtf')
    params:  
        gtf=GTF
    log:
        join(OUT_DIR, '{sample}', 'cufflinks.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'cufflinks.benchmark.tsv')
    message: 
        """--- Assembling "{wildcards.sample}" transcripts with cufflinks."""
    run: 
        shell('cufflinks'
              ' -g {params.gtf}'
              ' -p 2'
              ' -o ' +  join(OUT_DIR, '{wildcards.sample}', 'cufflinks_out') + '/'
              ' {input.bam}'
              ' &> {log}')

## Rule for merging Cufflinks assembled transcripts
rule cuffmerge:
    input:
        asmblys = expand(join(OUT_DIR, '{sample}', 'cufflinks_out', 'transcripts.gtf'), sample = SAMPLES)
    output: 
        merged = join(OUT_DIR, 'Cuffmerge', 'merged.gtf')
    params: 
        gtf =GTF, 
        fa = DNA
    log:
        join(OUT_DIR, 'Cuffmerge', 'cuffmerge.log')
    benchmark:
        join(OUT_DIR, 'Cuffmerge', 'cuffmerge.benchmark.tsv')
    message: 
        "--- Comparing transcripts to the reference and outputting merged gtf file."
    run: 
        # generate the assemblies text file
        shell('ls -1 ' + join(OUT_DIR) + '/*/cufflinks_out/transcripts.gtf > ' + join(OUT_DIR, 'assemblies.txt'))
        # run cuffmerge
        shell('cuffmerge'
              ' -o ' + join(OUT_DIR, 'Cuffmerge') +
              ' -g {params.gtf}'
              ' --keep-tmp'
              ' -s {params.fa}' 
              ' -p 2 ' + join(OUT_DIR, 'assemblies.txt') + 
              ' &> {log}')

## Rule for quantifying abundance with Cuffquant
rule cuffquant:
    input:
        bam = rules.tophat.output.bam,
        gtf = rules.cuffmerge.output.merged,
        dna = DNA
    output: 
        join(OUT_DIR, '{sample}', 'Cuffquant', 'abundances.cxb')
    log:
        join(OUT_DIR, '{sample}', 'Cuffquant', 'cffqnt.map.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'Cuffquant', 'cffqnt.benchmark.tsv')
    message: 
        """--- Quantifying abundances with Cuffquant for sample "{wildcards.sample}"."""
    run: 
        # run cuffmerge
        shell('cuffquant'
              ' -o ' + join(OUT_DIR, '{wildcards.sample}', 'Cuffquant') +
              # ' -p 8'
              ' -b {input.dna}'
              ' -u' 
              ' {input.gtf}'
              ' {input.bam}'
              ' &> {log}')

## Rule for quantifying abundance with Cuffquant
rule cuffnorm:
    input:
        cxb = expand(join(OUT_DIR, '{sample}', 'Cuffquant', 'abundances.cxb'), sample=SAMPLES),
        gtf = rules.cuffmerge.output.merged
    output: 
        join(OUT_DIR, 'Cuffnorm', 'expression_data', 'run.info')
    log:
        join(OUT_DIR, 'Cuffnorm', 'cffnrm.map.log')
    benchmark:
        join(OUT_DIR, 'Cuffnorm', 'cffnrm.benchmark.tsv')
    message: 
        """--- Merge Cuffquant abundances with Cuffnorm."""
    run:
        # create sample sheet
        shell('echo -e "sample_name\tgroup" > ' + join(OUT_DIR, 'Cuffnorm', 'sample_sheet.txt'))
        shell('ls -1 ' + join(OUT_DIR, '*', 'Cuffquant', 'abundances.cxb') + ' > ' + join(OUT_DIR, 'Cuffnorm', 'cq_abndces.txt'))
        shell('cat ' + join(OUT_DIR, 'Cuffnorm', 'cq_abndces.txt') + ' | sed "s/\/Cuffquant.*//g" | sed "s/.*\///g" | paste -d"\t" ' + join(OUT_DIR, 'Cuffnorm', 'cq_abndces.txt') +' - >> ' + join(OUT_DIR, 'Cuffnorm', 'sample_sheet.txt'))
        # run cuffnorm
        shell('cuffnorm'
              ' --use-sample-sheet'
              ' -o ' + join(OUT_DIR, 'Cuffnorm', 'expression_data') +
              # ' -p 8'
              ' {input.gtf} ' + join(OUT_DIR, 'Cuffnorm', 'sample_sheet.txt') +
              ' &> {log}')


## Rule for extracting novel transcripts not present in the parent GTF file
rule select_novel_transcripts:
    input:
        gtf = rules.cuffmerge.output.merged
    output:
        nvl_transcripts = join(OUT_DIR, 'Cuffmerge', 'nvl_transcripts.gtf')
    message: 
        "--- Extracting novel isoforms."
    run: 
        shell('grep -v "=" {input.gtf} > {output.nvl_transcripts}')

## Rule for adding XLOC gene IDs as gene names to novel transcripts
rule add_gene_name_to_unknown:
    input: 
        nvl_transcripts = rules.select_novel_transcripts.output.nvl_transcripts
    output: 
        nvlGeneName = join(OUT_DIR, 'Cuffmerge', 'nvl_transcripts_gn.gtf')
    message: 
        "--- Adding gene names to novel transcripts."
    run:
        import re 
        fh_in = open(input[0], "r")
        fh_out = open(output[0], "w")
        for line in fh_in:
          line = line.rstrip("\n")
          if not re.search("gene_name", line):
            gene_id = re.match('.*gene_id "(.*?)"', line).group(1)
            fh_out.write(line + ' gene_name "' + gene_id + '";\n')

## Rule for merging novel transcripts with known ones
rule merge_novel_and_known:
    input: 
        novel = rules.add_gene_name_to_unknown.output.nvlGeneName, 
        known = GTF
    output: 
        merged_all = join(OUT_DIR, 'Cuffmerge', 'merged_all.gtf')
    message: 
        "--- Merging known and novel transcripts."
    run:  
        shell('cat {input.novel} {input.known} > {output.merged_all}')

# Rule for making a cDNA file from the new GTF file and the DNA, then prep index for bowtie2.
rule make_cdna:
    input:
        dna = DNA,
        gtf = rules.merge_novel_and_known.output.merged_all
    output:
        cdna = join(OUT_DIR, 'transcriptome', 'gffread_transcripts.fa'),
        geneTrans = join(OUT_DIR, 'transcriptome', 'gffread_transcripts.gene_trans_map'),
        bt2_trans_indx = join(OUT_DIR, 'transcriptome', 'gffread_transcripts.fa.bowtie2.ok')
    log:
        gffread = join(OUT_DIR, 'transcriptome', 'logs', 'gffread.log'),
        trans_bt2 = join(OUT_DIR, 'transcriptome', 'logs', 'trans_bt2.log')
    benchmark:
        join(OUT_DIR, 'transcriptome', 'logs', 'gffread.benchmark.tsv')
    message: 
        "--- Building bowtie2 transcriptome index for gffread transcripts."
    run:
        # Extract a sequence for each transcript in the GTF file.
        shell('gffread -F -w {output.cdna} -g {input.dna} {input.gtf} > {log.gffread}')
        # Extract the FASTA header from the cDNA file and make into
        # trans_map file.
        shell('grep ">" {output.cdna} | sed "s/>//g" | sed "s/gene=//g" | awk \'{{print $2"\t"$1}}\' | sort -u > {output.geneTrans}')
        # And finally make the index files.
        shell('align_and_estimate_abundance.pl' 
              ' --transcripts {output.cdna}'
              ' --gene_trans_map {output.geneTrans}'
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --prep_reference'
              ' --output_dir ' + join(OUT_DIR, 'transcriptome') +
              ' > {log.trans_bt2} 2>&1')

# Rule for mapping reads to the new transcriptome file with bowtie2 and quantifying abundance with eXpress
rule express:
    input:
        reads = lambda wildcards: FILES[wildcards.sample]['R1'],
        cdna = rules.make_cdna.output.cdna,
        geneTrans = rules.make_cdna.output.geneTrans
    output:
        results = join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs')
    log:
        join(OUT_DIR, 'eXpress', '{sample}', 'logs', 'eXpress.log')
    benchmark:
        join(OUT_DIR, 'eXpress', '{sample}', 'logs', 'eXpress.benchmark.tsv')
    message: 
        """--- Mapping "{wildcards.sample}" reads to transcriptome with bowtie2 and quantifying abundance with eXpress."""
    run:
        shell('align_and_estimate_abundance.pl' 
              ' --transcripts {input.cdna}'
              ' --seqType fq'
              ' --single {input.reads}'
              ' --gene_trans_map {input.geneTrans}'
              ' --thread_count 8'  
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --output_dir ' + join(OUT_DIR, 'eXpress', '{wildcards.sample}') +  
              ' > {log} 2>&1')


rule merge_abundance:
    input:
        quants = expand(join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs'), sample = SAMPLES)
    output:
        abundances = join(OUT_DIR, 'eXpress', 'genes.TMM.EXPR.matrix'),
        samplesList = join(OUT_DIR, 'eXpress', 'genes.samples.list')
    log:
        join(OUT_DIR, 'eXpress', 'abnd_merge.log')
    benchmark:
        join(OUT_DIR, 'eXpress', 'abnd_merge.benchmark.tsv')        
    message: 
        "--- Merging eXpress outputs from all samples"
    run:
        shell('ls -1 ' + join(OUT_DIR, 'eXpress', '*', 'results.xprs.genes') + ' > ' + join(OUT_DIR, 'eXpress', 'genes.samples.list'))
        shell('ls -1 ' + join(OUT_DIR, 'eXpress', '*', 'results.xprs') + ' > ' + join(OUT_DIR, 'eXpress', 'isoforms.samples.list'))
        shell('cd ' + join(OUT_DIR, 'eXpress') +
                ' && abundance_estimates_to_matrix.pl'
                ' --est_method eXpress'
                ' --name_sample_by_basedir'
                ' --out_prefix genes'
                ' ' + join(OUT_DIR, 'eXpress', 'genes.samples.list') +
                ' > {log} 2>&1')
        shell('cd ' + join(OUT_DIR, 'eXpress') +
                ' && abundance_estimates_to_matrix.pl'
                ' --est_method eXpress'
                ' --name_sample_by_basedir'
                ' --out_prefix isoforms'
                ' ' + join(OUT_DIR, 'eXpress', 'isoforms.samples.list') +
                ' > {log} 2>&1')