# SCIPhI: Single-cell mutation identification via phylogenetic inference
#
# Copyright (C) 2018 ETH Zurich, Jochen Singer
#
# This file is part of SCIPhI.
#
# SCIPhI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SCIPhI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
# 
# @author: Jochen Singer

def finalTSVNames():
    out = []
    for srr, sample in SRR2SAMPLE.items():
        out.append('sra/' + sample + '/PAIREDEND/' + srr + '.tsv')
    return out
FINALTSVNAMES = finalTSVNames()

rule all:
    input: 
        expand(OUTDIR + '{file}', file = FINALFASTQNAMES),
        expand(OUTDIR + '{file}', file = FINALTSVNAMES),
        DATABASEDIR + 'ucsc.hg19.fasta.amb',
        DATABASEDIR + 'ucsc.hg19.fasta',
        DATABASEDIR + 'ucsc.hg19.fasta.fai',
        DATABASEDIR + 'ucsc.hg19.dict',
	DATABASEDIR + 'Mills_and_1000G_gold_standard.indels.hg19.sites.vcf',
	DATABASEDIR + 'Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx'

def getShortId(wildcards):
    return wildcards.id[0:6]

rule download:
    output: temp(OUTDIR + 'sra/{id}.sra')
    params: 
        shortId = getShortId,
        lsfoutfile = OUTDIR + 'sra/{id}.sra.out',
        lsferrfile = OUTDIR + 'sra/{id}.sra.err',
        outdir = OUTDIR,
        scratch = '10000',
        mem = '10000',
        time = '15'
    shell: 'wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{params.shortId}/{wildcards.id}/{wildcards.id}.sra -O {params.outdir}sra/{wildcards.id}.sra'

rule fastqDumpPaired:
    input: 
        sra = OUTDIR + 'sra/{id}.sra'
    output: 
        R1 = OUTDIR + 'sra/{id}/PAIREDEND/{id}_1.fastq.gz',
        R2 = OUTDIR + 'sra/{id}/PAIREDEND/{id}_2.fastq.gz'
    params: 
        outDir = OUTDIR + 'sra/{id}/PAIREDEND/',
        lsfoutfile = OUTDIR + 'sra/{id}.fastq.out',
        lsferrfile = OUTDIR + 'sra/{id}.fastq.err',
        scratch = '10000',
        mem = '10000',
        time = '120'
    conda:
        'sra.yaml'
    shell:
        'fastq-dump --split-3 --gzip --outdir {params.outDir} {input.sra}'

rule fastqDumpSingle:
input:
    sra = OUTDIR + 'sra/{id}.sra'
output:
    OUTDIR + 'sra/{id}/SINGLEEND/{id}.fastq.gz',
params:
    outDir = OUTDIR + 'sra/{id}/SINGLEEND/',
    lsfoutfile = OUTDIR + 'sra/{id}.fastq.out',
    lsferrfile = OUTDIR + 'sra/{id}.fastq.err'
shell:
    'fastq-dump --split-3 --gzip --outdir {params.outDir} {input.sra}'

localrules: linkFastqsPaired
rule linkFastqsPaired:
    input:
        fastq = OUTDIR + 'sra/{srr}/PAIREDEND/{srr}_{mate}.fastq.gz'
    output:
        fastq = OUTDIR + 'sra/{sample}/PAIREDEND/{srr}_R{mate}.fastq.gz'
    params:
        outdir = OUTDIR
    shell:
        'cd {params.outdir}sra/{wildcards.sample}/PAIREDEND/; ln -s ../../{wildcards.srr}/PAIREDEND/{wildcards.srr}_{wildcards.mate}.fastq.gz {wildcards.srr}_R{wildcards.mate}.fastq.gz; cd - ; touch -h {output.fastq}'

localrules: linkFastqsSingle
rule linkFastqsSingle:
    input:
        fastq = OUTDIR + 'sra/{srr}/SINGLEEND/{srr}.fastq.gz'
    output:
        fastq = OUTDIR + 'sra/{sample}/SINGLEEND/{srr}.fastq.gz'
    shell:
        'ln -s {input.fastq} {output.fastq} && touch -h {output.fastq}'

localrules: createTSVsPaired
rule createTSVsPaired:
    input:
        fastq = OUTDIR + 'sra/{srr}/PAIREDEND/{srr}_1.fastq.gz'
    output:
        fastq = OUTDIR + 'sra/{sample}/PAIREDEND/{srr}.tsv'
    run:
        f = open(output.fastq, 'w')
        f.write("RUN_NAME_FOLDER\t" + wildcards.sample + "\n")
        f.write("LANE_NUMBER\t" + wildcards.sample + "\n")
        f.write("SAMPLE_CODE\t" + wildcards.sample + "\n")
        f.write("SAMPLE_TYPE\tILLUMINA\n")
        f.close()

localrules: createTSVsSingle
rule createTSVsSingle:
    input:
        fastq = OUTDIR + 'sra/{srr}/SINGLEEND/{srr}.fastq.gz'
    output:
        fastq = OUTDIR + 'sra/{sample}/SINGLEEND/{srr}.tsv'
    run:
        f = open(output.fastq, 'w')
        f.write("RUN_NAME_FOLDER\t" + wildcards.sample + "\n")
        f.write("LANE_NUMBER\t" + wildcards.sample + "\n")
        f.write("SAMPLE_CODE\t" + wildcards.sample + "\n")
        f.write("SAMPLE_TYPE\tILLUMINA\n")
        f.close()

rule downloadReference:
    output:
        fa = 'databases/ucsc.hg19.fasta',
        fai = 'databases/ucsc.hg19.fasta.fai',
        dict = 'databases/ucsc.hg19.dict',
	mills = 'databases/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf',
	millsIdx = 'databases/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx'
    params:
        lsfoutfile = 'databases/ucsc.fasta.hg19.lsfout.log',
        lsferrfile = 'databases/ucsc.fasta.hg19.lsferr.log',
        scratch = '2000',
        mem = '2000',
        time = '120'
    shell:
        ('mkdir -p databases; cd databases; ' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz; gunzip ucsc.hg19.fasta.gz; ' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz; gunzip ucsc.hg19.fasta.fai.gz; ' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz; gunzip ucsc.hg19.dict.gz;' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz; gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz; ' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz; gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz; ' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz; gunzip 1000G_phase1.indels.hg19.sites.vcf.gz; ' +
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz; gunzip 1000G_phase1.indels.hg19.sites.vcf.idx.gz;')

rule createBwaIndex:
    input:
        ref = 'databases/ucsc.hg19.fasta'
    output:
        index1 = 'databases/ucsc.hg19.fasta.amb',
        index2 = 'databases/ucsc.hg19.fasta.ann',
        index3 = 'databases/ucsc.hg19.fasta.bwt',
        index4 = 'databases/ucsc.hg19.fasta.pac',
        index5 = 'databases/ucsc.hg19.fasta.sa'
    params:
        lsfoutfile = 'databases/ucsc.hg19.fasta.bwa.lsfout.log',
        lsferrfile = 'databases/ucsc.hg19.fasta.bwa.lsferr.log',
        scratch = '10000',
        mem = '10000',
        time = '240',
        index = 'databases/ucsc.hg19.fasta'
    conda:
        'bwa.yaml'
    shell:
        'bwa index {input.ref}'
