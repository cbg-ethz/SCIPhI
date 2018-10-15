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

#These rules and functions were copied from https://github.com/cbg-ethz/NGS-pipe/

import ntpath

REGIONSOUT = OUTDIR + 'regions/'
localrules: getRegionsOfChromosome
rule getRegionsOfChromosome:
    input:
        regions = config['resources'][ORGANISM]['regions']
    output:
        chrRegion = REGIONSOUT + config['resources'][ORGANISM]['regions'].strip().split("/")[-1].replace('.bed','') + '_{chrom}.bed'
    shell:
        'awk \'$1 == \"{wildcards.chrom}\"\' {input.regions} > {output.chrRegion}'

# This function returns the TSV file connected to a given FASTQ file
def getTSV(wildcards):
    localFastq =''
    if wildcards.fastq[-3:] == '_R1' or wildcards.fastq[-3:] == '-R1' or wildcards.fastq[-3:] == '_R2' or wildcards.fastq[-3:] == '-R2':
        localFastq = wildcards.fastq[0:-3]
    else:
        localFastq = wildcards.fastq
    return FASTQDIR + wildcards.sample + '/' + wildcards.type.replace('ORPAHN/','') + '/' + localFastq + '.tsv'

# This function extracts the information to create the read-group identifier
# given to the read mappers
def createReadGroup(wildcards):
    out = [] # list containing flowcellID, lane, library, platform in that order
    localFastq =''
    if wildcards.fastq[-3:] == '_R1' or wildcards.fastq[-3:] == '-R1' or wildcards.fastq[-3:] == '_R2' or wildcards.fastq[-3:] == '-R2':
        localFastq = wildcards.fastq[0:-3]
    else:
        localFastq = wildcards.fastq
    tsv = open(FASTQDIR + wildcards.sample + '/' + wildcards.type.replace('ORPAHN/','') + '/' + localFastq + '.tsv')
    flowcellID = ''
    lane = ''
    library = ''
    platform = ''
    for line in tsv:
        if line.strip().startswith('RUN_NAME_FOLDER'):
            if flowcellID == '':
                flowcellID = line.strip().split('\t')[1]
        elif line.strip().startswith('LANE_NUMBER'):
            if lane == '':
                lane = line.strip().split('\t')[1]
        elif line.strip().startswith('SAMPLE_CODE'):
            if library == '':
                library = line.strip().split('\t')[1]
        elif line.strip().startswith('SAMPLE_TYPE'):
            if platform == '':
                platform = line.strip().split('\t')[1]
    tsv.close()
    out.append(flowcellID)
    out.append(lane)
    out.append(library)
    out.append(platform)
    return out

# This function creates the read-group entry for bowtie2
def createReadGroupBowtie2(wildcards):
    values = createReadGroup(wildcards)
    return '--rg-id ' + wildcards.sample + '.' + values[0] + '.' + values[1] + ' --rg LB:' + values[2] + ' --rg PL:' + values[3] + ' --rg PU:' + values[0] + '.' + values[1] + '.' + values[2] + ' --rg SM:' + wildcards.sample

# This function creates the read-group entry for yara
def createReadGroupYara(wildcards):
    values = createReadGroup(wildcards)
    return '\'@RG\\tID:' + wildcards.sample + '.' + values[0] + '.' + values[1] + '\\tLB:' + values[2] + '\\tPL:' + values[3] + '\\tPU:' + values[0] + '.' + values[1] + '.' + values[2]  + '\\tSM:' + wildcards.sample + '\''

# This function creates the read-group entry for bwa
def createReadGroupBwa(wildcards):
    values = createReadGroup(wildcards)
    return '\'@RG\\tID:' + wildcards.sample + '.' + values[0] + '.' + values[1] + '\\tLB:' + values[2] + '\\tPL:' + values[3] + '\\tPU:' + values[0] + '.' + values[1] + '.' + values[2]  + '\\tSM:' + wildcards.sample + '\''

# This rule aligns unpaired reads using bwa-mem
rule bwa_mem_paired:
    input:
        fastqR1 = BWAIN + '{sample}/{type}/{fastq}_R1.fastq.gz',
        fastqR2 = BWAIN + '{sample}/{type}/{fastq}_R2.fastq.gz',
        index1 = config['resources'][ORGANISM]['bwaIndex'] + '.amb',
        index2 = config['resources'][ORGANISM]['bwaIndex'] + '.ann',
        index3 = config['resources'][ORGANISM]['bwaIndex'] + '.bwt',
        index4 = config['resources'][ORGANISM]['bwaIndex'] + '.pac',
        index5 = config['resources'][ORGANISM]['bwaIndex'] + '.sa',
        tsv = getTSV
    output:
        bam=temp(BWAOUT + '{sample}/{type}/{fastq}.bam')
    params:
        lsfoutfile = BWAOUT + '{sample}/{type}/{fastq}.bam.lsfout.log',
        lsferrfile = BWAOUT + '{sample}/{type}/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bwa']['mem']['scratch'],
        mem = config['tools']['bwa']['mem']['memory'],
        time = config['tools']['bwa']['mem']['time'],
        index = config['resources'][ORGANISM]['bwaIndex'],
        params = config['tools']['bwa']['mem']['params'],
        rg = createReadGroupBwa
    conda:
        'envs/bwa.yaml'
    benchmark:
        BWAOUT + '{sample}/{type}/{fastq}.bam.benchmark'
    threads:
        config['tools']['bwa']['mem']['threads']
    log:
        BWAOUT + '{sample}/{type}/{fastq}.bam.log'
    shell:
        ('bwa mem ' +
        '{params.params} ' +
        '-R {params.rg} ' +
        '-t {threads} ' +
        '{params.index} ' +
        '{input.fastqR1} ' +
        '{input.fastqR2} ' +
        '2>{log}| samtools view -bhS - > {output.bam}')

rule bwa_mem_single:
    input:
        fastq = BWAIN + '{sample}/{type}/{fastq}.fastq.gz',
        index1 = config['resources'][ORGANISM]['bwaIndex'] + '.amb',
        index2 = config['resources'][ORGANISM]['bwaIndex'] + '.ann',
        index3 = config['resources'][ORGANISM]['bwaIndex'] + '.bwt',
        index4 = config['resources'][ORGANISM]['bwaIndex'] + '.pac',
        index5 = config['resources'][ORGANISM]['bwaIndex'] + '.sa',
        tsv = getTSV
    output:
        bam=temp(BWAOUT + '{sample}/{type}/{fastq}.bam')
    params:
        lsfoutfile = BWAOUT + '{sample}/{type}/{fastq}.bam.lsfout.log',
        lsferrfile = BWAOUT + '{sample}/{type}/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bwa']['mem']['scratch'],
        mem = config['tools']['bwa']['mem']['memory'],
        time = config['tools']['bwa']['mem']['time'],
        index = config['resources'][ORGANISM]['bwaIndex'],
        params = config['tools']['bwa']['mem']['params'],
        rg = createReadGroupBwa
    conda:
       'envs/bwa.yaml'
    benchmark:
        BWAOUT + '{sample}/{type}/{fastq}.bam.benchmark'
    threads:
        config['tools']['bwa']['mem']['threads']
    log:
        BWAOUT + '{sample}/{type}/{fastq}.bam.log'
    shell:
        ('bwa mem ' +
        '{params.params} ' +
        '-R {params.rg} ' +
        '-t {threads} ' +
        '{params.index} ' +
        '{input.fastq} ' + 
        '2>{log} | samtools view -bhS - > {output.bam}')

# This rule creates an indes of a BAM file
rule samtools_create_index:
    input:
        bam = '{sample}.bam',
    output:
        idx = '{sample}.bai',
    params:
        lsfoutfile = '{sample}.bai.lsfout.log',
        lsferrfile = '{sample}.bai.lsferr.log',
        scratch = config['tools']['samtools']['index']['scratch'],
        mem = config['tools']['samtools']['index']['mem'],
        time = config['tools']['samtools']['index']['time']
    conda:
        'envs/samtools.yaml'
    threads:
        config['tools']['samtools']['index']['threads']
    benchmark:
        '{sample}.bai.benchmark'
    shell:
        'samtools index {input.bam} && mv {input.bam}.bai {output.idx}'

# This rule create an symbolic link to an existing index
rule linkIndex:
    input:
        bai = '{sample}.bai',
    output:
        bai = '{sample}.bam.bai',
    params:
        lsfoutfile = '{sample}.bam.bai.lsfout.log',
        lsferrfile = '{sample}.bam.bai.lsferr.log',
        scratch = '1000', 
        mem = '1000', 
        time = '1' 
    threads:
        1
    benchmark:
        '{sample}.bam.bai.benchmark'
    shell:
        'dirName=$(dirname "{input.bai}"); inBai=$(basename "{input.bai}"); outBai=$(basename "{output.bai}"); cd "$dirName"; ln -s "$inBai" "$outBai"'

# This rule sorts a BAM file and fixes mate pair information if necessary
rule picards_fix_mate_pair_and_sort:
    input:
        bam=FIXMATEANDSORTIN + '{sample}.bam'
    output:
        bam=temp(FIXMATEANDSORTOUT + '{sample}.bam')
    params:
        lsfoutfile = FIXMATEANDSORTOUT + '{sample}.bam.lsfout.log',
        lsferrfile = FIXMATEANDSORTOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['picard']['fixMateInformation']['scratch'],
        mem = config['tools']['picard']['fixMateInformation']['mem'],
        time = config['tools']['picard']['fixMateInformation']['time'],
        sortOrder = config['tools']['picard']['fixMateInformation']['sortOrder'],
        assume_sorted = config['tools']['picard']['fixMateInformation']['assume_sorted'],
        params = config['tools']['picard']['fixMateInformation']['params']
    conda:
        'envs/picard.yaml'
    benchmark:
        FIXMATEANDSORTOUT + '{sample}.bam.benchmark'
    threads:
        config['tools']['picard']['fixMateInformation']['threads']
    log:
        FIXMATEANDSORTOUT + '{sample}.bam.log'
    shell:
        ('picard FixMateInformation ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam} ' +
        'SORT_ORDER={params.sortOrder} ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        '{params.params} ' +
        'TMP_DIR={TMPDIR} ' +
        '2> {log}')

# This functiom creates a list of BAM files created by the read mapper
def getAlignerBams():
    out = []
    if config['tools']['picard']['mergeBams']['useOrphans'] != "Y":
        for f in PAIREDFASTQFILESWITHOUTR:
            out.append( f + '.bam')
        for f in SINGLEFASTQFILES:
            out.append( f + '.bam')
    else:
        for f in PAIREDFASTQFILES:
            out.append(os.path.dirname(f) +'/ORPHAN/' + ntpath.basename(f) + '.bam')
        for f in PAIREDFASTQFILESWITHOUTR:
            out.append( f + '.bam')
    return out

# This function is a helper function to get the BAM file names which are then merged.
def getBamsToMerge(wildcards):
    out = []
    allBams = getAlignerBams()
    for bam in allBams:
        if wildcards.sample == bam.split("/")[0]: 
            out.append(MERGEBAMSIN + bam)
    if not out:
        return ['ERROR']
    return out

# This function prepends 'INPUT=' in front of every BAM file that is to be merged.
def prependBamsToMerge(wildcards):
    bamsToMerge = getBamsToMerge(wildcards)
    return ''.join(['INPUT='+bam+' ' for bam in bamsToMerge])

# This rule merges different BAM files.
rule picard_merge_bams:
    input:
        bams = getBamsToMerge
    output:
        bam = temp(MERGEBAMSOUT + '{sample}.bam')
    params:
        lsfoutfile = MERGEBAMSOUT + '{sample}.bam.lsfout.log',
        lsferrfile = MERGEBAMSOUT + '{sample}.bam.lsferr.log',
        assume_sorted = config['tools']['picard']['mergeBams']['assume_sorted'],
        scratch = config['tools']['picard']['mergeBams']['scratch'],
        mem = config['tools']['picard']['mergeBams']['mem'],
        time = config['tools']['picard']['mergeBams']['time'],
        params = config['tools']['picard']['mergeBams']['params'],
        input = prependBamsToMerge
    conda:
        'envs/picard.yaml'
    threads:
        config['tools']['picard']['mergeBams']['threads']
    benchmark:
        MERGEBAMSOUT + '{sample}.bam.benchmark'
    log:
        log = MERGEBAMSOUT + '{sample}.bam.log',
        metrics = MERGEBAMSOUT + '{sample}.bam.metrics'
    shell:
        ('picard ' +
        'MergeSamFiles ' +
        '{params.input} ' +
        'OUTPUT={output.bam} ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        '{params.params} ' +
        '2> {log.log}')

# This rule markes PCR duplicates
if not 'MARKPCRDUBLICATESIN' in globals():
    MARKPCRDUBLICATESIN = MERGEBAMSOUT
if not 'MARKPCRDUBLICATESOUT' in globals():
    MARKPCRDUBLICATESOUT = OUTDIR + 'markedDuplicates/'
rule picards_mark_PCR_duplicates:
    input:
        bam=MARKPCRDUBLICATESIN + '{sample}.bam',
    output:
        bam=temp(MARKPCRDUBLICATESOUT + '{sample}.bam'),
    params:
        lsfoutfile = MARKPCRDUBLICATESOUT + '{sample}.bam.lsfout.log',
        lsferrfile = MARKPCRDUBLICATESOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['picard']['markduplicates']['scratch'],
        mem = config['tools']['picard']['markduplicates']['mem'],
        time = config['tools']['picard']['markduplicates']['time'],
        params = config['tools']['picard']['markduplicates']['params']
    conda:
        'envs/picard.yaml'
    threads:
        config['tools']['picard']['markduplicates']['threads']
    benchmark:
        MARKPCRDUBLICATESOUT + '{sample}.bam.benchmark'
    log:
        log = MARKPCRDUBLICATESOUT + '{sample}.bam.log',
        metrics = MARKPCRDUBLICATESOUT + '{sample}.bam.metrics'
    shell:
        ('picard ' +
        'MarkDuplicates ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam}  ' +
        '{params.params} ' +
        'TMP_DIR={TMPDIR} ' +
        #'--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ' +
        'METRICS_FILE={log.metrics} ' +
        '2> {log.log}')

def getSamplesFromExperimentId(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    expMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == 'Y':
                    if exp not in expMap.keys():
                        expMap[exp] = []
                    expMap[exp].append(sample)
                elif config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "N":
                    if sample in expMap.keys():
                        raise ValueError(sample = " is not uniq in the sample mapping file.")
                    expMap[sample] = []
                    expMap[sample].append(sample)
                else:
                    return "Unknown parameter " + config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] + " to specify whether all bams of one experiment should be realiged together."

    if config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "Y":
        if wildcards.experiment not in expMap.keys():
            return "UnknownExperiment"
    elif config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "N":
        if wildcards.experiment not in expMap.keys():
            #raise ValueError(wildcards.experiment + " is not a valid sample name!")
            return "UnknownSample"
    return expMap[wildcards.experiment]

def getBamsFromExperimentId(wildcards):
    return expand('{sample}.bam', sample = getSamplesFromExperimentId(wildcards))

def getBaisFromExperimentId(wildcards):
    return expand('{sample}.bai', sample = getSamplesFromExperimentId(wildcards))

def getBamsToRealingFromExperimentId(wildcards):
    return expand(REALIGNINDELSIN + '{bam}', bam = getBamsFromExperimentId(wildcards))

def getBaisToRealingFromExperimentId(wildcards):
    return expand(REALIGNINDELSIN + '{bai}', bai = getBaisFromExperimentId(wildcards))

def prependBamsToRealign(wildcards):
    bamsToRealign = getBamsToRealingFromExperimentId(wildcards)
    return ''.join(['-I '+bam+' ' for bam in bamsToRealign])

def getDataBasisForRealign():
    out = []
    out.append(config['resources'][ORGANISM]['reference']) # this is a dummy such that something is retured
    if config['tools']['GATK']['realign']['Mills_indels'] == "Y":
        out.append(config['resources'][ORGANISM]['Mills_indels'])
    if config['tools']['GATK']['realign']['1000G_indels'] == "Y":
        out.append(config['resources'][ORGANISM]['1000G_indels'])
    return out

def prependDataBasisForRealignTargetCreator():
    out = ""
    if config['tools']['GATK']['realign']['Mills_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['Mills_indels'], Error):
            print("You have not specified config[resources][ORGANISM][Mills_indels]")
        out += " --known " + config['resources'][ORGANISM]['Mills_indels']
    if config['tools']['GATK']['realign']['1000G_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['1000G_indels'], Error):
            print("You have not specified config[resources][ORGANISM][1000G_indels]")
        out += " --known " + config['resources'][ORGANISM]['1000G_indels']
    return out

# Rule to create a file containing the regions to perfrom the indel realignment on
if not 'REALIGNINDELSIN' in globals():
    REALIGNINDELSIN = REMOVEPCRDUBLICATESOUT
if not 'REALIGNINDELSOUT' in globals():
    REALIGNINDELSOUT = OUTDIR + 'realignedIndels/'
rule gatk_realign_target_creation:
    input:
        bam = getBamsToRealingFromExperimentId,
        bai = getBaisToRealingFromExperimentId,
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions'],
        databasis = getDataBasisForRealign()
    output:
        intervals = temp(REALIGNINDELSOUT + '{experiment}.intervals'),
    params:
        lsfoutfile = REALIGNINDELSOUT + '{experiment}.intervals.lsfout.log',
        lsferrfile = REALIGNINDELSOUT + '{experiment}.intervals.lsferr.log',
        scratch = config['tools']['GATK']['realign']['targetCreator']['scratch'],
        mem = config['tools']['GATK']['realign']['targetCreator']['mem'],
        time = config['tools']['GATK']['realign']['targetCreator']['time'],
        params = config['tools']['GATK']['realign']['targetCreator']['params'],
        input = prependBamsToRealign,
        known = prependDataBasisForRealignTargetCreator(),
    benchmark:
        REALIGNINDELSOUT + '{experiment}.intervals.benchmark'
    threads:
        config['tools']['GATK']['realign']['targetCreator']['threads']
    shell:
        ('{config[tools][GATK][call]} -T RealignerTargetCreator ' +
        '-R {input.reference} ' +
        '-L {input.regions} ' +
        '{params.input} ' +
        '{params.known} ' +
        '-o {output.intervals} ' +
        '-nt {threads} ' +
        '{params.params}')

localrules: createRealingIndelsInOutMapping
rule createRealingIndelsInOutMapping:
    output:
        exp = REALIGNINDELSOUT + '{experiment}.map'
    run:
        import sys
        try:
            outFile = open(output.exp, 'x')
            bams = getBamsFromExperimentId(wildcards)
            for bam in bams:
                outFile.write(bam + "\t" + REALIGNINDELSOUT + "ORIGINAL_" + bam + "\n")
        except IOError:
            print("Could not open file: ", output.exp)

def prependDataBasisForTargetRealigner():
    out = ""
    if config['tools']['GATK']['realign']['Mills_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['Mills_indels'], Error):
            print("You have not specified config[resources][ORGANISM][Mills_indels]")
        out += " -known " + config['resources'][ORGANISM]['Mills_indels']
    if config['tools']['GATK']['realign']['1000G_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['1000G_indels'], Error):
            print("You have not specified config[resources][ORGANISM][1000G_indels]")
        out += " -known " + config['resources'][ORGANISM]['1000G_indels']
    return out

# Rule to perform the indel realignment
# This is a GATK tool
#ruleorder: realignIndels > createIndex
rule gatk_realign_indels:
    input:
        bam = getBamsToRealingFromExperimentId,
        bai = getBaisToRealingFromExperimentId,
        map = REALIGNINDELSOUT + '{experiment}.map',
        regions = config['resources'][ORGANISM]['regions'],
        intervals = REALIGNINDELSOUT + '{experiment}.intervals',
        reference = config['resources'][ORGANISM]['reference'],
        databasis = getDataBasisForRealign()
    output:
        txt = REALIGNINDELSOUT + '{experiment}.realigned.txt',
    params:
        lsfoutfile = REALIGNINDELSOUT + '{experiment}.realigned.txt.lsfout.log',
        lsferrfile = REALIGNINDELSOUT + '{experiment}.realigned.txt.lsferr.log',
        scratch = config['tools']['GATK']['realign']['realignIndels']['scratch'],
        mem = config['tools']['GATK']['realign']['realignIndels']['mem'],
        time = config['tools']['GATK']['realign']['realignIndels']['time'],
        params = config['tools']['GATK']['realign']['realignIndels']['params'],
        input = prependBamsToRealign,
        known = prependDataBasisForTargetRealigner(),
    benchmark:
        REALIGNINDELSOUT + '{experiment}.realigned.txt.benchmark'
    threads:
        config['tools']['GATK']['realign']['realignIndels']['threads']
    shell:
        ('{config[tools][GATK][call]} -T IndelRealigner ' +
        '-R {input.reference} ' +
        '-L {input.regions} ' +
        '{params.input} ' +
        '{params.known} ' +
        '--nWayOut {input.map} ' +
        '-targetIntervals {input.intervals} ' +
        '{params.params} ' +
        '&& touch {output.txt}')

def getExperimentIdFromBam(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "Y":
                    sampleMap[sample] = exp
                elif config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "N":
                    sampleMap[sample] = sample
                else:
                    return "Unknown parameter " + config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] + " to specify whether all bams of one experiment should be realiged together."
    
    if wildcards.sample not in sampleMap.keys():
        #raise ValueError(wildcards.sample + " is not a sample ID")
        return "UnknownSample"
    return sampleMap[wildcards.sample]

def getRealignedExperiment(wildcards):
    return expand(REALIGNINDELSOUT + '{experiment}.realigned.txt', experiment = getExperimentIdFromBam(wildcards))

# This is a "dummy" rule
localrules: getRealignedBam 
rule getRealignedBam:
    input:
        txt = getRealignedExperiment
    output:
        bam = REALIGNINDELSOUT + '{sample}.bam'
    params:
        dirName = REALIGNINDELSOUT,
        originalBam = REALIGNINDELSOUT + 'ORIGINAL_{sample}.bam'
    shell:
        'cd {params.dirName}; ln -s ORIGINAL_{wildcards.sample}.bam {wildcards.sample}.bam'

rule mpileupMpileup:
    input:
        bam = MPILEUPIN + '{sample}.bam',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        mpileup = temp(MPILEUPOUT + '{sample}.mpileup')
    params:
        lsfoutfile = MPILEUPOUT + '{sample}.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{sample}.mpileup.lsferr.log',
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
        params = config['tools']['samtools']['mpileup']['params']
    conda:
        'envs/samtools.yaml'
    threads:
        config['tools']['samtools']['mpileup']['threads']
    benchmark:
        MPILEUPOUT + '{sample}.mpileup.benchmark'
    log:
        MPILEUPOUT + '{sample}.mpileup.log'
    shell:
        ('samtools mpileup ' +
        '{params.params} ' + 
        '-f {input.reference} ' + 
        '-o {output.mpileup} ' + 
        '-l {input.regions} ' +
        '{input.bam}')

localrules: sciphiToHeatMap
rule sciphiToHeatMap:
    input:
        probs = '{sciphi}.probs',
        gv = '{sciphi}.gv'
    output:
        tsv = '{sciphi}_heatTree.tsv',
        pdf = '{sciphi}_heatTree.pdf'
    shell:
        '{config[tools][sciphiToHeatMap][call1]} -i {input.gv} -p {input.probs} -o {output.tsv}; {config[tools][sciphiToHeatMap][call2]} {output.tsv} {output.pdf}'

localrules: monovarToHeatMap
rule monovarToHeatMap:
    input:
        tsv = '{dir}/sciphi/1/{sample}_heatTree.tsv',
        vcf = '{dir}/monovar/{sample}.vcf',
        gv = '{dir}/sciphi/1/{sample}.gv'
    output:
        tsv = '{dir}/monovar/{sample}.probs',
        pdf = '{dir}/monovar/{sample}_heatTree.pdf',
    shell:
        '{config[tools][monovarToHeatMap][call1]} {input.tsv} {input.vcf} {output.tsv}; {config[tools][monovarToHeatMap][call2]} {output.tsv} {output.pdf}'

localrules: sciphiOverviewGraph
rule sciphiOverviewGraph:
    input:
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_all_complete.mpileup',
        fileNames = SCIPHIOUT + '{experiment}_all_bamFileNames.txt',
        bTree = SCIPHIOUT + '{run}/{experiment}/best_index/tree.gv',
        sTree = SCIPHIOUT + '{run}/{experiment}.gv'
    output:
        smt = SCIPHIOUT + '{run}/{experiment}_smt.tsv',
        gv = SCIPHIOUT + '{run}/{experiment}_overview.gv',
        pdf = SCIPHIOUT + '{run}/{experiment}_overview.pdf'
    params:
        index = SCIPHIOUT + '{run}/{experiment}/best_index/'
    shell:
        ('{config[tools][sciphi][call]} --smt {output.smt} --il {params.index} {input.mpileup}; ' +
        '{config[tools][sciphiOverviewGraph][call1]} {output.smt} {input.bTree} {input.fileNames} {input.sTree} > {output.gv}; {config[tools][sciphiOverviewGraph][call2]} -Tpdf {output.gv} > {output.pdf}')

localrules: monovarCluster
rule monovarCluster:
    input:
        probs = '{sample}.probs'
    output:
        probs = '{sample}_cluster.probs',
        pdf = '{sample}_cluster.pdf',
    shell:
        ('cp {input.probs} {output.probs}; sed -i \'s/.bam//g\' {output.probs}; ' +
        '{config[tools][monovarCluster][call]} {output.probs} {output.pdf}')
