rule samtools_mpileup_single_cell:
    input:
        bams = getFinalBams,
        fileNames = FINALBAMOUTDIR + '{experiment}_{type}_bamFileNames.txt',
        ref = config['resources'][ORGANISM]['reference'],
        chrRegion = REGIONSOUT + config['resources'][ORGANISM]['regions'].strip().split("/")[-1].replace('.bed','') + '_{chrom}.bed'
    output:
        mpileup = temp(MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup')
    params:
        lsfoutfile = MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup.lsferr.log',
        params = config['tools']['samtools']['mpileup']['params'],
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
    conda:
        'envs/samtools.yaml'
    benchmark:
        MPILEUPOUT + '{experiment}_{type}-{chrom}.mpileup.benchmark'
    shell:
        'samtools mpileup -f {input.ref} {params.params} -b {input.fileNames} -l {input.chrRegion} > {output.mpileup}'

rule samtoolsCombineMpileup:
    input:
        mpileup = expand(MPILEUPOUT + '{{experiment}}_{{type}}-{chrom}.mpileup', chrom=getContigNames())
    output:
        mpileup = temp(MPILEUPOUT + '{experiment}_{type}_complete.mpileup'),
        gz = MPILEUPOUT + '{experiment}_{type}.mpileup.gz'
    params:
        lsfoutfile = MPILEUPOUT + '{experiment}_{type}_complete.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{experiment}_{type}_complete.mpileup.lsferr.log',
        params = config['tools']['samtools']['mpileup']['params'],
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
    conda:
        'envs/samtools.yaml'
    benchmark:
        MPILEUPOUT + '{experiment}_{type}_complete.mpileup.benchmark'
    shell:
        'cat {input.mpileup} > {output.mpileup}; gzip < {output.mpileup} > {output.gz}'

ruleorder: createBamFileSummaryScate > createBamFileSummary
localrules: createBamFileSummaryScate
rule createBamFileSummaryScate:
    input:
        bams = getFinalMpileupBams,
    output:
        SCIPHIOUT + '{experiment}_all_bamFileNames.txt'
    run:
        sampleMappingFile = open(SAMPLEMAPPING, 'r')
        sampleMapping = {}
        for line in sampleMappingFile:
            sampleMapping[line.strip().split('\t')[1]] = line.strip().split('\t')[2]

        outfile = open(str(output), "w")
        for entry in input.bams:
            sample = entry.split('/')[-1].replace('.bam','')
            outfile.write(entry + '\t' + sampleMapping[sample] + '\n')
        outfile.close()

rule sciphi:
    input:
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_all_complete.mpileup',
        fileNames = SCIPHIOUT + '{experiment}_all_bamFileNames.txt'
    output:
        tsv = SCIPHIOUT + '{run}/{experiment}_mut2Sample.tsv',
        probs = SCIPHIOUT + '{run}/{experiment}.probs',
        gv = SCIPHIOUT + '{run}/{experiment}.gv',
        params = SCIPHIOUT + '{run}/{experiment}.params.txt',
        vcf = SCIPHIOUT + '{run}/{experiment}.vcf',
        bParams = SCIPHIOUT + '{run}/{experiment}/best_index/nuc.tsv',
        bTree = SCIPHIOUT + '{run}/{experiment}/best_index/tree.gv' 
    params:
        lsfoutfile = SCIPHIOUT + '{run}/{experiment}.lsfout.log',
        lsferrfile = SCIPHIOUT + '{run}/{experiment}.lsferr.log',
        scratch = config['tools']['sciphi']['scratch'],
        mem = config['tools']['sciphi']['mem'],
        time = config['tools']['sciphi']['time'],
        out = SCIPHIOUT + '{run}/{experiment}',
        params = config['tools']['sciphi']['params'],
        outIndex = SCIPHIOUT + '{run}/{experiment}/index'
    benchmark:
        SCIPHIOUT + '{run}/{experiment}.benchmark'
    threads:
        1
    log:
        SCIPHIOUT + '{run}/{experiment}.log'
    shell:
        ('{config[tools][sciphi][call]} ' +
        '-o {params.out} ' +
        '--ol {params.outIndex} ' + 
        '-i {input.ref} ' +
        '--in {input.fileNames} ' +
        #'--cwm 2 ' +
        #'--ms 3 ' +
        #'--nmc 2 ' +
        '--lz 1 ' +
        '--seed {wildcards.run} ' +
        '{params.params} ' +
        '{input.mpileup}')

if not 'HAPLOTYPECALLERIN' in globals():
    HAPLOTYPECALLERIN = 'placeholder'
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = 'placeholder'
rule gatkHaplotypeCaller:
    input:
        bam = HAPLOTYPECALLERIN + '{sample}.bam',
        bai = HAPLOTYPECALLERIN + '{sample}.bai',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        vcf = HAPLOTYPECALLEROUT + '{sample}.g.vcf',
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + '{sample}.g.vcf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + '{sample}.g.vcf.lsferr.log',
        scratch = config['tools']['GATK']['haplotypeCaller']['scratch'],
        mem = config['tools']['GATK']['haplotypeCaller']['mem'],
        time = config['tools']['GATK']['haplotypeCaller']['time'],
        params = config['tools']['GATK']['haplotypeCaller']['params'],
    threads:
        config['tools']['GATK']['haplotypeCaller']['threads']
    benchmark:
        HAPLOTYPECALLEROUT + '{sample}.g.vcf.benchmark'
    log:
        HAPLOTYPECALLEROUT + '{sample}.g.vcf.log'
    conda:
        'envs/gatk.yaml'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T HaplotypeCaller ' +
        '{params.params} ' + 
        '-R {input.reference} ' + 
        '-I {input.bam} ' +
        '--emitRefConfidence GVCF ' +
        '-L {input.regions} ' +
        '-mmq 40 ' +
        '-mbq 30 ' + 
        '-o {output.vcf}')

rule gatkGenotypeGVCFs:
    input:
        vcf = HAPLOTYPECALLEROUT + '{sample}.g.vcf',
        reference = config['resources'][ORGANISM]['reference']
    output:
        vcf = HAPLOTYPECALLEROUT + '{sample}.vcf'
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + '{sample}.vcf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['GATK']['genotypeGVCFs']['scratch'],
        mem = config['tools']['GATK']['genotypeGVCFs']['mem'],
        time = config['tools']['GATK']['genotypeGVCFs']['time'],
        reference = config['resources'][ORGANISM]['reference'],
        params = config['tools']['GATK']['genotypeGVCFs']['params']
    threads:
        config['tools']['GATK']['genotypeGVCFs']['threads']
    benchmark:
        HAPLOTYPECALLEROUT + '{sample}.vcf.benchmark'
    log:
        HAPLOTYPECALLEROUT + '{sample}.vcf.log'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T GenotypeGVCFs ' +
        '{params.params}' +
        '-R {input.reference} ' +
        '--variant {input.vcf} ' +
        '-o {output.vcf}')
