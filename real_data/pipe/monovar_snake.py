rule monovar:
    input:
        bams = getFinalTumorBams,
        ref = config['resources'][ORGANISM]['reference'],
        mpileup = MPILEUPOUT + '{experiment}_tumor_complete.mpileup',
        fileNames = FINALBAMOUTDIR + '{experiment}_tumor_bamFileNames.txt'
    output:
        MONOVAROUT + '{experiment}.vcf'
    params:
        lsfoutfile = MONOVAROUT + 'monovar.vcf.lsfout.log',
        lsferrfile = MONOVAROUT + 'monovar.vcf.lsferr.log',
        scratch = config['tools']['monovar']['scratch'],
        mem = config['tools']['monovar']['mem'],
        time = config['tools']['monovar']['time']
    benchmark:
        MONOVAROUT + 'monovar.vcf.benchmark'
    threads:
        config['tools']['monovar']['threads']
    shell:
        ('{config[tools][monovar][activate]}; ' +
        'cat {input.mpileup} | ' +
        '{config[tools][monovar][call]} ' +
        '-p 0.002 -a 0.2 -t 0.05 ' +
        '-f {input.ref} -b {input.fileNames} -o {output}')
