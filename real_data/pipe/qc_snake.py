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

# This rule uses qualimap to generate a report about several statistics of a BAM file.
localrules: createBedForQualimap
rule createBedForQualimap:
    input:
        regions = config['resources'][ORGANISM]['regions']
    output:
        regions = config['resources'][ORGANISM]['regions'] + '_qual.bed'
    shell:
        'awk \'{{sub(/\r/,\"\"); if(NF > 1){{printf $0; for(i = NF; i < 6; ++i){{printf \"\\t*\"}}; printf \"\\n\"}}}}\' {input.regions} > {output.regions}'

rule qualimap_HTML:
    input:
        bam = '{sample}.bam',
        regions = config['resources'][ORGANISM]['regions'] + '_qual.bed'
    output:
        dir = '{sample}.bam_stats',
        file = '{sample}.bam_stats/qualimapReport.html'
    params:
        lsfoutfile = '{sample}.bam_stats/qualimapReport.html.lsfout.log',
        lsferrfile = '{sample}.bam_stats/qualimapReport.html.lsferr.log',
        scratch = config['tools']['qualimap']['scratch'],
        mem = config['tools']['qualimap']['mem'],
        time = config['tools']['qualimap']['time'],
        params = config['tools']['qualimap']['params']
    benchmark:
        '{sample}.bam_stats/qualimapReport.html.benchmark'
    threads:
        config['tools']['qualimap']['threads']
    conda:
        'envs/qualimap.yaml'
    shell:
        'if [[ ! -n $(samtools view {input.bam} | head -n 1) ]]; then touch {output.file}; else qualimap bamqc -bam {input.bam} -outdir {output.dir} -outformat HTML -os -feature-file {input.regions} --java-mem-size={config[tools][qualimap][mem]}M {params.params}; fi'

rule samtools_flagstat:
    input:
        bam = '{sample}.bam',
    output:
        flagstat = '{sample}.bam.flagstat',
    params:
        mem = config['tools']['samtools']['flagstat']['mem'],
        scratch = config['tools']['samtools']['flagstat']['scratch'],
        time = config['tools']['samtools']['flagstat']['time'],
        lsfoutfile = '{sample}.bam.flagstat.lsfout.log',
        lsferrfile = '{sample}.bam.flagstat.lsferr.log',
        params = config['tools']['samtools']['flagstat']['params']
    threads:
        config['tools']['samtools']['flagstat']['threads']
    conda:
        'envs/samtools.yaml'
    shell:
        'samtools flagstat {params.params} {input.bam} > {output.flagstat}'

rule multiqcBam:
    input:
        qualimap = expand('{{analysisStep}}/{sample}.bam_stats/qualimapReport.html', sample = SAMPLENAMES),
        samtools = expand('{{analysisStep}}/{sample}.bam.flagstat', sample = SAMPLENAMES),
    output:
        '{analysisStep}/multiqc_report.html'
    params:
        lsfoutfile = '{analysisStep}/multiqc_report.html.lsfout.log',
        lsferrfile = '{analysisStep}/multiqc_report.html.lsferr.log',
        scratch = config['tools']['multiqc']['scratch'],
        mem = config['tools']['multiqc']['mem'],
        time = config['tools']['multiqc']['time'],
        params = config['tools']['multiqc']['params']
    benchmark:
        '{analysisStep}/multiqc_report.html.benchmark'
    threads:
        config['tools']['multiqc']['threads']
    conda:
        'envs/multiqc.yaml'
    shell:
        '{config[tools][multiqc][call]} --outdir {wildcards.analysisStep} --filename multiqc_report.html {params.params} {wildcards.analysisStep}'
