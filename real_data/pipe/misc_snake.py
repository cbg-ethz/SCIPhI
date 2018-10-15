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


#This file is copied from https://github.com/cbg-ethz/NGS-pipe/blob/master/snake/common/misc/misc_snake.py

import os.path
import sys
import inspect
import copy

fail_instantly = False

class Error(object):
    def __init__(self, key, name):
        self.__key = key
        self.__name = name

    def __add__(self, other):
        return self

    def __call__(self, wildcards=None):
        sys.exit(
            """
            ===============================================
            You have not specified '{}' for '{}'
            ===============================================
            """.format(self.__key, self.__name))

    def __getitem__(self, value):
        return Error(key=self.__key, name=self.__name)

class Config(object):
    def __init__(self, kwargs, name='Config'):
        self.__name = name
        self.__members = {}
        for (key, value) in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=key)
            else:
                self.__members[key] = value
    
    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            if fail_instantly:
                sys.exit(
                    """
                    ===============================================
                    You have not specified '{}' for '{}'
                    ===============================================
                    """.format(key, self.__name))
            else:
                return Error(key=key, name=self.__name)

config = Config(config)

def getTumorSampleNames(type):
    output = []
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            lineSplit = line.strip().split()
            if lineSplit[2] == type:
                sample = lineSplit[1]
                output.append(sample)
    return output

def getSampleNames():
    output = [] #[samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
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
                    sample = lineSplit[1]
                    if not (sample in output):
                        output.append(sample)
    return output

def getExperimentNames():
    output = [] #[samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
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
                    if not (exp in output):
                        output.append(exp)
    return output

def getSampleNamesFromExperimentNames(wildcards):
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
                if exp not in expMap.keys():
                    expMap[exp] = []
                expMap[exp].append(sample)
    if wildcards.experiment in expMap:
        return expMap[wildcards.experiment]
    return ['ERROR']

def checkFilesAgainstSampleNames(files, sampleNames):
    finalFiles = []
    for f in files:
        for name in sampleNames:
            if name + "/" == f[0:len(name+"/")]:
                finalFiles.append(f)

    return finalFiles

def getSingleFastqFiles(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

    #return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
    #return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]

def getPairedFastqFiles(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq')]
   
    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

    #return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq.gz')]
    #return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq')]

def getPairedFastqFilesWithoutR(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

def getAllPairedFastqFilesWithoutR():
    files = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    files = files + [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

    return files
    #return [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    #return [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

def getNormalTumorFiles():
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    output = []
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if exp not in sampleMap.keys():
                    sampleMap[exp] = dict()
                if tpoint not in sampleMap[exp].keys():
                    sampleMap[exp][tpoint] = dict()
                if sampleType not in sampleMap[exp][tpoint].keys():
                    sampleMap[exp][tpoint][sampleType] = []
                sampleMap[exp][tpoint][sampleType].append(sample)
    for expKey, expValue in sampleMap.items():
        for tpointKey, tpointValue in expValue.items():
            if 'T' in tpointValue and 'N' in tpointValue:
                for sampleTumor in tpointValue['T']:
                    for sampleNormal in tpointValue['N']:
                        output.append(sampleTumor + '_vs_' + sampleNormal)
    return output

def getFinalTumorBams(wildcards):
    output = []
    samplesFromExp = getSampleNamesFromExperimentNames(wildcards)
    for name in TUMORSAMPLENAMES:
        if name in samplesFromExp:
            output.append(FINALBAMOUTDIR + name + ".bam")
    return output

def getFinalMpileupBams(wildcards):
    output = []
    samplesFromExp = getSampleNamesFromExperimentNames(wildcards)
    for name in SAMPLENAMES:
        if name in samplesFromExp:
            output.append(FINALBAMOUTDIR + name + ".bam")
    return output

localrules: createTumorBamFileSummary
rule createTumorBamFileSummary:
    input:
        bams = getFinalTumorBams,
    output:
        OUTDIR + '{location}/{experiment}_tumor_bamFileNames.txt'
    run:
        outfile = open(str(output), "w")
        for entry in input.bams:
            outfile.write(entry + '\n')
        outfile.close()

localrules: createBamFileSummary
rule createBamFileSummary:
    input:
        bams = getFinalMpileupBams,
    output:
        OUTDIR + '{location}/{experiment}_all_bamFileNames.txt'
    run:
        outfile = open(str(output), "w")
        for entry in input.bams:
            outfile.write(entry + '\n')
        outfile.close()

def getFinalBams(wildcards):
    if wildcards.type == "all":
        return getFinalMpileupBams(wildcards)
    if wildcards.type == "tumor":
        return getFinalTumorBams(wildcards)
    return "ERROR"

def getContigNames():
    regionsFile = open(config['resources'][ORGANISM]['regions'], 'r')
    contigs = set()
    for line in regionsFile:
        contigs.add(line.split("\t")[0])
    print(sorted(contigs))
    return sorted(contigs)


