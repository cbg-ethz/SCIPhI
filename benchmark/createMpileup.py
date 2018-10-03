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
import sys
import random
import scipy.stats
import numpy
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('-o', '--outFile', type=str)
parser.add_argument('-i', '--inFile', type=str)
parser.add_argument('-f', '--inMutFile', type=str)
parser.add_argument('-n', '--numPos', type=int)
parser.add_argument('-c', '--covMean', type=float)
parser.add_argument('-v', '--covVar', type=float)
parser.add_argument('-l', '--averageRegionLength', type=int)
parser.add_argument('-N', '--numCells', type=int)
parser.add_argument('-m', '--mdaErrorRate', type=float)
parser.add_argument('-e', '--sequencingErrorRate', type=float)
parser.add_argument('-d', '--mdaDupReplacement', type=float, default=0.0)
parser.add_argument('-z', '--cygocityCoEff', type=float, default=0.1)
parser.add_argument('-s', '--seed', type=int, default = 42)
parser.add_argument('-M', '--missingInfo', type=float, default = 0.0)
parser.add_argument('-a', '--alpha', type=int)
parser.add_argument('-b', '--beta', type=int)
args = parser.parse_args()

outFileName = args.outFile
refFileName = args.inFile
numPos = args.numPos
covMean = args.covMean
covVar = args.covVar
averageRegionLength = args.averageRegionLength
numCells = args.numCells
mdaEr = args.mdaErrorRate
er = args.sequencingErrorRate
mutPosFileName = args.inMutFile
alpha = args.alpha
beta = args.beta

random.seed(args.seed)
numpy.random.seed(args.seed)

def indexToChar(index):
    if (index == 0):
        return 'A'
    elif (index == 1):
        return 'C'
    elif (index == 2):
        return 'G'
    elif (index == 3):
        return 'T'

def charToIndex(char):
    if (char == 'A'):
        return 0
    elif (char == 'C'):
        return 1
    elif (char == 'G'):
        return 2
    elif (char == 'T'):
        return 3

def getStrand(nuc):
    if nuc == '.':
        if random.random() < 0.5:
            return '.'
        else:
            return ','
    else:
        if random.random() < 0.5:
            return nuc
        else:
            return nuc.lower()


def createFrquencies(cov, counts):
    finalCounts = [0, 0, 0, 0]
    countsCopy = counts.copy()
    sumCounts = sum(counts);
    i = 0
    while i <  cov:
        currentCov = sum(counts)
        
        rand = random.randint(1, currentCov)
        currentAltSum = counts[0]
        pos = 1
        while pos < 4:
            if (rand > currentAltSum):
               currentAltSum += counts[pos]
            else:
                break
            pos += 1
        pos -= 1
    
        if random.random() > mdaEr:
            counts[pos] += 1
            while random.random() < args.mdaDupReplacement and i < cov - 1:
                counts[pos] += 1
                i += 1
        else:
            newPos = int(random.randint(0,3))
            while newPos == pos: 
                newPos = int(random.randint(0,3))
            counts[newPos] += 1
            while random.random() < args.mdaDupReplacement and i < cov - 1:
                counts[newPos] += 1
                i += 1
        i += 1

    for i in range(0, 4):
        counts[i] -= countsCopy[i]
        for j in range(0, counts[i]):
            if random.random() > er:
                finalCounts[i] += 1
            else:
                newPos = int(random.randint(0,3))
                while newPos == i: 
                    newPos = int(random.randint(0,3))
                finalCounts[newPos] += 1

    return finalCounts

    
def writeMutType(outFile, refNuc, mutNuc, alleleAffacted, cov):
    nucs = ""
    counts = [0, 0, 0, 0]

    if alleleAffacted == '0':               # 2 chrom reference
        counts[charToIndex(refNuc)] += alpha + beta
    elif alleleAffacted == '1':             # 2 chrom: 1 chrom ref, 1 chrom mut
        counts[charToIndex(refNuc)] += alpha
        counts[charToIndex(mutNuc)] += beta
    elif alleleAffacted == '2' or alleleAffacted == '4': # 1 chrom reference
        counts[charToIndex(refNuc)] += alpha
    elif alleleAffacted == '3' or alleleAffacted == '5': # 1 chrom mut
        counts[charToIndex(mutNuc)] += beta
    elif int(alleleAffacted) >= 6: # copy number
        counts[charToIndex(refNuc)] += alpha
        counts[charToIndex(mutNuc)] += beta + beta * (int(alleleAffacted) - 5)
    elif int(alleleAffacted) < 0: # copy number
        counts[charToIndex(refNuc)] += alpha + beta * (-1 * int(alleleAffacted) - 5)
        counts[charToIndex(mutNuc)] += beta
    else:
        sys.exit(1)



    counts = createFrquencies(cov, counts)

    for i in range(0, 4):
        for j in range(0, counts[i]):
            if i == charToIndex(refNuc):
                nucs += getStrand('.')
            else:
                nucs += getStrand(indexToChar(i))
    outFile.write("\t" + nucs)

def writeWildType(outFile, refNuc, cov):
    writeMutType(outFile, refNuc, refNuc, '0', cov)

def writeQuals(outFile, localCov):
    outFile.write("\t" + "I" * localCov)

def createCoverageRegions():
    u = covMean / covVar
    r = covMean ** 2 / (covVar - covMean)
    pOpenRegion = 1.0 / averageRegionLength;

    covRegions = [[] for i in range(numCells)]
    for cell in range(0, numCells):
        start = 0
        pos = 1
        covRegions.append([])
        while pos < numPos - 1:
            if random.uniform(0, 1) < pOpenRegion:
                end = pos - 1
                if args.missingInfo > 0.0 and random.uniform(0, 1) < args.missingInfo:
                    covRegions[cell].append((start, end, 0))
                else:
                    covRegions[cell].append((start, end, numpy.random.negative_binomial(r, u)))
                start = pos
            pos += 1
        end = numPos
        covRegions[cell].append((start, end, numpy.random.negative_binomial(r, u)))

    return covRegions

posMap = {}
posFile = open(mutPosFileName, 'r')
for line in posFile:
    if not line.startswith("#"):
        lineSplit = line.strip().split("\t")
        posMap[lineSplit[0]] = lineSplit[1:]

covRegions = createCoverageRegions()

outFile = open(outFileName, 'w')
fastaFile = open(refFileName, 'w')
fastaFile.write(">chr1\n")
posInCovRange = [0] * numCells
for pos in range (1, numPos + 1):
    outFile.write("chr1\t" + str(pos))
    ref = indexToChar(random.randint(0,3))
    outFile.write("\t" + ref)
    fastaFile.write(ref)
    if pos % 80 == 0:
        fastaFile.write("\n")
    alt = indexToChar(random.randint(0,3))
    while alt == ref:
        alt = indexToChar(random.randint(0,3))
    for cell in range(0, numCells):
        localCov = covRegions[cell][posInCovRange[cell]][2]
        if localCov == 0:
            #print("test: " + str(localCov))
            outFile.write("\t0\t*\t*")
        else:
            localCov = round(numpy.random.normal(localCov, localCov / 10.0))
            #print("test: " + str(localCov))
            outFile.write("\t" + str(localCov))
            if str(pos) in posMap:
                writeMutType(outFile, ref, alt, posMap[str(pos)][cell], localCov)
            else:
                writeWildType(outFile, ref, localCov)
            writeQuals(outFile, localCov)
        if pos == covRegions[cell][posInCovRange[cell]][1]:
            posInCovRange[cell] += 1
    outFile.write("\n")



