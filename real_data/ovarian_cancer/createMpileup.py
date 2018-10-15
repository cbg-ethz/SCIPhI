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

class Locus:
    def __init__(self, chrom, pos, ref, alt):
        self.c = chrom
        self.p = pos
        self.r = ref
        self.a = alt
    def str(self):
        return self.c + "\t" + self.p + "\t" + self.r


class MpileupEntry:
    def __init__(self, ref, alt, altChar):
        self.r = ref
        self.a = alt
        self.ac = altChar
    def str(self):
        cov = self.r + self.a
        if cov == 0:
            return '\t' + str(cov) + "\t*\t*"
        seq = []
        for i in range(self.a):
            if random.random() < 0.5:
                seq.append(self.ac.upper())
            else:
                seq.append(self.ac.lower())
        for i in range(self.a, cov):
            if random.random() < 0.5:
                seq.append('.')
            else:
                seq.append(',')
        return '\t' + str(cov) + "\t" + ''.join(seq) + '\t' + ('I' * cov)

inFileName = sys.argv[1]
outFileName = sys.argv[2]
nameFileName = sys.argv[3]
bedFileName = sys.argv[4]

inFile = open(inFileName, 'r')
line = inFile.readline()
lineSplit = line.strip().split(',')

pileupDict = []

bedFile = open(bedFileName, 'w')
for pos in lineSplit[1::2]:
    posSplit = pos.split("_")
    bedFile.write(posSplit[0] + "\t" + posSplit[1] + "\t" + posSplit[1] + "\n")
    pileupDict.append([Locus(posSplit[0],posSplit[1],posSplit[2],posSplit[3])])

nameFile = open(nameFileName, 'w')
for line in inFile:
    lineSplit = line.strip().split(',')
    nameFile.write(lineSplit[0] + "\tCT\n")
    counter = 0
    for i in range(1,len(lineSplit),2):
        refCount = int(lineSplit[i])
        altCount = int(lineSplit[i+1])
        idx = int((i-1)/2)
        pileupDict[idx].append(MpileupEntry(refCount, altCount, pileupDict[idx][0].a))

outFile = open(outFileName, 'w')
for value in pileupDict:
    for entry in value:
        outFile.write(entry.str())
    outFile.write('\n')


