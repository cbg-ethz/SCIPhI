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

inFileArtName = sys.argv[1]
inFileMpileupName = sys.argv[2]
inFileToolName = sys.argv[3]
outFileTxtName = sys.argv[4]
outFileFPName = sys.argv[5]
outFileFNName = sys.argv[6]

infileArt = open(inFileArtName, 'r')
infileArt.readline() # skip the header line

posDictTruth = {}
for line in infileArt:
    lineSplit = line.strip().split('\t')
    pos = lineSplit[0]
    posDictTruth[pos] = lineSplit[1:]
infileArt.close()

infileMpileup = open(inFileMpileupName, 'r')
posDictMpileup = {}
for line in infileMpileup:
    lineSplit = line.strip().split('\t')
    pos = lineSplit[1]
    posDictMpileup[pos] = lineSplit[3::3]
infileMpileup.close()
    
outfileFP = open(outFileFPName, 'w')
outfileFN = open(outFileFNName, 'w')

numRows = 0
truePositives = 0
trueNegatives = 0
falsePositives = 0
falseNegatives = 0
posDictInfered = {}
infileTool = open(inFileToolName, 'r')

ignoreNAs = False
if "noNAs" in inFileToolName:
    ignoreNAs = True

for line in infileTool:
    lineSplit = line.strip().split('\t')
    if not line.startswith("#") and not line.startswith("chrom"):
        numRows += 1
        pos = lineSplit[1]
        posDictInfered[pos] = lineSplit[1:]
        if pos in posDictTruth:
            cellId = 0
            for value in lineSplit[2:]:
                # ignore missing value fields if requires
                if ignoreNAs or int(posDictMpileup[pos][cellId]) > 0:

                    if value == '1':
                        if posDictTruth[pos][cellId] != '0' and posDictTruth[pos][cellId] != '4':
                            truePositives += 1
                        else:
                            falsePositives += 1
                            outfileFP.write(pos + "\t" + str(cellId) + "\n")
                    elif value == '0':
                        if posDictTruth[pos][cellId] != '0' and posDictTruth[pos][cellId] != '4':
                            falseNegatives += 1
                            outfileFN.write(pos + "\t" + str(cellId) + "\n")
                        else:
                            trueNegatives += 1
                    elif value == 'NA':
                        pass
                    else:
                        print("Unknown value: ", value)
                else:
                    if value != 'NA':
                        print("PROBLEM: 1:" + pos + " " + value)
                cellId += 1
        else:
            cellId = 0
            for value in lineSplit[2:]:
                # ignore missing value fields
                if ignoreNAs or int(posDictMpileup[pos][cellId]) > 0:

                    if value == '1':
                        falsePositives += 1
                        outfileFP.write(pos + "\t" + str(cellId) + "\n")
                    elif value == '0':
                        trueNegatives += 1
                    elif value == 'NA':
                        pass
                    else:
                        print("Unknown value: ", value)
                else:
                    if value != 'NA':
                        print("PROBLEM: 2")
                cellId += 1
infileTool.close()

for pos in posDictTruth:
    if not (pos in posDictInfered):
        numRows += 1
        cellId = 0
        for value in posDictTruth[pos]:
            # ignore missing value fields
            if ignoreNAs or int(posDictMpileup[pos][cellId]) > 0:

                if value == '0' or value == '4':
                    trueNegatives += 1
                else:
                    falseNegatives += 1
                    outfileFN.write(pos + "\t" + str(cellId) + "\n")
            cellId += 1

trueNegatives = "NA"
                
outfile = open(outFileTxtName, 'w')
outfile.write("truePositives\ttrueNegatives\tfalsePositives\tfalseNegatives\n")
outfile.write(str(truePositives) + "\t" + str(trueNegatives) + "\t" + str(falsePositives) + "\t" + str(falseNegatives) + "\n")
outfile.close()

