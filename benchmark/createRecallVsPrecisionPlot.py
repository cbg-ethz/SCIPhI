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
import os, sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

legends = []
inFileNames = []
pos = 1
while pos < len(sys.argv):
    if sys.argv[pos] == "-t":
        pos += 1;
        truthFileName = sys.argv[pos]
        pos += 1;
    elif sys.argv[pos] == "-o":
        pos += 1;
        outFileName = sys.argv[pos]
        pos += 1;
    elif sys.argv[pos] == "-i":
        pos += 1;
        legends.append(sys.argv[pos])
        pos += 1
        inFileNames.append(sys.argv[pos])
        pos += 1;
    else:
        print("Did not understand " + sys.argv[pos] + "!")
        pos += 1;

truthMap = {}
truthFile = open(truthFileName, 'r')
for line in truthFile:
    if not line.startswith("#"):
        lineSplit = line.strip().split("\t")
        truthMap[lineSplit[0]] = lineSplit[1:]
truthFile.close()

def getCategory(truthValue, testGenotype):
    if testGenotype == "./.":
        return "na"
    # case 1:1
    if (testGenotype == "0/1" or testGenotype == "1/1") and truthValue == "1":
        return 'TP'
    # case 0:0
    if testGenotype == "0/0" and truthValue == "0":
        return 'TN'
    # case 0:1
    if testGenotype == "0/0" and truthValue == "1":
        return 'FN'
    # case 1:0
    if (testGenotype == "0/1" or testGenotype == "1/1") and truthValue == "0":
        return 'FP'

    print(truthValue, testGenotype)
    sys.exit()

# get number of TP
tpTotal = 0
for key in truthMap:
    for value in truthMap[key]:
        tpTotal += int(value)

for idx, inFileName in enumerate(inFileNames):
    mut = []
    realFile = open(inFileName, 'r')
    for line in realFile:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                numCells = len(line.strip().split("\t")) - 9
        else:
            lineSplit = line.strip().split("\t")
            pos = lineSplit[1]
            tags = lineSplit[8].split(":")
            gq = -1 
            for tag in range(0, len(tags)):
                if tags[tag] == "GQ":
                    gq = tag
                    break
            for cell in range(9, 9 + numCells):
                splitCell = lineSplit[cell].split(":")
                gt = splitCell[0]
                if pos in truthMap:
                    category = getCategory(truthMap[pos][cell - 9], gt)
                    if category != "na":
                        mut.append((category, int(splitCell[gq]) if len(splitCell) > 1 else 0))
                else:
                    category = getCategory("0", gt)
                    if category != "na":
                        mut.append((category, int(splitCell[gq]) if len(splitCell) > 1 else 0))

    mut.sort(key=lambda tup: tup[1])
    mut = mut[::-1]

    print(mut)

    tpCounter = float(0)
    fpCounter = float(0)
    recalls = []
    precisions = []
    for value in mut:
        if value[0] == "TP":
            tpCounter += 1.0
        if value[0] == "FP":
            fpCounter += 1.0
        try:
            recall = tpCounter / tpTotal
        except ZeroDivisionError:
            recall = 0;
        recalls.append(recall)
        try:
            precision = tpCounter / (tpCounter + fpCounter)
        except ZeroDivisionError:
            precision = 0
        precisions.append(precision)

    print(precisions)

    plt.plot(recalls, precisions, label = legends[idx])
    
plt.axis([-0.05, 1.05, -0.05, 1.05])
lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylabel('Precision')
plt.xlabel('Recall')
plt.savefig(outFileName, bbox_extra_artists=(lgd,), bbox_inches='tight')




