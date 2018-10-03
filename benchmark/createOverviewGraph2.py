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
import math

inFileMutsName = sys.argv[1]
inFileGvName = sys.argv[2]
inFileNamesName = sys.argv[3]
inFileNameSimpleTree = sys.argv[4]

nameOrder = []
inFileSimpleTree = open(inFileNameSimpleTree, 'r')
for line in inFileSimpleTree:
    if "box" in line:
        nameOrder.append(line.strip().split(".bam")[0].split("\"")[-1])

names = []
inFileNames = open(inFileNamesName, 'r')
for line in inFileNames:
    lineSplit = line.strip().split("\t")
    if lineSplit[-1].startswith("C"):
        names.append(lineSplit[0].split("/")[-1].split(".bam")[0])


inFileMuts = open(inFileMutsName, 'r')
lineSplit = inFileMuts.readline().strip().split("\t")
muts = [0] * len(lineSplit)

for i in range(len(lineSplit)):
    muts[i] = float(lineSplit[i])

for line in inFileMuts:
    lineSplit = line.strip().split("\t")
    for i in range(len(lineSplit)):
        muts[i] += float(lineSplit[i])

for i in range(len(lineSplit)):
    muts[i] = int(round(muts[i]))

def getSameRank(nodeId2Label, nameOrder):
    result = []
    for value in nameOrder:
        result.append(str(nodeId2Label[value]))

    return result

sameRank = []
inFile = open(inFileGvName, 'r')
nodeId2Label = {}
for line in inFile:

    if "digraph" in line:
        print(line.strip())
        print("rankdir = TB;\nsubgraph {")
        continue

    if line.strip() == "}":
        sameRank = getSameRank(nodeId2Label, nameOrder)
        print("{rank = same; " + "->".join(sameRank) + "[style=invis];}\n}\n}")
        continue



    if not "label" in line:
        if "->" in line:
            if int(line.split("-")[0]) == 2* len(names) - 1:
                continue
        print(line.strip())
        continue

    if "shape=box," in line:
        lineSplit = line.strip().split("\"")
        index = int(lineSplit[1])
        index2 = line.strip().split("[")[0]
        index = index - len(names) + 1
        nodeId2Label[names[index]] = index2
        #sameRank.append(lineSplit[0].split("[")[0])
        print(lineSplit[0].split("[")[0] + "[shape=box,style=filled, fillcolor=white,label=\"" + names[index] + "\";fontsize=30;];")

    elif "label=\"" in line:
        lineSplit = line.strip().split("\"")
        index = int(lineSplit[1])
        if(index) < len(muts):
            num = muts[index]
            if (math.sqrt(num)/10) > 0.5:
                print(line.split("[")[0] + "[label=\"" + str(num) + "\";fontcolor=\"white\";shape=\"circle\"; width=" + str(math.sqrt(num)/8) + ";fillcolor=\"#4682B4\";style=filled;fontsize=30;];")
            else:
                print(line.split("[")[0] + "[label=\"\";shape=\"circle\"; width=" + str(math.sqrt(num)/8) + ";fillcolor=\"#4682B4\";style=filled;];")

    else:
        print("SHOOOT")

    

