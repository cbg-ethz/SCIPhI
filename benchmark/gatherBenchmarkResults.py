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

outFileName = sys.argv[1]
files = sys.argv[2:]

outFile = open(outFileName, 'w')
outFile.write("cells\tmutations\tdrop\tzyg\tcpn\tvio\trun\tprior\ttool\tTP\tTN\tFP\tFN\trecall\tprecision\tf1\n") 
if 'sciphi' in files[0]:
    baseDir = len(files[0].strip().split("/")) - 11
else:
    baseDir = len(files[0].strip().split("/")) - 9

#print(files)
#print(baseDir)
for file in files:
    out = []
    fileSplit = file.strip().split("/")
    numCells = fileSplit[baseDir].split("_")[1]
    numMuts = fileSplit[baseDir + 1].split("_")[1]
    drop = fileSplit[baseDir + 2].split("_")[1]
    clbm = fileSplit[baseDir + 3].split("_")[1]
    cpn = fileSplit[baseDir + 4].split("_")[1]
    vio = fileSplit[baseDir + 5].split("_")[1]
    run = fileSplit[baseDir + 6].split("_")[1]
    prior = fileSplit[baseDir + 7].split("_")[1]
    tool = fileSplit[baseDir + 8]
    if tool == "sciphi":
        tool = "SCIPhI"
    if tool == "monovar":
        tool = "Monovar"
    out.append(numCells)
    out.append(numMuts)
    out.append(drop)
    out.append(clbm)
    out.append(cpn)
    out.append(vio)
    out.append(run)
    out.append(prior)
    out.append(tool)

    inFile = open(file, 'r')
    inFile.readline()
    line = inFile.readline()
    line = line.strip()
    out.append(line)

    lineSplit = line.split('\t')
    TP = float(lineSplit[0])
    #TN = float(lineSplit[1])
    FP = float(lineSplit[2])
    FN = float(lineSplit[3])
    recall = TP / (TP + FN)
    precision = TP / (TP + FP)
    f1 = 2.0 * precision * recall / (precision + recall)

    out.append(str(recall))
    out.append(str(precision))
    out.append(str(f1))

    outFile.write("\t".join(out) + '\n')

