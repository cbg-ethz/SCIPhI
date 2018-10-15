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

def charToIndex(char):
    if char == 'A' or char == 'a':
        return 0
    if char == 'C' or char == 'c':
        return 1
    if char == 'G' or char == 'g':
        return 2
    if char == 'T' or char == 't':
        return 3


inFileNameSciphi = sys.argv[1]
inFileNameMonovar = sys.argv[2]
inFileNameHaplo = sys.argv[3]
inFileNameMpileup = sys.argv[4]
inFileNameRefPos = sys.argv[5]
minCov = int(sys.argv[6])
outFileNameFP = sys.argv[7]

cellMpileupMap = {}
cellInMpileup = open(inFileNameMpileup, 'r')
for line in cellInMpileup:
    cov = []
    lineSplit = line.split("\t")
    pos = lineSplit[0] + '_' + lineSplit[1]
    for i in range(3, len(lineSplit), 3):
        cov.append(int(lineSplit[i]))
    cellMpileupMap[pos] = cov

posMap = {}
inPos = open(inFileNameRefPos, 'r')
for line in inPos:
    lineSplit = line.split("\t")
    if int(lineSplit[3]) >= minCov:
        posMap[lineSplit[0] + '_' + lineSplit[1]] = [int(lineSplit[4]),int(lineSplit[5]),int(lineSplit[6]),int(lineSplit[7])]

haploMap = {}
inHaplo = open(inFileNameHaplo, 'r')
for line in inHaplo:
    if line.startswith("#"):
        continue
    lineSplit = line.split("\t")
    pos = lineSplit[0] + '_' + lineSplit[1]
    index = lineSplit[0] + '_' + lineSplit[1] + '_' + lineSplit[3].upper() + '_' + lineSplit[4].upper()
    if pos in posMap and len(lineSplit[3]) == 1 and len(lineSplit[4]) == 1:
        haploMap[index] = float(lineSplit[5])
inHaplo.close()

sciphiCounts = [0, 0, 0]
sciphiCountsCell = [0, 0, 0]
sciphiCountsCellNM = [0, 0, 0]
inSciphi = open(inFileNameSciphi, 'r')
sciphiFP = {}
sciphiMap = {}
noCovSciphi = {}


for line in inSciphi:
    if line.startswith("#"):
        continue
    lineSplit = line.split("\t")
    pos = lineSplit[0] + '_' + lineSplit[1]
    index = lineSplit[0] + '_' + lineSplit[1] + '_' + lineSplit[3].upper() + '_' + lineSplit[4].upper()
    alt = lineSplit[4]
    sciphiMap[index] = 1
    if pos in posMap:
        if index in haploMap:
            sciphiCounts[0] += 1
            for i in range(9,len(lineSplit)-1):
                splitCell = lineSplit[i].split(":")
                if splitCell[0] != '0/0':
                    sciphiCountsCell[0] += 1
                    if int(splitCell[2]) > 0:
                        sciphiCountsCellNM[0] += 1
                    else:
                        noCovSciphi[pos + "_" + str(i-9)] = 1
                else:
                    sciphiCountsCell[2] += 1
                    if int(splitCell[2]) > 0:
                        sciphiCountsCellNM[2] += 1
                    else:
                        noCovSciphi[pos + "_" + str(i-9)] = 1

        else:
            sciphiCounts[1] += 1
            counter = 0
            covCounter = 0
            for i in range(9,len(lineSplit)-1):
                splitCell = lineSplit[i].split(":")
                if splitCell[0] != '0/0':
                    sciphiCountsCell[1] += 1
                    if int(splitCell[2]) > 0:
                        sciphiCountsCellNM[1] += 1
                    if int(splitCell[2]) >= 5:
                        covCounter += 1
                        if int(splitCell[1])/int(splitCell[2]) >= 0.5:
                            counter += 1
            sciphiFP[index] = [posMap[pos][charToIndex(alt)], counter, covCounter]
inSciphi.close()

#inHaplo = open(inFileNameHaplo, 'r')
for index in haploMap:
    if not (index in sciphiMap):
        sciphiCounts[2] += 1
        indexSplit = index.split("_")
        pos = indexSplit[0] + '_' + indexSplit[1]
        if pos in cellMpileupMap:
            cov = cellMpileupMap[pos]
            sciphiCountsCell[2] += len(cov)
            counter = 0
            for value in cov:
                if value > 0:
                    sciphiCountsCellNM[2] += 1
                else:
                    noCovSciphi[pos + "_" + str(counter)] = 1
                counter += 1
#inHaplo.close()

def getF1(values):
    r = values[0] / (values[0] + values[2])
    p = values[0] / (values[0] + values[1])
    return 2 * r * p / (r+p)

print("[TP, FP, FN, F1]", sciphiCounts, getF1(sciphiCounts))
print("[TP, FP, FN, F1]", sciphiCountsCell, getF1(sciphiCountsCell))
print("[TP, FP, FN, F1]", sciphiCountsCellNM,getF1(sciphiCountsCellNM))

def applySupFilter(lineSplit):
    counter = 0
    for entry in lineSplit[9:-1]:
        splitCell = entry.split(":")
        if entry != './.':
            altCell = int(splitCell[1].split(",")[1])
            if altCell >= 3:
                counter += 1
                if counter >= 2:
                    return True
    return False

monovarCounts = [0, 0, 0]
monovarCountsCell = [0, 0, 0]
monovarCountsCellNM = [0, 0, 0]
inMonovar = open(inFileNameMonovar, 'r')
monovarFP = {}
monovarMap = {}
covInMono = {}
for line in inMonovar:
    if line.startswith("#"):
        continue

    lineSplit = line.split("\t")
    if not applySupFilter(lineSplit):
        continue

    pos = lineSplit[0] + '_' + lineSplit[1]
    alt = lineSplit[4]
    index = lineSplit[0] + '_' + lineSplit[1] + '_' + lineSplit[3].upper() + '_' + lineSplit[4].upper()
    if pos in posMap:
        if index in haploMap:
            if lineSplit[6] == "PASS":
                monovarCounts[0] += 1
                monovarMap[index] = 1
                for i in range(9,len(lineSplit)-1):
                    splitCell = lineSplit[i].split(":")
                    if splitCell[0] != './.':
                        if splitCell[0] != '0/0':
                            if int(splitCell[2]) > 0:
                                monovarCountsCellNM[0] += 1
                        else:
                            if int(splitCell[2]) > 0:
                                monovarCountsCellNM[2] += 1
                    else:
                        covInMono[pos + "_" + str(i - 9)] = 1
                        if cellMpileupMap[pos][i - 9] != 0:
                            print("test:", pos, i -9)
        else:
            if lineSplit[6] == "PASS":
                monovarCounts[1] += 1
                counter = 0
                covCounter = 0
                for i in range(9,len(lineSplit)-1):
                    splitCell = lineSplit[i].split(":")
                    if len(splitCell) > 1 and splitCell[0] != '0/0':
                        altCell = int(splitCell[1].split(",")[1])
                        covCell = int(splitCell[2])
                        if covCell > 0:
                            monovarCountsCellNM[1] += 1
                        if covCell >= 5:
                            covCounter += 1
                            if altCell/covCell >= 0.5:
                                counter += 1
                monovarFP[index] = [posMap[pos][charToIndex(alt)], counter, covCounter]
inMonovar.close()

for index in haploMap:
    if not (index in monovarMap):
        monovarCounts[2] += 1
        indexSplit = index.split("_")
        pos = indexSplit[0] + '_' + indexSplit[1]
        if pos in cellMpileupMap:
            cov = cellMpileupMap[pos]
            monovarCountsCell[2] += len(cov)
            counter = 0
            for value in cov:
                if value > 0:
                    monovarCountsCellNM[2] += 1
                else:
                    covInMono[pos + "_" + str(counter)] = 1
                counter += 1
inHaplo.close()

print("[TP, FP, FN, F1]", monovarCounts, getF1(monovarCounts))
print("[TP, FP, FN, F1]", monovarCountsCellNM, getF1(monovarCountsCellNM))

outFP = open(outFileNameFP, 'w')
for key, value in sciphiFP.items():
    outFP.write("sciphi\t" + key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\n")

for key, value in monovarFP.items():
    outFP.write("monovar\t" + key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\n")
outFP.close()
