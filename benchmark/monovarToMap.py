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

inSCIPhiName = sys.argv[1]
inMonovarName = sys.argv[2]
outFileName =  sys.argv[3]

def getProb(g0, g01, g11):
    g0Prob = pow(10, -1 * g0/10)
    g01Prob = pow(10, -1 * g01/10)
    g11Prob = pow(10, -1 * g11/10)

    mutProb = g01Prob + g11Prob
    sumProb = g0Prob + g01Prob + g11Prob
    return mutProb / sumProb

monovarMap = {}
monovarIn = open(inMonovarName, 'r')
probs = []
monovarCellOrder = {}
for line in monovarIn:
    if line.startswith("##"):
        continue

    lineSplit = line.strip().split("\t")
    if line.startswith("#CHROM"):
        counter = 0
        names = lineSplit[9:]
        for counter in range(len(names)):
            monovarCellOrder[names[counter].split(".")[0]] = counter  
        continue
    
    probs = lineSplit[9:-1]
    counter = 0
    for counter in range(len(probs)):
        if probs[counter] == './.':
            probs[counter] = 'NaN'
        else:
            #print(probs[counter])
            probs[counter] = getProb(float(probs[counter].split(',')[-3].split(":")[-1]), float(probs[counter].split(',')[-2]), float(probs[counter].split(',')[-1]))
        counter += 1

    monovarMap[lineSplit[0] + "_" + lineSplit[1]] = probs



posMap = {}
inSCIPhi = open(inSCIPhiName, 'r')
line = inSCIPhi.readline()
sciphiCellOrder = []
counter = 0
sciphiNames = line.strip().split("\t")[1:]
for value in line.strip().split("\t")[1:]:
    sciphiCellOrder.append(value.split(".")[0])
    #print("value: ", value)
    #cellMap[value] = counter
    #counter+=1

#print(monovarCellOrder)
#print(sciphiCellOrder)

outFile = open(outFileName, 'w')
outFile.write("cellName\t")
outFile.write("\t".join(sciphiNames))
outFile.write("\n")

for line in inSCIPhi:
    counter = 0
    lineSplit = line.strip().split("\t")

    if lineSplit[0] in monovarMap:
        tempMono = monovarMap[lineSplit[0]]
        outFile.write(lineSplit[0] + "\t")
        temp = len(tempMono) * [None]
        for counter in range(len(tempMono)):
            posInMono = monovarCellOrder[sciphiCellOrder[counter]]
            temp[counter] = str(tempMono[posInMono])
            counter += 1
        outFile.write("\t".join(temp))
        outFile.write("\n")



