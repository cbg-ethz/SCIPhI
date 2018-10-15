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

probsFileName = sys.argv[1]
vcfFileName = sys.argv[2]

namesMap = {}
probsFile = open(probsFileName, 'r')
line = probsFile.readline()
lineSplit = line.strip().split("\t")
for i in range(1, len(lineSplit)):
    namesMap[lineSplit[i]] = i

posMap = {}
counter = 0
for line in probsFile:
    lineSplit = line.strip().split("\t")
    posMap[lineSplit[0]] = counter
    counter += 1

#print(posMap)

vcfFile = open(vcfFileName, 'r')
vcfNames = []
for line in vcfFile:
    if line.startswith("##"):
        continue

    lineSplit = line.strip().split("\t")

    if line.startswith("#"):
        vcfNames = lineSplit
        continue

    #if posMap[lineSplit[0] + "_" + lineSplit[1]] in [7,24]:
    if posMap[lineSplit[0] + "_" + lineSplit[1]] in [0]:
        for i in range(9, len(lineSplit)):
            if namesMap[vcfNames[i]] >= 340: #and namesMap[vcfNames[i]] < 30:
                cellSplit = lineSplit[i].split(":")
                if (int(cellSplit[2])>0):
                    print(cellSplit[1] + "\t" + cellSplit[2] + "\t" + str(int(cellSplit[1])/int(cellSplit[2])) + "\t" + lineSplit[0] + "_" + lineSplit[1] + "\t" + vcfNames[i])

