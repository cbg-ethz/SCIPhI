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
import csv

def column(matrix, i):
    return [row[i] for row in matrix]

def distance(cell1, cell2):
    dist = 0
    for pos in range(len(cell1)):
        if (cell1[pos] != cell2[pos]):
            dist = dist + 1;
    return dist

def pairwiseCellDistance(cell):
    cellDistMatrix = [[0 for x in range(len(cell[0]))] for y in range(len(cell[0]))]
    for x in range(len(cell[0])):
        for y in range(x + 1, len(cell[0])):
            cellDistMatrix[x][y] = distance(column(cell, x), column(cell ,y))
    return cellDistMatrix

def matrixDistance(cell1DistMatrix, cell2DistMatrix):
    dist = 0
    for x in range(len(cell1DistMatrix)):
        for y in range(x, len(cell1DistMatrix)):
            dist += (abs(cell1DistMatrix[x][y] - cell2DistMatrix[x][y]))
    return dist

def columnWiseDist(cell1, cell2):
    dist = 0
    for col in range(len(cell1[0])):
        dist += distance(column(cell1, col), column(cell2, col))
    return dist

cell1FileName = sys.argv[1]
cell2FileName = sys.argv[2]

file = open(cell1FileName, 'r')
file.readline()
cell1 = []
for line in file:
    cell1.append(line.strip().split("\t"))
cell1.sort(key=lambda x: x[0])

file = open(cell2FileName, 'r')
file.readline()
cell2 = []
for line in file:
    cell2.append(line.strip().split("\t"))
cell2.sort(key=lambda x: int(x[0]))

# insert rows missing in cell1 or cell2
cell1Extended = []
cell2Extended = []
posCell1 = 0; posCell2 = 0;

if len(cell1) == 0 and len(cell2) == 0:
    print("0\t0")
    sys.exit(0)

if len(cell1) > 0 and len(cell2) > 0:
    while (posCell1 < len(cell1)) or (posCell2 < len(cell2)):
        if cell1[posCell1][0] == cell2[posCell2][0]:
            cell1Extended.append(cell1[posCell1][1:])
            cell2Extended.append(cell2[posCell2][1:])
            posCell1 += 1
            posCell2 += 1
        elif cell1[posCell1][0] < cell2[posCell2][0]:
            cell1Extended.append(cell1[posCell1][1:])
            cell2Extended.append(['0'] * (len(cell1[posCell1]) - 1))
            posCell1 += 1
        else:
            cell1Extended.append(['0'] * (len(cell2[posCell2]) - 1))
            cell2Extended.append(cell2[posCell2][1:])
            posCell2 += 1

elif len(cell1) == 0:
    while posCell2 < len(cell2):
        cell1Extended.append(['0'] * (len(cell2[posCell2]) - 1))
        cell2Extended.append(cell2[posCell2][1:])
        posCell2 += 1

elif len(cell2) == 0:
    while posCell1 < len(cell1):
        cell2Extended.append(['0'] * (len(cell1[posCell1]) - 1))
        cell1Extended.append(cell1[posCell1][1:])
        posCell1 += 1

else:
    print("Housten, we have a problem")

cell1DistMatrix = pairwiseCellDistance(cell1Extended)
cell2DistMatrix = pairwiseCellDistance(cell2Extended)
averageDis = matrixDistance(cell1DistMatrix, cell2DistMatrix)
columnDistance = columnWiseDist(cell1Extended, cell2Extended)
print(averageDis, "\t", columnDistance, sep = '')
