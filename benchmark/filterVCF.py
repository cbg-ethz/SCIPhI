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

inFileName = sys.argv[1]
outFileName = sys.argv[2]

inFile = open(inFileName, 'r')
outFile = open(outFileName, 'w')

for line in inFile:
    if line.startswith("#"):
        outFile.write(line)
    else:
        lineSplit = line.strip().split("\t")
        numMutCells = 0
        for pos in range(9, len(lineSplit)):
            gt = lineSplit[pos].split(":")[0]
            if gt == "0/1" or gt == "1/1":
                numMutCells += 1
                if numMutCells > 2:
                    outFile.write(line)
                    break
