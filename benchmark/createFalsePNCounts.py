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

positionsFileName = sys.argv[1]
mpileupFileName = sys.argv[2]
falsePNNucleotidesName = sys.argv[3]

positionsFile = open(positionsFileName, 'r')
positions = {}
for line in positionsFile:
    lineSplit = line.strip().split('\t')
    if lineSplit:
        positions[lineSplit[0]] = int(lineSplit[1])

mpileupFile = open(mpileupFileName, 'r')
falsePN = open(falsePNNucleotidesName, 'w')
for line in mpileupFile:
    lineSplit = line.strip().split('\t')
    pos = lineSplit[1]
    if pos in positions:
        falsePN.write(pos + "\t" + str(positions[pos]) + "\t" + lineSplit[4 + positions[pos] * 3] + "\n")
    
