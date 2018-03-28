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

inFileNameTool = sys.argv[1]
inFileNameMpileup = sys.argv[2]

inFileTool = open(inFileNameTool, 'r')
inTool = {}
for line in inFileTool:
    lineSplit = line.strip().split("\t")
    if lineSplit[0] in inTool:
        inTool[lineSplit[0]].append(lineSplit[1])
    else: 
        inTool[lineSplit[0]] = [lineSplit[1]]

inMpileup = open(inFileNameMpileup, 'r')
for line in inMpileup:
    lineSplit = line.strip().split("\t")
    if lineSplit[1] in inTool:
        for value in inTool[lineSplit[1]]:
            print(lineSplit[1] + "\t" + value + "\t" + lineSplit[4 + 3 * int(value)] + "\n")
