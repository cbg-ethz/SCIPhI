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

inFileName = sys.argv[1]

sameRank = []
inFile = open(inFileName, 'r')
for line in inFile:

    if "digraph" in line:
        print(line.strip())
        print("rankdir = TB;\nsubgraph {")
        continue

    if line.strip() == "}":
        print("{rank = same; " + ";".join(sameRank) + ";}\n}\n}")
        continue

    if not "label" in line:
        print(line.strip())
        continue

    if "shape=box," in line:
        lineSplit = line.strip().split('.bam')
        sameRank.append(lineSplit[0].split("[")[0])
        print(lineSplit[0] + "\";fontsize=30;];")

    elif "label=\"" in line:
        lineSplit = line.strip().split("\\n")
        if (math.sqrt(len(lineSplit))/10) > 0.5:
            print(line.split("[")[0] + "[label=\"" + str(len(lineSplit)-1) + "\";fontcolor=\"white\";shape=\"circle\"; width=" + str(math.sqrt(len(lineSplit)-1)/8) + ";fillcolor=\"#4682B4\";style=filled;fontsize=30;];")
        else:
            print(line.split("[")[0] + "[label=\"\";shape=\"circle\"; width=" + str(math.sqrt(len(lineSplit)-1)/8) + ";fillcolor=\"#4682B4\";style=filled;];")



    else:
        print("SHOOOT")

    

