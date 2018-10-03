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

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
outfileNoNAs = open(sys.argv[2].replace('.bin', '_noNAs.bin'), 'w')

outfile.write("#chrom\tpos\t")
outfileNoNAs.write("#chrom\tpos\t")
for line in infile:
    if line.startswith("#CHROM"):
        numFields = len(line.strip().split("\t"))
        for i in range(0, len(line.strip().split("\t")) - 9):
            outfile.write("cell" + str(i) + "\t")
            outfileNoNAs.write("cell" + str(i) + "\t")
        outfile.write("\n")
        outfileNoNAs.write("\n")
    if not line.startswith("#"): # skip the header lines
        splitLine = line.strip().split("\t")

        if splitLine[6] != "PASS":
            continue

        outfile.write(splitLine[0] + "\t" + splitLine[1] + "\t")
        outfileNoNAs.write(splitLine[0] + "\t" + splitLine[1] + "\t")
        for value in splitLine[9:numFields]:
            gt = value.split(":")[0]
            cov = value.split(":")[2]
            # normal distance
            if cov == '0':
                outfile.write('NA' + "\t")
            else:
                if gt == "0/1" or gt == "1/1":
                    outfile.write('1' + "\t")
                elif gt == "0/0":
                    outfile.write('0' + "\t")
                else:
                    print("Unknown genotype ", gt)

            # no NAs distance
            if gt == "0/1" or gt == "1/1":
                outfileNoNAs.write('1' + "\t")
            elif gt == "0/0":
                outfileNoNAs.write('0' + "\t")
            else:
                print("Unknown genotype ", gt)
        outfile.write("\n")
        outfileNoNAs.write("\n")
outfile.close()
outfileNoNAs.close()
