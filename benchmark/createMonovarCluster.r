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
args <- commandArgs(trailingOnly = TRUE)
inputName <- args[1]
outputName <- args[2]

library(ComplexHeatmap)
library(circlize)
df <- read.table(inputName, header = TRUE)
monovarDist <- data.matrix(df[,-1])
ht2 = Heatmap(monovarDist,
              name = "Pm",
              column_title = "",
              col = colorRampPalette(c("white", "steelblue"))(100),
              show_row_dend = FALSE,
              #show_column_dend = FALSE,
              #show_column_names = FALSE,
              column_names_gp = gpar(fontsize = 20),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 20)))
pdf(outputName)
draw(ht2)
dev.off()
