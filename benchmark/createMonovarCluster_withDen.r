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
library(dplyr)
library(gplots)
library(ggplot2)
library(reshape)
library(ggdendro)

df <- read.table(inputName, header = TRUE)
df <- df %>% dplyr::select_if(function(X) sum(!is.na(X))>=1)
df = data.matrix(df[,-1])
df[df=="NaN"]<-0.5
monovarDist <- df

ht2 <-heatmap.2(monovarDist,
          labRow = "",
          dendrogram = "both",
          trace = "none",
          col = colorRampPalette(c("white", "steelblue")))

mono = data.frame(monovarDist[rev(ht2$rowInd),ht2$colInd])
header <- colnames(mono)
header <- unlist(strsplit(header, '\t'))
header <- gsub(header, pattern = '.bam', replacement = '')
meltMono = melt(t(mono))
plot.dd <- ggplot(meltMono,  aes(X1,X2, fill=value)) +
  scale_fill_gradient(name = expression("P"[m]),
                    low = "white",
                    high = "steelblue",
                    limits=c(0,1)) +
  geom_raster() +
  scale_y_reverse() +
  ylab("Mutations") +
  xlab("Cells") +
  scale_x_discrete(limits=header, labels=header, breaks=header) +
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="right",
        legend.direction="vertical",
        legend.box='horizontal',
        #plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text("",margin=margin(-20,0,0,0)),
        #axis.title.x=element_text("",margin=margin(-3,0,0,0)),
        #axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

dendro.plot <- ggdendrogram(data = ht2$colDendrogram) + theme_void()
pdf(outputName)
grid.newpage()
print(plot.dd, vp = viewport(x = 0.5, y = .45, width = 1.0, height = 0.9))
print(dendro.plot, vp = viewport(x = 0.48, y = .918, width = 0.712, height = 0.15))
dev.off()
