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
library(ggplot2)

args <- commandArgs(TRUE)

inputName <- args[1]
df <- read.table(inputName, header = TRUE)
#df$dist =df$dist + 1
df$dist = as.numeric(df$dist)
df$cells <- as.factor(df$cells)
ggplot(data = df, aes(x = cells, y = dist, fill = exp)) +
  geom_boxplot(alpha = 0.5) +
  xlab("Number of cells") +
  ylab("Distance") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = c("firebrick3", "steelblue"),labels=c("Random", expression(paste("SCI", Phi)))) +
  theme(legend.title=element_blank(),
        #legend.position = c(0.8, 0.2),
        text = element_text(size=25),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines'))
ggsave(paste(gsub(".txt","",inputName), ".pdf", sep=""), width = 10, height = 7)

