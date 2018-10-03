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
df$cells <- as.factor(df$cells)

df$tool <- factor(df$tool, levels = c("Monovar_1", "Monovar_2","Monovar_10","Monovar_50","Monovar_100", "SCIPhI_1", "SCIPhI_2", "SCIPhI_10", "SCIPhI_50", "SCIPhI_100"))

ggplot(data = df, aes(x = cells, y = recall, fill = tool)) + 
  #geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Number of cells") +
  ylab("Recall") +
  scale_y_continuous(limits = c(0.6, 1)) +
  scale_fill_manual(values=c("#008D91", "#00A6AB", "#00BFC4", "#00D8DE", "#00F1F7", "#F6483C", "#F75F55", "#F8766D", "#F98D85", "#FAA49E")) +
  theme(legend.title=element_blank(),
        text = element_text(size=25))
  #scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_rec.pdf", sep=""), width = 14, height = 7)
  
ggplot(data = df, aes(x = cells, y = precision, fill = tool)) + 
  #geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Number of cells") +
  ylab("Precision") +
  scale_y_continuous(limits = c(0.95, 1)) +
  scale_fill_manual(values=c("#008D91", "#00A6AB", "#00BFC4", "#00D8DE", "#00F1F7", "#F6483C", "#F75F55", "#F8766D", "#F98D85", "#FAA49E")) +
  theme(legend.title=element_blank(),
        text = element_text(size=25))
  #scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_pre.pdf", sep=""), width = 14, height = 7)

ggplot(data = df, aes(x = cells, y = f1, fill = tool)) + 
  #geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Number of cells") +
  ylab("F1 score") +
  scale_y_continuous(limits = c(0.75, 1)) +
  scale_fill_manual(values=c("#008D91", "#00A6AB", "#00BFC4", "#00D8DE", "#00F1F7", "#F6483C", "#F75F55", "#F8766D", "#F98D85", "#FAA49E")) +
  theme(legend.title=element_blank(),
        text = element_text(size=25))
  #scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_f1.pdf", sep=""), width = 14, height = 7)

