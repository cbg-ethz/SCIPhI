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
df$zyg <- as.factor(df$zyg)

ggplot(data = df, aes(x = zyg, y = recall, fill = tool)) + 
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Homozygosity rate") +
  ylab("Recall") +
  expand_limits(y=0.55) +
  expand_limits(y=1) +
  theme(legend.position = c(0.8, 0.2),
        legend.title=element_blank(),
        legend.text.align = 0,
        text = element_text(size=25),
        legend.key.size = unit(3., 'lines'))+
  scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_rec.pdf", sep=""))

ggplot(data = df, aes(x = zyg, y = precision, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Homozygosity rate") +
  ylab("Precision") +
  expand_limits(y=0.9) +
  expand_limits(y=1) +
  theme(legend.position = c(0.8, 0.2),
        legend.title=element_blank(),
        legend.text.align = 0,
        text = element_text(size=25),
        legend.key.size = unit(3., 'lines')) +
  scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_pre.pdf", sep=""))

ggplot(data = df, aes(x = zyg, y = f1, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Homozygosity rate") +
  ylab("F1 score") +
  expand_limits(y=0.6) +
  expand_limits(y=1) +
  theme(legend.position = c(0.8, 0.2),
        legend.title=element_blank(),
        legend.text.align = 0,
        text = element_text(size=25),
        legend.key.size = unit(3., 'lines')) +
  scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_f1.pdf", sep=""))
