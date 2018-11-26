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
df$alpha <- as.factor(df$alpha)
df$tool <- factor(df$tool, levels = c("Monovar", "SCIPhI"))

ggplot(data = df, aes(x = cells, y = recall, fill = tool, alpha = alpha)) + 
  #geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot() +
  xlab("Number of cells") +
  ylab("Recall") +
  #scale_y_continuous(limits = c(0.75, 1)) +
  theme_bw() +
  theme(text = element_text(size=25)) +
  labs(fill='Method', alpha="Initial copies") +
  #scale_color_manual(values = c("firebrick3", "steelblue")) +
  scale_fill_manual(values = c("firebrick3", "steelblue"), labels=c("Monovar", expression(paste("SCI", Phi)))) +
  scale_alpha_manual(values=c(0.1, 0.25, 0.5, 0.75, 1)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1, 0.25, 0.5, 0.75, 1)),
                                              colour=NA))) 

ggsave(paste(gsub(".txt","",inputName), "_rec.pdf", sep=""), width = 14, height = 7)

ggplot(data = df, aes(x = cells, y = precision, fill = tool, alpha = alpha)) + 
  #geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot() +
  xlab("Number of cells") +
  ylab("Precision") +
  #scale_y_continuous(limits = c(0.7, 1)) +
  theme(text = element_text(size=25)) +
  labs(fill='Method', alpha="Initial copies") +
  #scale_color_manual(values = c("firebrick3", "steelblue")) +
  scale_fill_manual(values = c("firebrick3", "steelblue"), labels=c("Monovar", expression(paste("SCI", Phi)))) +
  scale_alpha_manual(values=c(0.1, 0.25, 0.5, 0.75, 1)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1, 0.25, 0.5, 0.75, 1)),
                                              colour=NA))) 
ggsave(paste(gsub(".txt","",inputName), "_pre.pdf", sep=""), width = 14, height = 7)

ggplot(data = df, aes(x = alpha, y = f1, fill = tool, alpha = cells)) + 
  geom_boxplot() +
  xlab("Initial chromosome copies") +
  ylab("F1 score") +
  scale_y_continuous(limits = c(min(0.7, df$f1 - 0.01), 1)) +
  #theme_bw() +
  theme(text = element_text(size=25)) +
  labs(fill='Method', alpha="Number of cells") +
  #scale_color_manual(values = c("firebrick3", "steelblue")) +
  scale_fill_manual(values = c("firebrick3", "steelblue"), labels=c("Monovar", expression(paste("SCI", Phi)))) +
  scale_alpha_manual(values=c(0.2,0.55, 1)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1, 0.55, 1)),
                                              colour=NA))) +
  theme(legend.text.align = 0)

ggsave(paste(gsub(".txt","",inputName), "_f1.pdf", sep=""), width = 14, height = 7)
