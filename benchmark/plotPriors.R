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
inputName = args[1]

options(scipen=999)
df <- read.table(inputName, header = TRUE)
df$prior <- as.factor(df$prior)
df$tool <- as.character(df$tool)

tool_names <- list(
  'Monovar'="Monovar",
  'SCIPhI'=expression(paste("SCI", Phi))
)
tool_labeller <- function(variable,value){
  return(tool_names[value])
}

ggplot(data = df, aes(x = prior, y = recall, fill = tool)) + 
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Prior") +
  ylab("Recall") +
  facet_grid(.~tool, scales = "free_x", labeller=tool_labeller) +
  guides(fill = "none") +
  scale_color_manual(values = c("firebrick3", "steelblue")) +
  scale_fill_manual(values = c("firebrick3", "steelblue")) +
  theme(text = element_text(size=25),
        axis.text.x=element_text(angle=45,hjust=1))
  #scale_fill_discrete("", labels=c("Monovar", expression(paste("SCI", Phi))))
ggsave(paste(gsub(".txt","",inputName), "_rec.pdf", sep=""))

ggplot(data = df, aes(x = prior, y = precision, fill = tool)) + 
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Prior") +
  ylab("Precision") +
  expand_limits(y=0.9) +
  facet_grid(.~tool, scales = "free_x", labeller=tool_labeller) +
  guides(fill = "none") +
  scale_color_manual(values = c("firebrick3", "steelblue")) +
  scale_fill_manual(values = c("firebrick3", "steelblue")) +
  theme(text = element_text(size=25),
        axis.text.x=element_text(angle=45,hjust=1))
ggsave(paste(gsub(".txt","",inputName), "_pre.pdf", sep=""))

ggplot(data = df, aes(x = prior, y = f1, fill = tool)) + 
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Prior") +
  ylab("F1 score") +
  facet_grid(.~tool, scales = "free_x", labeller=tool_labeller) +
  guides(fill = "none") +
  scale_color_manual(values = c("firebrick3", "steelblue")) +
  scale_fill_manual(values = c("firebrick3", "steelblue")) +
  theme(text = element_text(size=25),
        axis.text.x=element_text(angle=45,hjust=1))
ggsave(paste(gsub(".txt","",inputName), "_f1.pdf", sep=""))
