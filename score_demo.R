rm(list=ls())
setwd('/xxx/')
getwd()
#######
library('dplyr')
library('ggplot2')
library('edgeR')
library('gplots')
library('ggpubr')
library(ggExtra)
library(aplot)
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
#############
scatter_data <- as.data.frame(read.csv('demo_expression_file.csv'
                                       ,header=T,row.names= 1))
head(scatter_data)
#############
names(scatter_data)
corr_fig <- ggplot(non_zero_data,aes(x=NINJ1,y=inflammatory_score))+
  geom_jitter(aes(color=group),alpha=0.5,
              size=1,width = .05,height = .05)+##
  geom_smooth(method='lm',
              fullrange=TRUE)+
  stat_cor(method = 'pearson') +#pearson,spearman
  scale_color_manual(values=c('#003049','#f77f00','#d62828'))+
  theme_classic()
print(corr_fig)
ggsave(corr_fig, file='NINJ1_vs_inflammatory.pdf', width = 6, height = 4)#