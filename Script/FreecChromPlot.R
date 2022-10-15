#!/usr/bin/env Rscript
# 20220929
# pzw

library(ggplot2)
library(dplyr)


args <- commandArgs(TRUE)


# 设置非科学技术法
options(scipen=200)

# 传参
inputFile <- args[1]
outputDir <- args[2]

data <- read.table(inputFile, header=TRUE)
if (!file.exists(outputDir)) {
  dir.create(outputDir)
}
sample <- strsplit(basename(inputFile), ".ratio.")[[1]][1]

# 计算log值
data$log2 <- log2(data$Ratio)
data$CNTrue <- ifelse(data$Ratio < 0, NaN, 2 * data$Ratio)
data$color <- ifelse(data$CopyNumber < 2, colors()[490], ifelse(data$CopyNumber > 2, colors()[556], colors()[200]))


# 画图
for (i in c(1:22, "X", "Y")) {
  chromSelect = i
  png(paste(outputDir, "/", sample, ".", i, ".png", sep = ""), width=1080, height=540)
  chrom <- data %>% filter(Chromosome == chromSelect)
  p <- ggplot(chrom)
  p <- p + geom_point(aes(Start, CNTrue, colour=color), colour=chrom$color, stat='identity', alpha=0.5)
  p <- p + scale_y_continuous(limits = c(0, 5))
  p <- p + xlab(paste("Position, ", as.character(chromSelect))) + ylab("Copy Number Profile")
  print(p)
  dev.off()
}

# end
