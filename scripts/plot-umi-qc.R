#!/bin/bash

# load the libraries
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

# to plot family sizes
families = read.table("family-sizes.txt", sep = "\t", header = TRUE, colClasses = c("numeric", "character", "numeric", "character"))
families$Sample = factor(families$Sample)

ggplot(filter(families, FamilyType=="All"), aes(FamilySize, Frequency, color=Sample)) + geom_point() + geom_line() + xlab("Family Size") + ylab("Frequency") + coord_cartesian(xlim = c(0, 40))
ggsave("family-sizes-all.pdf", width=8, height=5)

ggplot(filter(families, FamilyType=="Simplex"), aes(FamilySize, Frequency, color=Sample)) + geom_point() + geom_line() + xlab("Family Size") + ylab("Frequency") + coord_cartesian(xlim = c(0, 40))
ggsave("family-sizes-simplex.pdf", width=8, height=5)

ggplot(filter(families, FamilyType=="Duplex"), aes(FamilySize, Frequency, color=Sample)) + geom_point() + geom_line() + xlab("Family Size") + ylab("Frequency") + coord_cartesian(xlim = c(0, 40))
ggsave("family-sizes-duplex.pdf", width=8, height=5)


# to plot clusters per position
clusters = read.table("clusters-per-position.txt", sep = "\t", header = TRUE, colClasses = c("numeric", "character"))
clusters$Sample = factor(clusters$Sample)
ggplot(clusters, aes(Clusters)) + geom_freqpoly(binwidth=1, aes(color=Sample)) + ylab("Frequency") + coord_cartesian(xlim = c(0, 30))
ggsave("clusters-per-position.pdf", width=8, height=5)

clusters = read.table("clusters-per-position-post-filtering.txt", sep = "\t", header = TRUE, colClasses = c("numeric", "character"))
clusters$Sample = factor(clusters$Sample)
ggplot(clusters, aes(Clusters)) + geom_freqpoly(binwidth=1, aes(color=Sample)) + ylab("Frequency") + coord_cartesian(xlim = c(0, 30))
ggsave("clusters-per-position-post-filtering.pdf", width=8, height=5)







