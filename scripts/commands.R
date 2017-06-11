
library(ggplot2)
library(dplyr)
library(grid)



# these 2 things need to be parameters. also the custom capture name, eg. ERBB2
setwd("/Users/patelj1/workspace/Marianas/custom-panel-5500-AQ")


# First get the needed fields and convert to fractions (for TOTAL):
# somehow collated.txt seems to have an extra empty column, hence t.
JuberTable = read.table("t", sep = "\t", header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 7)))
JuberTable$FractionCustom_Total <- JuberTable$total.AL.Switch.reads / JuberTable$TotalMappedReads
JuberTable$FractionOffTarget_Total <- (JuberTable$TotalMappedReads - JuberTable$total.AL.Switch.reads)/ JuberTable$TotalMappedReads
JuberTable$FractionCustom_Unique <- JuberTable$unique.AL.Switch.reads / (JuberTable$TotalMappedReads*(1-JuberTable$DuplicateFraction))
JuberTable$FractionOffTarget_Unique <- (JuberTable$UniqueMappedReads - JuberTable$unique.AL.Switch.reads)/ (JuberTable$TotalMappedReads*(1-JuberTable$DuplicateFraction))



# plot on-target and off-target read numbers
title <- "Read Distribution: TOTAL reads"
pdf("Fractions_totalReads_5500-AL.pdf",width=15)
ggplot() + geom_boxplot(data=JuberTable, aes(factor(ID), FractionCustom_Total, colour='On Target Fraction'),size=2)  + geom_boxplot(data=JuberTable, aes(factor(ID), FractionOffTarget_Total, colour='Off Target Fraction'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

title <- "Read Distribution: UNIQUE reads"
pdf("Fractions_uniqueReads_5500-AL.pdf",width=15)
ggplot() + geom_boxplot(data=JuberTable, aes(factor(ID), FractionCustom_Unique, colour='On Target Fraction'),size=2) + geom_boxplot(data=JuberTable, aes(factor(ID), FractionOffTarget_Unique, colour='Off Target Fraction'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#Plot mean coverage from LocusPocus Metrics output (locuspocus-coverage.txt)
Coverage = read.table("locuspocus-coverage.txt", sep = "\t", header = TRUE, colClasses = c("character", rep("numeric", 2)))

# if you want to order X axis by Library Yield
Coverage <- transform(Coverage, Sample=reorder(Sample, LibraryYield))

title <- "Mean On Target Coverage: TOTAL reads"
pdf("Coverage_totalReads_5500_AL.pdf",width=15)
ggplot() + geom_boxplot(data=Coverage, aes(factor(Sample), TotalCoverage, colour='Coverage (TOTAL Reads)'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

title <- "Mean On Target Coverage: UNIQUE reads"
pdf("Coverage_uniqueReads_5500_AL.pdf",width=15)
ggplot() + geom_boxplot(data=Coverage, aes(factor(Sample), UniqueCoverage, colour='Coverage (UNIQUE Reads)'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + 	theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




# plot mean coverage from collated.txt (or t)
title <- "Mean Coverage Among the Targets: TOTAL reads"
pdf("Coverage_totalReads_5500_AR.pdf",width=15)
ggplot() + geom_boxplot(data=JuberTable, aes(factor(ID), covCustom_Total, colour='Coverage Custom'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


title <- "Mean Coverage Among the Targets: UNIQUE reads"
pdf("Coverage_uniqueReads_5500_AR.pdf",width=15)
ggplot() + geom_boxplot(data=JuberTable, aes(factor(ID), covCustom_Unique, colour='Coverage Custom'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()







ggplot(altRatios, aes(x=LocusPocus, y=Pipeline)) + geom_point() + geom_smooth(method=lm)




# LocusPocus Genotyping comparison with pipeline
altRatios = read.table("locuspocus-pipeline-alt-ratios.txt", sep = "\t", header = TRUE, colClasses = c("character", "numeric", rep("character", 5), rep("numeric", 2)))
ggplot(altRatios, aes(x=LocusPocus, y=Pipeline, color=EventType)) + geom_point() + geom_smooth(method=lm) + ggtitle("LocusPocus-Pipeline Alt Allele Ratios")
ggsave("LocusPocus-Pipeline-Alt-Allele-Rations.pdf")


# make a heatmap
ggplot(t4, aes(Patient, EventName)) + geom_tile(aes(fill=Ratio), colour="white", interpolate = FALSE) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_gradient(low="white", high="navyblue") + ggtitle("ERBB2 Mutations Called in cfDNA")


# scatter plot with a box
ggplot(t5, aes(Total, AlleleFrequency, color=EventName)) + geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_rect(xmin=0, xmax=1000, ymin=0, ymax=0.025, fill="white", alpha=0.01, color="black", inherit.aes = FALSE) + ggtitle("ERBB2 Mutations with Freq in [0.005, 0.1] in cfDNA") + xlab("Coverage")





















































































###  Code below this needs to be verified!!!






#########################################################################################
# Plot theoretical max coverage basd on input mass along with observed unique coverage for each sample:
#########################################################################################

sampleInputMasses <- read.table("sample_inputMass.tsv",sep="\t",header=1)
JuberTableWithmasses <- merge(sampleInputMasses, JuberTable, by = "ID")
M <- 1e-12  # mass in 1 base pair in ng
B <- 6e9  # number of bp in a single genome
JuberTableWithmasses$theoreticalMaxCov <- (JuberTableWithmasses$LIBRARY.INPUT.ng/(M*B)) * 2  # theoretical max coverage (assuming 100% ligation efficiency and no DNA lost during prep)


ggplot(JuberTableWithmasses,aes(x = LIBRARY.INPUT.ng, y = covCustom_Unique)) + geom_point(size=10) +
  theme(legend.title = element_text(size = 35)) + theme(legend.text = element_text(size = 10)) + theme(legend.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size=30, face="bold", colour="black")) +  theme(axis.title.x = element_text(face="bold", colour="black", size=38)) +
  theme(axis.text.y = element_text(size=30, face="bold", colour="black")) +  theme(axis.title.y = element_text(face="bold", colour="black", size=35)) +
  theme(legend.text=element_text(size=35, colour="black", face="bold")) + theme(legend.key.height=unit(3,"line"))



##########################################################################################################################################################################################
# Use Juber's data for making probe-specific coverage distribution plots:
##########################################################################################################################################################################################

# First split the table into which condition got the 1-gene panel and which got the 11-gene panel:
library(dplyr)

TOTAL_Juber_probe_cov_CUSTOM <- read.table("ERBB2-coverage-with-duplicates.txt", sep="\t", header=1)
UNIQUE_Juber_probe_cov_CUSTOM <- read.table("ERBB2-coverage.txt", sep="\t", header=1)


## Now plot the distributions:

TOTAL_JuberCov <- gather(TOTAL_Juber_probe_cov_CUSTOM, sampleName, MeanCoverage, X:X.33)
UNIQUE_JuberCov <- gather(UNIQUE_Juber_probe_cov_CUSTOM, sampleName, MeanCoverage, X:X.33)

pdf("TOTAL_Coverage_byTarget_5500_AL.pdf", width=15)
ggplot(TOTAL_JuberCov, aes(factor(Probe), MeanCoverage)) + geom_boxplot() +
  facet_grid(. ~Gene, scales="free") +
  theme(strip.text.x = element_text(size = 20)) +
  xlab("Target intervals")+ theme(axis.text.x = element_text(size=15, face="bold", colour="black")) +  theme(axis.title.x = element_text(face="bold", colour="black", size=25)) +
  theme(axis.text.y = element_text(size=15, face="bold", colour="black")) +  theme(axis.title.y = element_text(face="bold", colour="black", size=25)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("UNIQUE_Coverage_byTarget_5500_AL.pdf", width=15)
ggplot(UNIQUE_JuberCov, aes(factor(Probe), MeanCoverage)) + geom_boxplot() +
  facet_grid(. ~Gene, scales="free") +
  theme(strip.text.x = element_text(size = 20)) +
  xlab("Target intervals")+ theme(axis.text.x = element_text(size=15, face="bold", colour="black")) +  theme(axis.title.x = element_text(face="bold", colour="black", size=25)) +
  theme(axis.text.y = element_text(size=15, face="bold", colour="black")) +  theme(axis.title.y = element_text(face="bold", colour="black", size=25)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()















####   !/usr/bin/Rscript

library(ggplot2)

setwd(".")

# load the file
con = read.table("info.txt", sep = "\t", header = TRUE, colClasses = c("character", "numeric", rep("character", 2), rep("numeric", 16)))

# plot the histogram

# using bin of 50, can be changed easily
qplot(data=con, Length, geom="histogram", binwidth=50)
ggsave("length-histogram.pdf", width=10, height=6)


# plot the cumulative density lines

# get the pseudo counts, normalize them
pseudoCounts = (con[, 5:20]+1)
normalized = pseudoCounts/rowSums(pseudoCounts)

# convert the normalized data to a vector
v = as.vector((as.matrix(normalized)))

# create the factor vector from dinucleotide column names
l = names(con)[5:20]
Dinucleotide = gl(length(l), length(v)/length(l), labels = l)

# create a data frame and plot using stat_ecdf()
df = data.frame(v, Dinucleotide)
ggplot(df, aes(v, colour = Dinucleotide)) + stat_ecdf()
ggsave("ecdf.pdf",  width=15, height=15)
