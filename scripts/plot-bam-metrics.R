#! /usr/bin/RScript

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(scales)


# read the read-counts.txt
onTarget = read.table("read-counts.txt", sep = "\t", header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 9)))
onTarget$ID = factor(onTarget$ID)
# use factors to order X axis
onTarget$ID = factor(onTarget$ID, levels=c("SK-PB-191-G-1-Loop-IGO-05500-CR-12", "SK-PB-191-G-3-Loop-IGO-05500-CR-11", "SK-PB-191-G-10-Loop-IGO-05500-CR-10", "SK-PB-191-G-30-Loop-IGO-05500-CR-9"))




# plot on-target and off-target numbers for both total and unique reads
ggplot() + geom_boxplot(data=onTarget, aes(ID, TotalOnTargetFraction, colour='On Target Fraction'), size=2)  + geom_boxplot(data=onTarget, aes(ID, 1-TotalOnTargetFraction, colour='Off Target Fraction'), size=2)  + ggtitle("Read Distribution: TOTAL reads")  + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in")) + xlab("Sample")
# ggplot() + geom_boxplot(data=onTarget, aes(ID, TotalOnTargetFraction, colour='On Target Fraction'), size=2)  + geom_boxplot(data=onTarget, aes(ID, 1-TotalOnTargetFraction, colour='Off Target Fraction'), size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(.1, .1, .1, 1), "in")) + xlab("Sample")
ggsave("On-Target-Total.pdf", width=8, height=5)



ggplot() + geom_boxplot(data=onTarget, aes(ID, UniqueOnTargetFraction, colour='On Target Fraction'), size=2)  + geom_boxplot(data=onTarget, aes(ID, 1-UniqueOnTargetFraction, colour='Off Target Fraction'), size=2)  + ggtitle("Read Distribution: UNIQUE reads")  + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in")) + xlab("Sample")
# ggplot() + geom_boxplot(data=onTarget, aes(ID, UniqueOnTargetFraction, colour='On Target Fraction'),size=2)  + geom_boxplot(data=onTarget, aes(ID, 1-UniqueOnTargetFraction, colour='Off Target Fraction'),size=2) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(.1, .1, .1, 1), "in")) + xlab("Sample")
ggsave("On-Target-Unique.pdf", width=8, height=5)



# plot duplicate fraction
ggplot(onTarget, aes(ID, DuplicateFraction)) + geom_bar(stat="identity") + ggtitle("Duplicate Fraction") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in"))
ggsave("Duplicate-Fraction.pdf", width=8, height=5)



#Plot mean coverage from Waltz Metrics output (locuspocus-coverage.txt)
Coverage = read.table("waltz-coverage.txt", sep = "\t", header = TRUE, colClasses = c("character", rep("numeric", 2)))
Coverage$Sample = factor(Coverage$Sample)

# if you want to order X axis by Library Yield
#Coverage <- transform(Coverage, Sample=reorder(Sample, LibraryYield))
# use factors to order X axis
Coverage$Sample = factor(Coverage$Sample, levels=c("SK-PB-191-G-1-Loop-IGO-05500-CR-12", "SK-PB-191-G-3-Loop-IGO-05500-CR-11", "SK-PB-191-G-10-Loop-IGO-05500-CR-10", "SK-PB-191-G-30-Loop-IGO-05500-CR-9"))

Coverage$Sample = factor(Coverage$Sample, levels=c("SK-PB-191-H-1-Loop-IGO-05500-CR-24", "SK-PB-191-H-3-Loop-IGO-05500-CR-23", "SK-PB-191-H-10-Loop-IGO-05500-CR-22", "SK-PB-191-H-30-Loop-IGO-05500-CR-21"))



# plot total coverage
ggplot(Coverage, aes(Sample, TotalCoverage)) + geom_bar(stat="identity") + ggtitle("Mean On Target Coverage: TOTAL reads") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in"))
ggsave("Coverage-Total.pdf", width=8, height=5)

# ggplot() + geom_boxplot(data=Coverage, aes(Sample, TotalCoverage, colour='Coverage (TOTAL Reads)'),size=2)   + ggtitle("Mean On Target Coverage: TOTAL reads") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in"))
# ggplot() + geom_boxplot(data=Coverage, aes(Sample, TotalCoverage, colour='Coverage (TOTAL Reads)'),size=2) + geom_text(data=Coverage, aes(label=Input, x=Sample, y=1, angle=45, fontface="bold")) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(.1, .1, .1, 1), "in"))


# plot unique coverage
ggplot(Coverage, aes(Sample, UniqueCoverage)) + geom_bar(stat="identity") + ggtitle("Mean On Target Coverage: UNIQUE reads") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in"))
ggsave("Coverage-Unique.pdf", width=8, height=5)

# ggplot() + geom_boxplot(data=Coverage, aes(Sample, UniqueCoverage, colour='Coverage (UNIQUE Reads)'),size=2)   + ggtitle("Mean On Target Coverage: UNIQUE reads") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in"))
# ggplot() + geom_boxplot(data=Coverage, aes(Sample, UniqueCoverage, colour='Coverage (UNIQUE Reads)'),size=2) + geom_text(data=Coverage, aes(label=Input, x=Sample, y=1, angle=45, fontface="bold")) + theme(strip.text.x=element_text(size=20, angle=75)) + ggtitle(title) + theme(plot.title = element_text(size=30, face="bold"))+ theme(legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(.1, .1, .1, 1), "in"))




# plot Intervals-Samples-TotalCoverage-GC
Matrix = read.table("t5", sep = "\t", header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 2)))
Matrix$Sample = factor(Matrix$Sample)

# use factor to order x axis
Matrix$Sample = factor(Matrix$Sample, levels=c("SK-PB-191-G-1-IDT-IGO-05500-CR-4", "SK-PB-191-G-1-Loop-IGO-05500-CR-12", "SK-PB-191-G-1-UA-IGO-05500-CR-8", "SK-PB-191-G-3-IDT-IGO-05500-CR-3", "SK-PB-191-G-3-Loop-IGO-05500-CR-11", "SK-PB-191-G-3-UA-IGO-05500-CR-7", "SK-PB-191-G-10-IDT-IGO-05500-CR-2", "SK-PB-191-G-10-Loop-IGO-05500-CR-10", "SK-PB-191-G-10-UA-IGO-05500-CR_6", "SK-PB-191-G-30-IDT-IGO-05500-CR-1", "SK-PB-191-G-30-Loop-IGO-05500-CR-9", "SK-PB-191-G-30-UA-IGO-05500-CR-5"))

ggplot(Matrix, aes(Sample, Interval, size = TotalCoverage, color = GC)) + geom_point() + ggtitle("Intervals-Samples-TotalCoverage-GC") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), axis.text.y = element_text(size = 5))
ggsave("Intervals-Samples-TotalCoverage-GC.pdf", width=10, height=10)



# plot Intervals-Samples-UniqueCoverage-GC
Matrix = read.table("t6", sep = "\t", header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 2)))
Matrix$Sample = factor(Matrix$Sample)

# use factor to order x axis
Matrix$Sample = factor(Matrix$Sample, levels=c("s_KAPA_1", "s_KAPA_2", "s_KAPA_5", "s_KAPA_10", "s_Runicon_1", "s_Runicon_2", "s_Runicon_5", "s_Runicon_10", "s_Swift_1", "s_Swift__2", "s_Swift_5", "s_Swift_10", "s_Normal_Pooled_UNK_1"))

ggplot(Matrix, aes(Sample, Interval, size = UniqueCoverage, color = GC)) + geom_point() + ggtitle("Intervals-Samples-UniqueCoverage-GC") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), axis.text.y = element_text(size = 5))
ggsave("Intervals-Samples-UniqueCoverage-GC.pdf", width=10, height=10)




# plot CN normalized covergae
Matrix = read.table("waltz-cn-total.txt", sep = "\t", header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 2)))
Matrix$Sample = factor(Matrix$Sample)
ggplot(Matrix, aes(Sample, Interval, size = TotalCoverage, color = GC)) + geom_point() + ggtitle("CN-probes-normalized Total Coverage") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face="bold"), axis.text.y = element_text(size = 8, face="bold"))
ggsave("waltz-cn-total.pdf", width=10, height=10)


Matrix = read.table("waltz-cn-unique.txt", sep = "\t", header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 2)))
Matrix$Sample = factor(Matrix$Sample)
ggplot(Matrix, aes(Sample, Interval, size = UniqueCoverage, color = GC)) + geom_point() + ggtitle("CN-probes-normalized Unique Coverage") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face="bold"), axis.text.y = element_text(size = 8, face="bold"))
ggsave("waltz-cn-unique.pdf", width=10, height=10)






# Grouped total and unique coverage
ggplot(Coverage, aes(Sample, Coverage)) + geom_bar(aes(fill=Type), stat="identity", position=position_dodge()) + ggtitle("Total and Unique Coverage") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in"))
ggsave("Coverage-Comparison.pdf", width=8, height=5)


# Factors


onTarget$ID = factor(onTarget$ID, levels=c("SK-PB-191-H-1-IDT-IGO-05500-CR-16", "SK-PB-191-H-1-Loop-IGO-05500-CR-24", "SK-PB-191-H-1-UA-IGO-05500-CR-20", "SK-PB-191-H-3-IDT-IGO-05500-CR-15", "SK-PB-191-H-3-Loop-IGO-05500-CR-23", "SK-PB-191-H-3-UA-IGO-05500-CR-19", "SK-PB-191-H-10-IDT-IGO-05500-CR-14", "SK-PB-191-H-10-Loop-IGO-05500-CR-22", "SK-PB-191-H-10-UA-IGO-05500-CR-18", "SK-PB-191-H-30-IDT-IGO-05500-CR-13", "SK-PB-191-H-30-Loop-IGO-05500-CR-21", "SK-PB-191-H-30-UA-IGO-05500-CR-17"))

onTarget$ID = factor(onTarget$ID, levels=c("SK-PB-191-G-1-IDT-IGO-05500-CR-4", "SK-PB-191-G-1-Loop-IGO-05500-CR-12", "SK-PB-191-G-1-UA-IGO-05500-CR-8", "SK-PB-191-G-3-IDT-IGO-05500-CR-3", "SK-PB-191-G-3-Loop-IGO-05500-CR-11", "SK-PB-191-G-3-UA-IGO-05500-CR-7", "SK-PB-191-G-10-IDT-IGO-05500-CR-2", "SK-PB-191-G-10-Loop-IGO-05500-CR-10", "SK-PB-191-G-10-UA-IGO-05500-CR-6", "SK-PB-191-G-30-IDT-IGO-05500-CR-1", "SK-PB-191-G-30-Loop-IGO-05500-CR-9", "SK-PB-191-G-30-UA-IGO-05500-CR-5"))


Coverage$Sample = factor(Coverage$Sample, levels=c("SK-PB-191-H-1-IDT-IGO-05500-CR-16", "SK-PB-191-H-1-Loop-IGO-05500-CR-24", "SK-PB-191-H-1-UA-IGO-05500-CR-20", "SK-PB-191-H-3-IDT-IGO-05500-CR-15", "SK-PB-191-H-3-Loop-IGO-05500-CR-23", "SK-PB-191-H-3-UA-IGO-05500-CR-19", "SK-PB-191-H-10-IDT-IGO-05500-CR-14", "SK-PB-191-H-10-Loop-IGO-05500-CR-22", "SK-PB-191-H-10-UA-IGO-05500-CR-18", "SK-PB-191-H-30-IDT-IGO-05500-CR-13", "SK-PB-191-H-30-Loop-IGO-05500-CR-21", "SK-PB-191-H-30-UA-IGO-05500-CR-17"))

Coverage$Sample = factor(Coverage$Sample, levels=c("SK-PB-191-G-1-IDT-IGO-05500-CR-4", "SK-PB-191-G-1-Loop-IGO-05500-CR-12", "SK-PB-191-G-1-UA-IGO-05500-CR-8", "SK-PB-191-G-3-IDT-IGO-05500-CR-3", "SK-PB-191-G-3-Loop-IGO-05500-CR-11", "SK-PB-191-G-3-UA-IGO-05500-CR-7", "SK-PB-191-G-10-IDT-IGO-05500-CR-2", "SK-PB-191-G-10-Loop-IGO-05500-CR-10", "SK-PB-191-G-10-UA-IGO-05500-CR-6", "SK-PB-191-G-30-IDT-IGO-05500-CR-1", "SK-PB-191-G-30-Loop-IGO-05500-CR-9", "SK-PB-191-G-30-UA-IGO-05500-CR-5"))




# plot viral presence
i = read.table("TD-hpv-P72_0_IGO_05500_EY_8_S8-intervals.txt", sep = "\t", header = TRUE, colClasses = c("character", "numeric", "numeric", "character", rep("numeric", 5)))
human=51325450
ggplot(filter(i, FragmentsMapped>0), aes(TargetName,FragmentsMapped/human)) + geom_point(size=5, color="red") + geom_text(aes(label=paste(FragmentsMapped, '/', human)), hjust=-0.1, vjust=0) + xlab("Target") + ylab("Viral Fragments / Human Fragments") + ggtitle("Viral Presence") + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), plot.margin = unit(c(.1, .1, .1, 1), "in")) + scale_y_continuous(labels=comma)
ggsave("P72-0.pdf", width=8, height=5)


















#





