#' ---
#' title: "Class 5: R Graphics"
#' author: "Enrique Sandoval"
#' date: "Thu Apr 18 2019"
#' ---

#Class 5 R graphics

#2A. Line plot
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)

plot(weight$Age, weight$Weight, xlab="Age(months)", ylab="Weight(kg)", pch=15,
     cex=1.5, lwd=2, ylim=c(2,10), main="Some title", type="b")


#2B. Bar plot
feat <- read.table("bimm143_05_rstats/feature_counts.txt",
           sep= "\t",header=TRUE)
par(mar=c(5,11,1,1))
barplot(feat$Count, names.arg =feat$Feature, horiz = TRUE,
        main="Some title", las=1)


#3A. color vectors
counts <- read.table("bimm143_05_rstats/male_female_counts.txt",
           sep="\t", header=TRUE)

#Could also use read.delim()function 
counts <- read.delim("bimm143_05_rstats/male_female_counts.txt",
                     sep="\t", header=TRUE)
barplot(counts$Count, names.arg = counts$Sample, 
        las=2, col=c("red", "blue"))
barplot(counts$Count, names.arg = counts$Sample, 
        las=2, col=rainbow(10))


