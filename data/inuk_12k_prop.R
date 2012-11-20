inuk_12k_prop <- read.table("inuk_12k_prop.txt", header=T, sep=',')
names(inuk_12k_prop) <- c("X1","X2","C1","T1")

inuk_12k_prop$S1 <- factor(paste(as.character(inuk_12k_prop$C1), as.character(inuk_12k_prop$T1)))

