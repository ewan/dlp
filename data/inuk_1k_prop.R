inuk_1k_prop <- read.table("inuk_1k_prop.txt", header=T, sep=',')
names(inuk_1k_prop) <- c("X1","X2","C1","T1")

inuk_1k_prop$S1 <- factor(paste(as.character(inuk_1k_prop$C1), as.character(inuk_1k_prop$T1)))
