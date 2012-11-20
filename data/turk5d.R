turk5d <- read.table("turk5d.txt", header=T)

# Fix classes
turk5d$S3 <- as.character(turk5d$S1)
turk5d$S3[turk5d$S3=='aa' & turk5d$T1==1] <- 'al'
turk5d$S2 <- turk5d$S3
turk5d$S2[turk5d$S2=='oe' & turk5d$T1==1] <- 'oel'
turk5d$S2 <- factor(turk5d$S2)
turk5d$S3 <- factor(turk5d$S3)

