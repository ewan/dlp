turk4 <- read.table("turk4.txt", header=T)
names(turk4) <- c("X1","X2","X3","C1","T1")

# Fix classes
turk4$C1 <- as.character(turk4$C1)
turk4[turk4$C1=='e',]$C1 <- 'ee'
turk4$S1 <- turk4$C1
unaffected <- turk4$X1 <= ((31/45)*(turk4$X2) - (802+2/9))
turk4[turk4$C1=='ee' & turk4$T1==1 & !(unaffected),]$S1 <- 'ae'
turk4$S2 <- turk4$S1
turk4[turk4$S2!='ae',]$S2 <- '!ae'
turk4$C1 <- factor(turk4$C1)
turk4$S1 <- factor(turk4$S1)
turk4$S2 <- factor(turk4$S2)
