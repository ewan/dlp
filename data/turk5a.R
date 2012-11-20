turk5a <- read.table("turk5a.txt", header=T)
names(turk5a) <- c("X1","X2","X3","C1","T1")

# Fix classes
turk5a$C1 <- as.character(turk5a$C1)
turk5a[turk5a$C1=='e',]$C1 <- 'ee'
turk5a[turk5a$C1=='a',]$C1 <- 'aa'
turk5a$S1 <- turk5a$C1
unaffected <- turk5a$X1 <= ((31/45)*(turk5a$X2) - (802+2/9))
turk5a[turk5a$C1=='ee' & turk5a$T1==1 & !(unaffected),]$S1 <- 'ae'
turk5a$C1 <- factor(turk5a$C1)
turk5a$S1 <- factor(turk5a$S1)
