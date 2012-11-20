turkx <- read.table("turkx.txt", header=T)
names(turkx) <- c("X1","X2","X3","C1","T1")

# Fix classes
turkx$C1 <- as.character(turkx$C1)
turkx[turkx$C1=='e',]$C1 <- 'ee'
turkx[turkx$C1=='a',]$C1 <- 'aa'
turkx$S1 <- turkx$C1
unaffected <- turkx$X1 <= ((31/45)*(turkx$X2) - (802+2/9))
turkx[turkx$C1=='ee' & turkx$T1==1 & !(unaffected),]$S1 <- 'ae'
turkx$C1 <- factor(turkx$C1)
turkx$S1 <- factor(turkx$S1)
