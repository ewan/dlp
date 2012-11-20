turk4s <- read.table("turk4s.txt", header=T)
names(turk4s) <- c("X1","X2","X3","C1", "T1", "T2", "T3", "T4", "T5", "T6", "T7",
"T8", "T9", "T10", "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18",
"T19", "T20", "T21", "T22", "T23", "T24", "T25", "T26", "T27", "T28")
# Fix classes
turk4s$C1 <- as.character(turk4s$C1)
turk4s[turk4s$C1=='e',]$C1 <- 'ee'
turk4s$C1 <- factor(turk4s$C1)

