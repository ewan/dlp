library(devtools)
library(rJava)

hack_load_data <- function(path) {
    files <- dir(paste(path,"/data",sep=''), pattern = "\\.R$", full.names = TRUE)
    for (i in files) {
      sys.source(i, envir = as.environment("package:dlp"), chdir=T)
    }
}

PKG_DIR <- paste(getwd(), "/../bin/dlp", sep='') 

load_all(PKG_DIR)
hack_load_data(PKG_DIR)

.jinit()
.jaddClassPath(paste(PKG_DIR, "/inst/java", sep=''))

m1 <- dpmlmb(c("X1","X2","X3"), c("T1","T2","T3"), "C", data=JavaTest, nburnin=0, nsamp=10)
plot(m1)
m2 <- dpmlmvb(c("X1","X2","X3"), c("T1","T2","T3"), "C", data=JavaTest, nburnin=0, nsamp=10)
plot(m2)
m3 <- dpmlmvb(c("X1","X2","X3"), c("T1","T2","T3"), "C", data=JavaTest, nburnin=2000, nsamp=100)
plot(m3)
