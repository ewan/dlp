#!/bin/sh

rm -rf bin

mkdir -p bin/dlp/R # FIXME
mkdir -p bin/dlp/man
mkdir -p bin/dlp/data
mkdir -p bin/dlp/inst/java

cp R/DESCRIPTION bin/dlp/
cp R/NAMESPACE bin/dlp/
cp R/src/* bin/dlp/R/

cp java/Gibbs/spanish_mfiau_f0.txt bin/dlp/data/JavaTest.txt
cat > bin/dlp/data/JavaTest.R << END
JavaTest <- read.table('JavaTest.txt',header=T,sep=',')
JavaTest\$T1 <- C(factor(paste('S',JavaTest\$T1,sep='')), contr.sum)
JavaTest\$T3 <- factor(ifelse(sin((1:536)^2) - (sin(1:536))^2 > -0.1, "A","B"))
END

cp -R data/*RData bin/dlp/data/

cp -R java/Gibbs/bin/org bin/dlp/inst/java/
# FIXME: do the right thing with the dependencies
for f in dependencies/java/*
do
  cp -R "$f" bin/dlp/inst/java/
done
