#Analyze metabolomics data
setwd("D:/Users/eduardk.NET/Dropbox/Postdoc/pnnl/metabolomics cn limitation")
data1 <- read.table(file="140326_F4-5_exometa2.txt", sep="\t", header=T)
data2 <- cbind (data1$F4C5,
                data1$F4C7,
                data1$F4C8,
                data1$F4N1,
                data1$F4N2,
                data1$F4N4,
                data1$F5C5,
                data1$F5C7,
                data1$F5C8,
                data1$F5N1,
                data1$F5N2,
                data1$F5N3)
colnames(data2)
pairs(data1)
arc.pca1 <- princomp(data1,scores=TRUE,cor=TRUE)
summary(arc.pca1)
plot(arc.pca1)
biplot(arc.pca1)
arc.pca1$loadings
arc.pca1$scores
