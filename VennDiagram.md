Znev <- c()
n <- 1
l <- nrow(ortgr)
while (n < l) {
  if (ortgr[n,7] > 0) {
    Znev <- append(Znev, ortgr[n,1])
  }
  n <- n+1
}

install.packages("sf")
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)
library(sf)
jpeg(file="4spec.jpg", width = 1200, height = 800)
orthologs <- list(Blatger,Cryptsec,Mnat,Ofo,Rspe,Znev)
ggVennDiagram(orthologs,
              category.names = c("Bger","Csec","Mnat","Ofor","Rspe","Znev"),
              label = "count") +
  scale_fill_gradient(low="#B07BAC", high="#5F7367")
dev.off()
