For creating Venn Diagram, you can use the maximum of 4 species for its legibility. Create the Venn Diagram from the Orthogroup_GeneCount.csv file.


```
Znev <- c(). ### You will need to change the name of the array for each species to create the four arrays of orthogroups
n <- 1 ### here we set the counter
l <- nrow(ortgr) ### here we tell R the maximum number of rows it will have to iterate over
while (n <= l) { ### here, the n will get bigger by one for every iteration until it reaches higher number than l and ends the loop
  if (ortgr[n,7] > 0) {. ### this is a condition, you have to change the number of the column for a species you want to work with, R starts with 1!!! not 0 like python. Bear that in mind. To translate: if the number in the specific row of a column is bigger than zero:
    Znev <- append(Znev, ortgr[n,1]) ### append the GOterm to our created blank array
  }
  n <- n+1 ### increase the counter by one and go back to the while loop
}

install.packages("sf") #you just need to install these packages and load them
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)
library(sf)
jpeg(file="4spec.jpg", width = 1200, height = 800) ### here we tell R to prepare a file
orthologs <- list(Blatger,Cryptsec,Mnat,Ofo,Rspe,Znev) ### here we merge all the arrays into one
ggVennDiagram(orthologs, ### we use the function ggVennDiagram and give it our merged arrays, we name the categories and set the label for count
              category.names = c("Bger","Csec","Mnat","Ofor","Rspe","Znev"),
              label = "count") +
  scale_fill_gradient(low="#B07BAC", high="#5F7367") ### here we chose colours for the diagram, you can choose yours, should be contrasting
dev.off() #here we close the creation of the image and save it.
```
