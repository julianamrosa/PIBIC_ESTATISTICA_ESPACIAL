if (!require(spgwr)){
  install.packages("spgwr")
  library(spgwr)
}

data(georgia)
georgia <- as.data.frame(gSRDF)

setwd("C:/Users/Juliana Rosa/OneDrive/Documents/PIBIC/pacote_final/mgwnbr/data")

save(georgia, file="georgia.RData")


setwd("C:/Users/Juliana Rosa/OneDrive/Documents/PIBIC/pacote_final")

library(roxygen2)
library(devtools)

document("mgwnbr")
install("mgwnbr")
