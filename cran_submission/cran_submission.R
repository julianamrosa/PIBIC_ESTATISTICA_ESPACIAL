library(roxygen2) # In-Line Documentation for R 
library(devtools) # Tools to Make Developing R Packages Easier
library(testthat) # Unit Testing for R
library(usethis)  # Automate Package and Project Setup
library(rhub)

setwd("C:/Users/Juliana Rosa/OneDrive/Documents/PIBIC/pacote_final")

# Check for CRAN specific requirements using rhub and save it in the results 
# objects
verify_mgwnbr <- rhub::check_for_cran("mgwnbr")
#check_rhub()
#rhub::check()
# Get the summary of your results
verify_mgwnbr$cran_summary()

# Generate your cran-comments.md, then you copy-paste the output from the function above
setwd("C:/Users/Juliana Rosa/OneDrive/Documents/PIBIC/pacote_final/mgwnbr")
usethis::use_cran_comments()
usethis::use_news_md()

# Submeter o pacote
devtools::release()
