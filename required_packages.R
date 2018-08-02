#This script installs all required packages for the 0.8 version of
#app.R
#app_db.R
#call as: source("required_packges.R")

testPkg <- function(x){
  if (!require(x, character.only=T)){
    install.packages(x)
    if (!require(x, character.only = T))
      stop(paste(x, " package could not be installed"))
  }
}

sapply(c("shiny","shinyjs","dplyr","dbplyr","readr","DBI","pool"),testPkg)
