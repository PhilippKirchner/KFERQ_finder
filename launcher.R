#For ease of use running the web app is wrapped in this laucher script
#This way shiny can be loaded automatically
#Only works if the wd is set to the directory containing the launcher script

#testing for the presence of dplyr and installing it if necessary

#testPkg <- function(x){
#if (!(x %in% rownames(installed.packages()))){
#  install.pakages(x)
#  if (!(x %in% rownames(installed.packages())))
#    stop(paste(x, " package could not be installed"))
#}
#}

testPkg <- function(x){
  if (!require(x, character.only=T)){
    install.packages(x)
    if (!require(x, character.only = T))
      stop(paste(x, " package could not be installed"))
  }
}

testPkg("shiny")
testPkg("dplyr")

#library(shiny)
#library(dplyr)

#running the web app
runApp("../KFERQ_finder")