library(shiny)

source("tree/combineplots.R")

runApp(appDir=getwd(),port=2102, launch.browser = FALSE,  host = getOption("shiny.host", "0.0.0.0"))