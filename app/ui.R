# Options for Spinner
library(shiny)
library(shinycssloaders)
library(shinyFiles)
library(rxode2)
library(tidyverse)
library(MASS)
library(units)
library(Hmisc)
library(glue)
library (shinyjs)

options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

# Define UI for dataset viewer application
ui <- fluidPage(
  shinyjs::useShinyjs(),  
  
  headerPanel(h3("Regimen simulator"),windowTitle = "Regimen simulator"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 h4("Day 1 Dose"),   
                 fluidRow(
                   column(6, numericInput("dose1a", "1st Dose (mg):", 
                                          300, min = 1, max = 500)),
                   column(6, numericInput("dose1b", "(Optional) 2nd Dose (mg):", 
                                          100, min = 0, max = 500),)
                 ),
                 h4("Day 2 Dose"), 
                 fluidRow(
                   column(6, numericInput("dose2", "Dose (mg):", 
                                          100, min = 1, max = 500)),
                   column(6, radioButtons("dfreq1", "Frequency:",
                                          c("qd" = 1,
                                            "bid" = 2)))
                 ),
                 h4("Maintenance Dose"),   
                 fluidRow(
                   column(5, numericInput("dose3", "Dose (mg), from Day 3:", 
                                          30, min = 1, max = 500)),
                   column(4, numericInput("ldoseDay", "until Day:", 
                                          10, min = 4, max = 20),),
                   column(3, radioButtons("dfreq2", "Frequency:",
                                          c("qd" = 1,
                                            "bid" = 2)))
                 ),
                 h4("BID Options"),   
                 fluidRow(
                   column(12, radioButtons("bidfreq", "Dosing Periods:",
                                          c("10 and 14 hours" = 10,
                                            "12 and 12 hours" = 12),
                                          inline=T))
                 ),
                 actionButton(
                   inputId = "submit_mod",
                   label = "Run Model"),
                 fluidRow(
                   tags$br()
                 ),
                 h4("Plot and Dataset Options:"), 
                 radioButtons("semilog", "Plotting Options:", 
                              choices=c("linear", "semi-log"),
                              selected = "linear",
                              inline=T),
                 fluidRow(
                   column(6, numericInput("day1", "PK parameters for Day:", 
                                          3, min = 1, max = 20)),
                   column(6, numericInput("day2", "and Day:", 
                                          10, min = 1, max = 20))
                 ),
                 fluidRow(
                   column(6, downloadButton(
                     outputId = "downloadCSV",
                     label = "Download Datasets",
                     icon = icon("file-download"))),
                   column(6, downloadButton(
                     outputId = "downloadPlot",
                     label = "Download Plot",
                     icon = icon("file-download")))
                 )
    ),
    
    
    mainPanel(width = 9,
              h4("Median, 5th and 95th pctl concentration-time profile"),
              withSpinner(plotOutput("CpPlot"), type = 1),
              fluidRow(
                # Input: Choose dataset ----
                selectInput("dataset", "View dataset:",
                            choices = c("Summary of PK parameters","Summary of PK",
                                        "PK parameters","PK (Sample of 10 Subjects)"))
              ),
              tableOutput("table"))
    
  )
)

