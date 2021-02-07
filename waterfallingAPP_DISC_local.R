#!/usr/bin/env Rscript
library(shiny)
library(shinyFiles)
library(htmltools)
library(wesanderson)
library(colourpicker)
library(shinydashboard)
# --- COLORS
color1 = wes_palette("FantasticFox1")
color2 = wes_palette("Royal1")
color3 = wes_palette("BottleRocket2")
color4 = wes_palette("Zissou1")
color5 = wes_palette("Chevalier1")
#---

ui <- dashboardPage(
  dashboardHeader(title = tags$code("BLINK")),
  dashboardSidebar(),
  dashboardBody(
    fluidRow(
      column(4, offset = 0,
             #tags$div(tags$h1(tags$b("Input Parameters"))),
             tags$hr(),
             tags$h3(tags$b("Experiment Set-Up")),
             "----------------------------------------------------------------",
             textInput(inputId = 'dataFolder', label = 'Data Folder Name', placeholder = 'Paste Data Folder Name Here'),
             "----------------------------------------------------------------",
             numericInput(inputId = 'baselineDur', label = 'Baseline Duration (in ms)', value = 200),
             numericInput(inputId = 'trialDur', label = 'Trial Duration (in ms)', value = 1200),
             "----------------------------------------------------------------",
             tags$b('Basic Analysis'),
             checkboxInput(inputId = 'doStartle', label = 'Analyze for Startle', value = 0),
             checkboxInput(inputId = 'doWaterfall', label = 'Analyze for Conditioned Responses', value = 0),
             tags$div("----------------------------------------------------------------"),
             tags$hr(),
             actionButton(inputId = "start",label = "Start")
      ), # End of column function
      column(4, offset = 0,
             tags$hr(tags$h3(tags$b("Advanced Options"))),
             tags$div("----------------------------------------------------------------"),
             tags$div(tags$em("The following options create detailed response visualizations but slow program 15-30x")),
             tags$br(tags$b(("Session Visualizaton"))),
             checkboxInput(inputId = 'include_session_waterfall', label = 'Create waterfall PNG (quick-view, faster)', value = 0),
             checkboxInput(inputId = 'postscript_flag', label = 'Create waterfall EPS (editable, slower)', value = 0),
             checkboxInput(inputId = 'include_timing_plots', label = 'Create timing plots', value = 0),
             tags$div("----------------------------------------------------------------"),
             tags$b("Visualze Average Traces"),
             checkboxInput(inputId = 'doTrajectories', label = 'Create average reponse trace for each animal', value = 0),
             numericInput(inputId = 'trajMin', label = 'Minimum responses required for inclusion', value = 10),
             tags$br(),
             checkboxInput(inputId = 'restrict', label = 'Control for learning?', value = 0),
             sliderInput("trajRestrict", "Response Range",
                         min = 0, max = 100, value = c(5,50))
      ), # End of column function
      
      column(4, offset = 0,
             tags$hr(tags$h3(tags$b("Plot Colors"))),
             tags$div("----------------------------------------------------------------"),
             h4("Analysis Plots"),
             colourInput(
               "group1_color", label = 'Group 1',
               palette = "limited",
               allowedCols = c('gray','red','blue','magenta','cyan','black',color1,color2,color3,color4,color5,'#44B655')),
             colourInput(
               "group2_color", label = 'Group 2',
               palette = "limited",
               allowedCols = c('blue','red','black','magenta','cyan','gray',color1,color2,color3,color4,color5,'#44B655')),
             colourInput(
               "onset_histogram", label = 'Onset Histograms',
               palette = "limited",
               allowedCols = c('black','blue','red','magenta','cyan','gray',color1,color2,color3,color4,color5,'#44B655')),
             colourInput(
               "peak_histogram", label = 'Peak Histograms',
               palette = "limited",
               allowedCols = c('black','blue','red','cyan','magenta','gray',color1,color2,color3,color4,color5,'#44B655')),
             tags$div("----------------------------------------------------------------"),
             h4("Session Waterfall Plots"),
             colourInput(
               "cue_color", label = 'Cue (CS)',
               palette = "limited",
               allowedCols = c(color5[2],color5[1],color5[3:4],'blue','red','black','magenta','cyan','gray',color1,color2,color3,color4,'#44B655')),
             colourInput(
               "airpuff_color", label = 'Air-puff (US)',
               palette = "limited",
               allowedCols = c(color5[3],color5[1:2],color5[4],'blue','red','black','magenta','cyan','gray',color1,color2,color3,color4,'#44B655'))
      ) # End of column function
    ) # End fluid Row
  ),# End dashboard body
  skin = 'black'
)

server = function(input,output){
  
  txt = observeEvent(input$start, {
    dataFolder = input$dataFolder
    trialDur = input$trialDur
    baselineDur = input$baselineDur
    doStartle = input$doWaterfall
    doWaterfall = input$doWaterfall
    include_session_waterfall = input$include_session_waterfall
    include_timing_plots = input$include_timing_plots
    group1_color = input$group1_color
    group2_color = input$group2_color
    cue_color = input$cue_color
    airpuff_color = input$airpuff_color
    onset_histogram = input$onset_histogram
    peak_histogram = input$peak_histogram
    postscript_flag = input$postscript_flag
    doTrajectories = input$doTrajectories
    trajRestrict = input$trajRestrict
    trajMin = input$trajMin
    restrict = input$restrict
    source('waterfalling_DISC_local.R', local = TRUE)
    
  })
  
} # End of Shiny Server function

shinyApp(ui = ui, server = server,options = list(port=7990))