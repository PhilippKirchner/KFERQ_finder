#testing layout options

library(shiny)
library(dplyr)

source("helpers.R")
#The current version only works with the human proteome
#Additional proteomes can be added using the "makeMotifList" function from the helpers.R script
all.protCMA.human <- readRDS("data/allMotifs_human.Rds")
all.protCMA.mouse <- readRDS("data/allMotifs_Mouse.Rds")
all.protCMA.rat <- readRDS("data/allMotifs_Rat.Rds")
#The list of all possible motifs is necessary to find motifs in unknown sequences
modRef <- readRDS("data/motif_modification_reference.Rds")

ui <- fluidPage(
  
  #This CSS is necessary to allow line breaks in the sequence which is one long word
  tags$head(
    tags$style(HTML("
                    div{word-wrap: break-word}"))
  ),
  
  navbarPage("KFERQ finder",
             tabPanel("motifs from database",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput(inputId = "organism", "select organism", choices = c("human","mouse","rat"), selected = c("human")),
                          #here selectizeInput is used to run the searches server side in r reducing crashes by 100%
                          selectizeInput(inputId = "selectEntries", "select proteins by UniProt ID", choices = NULL, multiple = T),
                          #as an alternative Uniprot IDs can be uploaded from file (this overrides any entry selection)
                          fileInput(inputId = "uploadFile", "upload UniProt IDs as csv file", 
                                    accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                          #To reduce the complexity only the canon. phos. and acetyl. motifs are searched for by default
                          #by klicking the checkbox additional choices like motifs with an asparagine instead of glutamine appear
                          checkboxInput(inputId = "expand_base", "limit motifs to canonical,\nphos. activated. and\nacetyl. activated", value = T),
                          checkboxGroupInput(inputId = "whichMotifs_base","select the kinds of motifs to include"),
                          actionButton(inputId = "go_base", label = "find motifs")
                        ),
                        mainPanel(
                          h4("motifs is selected proteins"),
                          tableOutput("motifList")
                        )
                      )
             ),
             tabPanel("motifs from sequence",
                      sidebarLayout(
                        sidebarPanel(
                          #Right now only one sequence can be analyzed at a time
                          textAreaInput(inputId = "sequence", "paste sequence as plain text", value = "", rows = 4),
                          checkboxInput(inputId = "expand_seq", "limit motifs to canonical,\nphos. activated. and\nacetyl. activated", value = T),
                          #same controls as the motif search from data base
                          checkboxGroupInput(inputId = "whichMotifs_seq","select the kinds of motifs to include"),
                          actionButton(inputId = "go_seq", label = "find motifs")
                        ),
                        mainPanel(
                          h4("list of motifs in the input sequence"),
                          wellPanel(
                            textOutput("header")),
                          wellPanel(
                            tableOutput("fromSeq")),
                          h4("motif location within the input sequence"),
                          wellPanel(
                            htmlOutput("colorSeq")))
                      )
             ),
             navbarMenu("More",
                        tabPanel("About",
                                 fluidRow(
                                   column(10, 
                                          h1("KFERQ finder app programmed in R using Shiny"),
                                          p(br(),"Version 0.64",br(),"23. 02. 2018",
                                            br(),"To do:",br(),"Change execution priority to selectizeInput > fileInput",
                                            br(),"use real motif names in the results",
                                            br(),"Export results as file",
                                            br(),"Old KFERQ finder output (colored motifs within protein sequence) also for data base searches"))
                                 )
                        )
             )
  )
)

server <- function(input, output, session){
  #server functions for finding motifs for a given UniProt entry
  #-----------------------
  #By default the motif kinds are limited to canon., phos. and acetyl.
  observe({
    if(input$expand_base == F)
      updateCheckboxGroupInput(session, "whichMotifs_base", choices = c("canonical" = "C","phos. act." = "P","acetyl.act" = "K",
                                                                        "acetyl. and phos. act." = "KP",
                                                                        "asparagine" = "N",
                                                                        "asparagine and phos. act" = "NP"),
                               selected = c("C","P","K"))
    if(input$expand_base == T)
      updateCheckboxGroupInput(session, "whichMotifs_base", choices = c("canonical" = "C","phos. act." = "P","acetyl.act" = "K"),
                               selected = c("C","P","K"))
  })
  
  dataSet <- reactive({
    return(
      switch(input$organism,
             "human" = all.protCMA.human,
             "mouse" = all.protCMA.mouse,
             "rat" = all.protCMA.rat
      ))
  })
  
  observe({
    updateSelectizeInput(session, "selectEntries", choices = dataSet()$Entry, server = T)
  })
  
  observe({
    if(!(is.null(input$uploadFile)))
      updateSelectizeInput(session, "selectEntries", selected = NULL  , server = T)
  })
  
  #There seems to be no straightforward function to revert fileInput to NULL if changes are made to
  #selectizeInput. As a result the page has to be reloaded if a file was uploaded and then Entries are
  #selected by hand
  
  motifsFromDataBase <- eventReactive(input$go_base, {
    #The script attempts to read a file
    #If no file is found the values from the selectizeInput are used
    file <- input$uploadFile
    if(is.null(file))
      Entries <- input$selectEntries
    else{
      table <- read.csv(file$datapath,  header = T)
      Entries <- table[,1]
    }
    Entries <- toupper(Entries)
    motifsFound <- findMotifs(Entries, dataSet = dataSet(), motifKinds = input$whichMotifs_base)
    invalidEntries <- Entries[!(Entries %in% motifsFound$Entry)]
    if (length(invalidEntries) > 0){
      invalidEntries <- data.frame(Entry = invalidEntries,
                                   motif = NA, type = NA, motifStart= NA,
                                   `Protein names` = "no matching Uniprot ID found",
                                   `Gene names` = "no matching Uniprot ID found")
      motifsFound <- rbind(motifsFound, invalidEntries)
    }
    return(motifsFound)
  })
  output$motifList <- renderTable(motifsFromDataBase())
  
  #server functions for finding motifs in an unknown sequence
  #--------------------
  #Like the motifs from data base the choice of motif types is limited by default
  observe({
    if(input$expand_seq == F)
      updateCheckboxGroupInput(session, "whichMotifs_seq", choices = c("canonical" = "C","phos. act." = "P","acetyl.act" = "K",
                                                                       "acetyl. and phos. act." = "KP",
                                                                       "asparagine" = "N",
                                                                       "asparagine and phos. act" = "NP"),
                               selected = c("C","P","K"))
    if(input$expand_seq == T)
      updateCheckboxGroupInput(session, "whichMotifs_seq", choices = c("canonical" = "C","phos. act." = "P","acetyl.act" = "K"),
                               selected = c("C","P","K"))
  })
  
  #The updating and rendering functions have to be kept separate to work
  #I don't know if this behaviour is intentional
  
  motifsFromSeq <- eventReactive(input$go_seq,{
    #entered text is automatically converted to upper case
    rawSeq <- unlist(strsplit(input$sequence,"\n"))
    if (substr(rawSeq[1],1,1) == ">"){
      #in this case the input is a FASTA format and the fist block is the header
      header <- substring(rawSeq[1],first = 2)
      #all other blocks contain the sequence, possible whitespace is removed
      body <- toupper(gsub("[[:space:]]","",paste(rawSeq[-1],collapse ="")))
    }
    else{
      header <- ""
      body <- toupper(gsub("[[:space:]]","",paste(rawSeq,collapse ="")))
    }  
    return(list(header = header, motifs = findMotifsFromSeq(body, input$whichMotifs_seq, modRef = modRef)))
  })
  
  
  output$header <- renderText({paste0("sequence name: ",motifsFromSeq()$header)})
  
  output$fromSeq <- renderTable({
    #if no motif was found the table contains an error message
    motifsFromSeq()$motifs[,c(1:4)]
  })
  
  output$colorSeq <- renderUI({
    HTML(paste0("<p>",paste0(motifsFromSeq()$motifs$html,collapse=""),"</p>",collapse=""))
  })
}

shinyApp(ui = ui, server = server)