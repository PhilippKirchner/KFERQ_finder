#V0.8
#20180730

#The following packages are required (for increased robustness the whole tidyverse could be loaded)
library(shiny)
library(dplyr)
library(readr)
#shinyjs allows updating widgets
library(shinyjs)
#functions to get the motif data are located in helpers.R 
source("helpers.R")

#For the time being only the human data set is loaded as part of the regular start up routine
#all.protCMA.human <- readRDS("data/allMotifs_human.Rds")
#all.prot.human <- readRDS("data/allProt_human.Rds")

#The list of all possible motifs is necessary to find motifs in unknown sequences
modRef <- readRDS("data/motif_reference.Rds")

ui = fluidPage(
  useShinyjs(),
  includeCSS("styles.css"),
  
  navbarPage("KFERQ finder V0.8",
             tabPanel("motifs from data base",
                      fluidPage(
                        fluidRow(
                          column(5,
                                 wellPanel(
                                   class = "headBox",
                                   h4(style = "margin-top:0;", "input options"),
                                   h5("select Entries from databse"),
                                   fluidRow(
                                     column(10, offset = 2,
                                            selectInput(
                                              inputId = "organism",
                                              "select organism",
                                              choices = c("human", "mouse", "rat"),
                                              selected = c("human")
                                            ),
                                            #here selectizeInput is used to run the searches server side in r reducing crashes by 100%
                                            selectizeInput(
                                              inputId = "selectEntries",
                                              "select proteins by UniProt ID",
                                              choices = NULL,
                                              multiple = T
                                            )
                                     )
                                   ),
                                   br(),
                                   h5("or upload UniProt IDs as csv file"),
                                   fluidRow(
                                     column(10, offset = 2,
                                            #as an alternative Uniprot IDs can be uploaded from file (this overrides any entry selection)
                                            fileInput(
                                              inputId = "uploadFile_db",
                                              label = NULL,
                                              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                                            )
                                     )
                                   )
                                 )),
                          column(5,
                                 wellPanel(
                                   class = "headBox",
                                   h4(style = "margin-top:0;","output options"),
                                   #In the initial version there was a checkbox whether to include more hypothetical motifs as well
                                   #I changed to this layout because I think it gives more info why some  motifs are not checked by default
                                   h5("select the kinds of motifs to include"),
                                   fluidRow(
                                     column(
                                       11,
                                       offset = 1,
                                       checkboxGroupInput(
                                         inputId = "stdMotifs_base",
                                         label = "standard motifs",
                                         choices = c("canonical" = "C", "phos. act." = "P", "acetyl.act." = "K"),
                                         selected = c("C", "P", "K")
                                       )
                                     )
                                   ),
                                   fluidRow(
                                     column(
                                       10,
                                       offset = 2,
                                       checkboxGroupInput(
                                         inputId = "advMotifs_base",
                                         label = "advanced motifs",
                                         choices = c(
                                           "acetyl. and phos. act." = "KP",
                                           "asparagine" = "N",
                                           "asparagine and phos. act." = "NP"
                                         ),
                                         selected = NULL
                                       )
                                     )
                                   ),
                                   checkboxInput(inputId = "showUb_db", "show inactivation through ubiquitination")
                                 )),
                          column(2)
                        ),
                        fluidRow(
                          column(3,
                                 actionButton(
                                   class = "goButton",
                                   inputId = "go_db",
                                   label = "find motifs"
                                 )
                          ),
                          column(9)
                        ),
                        fluidRow(
                          br()
                        ),
                        fluidRow(column(12,
                                        wellPanel(class = "bottomBox",
                                                  fluidRow(
                                                    column(10,
                                                           h4("motifs in the selected proteins:")),
                                                    column(2,
                                                           downloadButton("downloadTable_db", "download"))
                                                  ),
                                                  fluidRow(column(
                                                    12,
                                                    wellPanel(class = "resultsBox",
                                                              tableOutput("motifList"))
                                                  ))
                                                  
                                        ))
                        )
                      )
             ),
             tabPanel(
               "motifs from sequence",
               fluidPage(
                 fluidRow(
                   column(5,
                          wellPanel(
                            class = "headBox",
                            h4(class = "IOHeader", "input options"),
                            h5("paste amino acid sequences"),
                            textAreaInput(
                              inputId = "sequence",
                              label = NULL,
                              value = "",
                              rows = 6
                            ),
                            br(),
                            h5("or upload sequences as plain text (.txt)"),
                            fileInput(
                              inputId = "uploadFile_seq",
                              label = NULL,
                              accept = c("text/plain", ".txt")
                            )
                          )),
                   column(5,
                          wellPanel(
                            class = "headBox",
                            h4(class = "IOHeader", "output options"),
                            #In the initial version there was a checkbox whether to include more hypothetical motifs as well
                            #I changed to this layout because I think it gives more info why some  motifs are not checked by default
                            h5("select the kinds of motifs to include"),
                            fluidRow(
                              column(
                                11,
                                offset = 1,
                                checkboxGroupInput(
                                  inputId = "stdMotifs_seq",
                                  label = "standard motifs",
                                  choices = c(
                                    "canonical" = "C",
                                    "phos. act." = "P",
                                    "acetyl.act." = "K"
                                  ),
                                  selected = c("C", "P", "K")
                                )
                              )),
                            fluidRow(
                              column(
                                10,
                                offset = 2,
                                checkboxGroupInput(
                                  inputId = "advMotifs_seq",
                                  label = "advanced motifs",
                                  choices = c(
                                    "acetyl. and phos. act." = "KP",
                                    "asparagine" = "N",
                                    "asparagine and phos. act." = "NP"
                                  ),
                                  selected = NULL
                                )
                              )),
                            checkboxInput(
                              inputId = "showUb_seq", 
                              label = "show inactivation through ubiquitination")
                          )
                   ),
                   column(2)
                 ),
                 fluidRow(
                   column(3,
                          actionButton(
                            class = "goButton",
                            inputId = "go_seq",
                            label = "find motifs"
                          )
                   ),
                   column(9)
                 ),
                 fluidRow(
                   br()
                 ),   
                 fluidRow(
                   column(12,
                          tabsetPanel(
                            tabPanel(
                              "table output",
                              column(
                                12,
                                class = "tabPaneContent",
                                br(),
                                fluidRow(column(10,
                                                h4(
                                                  "motifs in the selected proteins:"
                                                )),
                                         column(
                                           2,
                                           downloadButton("downloadTable_seq", "download")
                                         )),
                                wellPanel(class = "resultsBox",
                                          tableOutput("fromSeq"))
                              )
                            ),
                            tabPanel(
                              "sequence output",
                              column(
                                12,
                                class = "tabPaneContent",
                                br(),
                                fluidRow(column(10,
                                                h4(
                                                  "motifs in the selected proteins:"
                                                )),
                                         column(2
                                                #This button does not yet work
                                                #I have to find out how to get the formatted text out of shiny
                                                #downloadButton("downloadColor_seq", "download"))),)),)),
                                         )),
                                wellPanel(class = "resultsBox",
                                          htmlOutput("colorSeq"))
                              )
                            )
                          )
                   )
                 )
               )
             ),
             tabPanel("about",
                      fluidRow(
                        column(
                          10,
                          h4("KFERQ finder app programmed in R using Shiny"),
                          br(),
                          column(
                            11,
                            offset = 1,
                            p(
                              "Version 0.8  2018 07 13",
                              br(),
                              br(),
                              "CMA-targeting (KFERQ-like) targeting motifs can be searched
                              in all human, mouse or rat proteins in UniProt",
                              br(),
                              "KFERQ-like motifs belong to different classes based on their amino acid composition.
                              For example, motifs with a serine, threonine or tyrosine can be grouped together as
                              phosphorylation-activated motifs and lysine instead of the terminal glutamine
                              belong to the group of acetylation-activated motifs.",
                              br(),
                              "There are more advanced groups of motifs combining different modification that are disabled
                              by default because less experimental evidence is available for their function in CMA.",
                              br(),
                              br()
                            )
                          ),
                          br(),
                          h4("1. motifs from database"),
                          column(
                            11,
                            offset = 1,
                            p(
                              h4("input options"),
                              "UniProt identifiers can be selected from a drop down list for the three available organisms
                              (human, mouse, rat). As an alternative, UniProt identifiers from these species can be uploaded
                              as comma-separated (.csv) tables for batch processing. 
                              The input file is assumed to contain a column header",
                              br(),
                              h4("output options"),
                              "Motif classes can be selected to be searched for in the input proteins. By default the most
                              well characterized classes are selected. If a protein contains none of the selected classes but
                              of other classes this will be mentioned in the results table.",
                              br(),
                              "KFERQ-like motifs may be inactivated through ubiquitination. Possible ubiquitination sites in
                              KFERQ-like motifs can be shown in the output table.",
                              br(),
                              h4("results"),
                              "The motifs identified and some information such as the protein names are displayed in the output.
                              The results table can be downloaded as comma-separated (.csv) file.",
                              br()
                            )
                          ),
                          br(),
                          h4("2. motifs from sequence"),
                          column(
                            11,
                            offset = 1,
                            p(
                              h4("input options"),
                              "To identify motifs in amino acid sequences not available in the UniProt data sets FASTA-style
                              sequences can be analyzed directly.",
                              br(),
                              "Lines starting with a \">\" are interpretated as headers containinf a sequence description.
                              Multiple sequences should be separated by a line break. Sequenes can be copy-pasted directly
                              into the text input field or uploaded as plain text (.txt) files.",
                              br(),
                              h4("output options"),
                              "The same output options are used as when identifying motifs from data base entries.",
                              br(),
                              h4("results"),
                              "In addition to the table output format motifs can also be displayed as colored text within the
                              amino acid sequence. This output format is limited to 100 sequences. Longer input 
                              should be split accordingly",
                              br(),
                              "RIGHT NOW THIS VISUALIZATION CANNOT BE DIRECTLY DOWNLOADED. HOWEVER IT IS POSSIBLE TO COPY-PASTE
                              IT INTO A TEXT EDITOR."
                            )
                          )
                          )
                        ))
             )
  
  )

server = function(input, output, session) {
  #This makes quits the app once the browser window is closed
  session$onSessionEnded(stopApp)

  #server functions for finding motifs for a given UniProt entry
  #-----------------------
  #The motif kinds are put together from the standard and advanced motif checkboxes
  whichMotifs_base <- reactive({
    return(
      c(input$stdMotifs_base,input$advMotifs_base))
  })
  
  #The current version works with the human mousa and rat proteome
  #Additional proteomes can be added using the "makeMotifList" function from the helpers.R script
  #In the future a relational db could help speeding things up a little
  #To save time data sets are only loaded when needed (the option "human" is enabled by default)
  all.protCMA <- reactive({
    return(
      switch(input$organism,
             "human" = readRDS("data/allMotifs_human.Rds"),#all.protCMA.human,
             "mouse" = readRDS("data/allMotifs_mouse.Rds"),#all.protCMA.mouse,
             "rat" = readRDS("data/allMotifs_rat.Rds")#all.protCMA.rat
      ))
  })
  
  all.prot <- reactive({
    return(
      switch(input$organism,
             "human" = readRDS("data/allProt_human.Rds"),#all.prot.human,
             "mouse" = readRDS("data/allProt_mouse.Rds"),#all.prot.mouse,
             "rat" = readRDS("data/allProt_rat.Rds")#all.prot.rat
      ))
  })
  
  #If the all.prot data set changes the available UniProt IDs are updated  
  observe({
    updateSelectizeInput(session, "selectEntries", choices = all.prot()$entry, server = T)
  })
  #Entries are saved in different places depending if they come from the drop down menu or
  #file upload (This is necessary to allow for both input methods to be mutually exclusive)
  Entries <- reactiveValues(
    select = NULL,
    file = NULL
  )
  #The two observe functions reset one input method, if the other is chosen (requires shinyjs)
  observe({
    Entries$select <- input$selectEntries
    if(!is.null(Entries$select)){
      reset("uploadFile_db")
      Entries$file <- NULL
    }
    
  })
  observe({
    req(input$uploadFile_db)
    Entries$file <- read.csv(input$uploadFile_db$datapath, header = T)[,1]
    reset("selectEntries")
  })
  
  
  motifsFromDataBase <- eventReactive(input$go_db, {
    if(is.null(Entries$file)){
      req(Entries$select)
      Entries <- input$selectEntries
    }
    else{
      #Uniprot identifiers have to be upper case
      Entries <- toupper(Entries$file)
    }
    motifsFound <- findMotifs(Entries, all.protCMA = all.protCMA(), all.prot = all.prot(),
                              motifKinds = whichMotifs_base())
    invalidEntries <- Entries[!(Entries %in% all.prot()$entry)]
    if (length(invalidEntries) > 0){
      invalidEntries <- data.frame(entry = invalidEntries, status = NA, protein_names = "no matching Uniprot ID found",
                                   gene_names = "no matching Uniprot ID found", length = NA,
                                   motif = NA, motif_start = NA, motif_type = NA)
      motifsFound <- rbind(motifsFound, invalidEntries)
    }
    
    
    if(input$showUb_db){
      motifsFound$ubi.inactivated <- ifelse(grepl("K",motifsFound$motif),TRUE,ifelse(is.na(motifsFound$motif),NA,FALSE))
    }
    return(motifsFound)
  })
  #output the results in the website
  output$motifList <- renderTable(
    motifsFromDataBase()
  )
  #download results to file
  output$downloadTable_db <- downloadHandler(
    filename = function(){
      paste0("KFERQ_finder_table_",format(Sys.time(), "%Y%m%d_%H%M"),".csv",sep ="")
    },
    content = function(file){
      write.csv(motifsFromDataBase(),file,row.names = F)
    }
  )
  
  #server functions for finding motifs in an unknown sequence
  #--------------------
  #The motif kinds are put together from the standard and advanced motif checkboxes
  whichMotifs_seq <- reactive({
    return(
      c(input$stdMotifs_seq,input$advMotifs_seq))
  })
  #sequences are saved in different places depending if they are entered manually or
  #uploades as a file (This is necessary to allow for both input methods to be mutually exclusive)
  FASTA <- reactiveValues(
    manual = NULL,
    file = NULL
  )
  #The two observe functions reset one input method, if the other is chosen (requires shinyjs)
  observe({
    FASTA$manual <- input$sequence
    if(nchar(FASTA$manual)!=0){
      reset("uploadFile_seq")
      FASTA$file <- NULL
    }
    
  })
  observe({
    req(input$uploadFile_seq)
    FASTA$file <- read_file(input$uploadFile_seq$datapath)
    reset("sequence")
  })
  #The unknown sequences are handed over to a function that splits the sequences 
  #and returns a list object that contains the header, body (the actual amino acid sequences) and identified motifs
  motifsFromSeq <- eventReactive(input$go_seq,{
    if(is.null(FASTA$file)){
      req(FASTA$manual)
      inputSeq <- FASTA$manual
    }
    else{
      inputSeq <- FASTA$file
    }
    findMotifsFromSeq(inputSeq, whichMotifs_seq(), modRef = modRef)
    })
  #dsplay the results file in the website
  output$fromSeq <- renderTable({
    #usually I would use tata.table::rbindlist but there is a similar function in deplyr::bind_rows
    bind_rows(lapply(seq_len(length(motifsFromSeq())),function(x){motifsFromSeq()[[x]]$motifs}))
  })
  #download the results table as file
  output$downloadTable_seq <- downloadHandler(
    filename = function(){
      paste0("KFERQ_finder_from_sequences_",format(Sys.time(), "%Y%m%d_%H%M"),".csv",sep ="")
    },
    content = function(file){
      write.csv(bind_rows(lapply(seq_len(length(motifsFromSeq())),function(x){motifsFromSeq()[[x]]$motifs})),
                file, row.names = F)
    }
  )
  #This function may take quite a while to execute for longer sequences
  #to avoid trouble with this kind of function maximally 100 sequences will be displayed
  output$colorSeq <- renderUI({
      HTML(paste0(makeColorText(motifsFromSeq())$html,collapse=""))
  })
}

shinyApp(ui = ui, server = server)