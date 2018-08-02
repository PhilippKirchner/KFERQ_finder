#This script provides analysis functions for the KFERQ_finder web app


#generating the data base for the web app
#the table proteins contains the protein information (status, names, length)
#the table sequence only contains the sequence information
#the table motivs contains the motif per protein
generateDB <- function(){
  library(RSQLite)
  library(DBI)
  setwd("/Users/kirchner/Einstein_Data_Kirchner/projects/06_KFERQ_analysis")
  #open connection to the data base
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "KFERQ_finder.db")
  #generate tables for the general protein information
  DBI::dbCreateTable(con, "proteins", c(entry = "text", status = "text",
                                        protein_names = "text", gene_names = "text", length = "integer"))
  #protein sequences
  DBI::dbCreateTable(con, "sequences", c(entry = "text", sequence = "text"))
  #motif information
  DBI::dbCreateTable(con, "motifs", c(entry = "text", motif = "text", motif_start = "integer", type = "text"))
  
  #loading the protein information
  loadProteins <- function(organism){
    path <- switch(organism,
                   "human" = "human/R/Human.proteome.annotated.Rds",
                   "mouse" = "mouse/R/Mouse.proteome.annotated.Rds",
                   "rat" = "rat/R/Rat.proteome.annotated.Rds"
    )
    prot <- readRDS(path)
    prot <- prot %>% select(entry = Entry, status = Status, protein_names = `Protein names`,
                            gene_names = `Gene names`, length = Length, sequence = Sequence)
    DBI::dbAppendTable(con, "proteins", prot[,c(1:5)], row.names = NULL)
    DBI::dbAppendTable(con, "sequences", prot[,c(1,6)], row.names = NULL)
  }
  sapply(c("human","mouse","rat"),loadProteins)
  
  #loading the motif information
  setwd("/Users/kirchner/Einstein_Data_Kirchner/projects/06_KFERQ_analysis/02_web_app/KFERQ_finder")
  loadMotifs <- function(organism){
    motifs <- readRDS(paste0("data/allMotifs_",organism,".Rds"))
    DBI::dbAppendTable(con, "motifs", motifs, row.names = NULL)
  }
  sapply(c("human","mouse","rat"),loadMotifs)
  
  DBI::dbDisconnect(con)
}

#similar to the motif list the reference list of motifs needs a type column
makeModRefType <- function(){
  modRef$type <- ifelse(modRef$canonical == 1, "C",
                        ifelse(modRef$phos == 1, "P",
                               ifelse(modRef$K == 1, "K",
                                      ifelse(modRef$K_phos == 1, "KP",
                                             ifelse(modRef$N == 1, "N",
                                                    ifelse(modRef$N_phos == 1, "NP", "no motif"))))))
  saveRDS(select(modRef,motif,type),"motif_reference.Rds")
}

#The following functions are called to find motifs in the data sets or in new sequences
#-------
#Motifs from a precomputed table
#-------
findMotifs <- function(IDs, con, motifKinds){
  #necessary protein and motif information is pulled from db
  protInfo <- con %>% tbl("proteins") %>% filter(entry %in% IDs) %>% collect()
  motifInfo <- con %>% tbl("motifs") %>% filter(entry %in% IDs) %>% collect()
  #Different from the initial approach proteins that contain motifs of a different class than the chosen one 
  #are now shon as "other types of motifs"
  #Only proteins that defenitely have no motif are shown as "no motif"
  withMotif <- motifInfo %>% filter(type %in% c(motifKinds,"no motif"))
  
  otherMotif <- motifInfo %>% filter(!entry %in% withMotif$entry) %>% filter(!duplicated(entry)) %>% 
                                  mutate(motif = NA, motif_start = NA, type = "other")
  
  result <- right_join(protInfo,rbind(withMotif,otherMotif),by="entry") %>% arrange(desc(entry),motif_start)
  #To make the output easier to read the types of motifs are reconverted to the spelling in the checkboxes
  origSpell <- data.frame(type = c("C","P","K","N","KP","NP","no motif","other"), 
                          motif_type = c("canonical","phos. act.","acetyl. act.","asparagine","acetyl. and phos. act.",
                                         "asparagine and phos. act.","no motif","other types of motifs"),
                          stringsAsFactors = F)
  result <- left_join(result, origSpell, by = "type") %>% select(-type)
  
  return(result)
}

makeColorText_db <- function(IDs, con, motifKinds){
  protInfo <- con %>% tbl("proteins") %>% filter(entry %in% IDs) %>% collect()
  seqInfo <- con %>% tbl("sequences") %>% filter(entry %in% IDs) %>% collect()
  motifInfo <- con %>% tbl("motifs") %>% filter(entry %in% IDs) %>% collect()
  motifInfo <- motifInfo %>% filter(type %in% motifKinds)
  
  colorizeMotifs <- function(ID){
    prot_name <- protInfo$protein_names[protInfo$entry == ID]
    #This catches faulty IDs from uploading a file
    if (length(prot_name) == 0){
      return(data.frame(html = paste0("<p>no UniProt entry was found for: ",ID,"</p>"),stringsAsFactors = F))
      break
    }
    sequence <- seqInfo$sequence[seqInfo$entry == ID]
    motifs <- motifInfo %>% filter(entry == ID)
    
    #The header comes before the amino acid sequence
    colorHTML <- list(html = paste0("<p><b>UniProt ID:</b> ",ID," <b>protein names:</b> ",prot_name,"<br>"))
    #If there is no motif in the sequence the motifPosition is NA
    if (is.na(motifs$motif_start[1])){
      colorHTML$html[2] <- sequence
    }
    else{
      index <- 1
      for(i in seq_len(nrow(motifs))){
        tag <- switch(motifs$type[i],
                      "C" = "#ffa500", #orange
                      "P" = "#0000ff", #blue
                      "K" = "#008000", #green
                      "N" = "#ff4500", #orangeRed
                      "KP" = "#800080", #purple
                      "NP" = "#800080") #purple
        if (motifs$motif_start[i] - index <0){
          motif_start <- index
          endMotif = 4-(index - motifs$motif_start[i])
        }
        else{
          motif_start <- motifs$motif_start[i]
          endMotif = 4
        }
        colorHTML$html[i+1] <- paste0(substr(sequence,index,motif_start-1),
                                      "<font color =\"",tag,"\">",
                                      substr(sequence,motif_start,motif_start+endMotif),
                                      "</font>", collapse = "")
        index <- motifs$motif_start[i]+5
      }
      colorHTML$html[nrow(motifs)+1] <- paste0(colorHTML$html[nrow(motifs)+1],
                                               substring(sequence,first = index),"</p>",collapse = "")
    }
    return(as.data.frame(colorHTML,stringsAsFactors = F ))
  }
  htmlOutput <- list()
  for (i in seq_len(min(100,length(IDs)))){
    htmlOutput[[i]] <- colorizeMotifs(IDs[i])
  }
  return(bind_rows(htmlOutput))
}

#-------
#Motifs from a list of custom sequences (FASTA)
#-------
findMotifsFromSeq <- function(rawSeq, motifKinds, modRef){
  #The raw string is split into individual FASTA formatted sequences where possible
  #The output contains is a list of lists containing the header, sequence and identified motifs
  findMotifs <- function(x){
    output <- list()
    rawSeq <- seqList[x]
    #FASTA files are characterized by a header starting with ">" followed by a line break and the
    #text between ">" and the first "\n" is considered the header
    if (grepl("^>",rawSeq)){
      rawSeq <- unlist(strsplit(rawSeq,"\n"))
      output$header <- substring(rawSeq[1],first = 2)
      output$body <- toupper(gsub("[[:space:]]","",paste(rawSeq[-1],collapse ="")))
    }
    else{
      output$header <- ""
      output$body <- toupper(gsub("[[:space:]]","",rawSeq))
    }
    
    #Motifs are identified in the processed input sequences
    leftLimit <- 1:(nchar(output$body)-(5-1))
    rightLimit  <- 5:nchar(output$body)
    pentas <- mapply(substr,output$body,leftLimit,rightLimit,USE.NAMES=F)
    filterModRef <- modRef %>% filter(type %in% motifKinds)
    motif_start <- leftLimit[pentas %in% filterModRef$motif]
    
    if (length(motif_start) > 0){
      motifs <- data.frame(sequence_name = output$header, sequence_length = nchar(output$body), motif = pentas[motif_start], motif_start, stringsAsFactors = F)
      motifs <- left_join(motifs, modRef, by = "motif")
      origSpell <- data.frame(type = c("C","P","K","N","KP","NP"), 
                              motif_type = c("canonical","phos. act.","acetyl. act.","asparagine","acetyl. and phos. act.",
                                             "asparagine and phos. act."),
                              stringsAsFactors = F)
      output$motifs <- left_join(motifs, origSpell, by = "type") %>% select(-type)
    }
    else{
      output$motifs <- data.frame(sequence_name = output$header, sequence_length = nchar(output$body), motif = NA, motif_start = NA,
                                  motif_type = "none of the specified motifs found",  
                           stringsAsFactors = F)
    }
    
    return(output)
  }
  
  #Possible leading whitespace before the ">" needs to be removed for strsplit to work properly
  rawSeq <- gsub("[[:space:]]*>",">",rawSeq)
  #Because the leading ">" has to be conserved to be able to tell if a string contains a header or not
  #a PERL regular expression (lookforward) is used that conserves the splitting string ahead of the split string
  seqList <- unlist(strsplit(rawSeq,"(?<=.)(?=[>])",perl=T))
  return(lapply(seq_len(length(seqList)),findMotifs))
}
   
makeColorText <- function(seqsToColor){
  #This function takes the list of lists with header, sequence and motifs as input and generates
  #a list of html formatting that can be displayed to genrate the visualization of motifs
  colorizeMotifs <- function(x){
    input <- seqsToColor[[x]]
    #The header comes before the amino acid sequence
    colorHTML <- list(html = paste0("<p><b>sequence name:</b> ",input$header,"<br>"))
    #If there is no motif in the sequence the motifPosition is NA
    if (is.na(input$motifs$motif_start[1])){
      colorHTML$html[2] <- input$body
    }
    else{
      index <- 1
      for(i in seq_len(nrow(input$motifs))){
        tag <- switch(input$motifs$motif_type[i],
                      "canonical" = "#ffa500", #orange
                      "phos. act." = "#0000ff", #blue
                      "acetyl. act." = "#008000", #green
                      "asparagine" = "#ff4500", #orangeRed
                      "acetyl. and phos. act." = "#800080", #purple
                      "asparagine and phos. act." = "#800080") #purple
        if (input$motifs$motif_start[i] - index <0){
          motif_start <- index
          endMotif = 4-(index - input$motifs$motif_start[i])
        }
        else{
          motif_start <- input$motifs$motif_start[i]
          endMotif = 4
        }
        colorHTML$html[i+1] <- paste0(substr(input$body,index,motif_start-1),
                                      "<font color =\"",tag,"\">",
                                      substr(input$body,motif_start,motif_start+endMotif),
                                      "</font>", collapse = "")
        index <- input$motifs$motif_start[i]+5
      }
      colorHTML$html[nrow(input$motifs)+1] <- paste0(colorHTML$html[nrow(input$motifs)+1],
                                                     substring(input$body,first = index),"</p>",collapse = "")
    }
    return(as.data.frame(colorHTML,stringsAsFactors = F ))
  }
  #To avoid gigantic output maximally 100 sequences will be displayed this way  
  return(bind_rows(lapply(seq_len(min(length(seqsToColor),100)),colorizeMotifs)))
}    
    
    
    