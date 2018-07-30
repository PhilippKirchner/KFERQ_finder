#This script provides analysis functions for the KFERQ_finder web app

#If the list of prteins with motifs is not already present it can be generated using the following function
#This requires tables formatted like all.prot, CMA.motifs of the human data set
#Because the shiny read.table function automatically fixes whitespace in variable names some columns have to be
#renamed to include underscores instead of whitespace
makeMotifList <- function(organism){
  #The main goal of this function is to reduce the amount of data ti make loading faster
  all.prot <- readRDS(paste0(organism,"_proteome.Rds"))
  saveRDS(select(all.prot, entry = Entry, status = Status, protein_names = `Protein names`,
                 gene_names = `Gene names`, length = Length), paste0("allProt_",organism,".Rds")) 
  #The list of all motifs is trimmed down as well
  all.motifs <- readRDS(paste0(organism,"_CMA.motifs.Rds"))
  #The motif kinds are coded into character classes
  all.motifs$type <- ifelse(all.motifs$canonical ==1, "C",
                        ifelse(all.motifs$phos ==1, "P", 
                               ifelse(all.motifs$K ==1, "K", 
                                      ifelse(all.motifs$K_phos ==1, "KP",
                                             ifelse(all.motifs$N ==1, "N", 
                                                    ifelse(all.motifs$N_phos ==1, "NP", "no motif"))))))
  all.motifs <- all.motifs %>% select(entry = Entry, motif, motif_start = motifStart, type)
  #To make downstream computations a little easier proteins with no motif at all are attached to this table
  no.motifs <- data.frame(entry = all.prot$Entry[! (all.prot$Entry %in% all.motifs$Entry)], motif = NA,
                          motif_start = NA, type = "no motif")
  saveRDS(rbind(all.motifs,no.motifs),paste0("allMotifs_",organism,".RDS"))
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
findMotifs <- function(IDs, all.protCMA, all.prot, motifKinds){
  #Different from the initial approach proteins that contain motifs of a different class than the chosen one 
  #are now shon as "other types of motifs"
  #Only proteins that defenitely have no motif are shown as "no motif"
  subSet <- all.protCMA %>% filter(entry %in% IDs) 
  withMotif <- subSet %>% filter(type %in% c(motifKinds,"no motif"))
  
  otherMotif <- subSet %>% filter(!entry %in% withMotif$entry) %>% filter(!duplicated(entry)) %>% 
                                  mutate(motif = NA, motif_start = NA, type = "other")
  
  result <- right_join(all.prot,rbind(withMotif,otherMotif),by="entry") %>% arrange(desc(entry),motif_start)
  #To make the output easier to read the types of motifs are reconverted to the spelling in the checkboxes
  origSpell <- data.frame(type = c("C","P","K","N","KP","NP","no motif","other"), 
                          motif_type = c("canonical","phos. act.","acetyl. act.","asparagine","acetyl. and phos. act.",
                                         "asparagine and phos. act.","no motif","other types of motifs"),
                          stringsAsFactors = F)
  result <- left_join(result, origSpell, by = "type") %>% select(-type)
  
  return(result)
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
    
    
    