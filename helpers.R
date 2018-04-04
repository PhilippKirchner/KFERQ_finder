#This script provides analysis functions for the KFERQ_finder web app

#If the list of prteins with motifs is not already present it can be generated using the following function
#This requires the tables all.prot, CMA.motifs, and CMA.protRel
makeMotifList <- function(){
  all.protCMA <- left_join(select(all.prot, Entry, `Entry name`, Status, `Protein names`, `Gene names`, Length),
                  select(CMA.motifs, Entry, motif, motifStart, canonical, phos, K, K_phos, N, N_phos), by = "Entry")
  all.protCMA$canonical[is.na(all.protCMA$canonical)] <- 0
  all.protCMA$phos[is.na(all.protCMA$phos)] <- 0
  all.protCMA$K[is.na(all.protCMA$K)] <- 0
  all.protCMA$K_phos[is.na(all.protCMA$K_phos)] <- 0
  all.protCMA$N[is.na(all.protCMA$N)] <- 0
  all.protCMA$N_phos[is.na(all.protCMA$N_phos)] <- 0
  all.protCMA$type <- ifelse(all.protCMA$canonical ==1, "C",
                      ifelse(all.protCMA$phos ==1, "P", 
                      ifelse(all.protCMA$K ==1, "K", 
                      ifelse(all.protCMA$K_phos ==1, "KP",
                      ifelse(all.protCMA$N ==1, "N", 
                      ifelse(all.protCMA$N_phos ==1, "NP", "no motif"))))))
  #To make it possible to sort proteins by their motif content all motifs are summed up
  motifContent <- all.protCMA %>% group_by(Entry) %>% summarise(motifContent = paste0(
                                      sum(canonical, na.rm = T),"C_",
                                      sum(phos, na.rm =T),"P_",
                                      sum(K, na.rm =T),"K_",
                                      sum(K_phos, na.rm = T), "KP_",
                                      sum(N, na.rm =T),"N_",
                                      sum(N_phos, na.rm = T), "NP"))
  all.protCMA <- left_join(all.protCMA, motifContent, by = "Entry")
  return(all.protCMA)
}

#similar to the motif list the reference list of motifs needs a type column
makeModRefType <- function(){
  modRef$type <- ifelse(modRef$canonical ==1, "C",
                 ifelse(modRef$phos ==1, "P", 
                 ifelse(modRef$K ==1, "K", 
                 ifelse(modRef$K_phos ==1, "KP",
                 ifelse(modRef$N ==1, "N", 
                 ifelse(modRef$N_phos ==1, "NP", "no motif"))))))
  return(modRef)
}

#Issues to deal with:
#If I choose to only display canonical motifs are all proteins with no canonical motif considered "no motif"?
#Which kinds of motifs should be accessible? If some are locked from the user can they be selected otherwise
findMotifs <- function(IDs, dataSet, motifKinds){
  subSet <- dataSet %>% filter(Entry %in% IDs)
  #The list has to be split in two because the proteins with none of the chosen motifs would still 
  #appear as "no motif" for all other motif instances
  withMotif <- subSet %>% filter(type %in% c(motifKinds))
  noMotif <- subSet %>% filter(!(Entry %in% withMotif$Entry)) %>% mutate(motif = NA, type = "no motif") %>% filter(!duplicated(Entry))
  #Because the shiny read.table function automatically fixes whitespace in variable names some columns have to be
  #renamed to include dots instead of whitespace
  result <- rbind(withMotif, noMotif) %>% 
            select(Entry, motif, type, motifStart, Protein.names = `Protein names`, Gene.names = `Gene names`) %>% 
            arrange(desc(Entry), motifStart)
  return(result)
}

findMotifsFromSeq <- function(sequence, motifKinds, modRef){
  leftLimit <- 1:(nchar(sequence)-(5-1))
  rightLimit  <- 5:nchar(sequence)
  pentas <- mapply(substr,sequence,leftLimit,rightLimit,USE.NAMES=F)
  filterModRef <- modRef %>% filter(type %in% motifKinds)
  motifStart <- leftLimit[pentas %in% filterModRef$motif]
  
  if (length(motifStart) > 0){
    motifs <- data.frame(motif = pentas[motifStart],motifStart, sequence.length = nchar(sequence), stringsAsFactors = F)
    motifs <- left_join(motifs, select(modRef, motif, type), by = "motif")
    
    motifs$html <- "x"
    
    index <- 1
    for(i in seq_len(nrow(motifs))){
      tag <- switch(motifs$type[i],
                    "C" = "gold",
                    "P" = "blue",
                    "K" = "green",
                    "N" = "orange",
                    "KP" = "purple",
                    "NP" = "purple")
      if (motifs$motifStart[i] - index <0){
        motifStart <- index
        endMotif = 4-(index - motifs$motifStart[i])
      }
      else{
        motifStart <- motifs$motifStart[i]
        endMotif = 4
      }
      motifs$html[i] <- paste0(substr(sequence,index,motifStart-1),
                                  "<font color =\"",tag,"\">",
                                  substr(sequence,motifStart,motifStart+endMotif),
                                  "</font>", collapse = "")
      index <- motifs$motifStart[i]+5
    }
    motifs$html[nrow(motifs)] <- paste0(motifs$html[nrow(motifs)],substring(sequence,first = index),collapse = "")
    
    return(motifs)
  }
  else{
    motifs <- data.frame(motif = paste0("no ",paste0(motifKinds, collapse = ", ")," motifs found"),
                         motifStart = NA, type = NA, sequence.length = nchar(sequence),html = sequence, 
                         stringsAsFactors = F)
    return(motifs)
  }
    
}