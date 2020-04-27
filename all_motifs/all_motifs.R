#' Import and merge three Pacbio outputs from the RS_Modification_and_Motif Analysis Protocol for data analysis.
#'
#'Requires tidyverse, Biostrings, reshape
#'
#' @param file1 the file path for the motifs.gff file in quotes
#' @param file2 the file path for the modifications.csv file in quotes
#' @param file3 the file path for the motifs.csv file in quotes
#' 
#' @return Returns a dataframe containing the full genome data detailing modification type (feature) and all possible motifs (motif).
#'
#' @example
#' all_motifs("Motifs.gff", "Modifications.csv", "Motifs.csv")
#'
#' @export

all_motifs <- function(file1, file2, file3){
  #Reads the motif file and removes spaces
  Rawfile1 <- read.table(file1, 
                         skip = 4,
                         sep = "\t",
                         na.strings = c("","NA"))
  Rawfile2 <- read.csv(file2, 
                       na.strings = c("","NA"))
  
  Rawfile3 <- read.csv(file3, 
                       stringsAsFactors = F,
                       na.strings = c("","NA"))
  
  #converts the motifs.gff dataframe into a tibble and selects important columns
  #Newfile1 <- tibble::as_tibble(select(Rawfile1, V3, V4, V7, V9)) #getting warnings V7
  Newfile1 <- select(Rawfile1, V3, V4, V7, V9)
  #converts strand column to characters
  Newfile1$V7 <- as.character(Newfile1$V7)
  #change forward/reverse denotation to 0, 1, respectively
  Newfile1$V7[Newfile1$V7 == "+"] <- "0"
  Newfile1$V7[Newfile1$V7 == "-"] <- "1"
  #converts strand column to integers
  Newfile1$V7 <- as.integer(Newfile1$V7)
  #change the column names in Newfile1
  colnames(Newfile1) <- c("feature","tpl","strand","attributes")
  #merge the edited motifs.gff file and modifications.csv dataframes by position, strand
  merged_file <- merge(Newfile1, Rawfile2, by=c("tpl", "strand"), all=T)
  #reorder the merged dataframe
  final_merge <- merged_file[,c(5,1,2,3,6,7,11,12,13,14,15,4)]
  final_merge <- as.data.frame(final_merge, stringsAsFactors=FALSE)
  
  #bin file by strand
  forward <- final_merge[final_merge[,3]!=1,]
  reverse <- final_merge[final_merge[,3]!=0,]
  
  #prompt for use input in deciding what motifs to use
  Rawfile3$to.keep <- NA # creates a new column to fill 
  for(i in 1:nrow(Rawfile3)) { # check out seq_along syntax - it is safer
    prompt <- cat("Do you want to analyze motif ", Rawfile3$motifString[i], "? (y/n) " )
    input <- readline(prompt=prompt)
    if(input == "y") { # THIS WOULD BE MORE ELEGANT WITH A SWITCH STATEMENT
      # Mark row as to be kept
      Rawfile3$to.keep[i] <- TRUE
    } else if(input == "n") {
      Rawfile3$to.keep[i] <- FALSE
    } else {
      warning(paste("The input at row", i, "was not y or n; I'll assume you do not want to analyze this one."))
      Rawfile3$to.keep[i] <- FALSE
    }
  }
  
  # Throw out motif strings we don't care about
  good_seqs <- Rawfile3[Rawfile3$to.keep == TRUE, ]
  
  
  #Create forward and reverse strings to search motifs for using biostrings xstring
  forward_string <- DNAString(paste((forward$base), collapse="")) #returns a single genome string "GGGAGAA..."
  reversed <- arrange(reverse, desc(tpl)) #rearranges the dataframe in descending order of tpl
  reverse_string <- DNAString(paste((reversed$base), collapse="")) #returns a single genome string "AAAAGTTT..."
  
  #initialize lists to store results for all iterations of the following loop  
  results_forward <- NULL
  results_reverse <- NULL
  
  # Identify all potential positions for motifs selected for analysis.. right now it is only keeping values of the last iteration.
  for(i in 1:nrow(good_seqs)) {
    # Process everything for the forward string
    motifs_f <- good_seqs$motifString #creates a vector for the motifString by element
    hits_for <- matchPattern(motifs_f[i], forward_string, fixed=FALSE) #Searches matches on the DNAString subject (test) with a column for start, end, width.  Default sets to match sequence.
    hits_for <- start(hits_for) #Extract the start column to yield a vector of integer positions
    hits_for_position <- hits_for + good_seqs$centerPos[i] #adjusts the list of position by centerPos data
    results_forward <- rbind(results_forward, data.frame(hits_for_position, motifs_f[i])) #appends the results of each iteration to results
    
    #Do the same for the reverse string
    motifs_r <- good_seqs$motifString
    hits_rev <- matchPattern(motifs_r[i], reverse_string, fixed=FALSE)
    hits_rev <- start(hits_rev)
    hits_rev_position <- hits_rev + good_seqs$centerPos[i]
    results_reverse <- rbind(results_reverse, data.frame(hits_rev_position, motifs_r[i])) # There is an issue here... I need to adjust for 
  }
  
  
  #Convert colnames for merging dataframes
  colnames(results_forward) <- c("tpl", "motif") #converts column name
  colnames(results_reverse) <- c("revtpl", "motif") #converts column name
  
  #Merge the forward file to get the full thing
  forward_final <- left_join(forward, results_forward, by=c("tpl"))
  
  #Merge the reverse files to get the correct tpl
  results_reverse$revtpl <- as.character(results_reverse$revtpl) #coerce to character for left_join merge
  reverse_final <- dplyr::left_join(reversed %>%
                                      mutate(revtpl = rownames(reversed)),
                                    results_reverse,
                                    by = "revtpl")
  reverse_final$revtpl <- NULL #removes revtpl row
  reverse_final <- arrange(reverse_final, tpl)
  
  #add the forward and reverse strands back together
  all_motifs <- rbind(forward_final, reverse_final)
  all_motifs <- arrange(all_motifs, tpl)
  all_motifs <- all_motifs[,c(1,2,3,4,5,6,12,7,8,9,10,11,13)]
  #all_motifs <- as.data.frame(all_motifs)
  
  #change the attributes info
  all_motifs$attributes <- gsub("(.*)(identificationQv=)(.*)", "\\3", all_motifs$attributes)
  all_motifs$attributes <- gsub("(co)(.*)", NA, all_motifs$attributes)
  colnames(all_motifs)[7] <- "score2"
  
  all_motifs$frac[is.na(all_motifs$frac)] <- 0
  
  return(all_motifs)
}
