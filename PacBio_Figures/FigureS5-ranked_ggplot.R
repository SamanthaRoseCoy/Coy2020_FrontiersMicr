###Combine pacbio files into one master file using all_motifs 
low30x_1C <- all_motifs("PBCV1-1C_CMD_30C_Motifs.gff", "PBCV1-1C_CMD_30C_Modifications.csv", "PBCV1-1C_CMD_30C_Motifs.csv")

#Function to make ranked modQV plot that is part of Figure S5
#' @param bulk_methyl percent used to set a threshold (e.g. percent of methylated adenines determined by LC-MS)
#' 
#' @return returns a rank-abundance curve of the modQV scores of your chosen nucleotide to analyze, color coded by motif association.  There is an option to set a threshold based on known values.
#'
#' @example
#' rankedmodQV_plot(low30x_1C, "A", 0.015)
#'
#' @export

##Must write base_type in quotes for some reason (e.g. "A")
rankedmodQV_plot <- function(df, base_type, bulk_methyl){
  rankedScore <- df %>%
    filter(df$base==base_type) %>%
    group_by(motif) %>%
    arrange(desc(score)) %>%
    mutate(id = row_number())
  
  #If else statement to decide how to graph the values
  if(missing(bulk_methyl)){
    #Graph this
    rankedplot_modQV <- ggplot(rankedScore, aes(y=score, x=id, color=motif)) +
      geom_point() +
      scale_y_continuous(breaks = seq(0,80, by =20)) +
      facet_grid(.~motif, scales="free") +
      labs(x="", y="ModificationQV") +
      theme(axis.text.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(1, "lines")) +
      scale_colour_manual(na.translate=TRUE, values=c("CATG"="#0000CC", "GATC"="#FF9900", "GNNNNVNH"="#FF6666", "CNNNNRNH"="#00CCCC"), na.value="#999999") +
      theme(axis.text=element_text(size=12, face="bold")) +
      theme(axis.title=element_text(size=12, face= "bold")) +
      theme(strip.text = element_text(face="bold", size=12))
  }
  else{
    
    #What is the name of your threshold
    prompt <- cat("Choose a threshold name: ")
    input <- readline(prompt=prompt)
    
    #Counts the number of Adenine sites...not sure about evaluation format
    total_basetype <-nrow(rankedScore)
    
    #Identifies the highest ipdR values based on batch methylation measurements of this base
    highIPDs <- total_basetype * bulk_methyl
    nhighIPDs <- as.integer(highIPDs)
    
    #Identifies the threshold
    value <- top_n(x=rankedScore, n=nhighIPDs, wt=rankedScore$score)
    threshold <- min(value$score)
    
    #Graph this
    rankedplot_modQV <- ggplot(rankedScore, aes(y=score, x=id, color=motif)) +
      geom_point() +
      scale_y_continuous(breaks = seq(0,80, by =20)) +
      facet_grid(.~motif, scales="free") +
      labs(x="", y="ModificationQV") +
      theme(axis.text.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(1, "lines")) +
      scale_colour_manual(na.translate=TRUE, values=c("CATG"="#0000CC", "GATC"="#FF9900", "GNNNNVNH"="#FF6666", "CNNNNRNH"="#00CCCC"), na.value="#999999") +
      theme(axis.text=element_text(size=12, face="bold")) +
      theme(axis.title=element_text(size=12, face= "bold")) +
      theme(strip.text = element_text(face="bold", size=12)) +
      geom_hline(aes(yintercept=threshold, linetype = input)) +
      scale_linetype_manual(name="threshold", values=threshold) 
    
  }
  
  return(rankedplot_modQV)
}

