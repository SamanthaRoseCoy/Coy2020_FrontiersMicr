#Code to make Figure 4

##Need three files from all_motifs.R to combine together. I named these low30x_1C, low30x_2A, and low30x_E1.  
##The naming should refelct the coverage and files you analyzed.

###Combine pacbio files into one master file using all_motifs 
low30x_1C <- all_motifs("PBCV1-1C_CMD_30C_Motifs.gff", "PBCV1-1C_CMD_30C_Modifications.csv", "PBCV1-1C_CMD_30C_Motifs.csv")
low30x_2A <- all_motifs("PBCV1-2A_CMD_30C_Motifs.gff", "PBCV1-2A_CMD_30C_Modifications.csv", "PBCV1-2A_CMD_30C_Motifs.csv")
low30x_E1 <- all_motifs("PBCV1-E1_CMD_30C_Motifs.gff", "PBCV1-E1_CMD_30C_Modifications.csv", "PBCV1-E1_CMD_30C_Motifs.csv")


###Extract the information of interest for methylFrac analysis
methylFrac_30X_R1 <- low30x_1C %>%
  filter(motif == "GATC"|motif == "CATG") %>%
  mutate(Sample = coverage) %>%
  select(tpl, Sample, frac, motif)

methylFrac_30X_R2 <- low30x_2A %>%
  filter(motif == "GATC"|motif == "CATG") %>%
  mutate(Sample = coverage) %>%
  select(tpl, Sample, frac, motif)

methylFrac_30X_R3 <- low30x_E1 %>%
  filter(motif == "GATC"|motif == "CATG") %>%
  mutate(Sample = coverage) %>%
  select(tpl, Sample, frac, motif)

###Label the replicates in the data (coverage column will be converted)
methylFrac_30X_R1$Sample[methylFrac_30X_R1$Sample == "30"] <- "Rep1"
methylFrac_30X_R2$Sample[methylFrac_30X_R2$Sample == "30"] <- "Rep2"
methylFrac_30X_R3$Sample[methylFrac_30X_R3$Sample == "30"] <- "Rep3"

###bind the three replicate data together and calculate mean, stdev
allreps_30X <- rbind(methylFrac_30X_R1,methylFrac_30X_R2,methylFrac_30X_R3) %>%
  group_by(tpl) %>%
  mutate(AverageFrac=mean(frac), SDFrac=sd(frac)) %>%
  filter(row_number() ==1) %>% #removes duplicate calculations
  arrange()

###graph it
plot_fig4 <- ggplot(allreps_30X, aes(x=AverageFrac, y=SDFrac, color=motif)) +
  geom_point() +
  labs(x="methylFrac (\U003BC, n=3)", y="methylFrac (\u03C3, n=3)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_colour_manual(na.translate=FALSE, values=c("CATG"="#0000CC", "GATC"="#FF9900")) +
  theme(panel.margin = unit(2, "lines")) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept=0.3) +
  geom_vline(xintercept=0.35)
ggMarginal(plot_fig5, type = "histogram", fill = "transparent") #Need ggExtra
