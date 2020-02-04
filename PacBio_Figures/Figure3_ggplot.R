#Code to make Figure 3
plot_fig3 <- ggplot(low30x_1C, aes(x=seq_along(low30x_1C$ipdRatio), y=ipdRatio, color=motif)) +
  geom_point() +
  facet_grid(.~base, scales="free_x") +
  labs(x="",y="ipdRatio") + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"),legend.text=element_text(size=12),legend.title=element_blank()) + 
  theme(panel.border = element_rect(fill=NA, colour = "black"), strip.background=element_blank(), axis.text.x=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), panel.background=element_blank()) +
  scale_colour_manual(na.translate=TRUE, values=c("CATG"="#0000CC", "GATC"="#FF9900", "GNNNNVNH"="#FF6666", "CNNNNRNH"="#00CCCC"), na.value="#999999") +
  theme(axis.text=element_text(size=12, face="bold")) +
  theme(axis.title=element_text(size=12, face= "bold")) +
  theme(strip.text = element_text(face="bold", size=12))
