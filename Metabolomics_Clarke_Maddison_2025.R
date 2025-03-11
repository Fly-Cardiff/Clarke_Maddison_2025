# Metabolomics visualizing code 

# Knockdown Groups
#C Control LacZ KD exp
#D KD LacZ KD exp
#G Con_AB KD exp
#H KD_AB KD exp 

# OE groups 
#I Con_LacZ OE exp
#J OE_LacZ OE exp
#K Con AB OE exp
#L OE AB OE exp 

library("ggplot2")
library("ggrepel")
################################# JUST KNOCKDOWN EXPERIMENT ##########################################

KD_All <- read.csv('/MetabolomicsData/2070_peaks_list_KD_dataset.csv')
KD_200 <- read.csv('/MetabolomicsData/Two_way_ANOVA_2720_KD_top_200.csv');

# Main effect of gene (C+G (control level) Vs. D+H) 
NoWwoxKD_Avg <- rowMeans(data.frame(KD_All$C1, KD_All$C2, KD_All$C3, KD_All$C4, KD_All$C5, KD_All$G1, KD_All$G2, KD_All$G3, KD_All$G4, KD_All$G5), na.rm=TRUE)
WwoxKD_Avg <- rowMeans(data.frame(KD_All$D1, KD_All$D2, KD_All$D3, KD_All$D4, KD_All$D5, KD_All$H1, KD_All$H2, KD_All$H3, KD_All$H4, KD_All$H5), na.rm=TRUE)
Foldchange <- WwoxKD_Avg[]/NoWwoxKD_Avg[];
LogFold <- log2(Foldchange)
index <- KD_200$X;
count = 0;
GeneMainEffect_inKD <- numeric(200)
  
for (x in index) {
  count = count+1;
  GeneMainEffect_inKD[count] <- LogFold[x]
}

# Main effect of AB (C+D (control level) Vs. G+H)
LacZ_Avg <- rowMeans(data.frame(KD_All$C1, KD_All$C2, KD_All$C3, KD_All$C4, KD_All$C5,KD_All$D1, KD_All$D2, KD_All$D3, KD_All$D4, KD_All$D5), na.rm=TRUE)
AB_Avg <- rowMeans(data.frame(KD_All$G1, KD_All$G2, KD_All$G3, KD_All$G4, KD_All$G5,KD_All$H1, KD_All$H2, KD_All$H3, KD_All$H4, KD_All$H5), na.rm=TRUE)
Foldchange2 <- AB_Avg[]/LacZ_Avg[];
LogFold2 <- log2(Foldchange2)
count2 = 0;
ABMainEffect_inKD <- numeric(200)

for (x in index) {
  count2 = count2+1;
  ABMainEffect_inKD[count2] <- LogFold2[x]
}

# Gene in AB (G (control level) Vs. H)
G_Avg <- rowMeans(data.frame(KD_All$G1, KD_All$G2, KD_All$G3, KD_All$G4, KD_All$G5), na.rm=TRUE)
H_Avg <- rowMeans(data.frame(KD_All$H1, KD_All$H2, KD_All$H3, KD_All$H4, KD_All$H5), na.rm=TRUE)
Foldchange3 <- H_Avg[]/G_Avg[];
LogFold3 <- log2(Foldchange3)
count3 = 0;
ABGeneInteraction_inKD <- numeric(200)

for (x in index) {
  count3 = count3+1;
  ABGeneInteraction_inKD[count3] <- LogFold3[x]
}


#VOLCANOES IN KNOCKDOWN 

#Plotting Volcano for gene main effect
plotdf <- data.frame(GeneMainEffect_inKD, KD_200$Gene.adj.p.)
plotdf$diffex<-"Non-Sig"
plotdf$diffex[GeneMainEffect_inKD > 0 & KD_200$Gene.adj.p. < 0.1] <- "Up"
plotdf$diffex[GeneMainEffect_inKD < 0 & KD_200$Gene.adj.p. < 0.1] <- "Down"

# Min and max for scaling 
max(plotdf$GeneMainEffect_inKD, na.rm=T) 
min(plotdf$GeneMainEffect_inKD, na.rm=T) 
max(-log10(KD_200$Gene.adj.p.), na.rm=T)
min(-log10(KD_200$Gene.adj.p.), na.rm=T) 


p <- ggplot(plotdf, aes(x=GeneMainEffect_inKD, y=-log10(KD_200$Gene.adj.p.), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-5,5) + ylim(0,10)
p2<- p + geom_hline(yintercept=-log10(0.1), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("#D81B60", "black", "#1E88E5"))
p4 <- p3 + labs(y="-Log10 Padj",x="Log2 Fold Change") + theme(axis.line=element_line(linewidth =1), axis.title = element_text(family = "Arial", color="black", face="bold", size=10)) 

#p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())


#### Adding labels 
plotdf$labels <- ""
plotdf$labels[9]="L-Lysine"
plotdf$labels[15]="L-Lysine"
plotdf$labels[27]="Orotate"
plotdf$labels[28]="D-Galactonate"
plotdf$labels[36]="4-Trimethylammoniobutanoate"
plotdf$labels[38]="2-Hydroxyglutarate"
plotdf$labels[44]="L-Lactate"
plotdf$labels[52]="Itaconate"

p5 <- p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p6 <- p5 + geom_text_repel(size=2.5,colour="black", fontface = "bold", label = plotdf$labels, nudge_x = 2, nudge_y=0.5) + theme(legend.position="none")

ggsave(plot = p6, units = "mm",
       filename = "/metab-genemaineffectkd.png",
       width=100,
       height = 90,
       dpi = 600
)

###############

#Plotting Volcano for AB main effect
plotdf <- data.frame(ABMainEffect_inKD, KD_200$lacZ_AB.adj.p.)
plotdf$diffex<-"Non-Sig"
plotdf$diffex[ABMainEffect_inKD > 0 & KD_200$lacZ_AB.adj.p < 0.1] <- "Up"
plotdf$diffex[ABMainEffect_inKD < 0 & KD_200$lacZ_AB.adj.p. < 0.1] <- "Down"

# Min and max for scaling 
max(plotdf$ABMainEffect_inKD, na.rm=T)
min(plotdf$ABMainEffect_inKD, na.rm=T) 
max(-log10(KD_200$lacZ_AB.adj.p), na.rm=T) 
min(-log10(KD_200$lacZ_AB.adj.p), na.rm=T) 

p <- ggplot(plotdf, aes(x=ABMainEffect_inKD, y=-log10(KD_200$lacZ_AB.adj.p.), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-5,5) + ylim(0,10)
p2<- p + geom_hline(yintercept=-log10(0.1), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("#D81B60", "black", "#1E88E5"))
p4 <- p3 + labs(y="-Log10 Padj",x="Log2 Fold Change") + theme(axis.line=element_line(linewidth =1), axis.title = element_text(family = "Arial", color="black", face="bold", size=10)) 

#p4 + theme(axis.text=element_text(family="Arial", fontface = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())

plotdf$labels <- ""
plotdf$labels[38]="2-Hydroxyglutarate"
plotdf$labels[44]="L-Lactate"
p5 <- p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p6 <- p5 + geom_text_repel(size=2.5,colour="black", fontface = "bold", label = plotdf$labels, nudge_x = 2, nudge_y = 1) + theme(legend.position="none")


ggsave(plot = p6, units = "mm",
       filename = "/metab-abmain.png",
       width=100,
       height = 90,
       dpi = 600
)


###############

#Plotting Volcano for interaction 
plotdf <- data.frame(ABGeneInteraction_inKD, KD_200$Interaction.adj.p.)
plotdf$diffex<-"Non-Sig"
plotdf$diffex[ABGeneInteraction_inKD> 0 & KD_200$Interaction.adj.p. < 0.1] <- "Up"
plotdf$diffex[ABGeneInteraction_inKD< 0 & KD_200$Interaction.adj.p. < 0.1] <- "Down"

# Min and max for scaling 
max(plotdf$ABGeneInteraction_inKD, na.rm=T) 
min(plotdf$ABGeneInteraction_inKD, na.rm=T) 
max(-log10(KD_200$Interaction.adj.p.), na.rm=T)
min(-log10(KD_200$Interaction.adj.p.), na.rm=T) 

p <- ggplot(plotdf, aes(x=ABGeneInteraction_inKD, y=-log10(KD_200$Interaction.adj.p.), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-5,5) + ylim(0,10)
p2<- p + geom_hline(yintercept=-log10(0.1), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("#D81B60", "black", "#1E88E5")) 
p4 <- p3 + labs(y="-Log10 Padj",x="Log2 Fold Change") + theme(axis.line=element_line(linewidth =1), axis.title = element_text(family = "Arial", color="black", face="bold", size=10)) 
#p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())

plotdf$labels <- ""
plotdf$labels[44]="L-Lactate"
p5 <- p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p6 <- p5 + geom_text_repel(size=2.5,colour="black",fontface='bold', label = plotdf$labels) + theme(legend.position="none")

ggsave(plot = p6, units = "mm",
       filename = "/metab-interaction-kd.png",
       width=100,
       height = 90,
       dpi = 600
)


KD_save <- data.frame(GeneMainEffect_inKD,ABMainEffect_inKD,ABGeneInteraction_inKD, KD_200$Standard.match, KD_200$X)
names(KD_save) <- c("Gene Main Effect", "AB Main Effect", "Interaction", "Standard Match","ID")
write.csv(KD_save, file="/Metab Lipid/KD_LFC.csv")


################################# OVEREXPRESSION EXPERIMENT ########################################## FIX the NA problems here 

OE_All <- read.csv('/Metab Lipid/MetabolomicsData/2070_peaks_list_OE_dataset.csv')
OE_200 <- read.csv('/Metab Lipid/MetabolomicsData/Two_way_ANOVA_2720_OE_top_200.csv');

# Main effect of gene (I K control level) Vs. J L )
NoWwoxOE_Avg <- rowMeans(data.frame(OE_All$I1, OE_All$I2, OE_All$I3, OE_All$I4, OE_All$I5, OE_All$K1, OE_All$K2, OE_All$K3, OE_All$K4, OE_All$K5), na.rm=TRUE)
WwoxOE_Avg <- rowMeans(data.frame(OE_All$J1, OE_All$J2, OE_All$J3, OE_All$J4, OE_All$J5, OE_All$L1, OE_All$L2, OE_All$L3, OE_All$L4, OE_All$L5), na.rm=TRUE)
Foldchange <- WwoxOE_Avg[]/NoWwoxOE_Avg[];
LogFold <- log2(Foldchange)

index <- OE_200$X;
count = 0;
GeneMainEffect_inOE <- numeric(200)

for (x in index) {
  count = count+1;
  GeneMainEffect_inOE[count] <- LogFold[x]
}

# Main effect of AB (IJ (control level) Vs. KL)
LacZ_Avg_OE <- rowMeans(data.frame(OE_All$I1, OE_All$I2, OE_All$I3, OE_All$I4, OE_All$I5,OE_All$J1, OE_All$J2, OE_All$J3, OE_All$J4, OE_All$J5), na.rm=TRUE)
AB_Avg_OE <- rowMeans(data.frame(OE_All$K1, OE_All$K2, OE_All$K3, OE_All$K4, OE_All$K5, OE_All$L1, OE_All$L2, OE_All$L3, OE_All$L4, OE_All$L5), na.rm=TRUE)
Foldchange2 <- AB_Avg_OE[]/LacZ_Avg_OE[];
LogFold2 <- log2(Foldchange2)
count2 = 0;
ABMainEffect_inOE <- numeric(200)

for (x in index) {
  count2 = count2+1;
  ABMainEffect_inOE[count2] <- LogFold2[x]
}

# Gene in AB (K (control level) Vs. L)
L_Avg <- rowMeans(data.frame(OE_All$L1, OE_All$L2, OE_All$L3, OE_All$L4, OE_All$L5), na.rm=TRUE)
K_Avg <- rowMeans(data.frame(OE_All$K1, OE_All$K2, OE_All$K3, OE_All$K4, OE_All$K5), na.rm=TRUE)
Foldchange3 <- L_Avg[]/K_Avg[];
LogFold3 <- log2(Foldchange3)
count3 = 0;
ABGeneInteraction_inOE <- numeric(200)

for (x in index) {
  count3 = count3+1;
  ABGeneInteraction_inOE[count3] <- LogFold3[x]
}


#VOLCANOES 
#Plotting Volcano for gene main effect
plotdf <- data.frame(GeneMainEffect_inOE, OE_200$OE.adj.p.)
plotdf$diffex<-"Non-Sig"
plotdf$diffex[GeneMainEffect_inOE > 0 & OE_200$OE.adj.p. < 0.1] <- "Up"
plotdf$diffex[GeneMainEffect_inOE < 0 & OE_200$OE.adj.p. < 0.1] <- "Down"

# Min and max for scaling
max(plotdf$GeneMainEffect_inOE, na.rm=T) # 4.3
min(plotdf$GeneMainEffect_inOE, na.rm=T) # -3
max(-log10(OE_200$OE.adj.p.), na.rm=T) # 9.2 
min(-log10(OE_200$OE.adj.p.), na.rm=T) # -0.00

p <- ggplot(plotdf, aes(x=GeneMainEffect_inOE, y=-log10(OE_200$OE.adj.p.), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-6,6) + ylim(0,10)
p2<- p + geom_hline(yintercept=-log10(0.1), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("#D81B60", "black", "#1E88E5")) 
p4 <- p3 + labs(y="-Log10 Padj",x="Log2 Fold Change") + theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=10)) 

#p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())

plotdf$labels <- ""
plotdf$labels[21]="L-Leucine"
plotdf$labels[23]="Sn glycerol-3-phosphocholine"
plotdf$labels[10]="L-Methionine"
plotdf$labels[22]="N(pi)-methyl-L-histidine"
plotdf$labels[16]="L-Methionine"

p5 <- p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p6 <- p5 + geom_text_repel(size=2.5,colour="black",label = plotdf$labels, nudge_y=1, nudge_x=-1,segment.linetype = 1, segment.size = 0.2) + theme(legend.position="none")


ggsave(plot = p6, units = "mm",
       filename = "/metab-genemaineffect-oe.png",
       width=100,
       height = 90,
       dpi = 600
) 


#Plotting Volcano for AB main effect
plotdf <- data.frame(ABMainEffect_inOE, OE_200$AB.adj.p.)
plotdf$diffex<-"Non-Sig"
plotdf$diffex[ABMainEffect_inOE > 0 & OE_200$AB.adj.p. < 0.1] <- "Up"
plotdf$diffex[ABMainEffect_inOE < 0 & OE_200$AB.adj.p. < 0.1] <- "Down"

# Min and max for scaling
max(plotdf$ABMainEffect_inOE, na.rm=T)
min(plotdf$ABMainEffect_inOE, na.rm=T) 
max(-log10(OE_200$AB.adj.p.), na.rm=T) 
min(-log10(OE_200$AB.adj.p.), na.rm=T) 

p <- ggplot(plotdf, aes(x=ABMainEffect_inOE, y=-log10(OE_200$AB.adj.p.), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-6,6) + ylim(0,10)
p2<- p + geom_hline(yintercept=-log10(0.1), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("#D81B60", "black", "#1E88E5")) 
p4 <- p3 + labs(y="-Log10 Padj",x="Log2 Fold Change") + theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=10)) 

#p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())

plotdf$labels <- ""
plotdf$labels[23]="Sn glycerol-3-phosphocholine"
plotdf$labels[63]="Betaine"
plotdf$labels[106]="Lactate"
plotdf$labels[66]="2-Hydroxyglutarate"
plotdf$labels[139]="Citrate"

p5 <- p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p6 <- p5 + geom_text_repel(size=2.5,colour="black",label = plotdf$labels,segment.linetype = 1, force=20, segment.size = 0.2) + theme(legend.position="none")


ggsave(plot = p6, units = "mm",
       filename = "/metab-abmain-oe.png",
       width=100,
       height = 90,
       dpi = 600
)



#Plotting Volcano for interaction
plotdf <- data.frame(ABGeneInteraction_inOE, OE_200$Interaction.adj.p.)
plotdf$diffex<-"Non-Sig"
plotdf$diffex[ABGeneInteraction_inOE> 0 & OE_200$Interaction.adj.p. < 0.1] <- "Up"
plotdf$diffex[ABGeneInteraction_inOE< 0 & OE_200$Interaction.adj.p. < 0.1] <- "Down"

# Min and max for scaling 
max(plotdf$ABGeneInteraction_inOE, na.rm=T) # 5.15
min(plotdf$ABGeneInteraction_inOE, na.rm=T) # -1.5
max(-log10(OE_200$Interaction.adj.p.), na.rm=T) # 9.2 - first but it is nan
min(-log10(OE_200$Interaction.adj.p.), na.rm=T) # -0.00

p <- ggplot(plotdf, aes(x=ABGeneInteraction_inOE, y=-log10(OE_200$Interaction.adj.p.), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-6,6) + ylim(0,10)
p2<- p + geom_hline(yintercept=-log10(0.1), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("#D81B60", "black", "#1E88E5")) 
p4 <- p3 + labs(y="-Log10 Padj",x="Log2 Fold Change") + theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=10)) 

#p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())

plotdf$labels <- ""
plotdf$labels[10]="L-Methionine"
plotdf$labels[16]="L-Methionine"


p5 <- p4 + theme(axis.text= element_text(family="Arial", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p6 <- p5 + geom_text_repel(size=2.5,colour="black",label = plotdf$labels,segment.linetype = 1, force=80,segment.linetype = 1, segment.size = 0.2,nudge_y=0.5) + theme(legend.position="none")

ggsave(plot = p6, units = "mm",
       filename = "/metab-interaction-oe.png",
       width=100,
       height = 90,
       dpi = 600
)



OE_save <- data.frame(GeneMainEffect_inOE,ABMainEffect_inOE,ABGeneInteraction_inOE, OE_200$X)
names(OE_save) <- c("Gene Main Effect", "AB Main Effect", "Interaction", "ID")
write.csv(OE_save, file="/Metab Lipid/OE_LFC.csv")
