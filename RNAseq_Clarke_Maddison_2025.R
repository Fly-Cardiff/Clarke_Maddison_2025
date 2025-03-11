rm(list=ls())
workDir <- "~/Data_Results/RNASeq/featureCounts"
setwd(workDir)

devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

library("DESeq2")
library("dplyr")
library("ggplot2")
library("SARTools")
library("RColorBrewer")
library("pheatmap")
library("ggrepel")
library("scales")

#####################################
# For Overexpression experiment
sample_ID <- c("J1", "J2", "J3", "J4", "K1", "K2", "K3", "K4", "L1", "L2", "L3", "L4", "M1", "M3", "M4")
gRNA <- factor(c("OE", "OE", "OE", "OE", "Con", "Con", "Con", "Con", "OE", "OE", "OE", "OE", "Con", "Con", "Con"),levels = c("Con","OE"))
Amyloid <- factor(c("LacZ", "LacZ", "LacZ", "LacZ","LacZ", "LacZ", "LacZ", "LacZ","AB","AB","AB","AB","AB","AB","AB"), levels=c("LacZ","AB"))
filenames <- c("J1_S37_merge.markdup.featurecount","J2_S38_merge.markdup.featurecount","J3_S39_merge.markdup.featurecount","J4_S40_merge.markdup.featurecount","K1_S41_merge.markdup.featurecount","K2_S42_merge.markdup.featurecount","K3_S43_merge.markdup.featurecount","K4_S44_merge.markdup.featurecount","L1_S45_merge.markdup.featurecount","L2_S46_merge.markdup.featurecount","L3_S47_merge.markdup.featurecount","L4_S48_merge.markdup.featurecount","M1_S49_merge.markdup.featurecount","M3_S50_merge.markdup.featurecount","M4_S51_merge.markdup.featurecount")
sampleTable = data.frame(sample_ID,filenames,gRNA,Amyloid)
colData_Exp = data.frame(gRNA,Amyloid)

#For Elav Knockdown experiment 
sample_ID <- c("G1", "G2", "G3", "G4", "F1", "F2", "F3", "F4", "H1", "H2", "H3", "H4", "I1", "I2", "I3", "I4")
RNAi <- factor(c("Con", "Con", "Con", "Con", "KD", "KD", "KD", "KD", "KD", "KD", "KD", "KD", "Con", "Con", "Con", "Con"), levels = c("Con","KD"))
Amyloid <- factor(c("LacZ", "LacZ", "LacZ", "LacZ","LacZ", "LacZ", "LacZ", "LacZ","AB","AB","AB","AB","AB","AB","AB","AB"), levels = c("LacZ","AB"))
filenames <- c("G1_S25_merge.markdup.featurecount","G2_S26_merge.markdup.featurecount","G3_S27_merge.markdup.featurecount","G4_S28_merge.markdup.featurecount","F1_S21_merge.markdup.featurecount","F2_S22_merge.markdup.featurecount","F3_S23_merge.markdup.featurecount","F4_S24_merge.markdup.featurecount","H1_S29_merge.markdup.featurecount","H2_S30_merge.markdup.featurecount","H3_S31_merge.markdup.featurecount","H4_S32_merge.markdup.featurecount","I1_S33_merge.markdup.featurecount","I2_S34_merge.markdup.featurecount","I3_S35_merge.markdup.featurecount","I4_S36_merge.markdup.featurecount")
sampleTable = data.frame(sample_ID,filenames,RNAi,Amyloid)
colData_Exp = data.frame(RNAi,Amyloid)

##################################################### 

#removing ncRNA and transposables 
featuresToRemove = NULL #If null counts is 191,000 elements long 
featuresToRemove <- as.character(read.table("ncRNAs.txt", stringsAsFactors = FALSE))
featuresToRemove <- as.character(read.table("~Data_Results/RNASeq/ncRNAs_insertions_transposable.txt", stringsAsFactors = FALSE))
counts <- loadCountData(target=sampleTable, rawDir=workDir, featuresToRemove=featuresToRemove)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData_Exp, design = ~RNAi + Amyloid +RNAi:Amyloid, tidy=FALSE, ignoreRank = FALSE)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData_Exp, design = ~gRNA + Amyloid + gRNA:Amyloid, tidy=FALSE, ignoreRank = FALSE)

# Run DESeq
dds <-DESeq(dds, test="Wald", fitType = "parametric")

#Add step here to save the dds value for plotting in future 
save(dds, file="/Data_Results/RNASeq/DiffferentialGene/NOTransposablesORinsertons/RMDUP/Elav_OE_WWOX/dds_Elav_OE_NoInsertions_NoTransposons")

##################################################################################

load("RNASeq/MARKDUP/Elav_OE_WWOX/dds_Elav_OE_NoInsertions_NoTransposons.RData")
load("/RNASeq/MARKDUP/Elav_KD_WWOX/dds_Elav_KD_NoInsertions_NoTransposons.RData")

######################### RESULTS #################
resultsNames(dds)

# effect of WWOX OE in LacZ <- Main effect 
res  <- results(dds, contrast=c("gRNA","OE", "Con"),independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH") 
# effect of amyloid in control group <- Main effect
res <- results(dds, contrast=c("Amyloid","AB","LacZ"),independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH") 
# effect of WWOX OE in AB <- Interaction effect 
res <- results(dds, list( c("gRNA_OE_vs_Con","gRNAOE.AmyloidAB") ), independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH")


#effect of KD in LacZ 
res  <- results(dds, contrast=c("RNAi","KD","Con"),independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH") 
#effect of amyloid in control 
res <- results(dds, contrast=c("Amyloid","AB","LacZ"),independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH") 
#effect of KD in AB 
res <- results(dds, list( c("RNAi_KD_vs_Con","RNAiKD.AmyloidAB") ), independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH")


#Summary of results 
summary(res)
show(dds)

#shrinkage using ashr 
res_ashr <- lfcShrink(dds=dds, res=res, type="ashr")
resOrdered <- res_ashr[order(res_ashr$padj),]

#Finding those that are significant using a cut off of 1.5 fold which is 0.58 log2fold change  
resSig <- subset(resOrdered, padj < 0.05)
Up <- subset(resSig, log2FoldChange > 0.58)
Down <- subset(resSig, log2FoldChange < -0.58)

# For graphing - min and max
max(resOrdered$log2FoldChange)
min(resOrdered$log2FoldChange) 

# Save DEGs
write.csv(as.data.frame(Up),file = "/RNASeq/MARKDUP/Elav_OE_WWOX/Up.csv")
write.csv(as.data.frame(Down),file = "/RNASeq/MARKDUP/Elav_OE_WWOX/Down.csv")

###################### PLOTTING ###################################################
#PCA and dendrograms 
VSTdds <-vst(dds, blind = TRUE) #Variance stabilisation transformation 

#Making cluster dendrogram 
sampleDists <- dist(t(assay(VSTdds)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(VSTdds))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#PCA
rv <- rowVars(assay(VSTdds))
select <- order(rv, decreasing=TRUE)
pca_results <- prcomp(t(assay(VSTdds)[select,]))
PC1=pca_results$x[,1]
PC2=pca_results$x[,2]
PC3=pca_results$x[,3]

Genotype <- dds$RNAi

Amyloid <- dds$Amyloid
Genotype <- dds$gRNA

PCA_three <- data.frame(PC1,PC2,PC3,Genotype,Amyloid)
factor3 <- interaction(Genotype, Amyloid, sep = "_")
PCA_three <- data.frame(PC1,PC2,PC3,factor3)

#Getting eigenvalues
eigenvalues = pca_results$sdev^2
eigenvalues = eigenvalues/sum(eigenvalues)
percentVar = round(eigenvalues*100)

#Scree plot
barplot(percentVar,xlab="Principal Component",ylab="Percent variance (%)",names.arg=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15"),ylim=c(0,50),las=2)


#PCA

p1<- ggplot(PCA_three, aes(x=PC1, y=PC2, col=factor3)) + xlim(-50,50) + ylim(-50,50) + geom_point(size=2)  +
  xlab(paste0("PC1: ",percentVar[1],"% Variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% Variance")) 
p2 <- p1 + theme_light() + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p3 <- p2 + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=12)) + theme(panel.grid = element_line(size = 0.75, linetype = 2))
p4 <- p3 + scale_color_manual(values=c("black","#E69F00","#56B4E9","#009E73")) 

p1<- ggplot(PCA_three, aes(x=PC1, y=PC3, col=factor3)) + xlim(-50,50) + ylim(-50,50) + geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% Variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% Variance")) 
p2 <- p1 + theme_light() + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p3 <- p2 + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=12)) + theme(panel.grid = element_line(size = 0.75, linetype = 2))
p4 <- p3 + scale_color_manual(values=c("black","#E69F00","#56B4E9","#009E73")) 

p1<- ggplot(PCA_three, aes(x=PC2, y=PC3, col=factor3)) + xlim(-50,50) + ylim(-50,50) + geom_point(size=2) +
  xlab(paste0("PC2: ",percentVar[2],"% Variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% Variance")) 
p2 <- p1 + theme_light() + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p3 <- p2 + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=12)) + theme(panel.grid = element_line(size = 0.75, linetype = 2))
p4 <- p3 + scale_color_manual(values=c("black","#E69F00","#56B4E9","#009E73")) 

ggsave(plot = p4, units = "mm",
       filename = "/PCA.png",
       width=150,
       height = 100,
       dpi = 600
)


# MA Plot 

# Create a data frame for plotting, ensuring we have baseMean, log2FoldChange, and adjusted p-value
ma_data <- as.data.frame(resOrdered)

# Create a new column for coloring based on log2FoldChange
ma_data$color <- "black"
ma_data$color[ma_data$log2FoldChange > 0.58 & ma_data$padj < 0.05] <- "#1E88E5"     # Upregulated genes
ma_data$color[ma_data$log2FoldChange < -0.58 & ma_data$padj < 0.05] <- "#D81B60" # Downregulated genes

# Create the MA plot
p1<- ggplot(ma_data, aes(x = baseMean, y = log2FoldChange, color = color)) + geom_point(size=0.8) +
  scale_color_identity() +  # Use the actual color values in the 'color' column
  scale_x_log10(labels = label_comma(drop0trailing = TRUE))  +
  labs(x = "Mean of normalized counts", y = "Log2 Fold Change") +
  theme_classic() + theme(axis.line=element_line(size=0.8),axis.title = element_text(family = "Arial", color="black", face="bold", size=13)) +
  ylim(-25, 25)  # Adjust the y-axis limit as needed
p2 <- p1 + geom_hline(yintercept=c(-0.58,0.58), col="darkgrey", lty="longdash")
p3 <- p2 + theme(axis.text = element_text(family="Arial", face = "bold", colour = "black", size = 13))

#Genes to label on plot
which(rownames(ma_data)=="FBgn0001258") # Code of LDH 
rownames(ma_data)[1] = "Ldh"
ma_data$labels <- ""
ma_data$labels <- ifelse(rownames(ma_data) == "Ldh",TRUE,FALSE)

p4 <- p3 + geom_text_repel(size=3.5,colour="black",fontface = "bold",label = ifelse(ma_data$labels, rownames(ma_data),""),nudge_x = 0, nudge_y = 10,max.overlaps = Inf)

# Save the plot 
ggsave(plot = p4, units = "mm",
       filename = "/Users/clarkeh/Dropbox/WWOX Manuscript 1",
       height = 90,
       dpi = 600
)

