


#install.packages("readxl")
library(readxl)
#install.packages("writexl")
library(writexl)

# for the cleaning tools and the object tools # 

source ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/olinkdata_objectTools_v1.1.R")

source ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/olinkdata_cleaningTools_v1.1.R")

outputDir = "M:/NL/6_Scientific Results/2_Students/Noha Salem/3 Output/20240212"

project_name = 'MAPT'


# defining the files that we will use # 
#1. olink file 

olinkDataFiles <- c("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/olink_npx_data.xlsx")
olinkDataSetNames <- c("data_MAPT")
OlinkData <- createOlinkDataObject(olinkDataFiles=olinkDataFiles,olinkDataSetNames=olinkDataSetNames,olinkLayoutFiles="")

#2. the organoids 
#2.a. to subsite the organoids 
sample_info <- read_excel("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/sample_info.xlsx")
Organoids <-sample_info[sample_info$group == "Organoid" , ]

excel_file <-("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/Organoids.xlsx")
write_xlsx(x=Organoids, path=excel_file)

olink_npx_data <- read_excel("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/olink_data_pilot_Mapt_contr_removed.xlsx")
olink_npx_data <- olink_npx_data[which(olink_npx_data$SampleID %in% Organoids$SampleID), ]

excel_file <-("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/olink_npx_data.xlsx")
write_xlsx(x=olink_npx_data, path=excel_file)

OlinkData@samples[['project']] = OlinkData@samples[['group']]


#annotation# 

Organoids <- ("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/Organoids.xlsx")
SampleIDcolumnsName <- ("SampleID")
OlinkData = annotateOlinkData(OlinkData,annotationFiles=Organoids,OlinkIDcolumnsName=SampleIDcolumnsName)

#set a specific color to the organoids# 

OlinkData@color.vector <- c (OlinkData@color.vector,setNames("#cc99a2" ,"Organoids"))


OlinkData@samples[['projectName']] = "Organoids"


OlinkData = set.design(OlinkData,
                       projectName = "Organoids",
                       samples.columnName.group = "Mutation",
                       samples.columnName.time = "time",
                       samples.columnName.time.label = "timepointLabel",
                       samples.columnName.sampleID = "SampleID",
                       samples.columnName.subjectID = "Donor")

designLayoutPlot = plot.design.layout(OlinkData, prct=F, colorbyProject = T)
plot (designLayoutPlot)

install.packages ("GGally2")
PCA_plot_scree_raw = plot.PCA.scree(OlinkData, slotName = "NPX.raw",nComponents=30, returnOut = "figure", verbose = F)
plot (PCA_plot_scree_raw)


PCA_plot_raw = plot.PCA(OlinkData,slot="NPX.raw",nComponents=3,samples.col.name.color="projectName")
PCA_plot_raw

# to detect the NPX below the LOD #
install.packages("diptest")
OlinkData = detectLowNPX(OlinkData,LOD.factor=1,numberCores=12)

#to create cut off plot and score with the default value 85% as we have more samples in Organoids

cutoffPlot_LOD = makeCutoffPlot(OlinkData, functionName = "detectLowNPX", designCol = "project")
plot (cutoffPlot_LOD) #almost 1600 proteins


#exclude the proteins that have NPX values below the LOD 
OlinkData = excludeLowNPX(OlinkData, prctPass = 85, FilterPerProject = T, numberCores = 12)


#to make sure how your cleaning step worked 

PCA_plot_LOD = plot.PCA(OlinkData, slot= "NPX", nComponents = 3, samples.col.name.color = "projectName")
PCA_plot_LOD

# check the irregularities #

#1. define the z-scores 
OlinkData <- zScores (OlinkData, method = "median")

cutoffPlot_capOutliers = makeCutoffPlot(OlinkData, functionName = "capOutliers", designCol = "project") #define the effect of cutoff, the default z-scores cut off - and + 5


#2. capping outlines 
OlinkData <- capOutliers(OlinkData, lower.zscore.cutoff = -5, upper.zscore.cutoff = 5, numberCores = 12, verbose = T)

PCA_plot_capped = plot.PCA(OlinkData, slot= "NPX.capped", nComponents= 3, samples.col.name.color = "projectName")# to check the effect of capping outliers 

PCA_plot_capped


OlinkData = interSampleCorrelation(OlinkList=OlinkData,panels=OlinkData@proteins$Panel,numberCores=12,corMethod="pearson",verbose=T) #


cutoffPlot_ISC = makeCutoffPlot(OlinkData, functionName="apply.ISC",designCol="project")#to check the panels excluded by the ISC looks like 
cutoffPlot_ISC

ISC_plot = plot.ISC(OlinkData,cutoff = 0.7,plotType = "ranked",order.column = "allPanels",verbose = F)
ISC_plot



#for the IVS10+16 mutation 
mut1_IVS10_2m <- OlinkData@NPX [ , c("9A","24H")]
wt1_IVS10_2m <- OlinkData@NPX [ , c("40H","20D")]
mut1_IVS10_4m <- OlinkData@NPX [ , c("32H")]
wt1_IVS10_4m <- OlinkData@NPX [ , c("39G")]
mut2_IVS10_2m <- OlinkData@NPX [ , c("7G","12D")]
wt2_IVS10_2m <- OlinkData@NPX [ , c("6F","13E")]
mut2_IVS10_4m <- OlinkData@NPX [ , c("38F")]
wt2_IVS10_4m <- OlinkData@NPX [ , c("26B")]

#to omit NA
mut1_IVS10_2m <- na.omit(mut1_IVS10_2m)
wt1_IVS10_2m <- na.omit(wt1_IVS10_2m)
mut1_IVS10_4m <- na.omit(mut1_IVS10_4m)
wt1_IVS10_4m <- na.omit(wt1_IVS10_4m)
mut2_IVS10_2m <- na.omit(mut2_IVS10_2m)
wt2_IVS10_2m <- na.omit(wt2_IVS10_2m)
mut2_IVS10_4m <- na.omit(mut2_IVS10_4m)
wt2_IVS10_4m <- na.omit(wt2_IVS10_4m)

#to calculate the mean 
mut1_IVS10_2m <- as.data.frame(mut1_IVS10_2m)
mut1_IVS10_2m$mut1_IVS10_2m <- rowMeans (mut1_IVS10_2m)
mut1_IVS10_2m$panelassay <- rownames(mut1_IVS10_2m)

wt1_IVS10_2m <- as.data.frame(wt1_IVS10_2m)
wt1_IVS10_2m$wt1_IVS10_2m <- rowMeans (wt1_IVS10_2m)
wt1_IVS10_2m$panelassay <- rownames(wt1_IVS10_2m)

mut1_IVS10_4m <- as.data.frame(mut1_IVS10_4m)
mut1_IVS10_4m$panelassay <- rownames (mut1_IVS10_4m)

wt1_IVS10_4m <- as.data.frame(wt1_IVS10_4m)
wt1_IVS10_4m$panelassay <- rownames(wt1_IVS10_4m)

mut2_IVS10_2m <- as.data.frame(mut2_IVS10_2m)
mut2_IVS10_2m$mut2_IVS10_2m <- rowMeans (mut2_IVS10_2m)
mut2_IVS10_2m$panelassay <- rownames(mut2_IVS10_2m)

wt2_IVS10_2m <- as.data.frame(wt2_IVS10_2m)
wt2_IVS10_2m$wt2_IVS10_2m <- rowMeans (wt2_IVS10_2m)
wt2_IVS10_2m$panelassay <- rownames(wt2_IVS10_2m)

mut2_IVS10_4m <- as.data.frame(mut2_IVS10_4m)
mut2_IVS10_4m$panelassay <- rownames (mut2_IVS10_4m)

wt2_IVS10_4m <- as.data.frame(wt2_IVS10_4m)
wt2_IVS10_4m$panelassay <- rownames(wt2_IVS10_4m)

#combine all the means for all samples
IVS10 <- as.data.frame(cbind (mut1_IVS10_2m$mut1_IVS10_2m,
                              wt1_IVS10_2m$wt1_IVS10_2m,
                              mut1_IVS10_4m$mut1_IVS10_4m,
                              wt1_IVS10_4m$wt1_IVS10_4m,
                              mut2_IVS10_2m$mut2_IVS10_2m,
                              wt2_IVS10_2m$wt2_IVS10_2m,
                              mut2_IVS10_4m$mut2_IVS10_4m,
                              wt2_IVS10_4m$wt2_IVS10_4m))


colnames(IVS10) <- c('mut1_IVS10_2m', 'wt1_IVS10_2m', 
                     'mut1_IVS10_4m', 'wt1_IVS10_4m',
                     'mut2_IVS10_2m', 'wt2_IVS10_2m',
                     'mut2_IVS10_4m', 'wt2_IVS10_4m')

rownames(IVS10) <- mut1_IVS10_2m$panelassay 
IVS10$PanelAssay <- mut1_IVS10_2m$panelassay 


#t-test for IVS10+16 

#1.adding the UNIPROT accession ID to the proteins 
uniprot <- OlinkData@proteins
uniprot$protein <- paste0(uniprot$panel , "_", uniprot$Assay)


#2.the new dataframe (uniprot) contains all values including those that are below the LOD and those that are NA (missingness)
missingness <- OlinkData@prct.NPX.below.LOD #assigning the prct.NPX.below.LOD to the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m') #separate the wt from the mut


#the thing that I did that here is that I can't combine based on the proteins as some proteins are in duplicates, I need UNIQUE value for each thing so I will combine using the PanelAssay
all_IVS10 <- merge(uniprot, IVS10, by="PanelAssay", all.y=TRUE)
#write.csv(all_IVS10, "all_IVS10.csv")

library ("FactoMineR") #for multivariate exploratory data analysis (MV-EDA), provides set of functions for PCA and clustering
library ("factoextra") #for extracting and visualizing the results of mulivariate data 


pca_data <- t(IVS10 [ ,1:8]) #the means for the 4 samples that I have 
PCA (pca_data, scale.unit = TRUE, ncp = 5, graph = TRUE) #scale unit is important to standardize all data (z-scores normalization)thus ensures that all data will contribute to the PCA making PCA more robust (decreases results skewness)
res.pca <- PCA(pca_data, graph = FALSE) #the results of PCA in a list

eig.val <- get_eigenvalue(res.pca) #eigenvalues represent the total amount of variance that can be explained by PCA. 
eig.val

fviz_eig (res.pca, addlabels = TRUE, ylim = c (0,50)) #to plot the eigenvalues/variance against the no. dimensions
var <- get_pca_var(res.pca)
var

#to identify 3 distinct groups (clusters) based on the similarity betweeen data points 
set.seed(123)
res.km <- kmeans (var$coord, centers = 3, nstart = 25) #to cluster the observations (var$coord) into 3 clusters based upon the nearest mean (centroid)
grp <- as.factor(res.km$cluster)

# coloring the clusters 
fviz_pca_var(res.pca, col.var = grp, palette = c ("#F2AFB4", "#9DC893", "#36454F"), legend.title = "Cluster")

# detailed description of the first 2 PCs  
res.desc <- dimdesc(res.pca, axes = c (1,2), proba = 0.05) #list for each dim contains the p-value and the correlation 

#to get info about the pca but for each dim (we have 2 dims and this will assign each sample with varaibles to specific location)
ind <- get_pca_ind(res.pca)
ind
grp <- as.factor (c ('mut1_IVS10_2m', 'wt1_IVS10_2m', 
                     'mut1_IVS10_4m', 'wt1_IVS10_4m',
                     'mut2_IVS10_2m', 'wt2_IVS10_2m',
                     'mut2_IVS10_4m', 'wt2_IVS10_4m'))

fviz_pca_ind(res.pca, fill.ind = grp, pointshape = 21, palette = c ("#F2AFB4", "#9DC893", "#36454F", "#005b96", "#265828", "#7293ec", "#790049", "#b21c0e"), repel = TRUE)

#t_test 
library(matrixTests)
library(plyr)
library(ggrepel)


all_IVS10 <- read.csv ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/all_IVS10.csv", row.names = 1)
#create NPX df that contains the means for the protein in each sample
NPX <-all_IVS10 [, (14:21)]
NPX$PanelAssay<-all_IVS10$PanelAssay
#define the con cloumns and mut columns separtely as it is paired t-tes 
#this arrangement is critical because we need to define the mean.diff (mean.diff= mean_mut - mean_con)

#to take in ur considerations the background 

mut2m_IVS10 <- c ('mut1_IVS10_2m', 'mut2_IVS10_2m')
wt2m_IVS10 <- c ('wt1_IVS10_2m', 'wt2_IVS10_2m')
mut4m_IVS10 <- c ('mut1_IVS10_4m', 'mut2_IVS10_4m')
wt4m_IVS10 <- c ('wt1_IVS10_4m', 'wt2_IVS10_4m')

all_IVS10[, mut2m_IVS10] 

#I will separate t-test according to the time point
#for 2m time point
t_test_protein_IVS10_2m <- row_t_paired(NPX [,colnames(NPX) %in% mut2m_IVS10],
                                        NPX [,colnames(NPX) %in% wt2m_IVS10])

t_test_protein_IVS10_2m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_IVS10_2m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_IVS10_2m$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_IVS10_2m <- merge(uniprot, t_test_protein_IVS10_2m, by = "PanelAssay", all.y = T)

t_test_protein_IVS10_2m$diffexpressed <- "Not-significant"  
t_test_protein_IVS10_2m$diffexpressed[t_test_protein_IVS10_2m$mean.diff>0.0 & t_test_protein_IVS10_2m$pvalue<0.05] <- "UP"
t_test_protein_IVS10_2m$diffexpressed[t_test_protein_IVS10_2m$mean.diff<0.0 & t_test_protein_IVS10_2m$pvalue<0.05] <- "DOWN"
t_test_protein_IVS10_2m$protein2 <- gsub(".*_", "", t_test_protein_IVS10_2m$PanelAssay)
t_test_protein_IVS10_2m$delabel <- NA
t_test_protein_IVS10_2m$delabel[t_test_protein_IVS10_2m$diffexpressed != "Not-significant"] <- t_test_protein_IVS10_2m$protein2[t_test_protein_IVS10_2m$diffexpressed != "Not-significant"]
t_test_protein_IVS10_2m$diffexpressed <- factor(t_test_protein_IVS10_2m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 


#enhanced vol plot
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_IVS10_2m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "IVS10+16_2m",
                     labels = c("Up (20)", "Down (8)", "Not-sig (1309)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)


ggplot(data = t_test_protein_IVS10_2m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression_IVS10_2m" , labels= c("Upregulated (n=20)", "Downregulated (n=8)", "Not significant (n=1309)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#for 4m time point 

t_test_protein_IVS10_4m <- row_t_paired(NPX [,colnames(NPX) %in% mut4m_IVS10],
                                        NPX [,colnames(NPX) %in% wt4m_IVS10])

t_test_protein_IVS10_4m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_IVS10_4m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_IVS10_4m$adjusted_p_value <- adjusted_p_value

uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_IVS10_4m <- merge(uniprot, t_test_protein_IVS10_4m, by = "PanelAssay", all.y = T)

t_test_protein_IVS10_4m$diffexpressed <- "Not-significant"  
t_test_protein_IVS10_4m$diffexpressed[t_test_protein_IVS10_4m$mean.diff>0.0 & t_test_protein_IVS10_4m$pvalue<0.05] <- "UP"
t_test_protein_IVS10_4m$diffexpressed[t_test_protein_IVS10_4m$mean.diff<0.0 & t_test_protein_IVS10_4m$pvalue<0.05] <- "DOWN"
t_test_protein_IVS10_4m$protein2 <- gsub(".*_", "", t_test_protein_IVS10_4m$PanelAssay)
t_test_protein_IVS10_4m$delabel <- NA
t_test_protein_IVS10_4m$delabel[t_test_protein_IVS10_4m$diffexpressed != "Not-significant"] <- t_test_protein_IVS10_4m$protein2[t_test_protein_IVS10_4m$diffexpressed != "Not-significant"]
t_test_protein_IVS10_4m$diffexpressed <- factor(t_test_protein_IVS10_4m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 


#enhanced vol plot
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_IVS10_4m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "IVS10+16_4m",
                     labels = c("Up (8)", "Down (33)", "Not-sig (1296)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)


ggplot(data = t_test_protein_IVS10_4m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 20, show.legend = F)+  
  scale_color_manual(name = "Protein expression_IVS10_4m" , labels= c("Upregulated (n=8)", "Downregulated (n=33)", "Not significant (n=1296)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#for the S305I mutation 
mut1_S305I_2m <- OlinkData@NPX [ , c("11C","22F")]
wt1_S305I_2m <- OlinkData@NPX [ , c("21E","18B")]
mut1_S305I_4m <- OlinkData@NPX [ , c("37E")]
wt1_S305I_4m <- OlinkData@NPX [ , c("28D")]
mut2_S305I_2m <- OlinkData@NPX [ , c("14F","19C")]
wt2_S305I_2m <- OlinkData@NPX [ , c("8H","16H")]
mut2_S305I_4m <- OlinkData@NPX [ , c("31G")]
wt2_S305I_4m <- OlinkData@NPX [ , c("30F")]

#to omit NA
mut1_S305I_2m <- na.omit(mut1_S305I_2m)
wt1_S305I_2m <- na.omit(wt1_S305I_2m)
mut1_S305I_4m <- na.omit(mut1_S305I_4m)
wt1_S305I_4m <- na.omit(wt1_S305I_4m)
mut2_S305I_2m <- na.omit(mut2_S305I_2m)
wt2_S305I_2m <- na.omit(wt2_S305I_2m)
mut2_S305I_4m <- na.omit(mut2_S305I_4m)
wt2_S305I_4m <- na.omit(wt2_S305I_4m)


#to calculate the mean 
mut1_S305I_2m <- as.data.frame(mut1_S305I_2m)
mut1_S305I_2m$mut1_S305I_2m <- rowMeans (mut1_S305I_2m)
mut1_S305I_2m$panelassay <- rownames(mut1_S305I_2m)

wt1_S305I_2m <- as.data.frame(wt1_S305I_2m)
wt1_S305I_2m$wt1_S305I_2m <- rowMeans (wt1_S305I_2m)
wt1_S305I_2m$panelassay <- rownames(wt1_S305I_2m)

mut1_S305I_4m <- as.data.frame(mut1_S305I_4m)
mut1_S305I_4m$panelassay <- rownames (mut1_S305I_4m)

wt1_S305I_4m <- as.data.frame(wt1_S305I_4m)
wt1_S305I_4m$panelassay <- rownames(wt1_S305I_4m)

mut2_S305I_2m <- as.data.frame(mut2_S305I_2m)
mut2_S305I_2m$mut2_S305I_2m <- rowMeans (mut2_S305I_2m)
mut2_S305I_2m$panelassay <- rownames(mut2_S305I_2m)

wt2_S305I_2m <- as.data.frame(wt2_S305I_2m)
wt2_S305I_2m$wt2_S305I_2m <- rowMeans (wt2_S305I_2m)
wt2_S305I_2m$panelassay <- rownames(wt2_S305I_2m)

mut2_S305I_4m <- as.data.frame(mut2_S305I_4m)
mut2_S305I_4m$panelassay <- rownames (mut2_S305I_4m)

wt2_S305I_4m <- as.data.frame(wt2_S305I_4m)
wt2_S305I_4m$panelassay <- rownames(wt2_S305I_4m)

#combine all the means for all samples
S305I <- as.data.frame(cbind (mut1_S305I_2m$mut1_S305I_2m,
                              wt1_S305I_2m$wt1_S305I_2m,
                              mut1_S305I_4m$mut1_S305I_4m,
                              wt1_S305I_4m$wt1_S305I_4m,
                              mut2_S305I_2m$mut2_S305I_2m,
                              wt2_S305I_2m$wt2_S305I_2m,
                              mut2_S305I_4m$mut2_S305I_4m,
                              wt2_S305I_4m$wt2_S305I_4m))


colnames(S305I) <- c('mut1_S305I_2m', 'wt1_S305I_2m', 
                     'mut1_S305I_4m', 'wt1_S305I_4m',
                     'mut2_S305I_2m', 'wt2_S305I_2m',
                     'mut2_S305I_4m', 'wt2_S305I_4m')

rownames(S305I) <- mut1_S305I_2m$panelassay 
S305I$PanelAssay <- mut1_S305I_2m$panelassay 


#t-test for IVS10+16 

#1.adding the UNIPROT accession ID to the proteins 
uniprot <- OlinkData@proteins
uniprot$protein <- paste0(uniprot$panel , "_", uniprot$Assay)


#2.the new dataframe (uniprot) contains all values including those that are below the LOD and those that are NA (missingness)
missingness <- OlinkData@prct.NPX.below.LOD #assigning the prct.NPX.below.LOD to the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m') #separate the wt from the mut


#the thing that I did that here is that I can't combine based on the proteins as some proteins are in duplicates, I need UNIQUE value for each thing so I will combine using the PanelAssay
all_S305I <- merge(uniprot, S305I, by="PanelAssay", all.y=TRUE)
write.csv(all_S305I, "all_S305I.csv")

library ("FactoMineR") #for multivariate exploratory data analysis (MV-EDA), provides set of functions for PCA and clustering
library ("factoextra") #for extracting and visualizing the results of mulivariate data 


pca_data <- t(all_S305I [ ,14:21]) #the means for the 4 samples that I have 
PCA (pca_data, scale.unit = TRUE, ncp = 5, graph = TRUE) #scale unit is important to standardize all data (z-scores normalization)thus ensures that all data will contribute to the PCA making PCA more robust (decreases results skewness)
res.pca <- PCA(pca_data, graph = FALSE) #the results of PCA in a list

eig.val <- get_eigenvalue(res.pca) #eigenvalues represent the total amount of variance that can be explained by PCA. 
eig.val

fviz_eig (res.pca, addlabels = TRUE, ylim = c (0,50)) #to plot the eigenvalues/variance against the no. dimensions
var <- get_pca_var(res.pca)
var

#to identify 3 distinct groups (clusters) based on the similarity betweeen data points 
set.seed(123)
res.km <- kmeans (var$coord, centers = 3, nstart = 25) #to cluster the observations (var$coord) into 3 clusters based upon the nearest mean (centroid)
grp <- as.factor(res.km$cluster)

# coloring the clusters 
fviz_pca_var(res.pca, col.var = grp, palette = c ("#F2AFB4", "#9DC893", "#36454F"), legend.title = "Cluster")

# detailed description of the first 2 PCs  
res.desc <- dimdesc(res.pca, axes = c (1,2), proba = 0.05) #list for each dim contains the p-value and the correlation 

#to get info about the pca but for each dim (we have 2 dims and this will assign each sample with varaibles to specific location)
ind <- get_pca_ind(res.pca)
ind
grp <- as.factor (c ('mut1_S305I_2m', 'wt1_S305I_2m', 
                     'mut1_S305I_4m', 'wt1_S305I_4m',
                     'mut2_S305I_2m', 'wt2_S305I_2m',
                     'mut2_S305I_4m', 'wt2_S305I_4m'))

fviz_pca_ind(res.pca, fill.ind = grp, pointshape = 21, palette = c ("#F2AFB4", "#9DC893", "#36454F", "#005b96", "#265828", "#7293ec", "#790049", "#b21c0e"), repel = TRUE)


#t_test 
library(matrixTests)
library(plyr)
library(ggrepel)


all_S305I <- read.csv ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/all_S305I.csv", row.names = 1)
#create NPX df that contains the means for the protein in each sample
NPX <- all_S305I[, (14:21)]
NPX$PanelAssay<-all_S305I$PanelAssay
#define the con cloumns and mut columns separtely as it is paired t-tes 
#this arrangement is critical because we need to define the mean.diff (mean.diff= mean_mut - mean_con)

#to take in ur considerations the background 

mut2m_S305I <- c ('mut1_S305I_2m', 'mut2_S305I_2m')
wt2m_S305I <- c ('wt1_S305I_2m', 'wt2_S305I_2m')
mut4m_S305I <- c ('mut1_S305I_4m', 'mut2_S305I_4m')
wt4m_S305I <- c ('wt1_S305I_4m', 'wt2_S305I_4m')

all_S305I[, mut2m_S305I]


#I will separate t-test according to the time point
#for 2m time point
t_test_protein_S305I_2m <- row_t_paired(NPX [,colnames(NPX) %in% mut2m_S305I],
                                        NPX [,colnames(NPX) %in% wt2m_S305I])

t_test_protein_S305I_2m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_S305I_2m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_S305I_2m$adjusted_p_value <- adjusted_p_value



uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_S305I_2m <- merge(uniprot, t_test_protein_S305I_2m, by = "PanelAssay", all.y = T)

t_test_protein_S305I_2m$diffexpressed <- "Not-significant"  
t_test_protein_S305I_2m$diffexpressed[t_test_protein_S305I_2m$mean.diff>0.0 & t_test_protein_S305I_2m$pvalue<0.05] <- "UP"
t_test_protein_S305I_2m$diffexpressed[t_test_protein_S305I_2m$mean.diff<0.0 & t_test_protein_S305I_2m$pvalue<0.05] <- "DOWN"
t_test_protein_S305I_2m$protein2 <- gsub(".*_", "", t_test_protein_S305I_2m$PanelAssay)
t_test_protein_S305I_2m$delabel <- NA
t_test_protein_S305I_2m$delabel[t_test_protein_S305I_2m$diffexpressed != "Not-significant"] <- t_test_protein_S305I_2m$protein2[t_test_protein_S305I_2m$diffexpressed != "Not-significant"]
t_test_protein_S305I_2m$diffexpressed <- factor(t_test_protein_S305I_2m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 


#enhanced vol plot
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_S305I_2m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "S305I_2m",
                     labels = c("Up (28)", "Down (30)", "Not-sig (1279)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)




ggplot(data = t_test_protein_S305I_2m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 12, show.legend = F)+  
  scale_color_manual(name = "Protein expression_S305I_2m" , labels= c("Upregulated (n=28)", "Downregulated (n=30)", "Not significant (n=1279)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()

#for 4m time point 
t_test_protein_S305I_4m <- row_t_paired(NPX [,colnames(NPX) %in% mut4m_S305I],
                                        NPX [,colnames(NPX) %in% wt4m_S305I])

t_test_protein_S305I_4m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_S305I_4m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_S305I_4m$adjusted_p_value <- adjusted_p_value

uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_S305I_4m <- merge(uniprot, t_test_protein_S305I_4m, by = "PanelAssay", all.y = T)

t_test_protein_S305I_4m$diffexpressed <- "Not-significant"  
t_test_protein_S305I_4m$diffexpressed[t_test_protein_S305I_4m$mean.diff>0.0 & t_test_protein_S305I_4m$pvalue<0.05] <- "UP"
t_test_protein_S305I_4m$diffexpressed[t_test_protein_S305I_4m$mean.diff<0.0 & t_test_protein_S305I_4m$pvalue<0.05] <- "DOWN"
t_test_protein_S305I_4m$protein2 <- gsub(".*_", "", t_test_protein_S305I_4m$PanelAssay)
t_test_protein_S305I_4m$delabel <- NA
t_test_protein_S305I_4m$delabel[t_test_protein_S305I_4m$diffexpressed != "Not-significant"] <- t_test_protein_S305I_4m$protein2[t_test_protein_S305I_4m$diffexpressed != "Not-significant"]
t_test_protein_S305I_4m$diffexpressed <- factor(t_test_protein_S305I_4m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 

#enhanced vol plot
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_S305I_4m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "S305I_4m",
                     labels = c("Up (30)", "Down (29)", "Not-sig (1278)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)








ggplot(data = t_test_protein_S305I_4m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression_S305I_4m" , labels= c("Upregulated (n=30)", "Downregulated (n=29)", "Not significant (n=1278)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#for the R406W
mut1_R406W_2m <- OlinkData@NPX [ , c("27C","23G")]
wt1_R406W_2m <- OlinkData@NPX [ , c("15G","29E")]
mut1_R406W_4m <-OlinkData@NPX [ , c("33A")]
wt1_R406W_4m <-OlinkData@NPX [ , c("25A")]
mut2_R406W_2m <- OlinkData@NPX [ , c("10B","17A")]
wt2_R406W_2m <-OlinkData@NPX [ , c("5E","34B")]
mut2_R406W_4m <- OlinkData@NPX [ , c("35C")]
wt2_R406W_4m <- OlinkData@NPX [ , c("36D")]


#to omit NA
mut1_R406W_2m <- na.omit(mut1_R406W_2m)
wt1_R406W_2m <- na.omit(wt1_R406W_2m)
mut1_R406W_4m <- na.omit(mut1_R406W_4m)
wt1_R406W_4m <- na.omit(wt1_R406W_4m)
mut2_R406W_2m <- na.omit(mut2_R406W_2m)
wt2_R406W_2m <- na.omit(wt2_R406W_2m)
mut2_R406W_4m <- na.omit(mut2_R406W_4m)
wt2_R406W_4m <- na.omit(wt2_R406W_4m)


#to calculate the mean
mut1_R406W_2m <- as.data.frame(mut1_R406W_2m)
mut1_R406W_2m$mut1_R406W_2m <- rowMeans (mut1_R406W_2m)
mut1_R406W_2m$panelassay <- rownames(mut1_R406W_2m)

wt1_R406W_2m <- as.data.frame(wt1_R406W_2m)
wt1_R406W_2m$wt1_R406W_2m <- rowMeans (wt1_R406W_2m)
wt1_R406W_2m$panelassay <- rownames(wt1_R406W_2m)

mut1_R406W_4m <- as.data.frame(mut1_R406W_4m)
mut1_R406W_4m$panelassay <- rownames(mut1_R406W_4m)

wt1_R406W_4m <- as.data.frame(wt1_R406W_4m)
wt1_R406W_4m$panelassay <- rownames(wt1_R406W_4m)


mut2_R406W_2m <- as.data.frame(mut2_R406W_2m)
mut2_R406W_2m$mut2_R406W_2m <- rowMeans (mut2_R406W_2m)
mut2_R406W_2m$panelassay <- rownames(mut2_R406W_2m)

wt2_R406W_2m <- as.data.frame(wt2_R406W_2m)
wt2_R406W_2m$wt2_R406W_2m <- rowMeans (wt2_R406W_2m)
wt2_R406W_2m$panelassay <- rownames(wt2_R406W_2m)


mut2_R406W_4m <- as.data.frame(mut2_R406W_4m)
mut2_R406W_4m$panelassay <- rownames(mut2_R406W_4m)


wt2_R406W_4m <- as.data.frame(wt2_R406W_4m)
wt2_R406W_4m$panelassay <- rownames(wt2_R406W_4m)

R406W <- as.data.frame(cbind (mut1_R406W_2m$mut1_R406W_2m,
                              wt1_R406W_2m$wt1_R406W_2m,
                              mut1_R406W_4m$mut1_R406W_4m,
                              wt1_R406W_4m$wt1_R406W_4m,
                              mut2_R406W_2m$mut2_R406W_2m,
                              wt2_R406W_2m$wt2_R406W_2m,
                              mut2_R406W_4m$mut2_R406W_4m,
                              wt2_R406W_4m$wt2_R406W_4m))


colnames(R406W) <- c('mut1_R406W_2m', 'wt1_R406W_2m', 
                     'mut1_R406W_4m', 'wt1_R406W_4m',
                     'mut2_R406W_2m', 'wt2_R406W_2m',
                     'mut2_R406W_4m', 'wt2_R406W_4m')

rownames(R406W) <- mut1_R406W_2m$panelassay 
R406W$PanelAssay <- mut1_R406W_2m$panelassay 

#t-test for IVS10+16 

#1.adding the UNIPROT accession ID to the proteins 
uniprot <- OlinkData@proteins
uniprot$protein <- paste0(uniprot$panel , "_", uniprot$Assay)


#2.the new dataframe (uniprot) contains all values including those that are below the LOD and those that are NA (missingness)
missingness <- OlinkData@prct.NPX.below.LOD #assigning the prct.NPX.below.LOD to the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m') #separate the wt from the mut


#the thing that I did that here is that I can't combine based on the proteins as some proteins are in duplicates, I need UNIQUE value for each thing so I will combine using the PanelAssay
all_R406W <- merge(uniprot, R406W, by="PanelAssay", all.y=TRUE)
write.csv(all_R406W, "all_S305I.csv")

library ("FactoMineR") #for multivariate exploratory data analysis (MV-EDA), provides set of functions for PCA and clustering
library ("factoextra") #for extracting and visualizing the results of mulivariate data 


pca_data <- t(all_R406W [ ,14:21]) #the means for the 4 samples that I have 
PCA (pca_data, scale.unit = TRUE, ncp = 5, graph = TRUE) #scale unit is important to standardize all data (z-scores normalization)thus ensures that all data will contribute to the PCA making PCA more robust (decreases results skewness)
res.pca <- PCA(pca_data, graph = FALSE) #the results of PCA in a list

eig.val <- get_eigenvalue(res.pca) #eigenvalues represent the total amount of variance that can be explained by PCA. 
eig.val

fviz_eig (res.pca, addlabels = TRUE, ylim = c (0,50)) #to plot the eigenvalues/variance against the no. dimensions
var <- get_pca_var(res.pca)
var

#to identify 3 distinct groups (clusters) based on the similarity betweeen data points 
set.seed(123)
res.km <- kmeans (var$coord, centers = 3, nstart = 25) #to cluster the observations (var$coord) into 3 clusters based upon the nearest mean (centroid)
grp <- as.factor(res.km$cluster)

# coloring the clusters 
fviz_pca_var(res.pca, col.var = grp, palette = c ("#F2AFB4", "#9DC893", "#36454F"), legend.title = "Cluster")

# detailed description of the first 2 PCs  
res.desc <- dimdesc(res.pca, axes = c (1,2), proba = 0.05) #list for each dim contains the p-value and the correlation 

#to get info about the pca but for each dim (we have 2 dims and this will assign each sample with varaibles to specific location)
ind <- get_pca_ind(res.pca)
ind
grp <- as.factor (c ('mut1_R406W_2m', 'wt1_R406W_2m', 
                     'mut1_R406W_4m', 'wt1_R406W_4m',
                     'mut2_R406W_2m', 'wt2_R406W_2m',
                     'mut2_R406W_4m', 'wt2_R406W_4m'))

fviz_pca_ind(res.pca, fill.ind = grp, pointshape = 21, palette = c ("#F2AFB4", "#9DC893", "#36454F", "#005b96", "#265828", "#7293ec", "#790049", "#b21c0e"), repel = TRUE)


#t_test 
library(matrixTests)
library(plyr)
library(ggrepel)


all_R406W <- read.csv ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/all_R406W.csv", row.names = 1)
#create NPX df that contains the means for the protein in each sample
NPX <-all_R406W [, (14:21)]
NPX$PanelAssay<-all_R406W$PanelAssay
#define the con cloumns and mut columns separtely as it is paired t-tes 
#this arrangement is critical because we need to define the mean.diff (mean.diff= mean_mut - mean_con)

#to take in ur considerations the background 

mut2m_R406W <- c ('mut1_R406W_2m', 'mut2_R406W_2m')
wt2m_R406W <- c ('wt1_R406W_2m', 'wt2_R406W_2m')
mut4m_R406W <- c ('mut1_R406W_4m', 'mut2_R406W_4m')
wt4m_R406W <- c ('wt1_R406W_4m', 'wt2_R406W_4m')


#I will separate t-test according to the time point
#for 2m time point
t_test_protein_R406W_2m <- row_t_paired(NPX [,colnames(NPX) %in% mut2m_R406W],
                                        NPX [,colnames(NPX) %in% wt2m_R406W])

t_test_protein_R406W_2m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_R406W_2m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_R406W_2m$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_R406W_2m <- merge(uniprot, t_test_protein_R406W_2m, by = "PanelAssay", all.y = T)

t_test_uniprot_R406W_2m$diffexpressed <- "Not-significant"  
t_test_uniprot_R406W_2m$diffexpressed[t_test_uniprot_R406W_2m$mean.diff>0.0 & t_test_uniprot_R406W_2m$pvalue<0.05] <- "UP"
t_test_uniprot_R406W_2m$diffexpressed[t_test_uniprot_R406W_2m$mean.diff<0.0 & t_test_uniprot_R406W_2m$pvalue<0.05] <- "DOWN"
t_test_uniprot_R406W_2m$protein2 <- gsub(".*_", "", t_test_uniprot_R406W_2m$PanelAssay)
t_test_uniprot_R406W_2m$delabel <- NA
t_test_uniprot_R406W_2m$delabel[t_test_uniprot_R406W_2m$diffexpressed != "Not-significant"] <- t_test_uniprot_R406W_2m$protein2[t_test_uniprot_R406W_2m$diffexpressed != "Not-significant"]
t_test_uniprot_R406W_2m$diffexpressed <- factor(t_test_uniprot_R406W_2m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 

#enhanced vol plot
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_uniprot_R406W_2m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "R406W_2m",
                     labels = c("Up (24)", "Down (20)", "Not-sig (1293)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)





ggplot(data = t_test_uniprot_R406W_2m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 12, show.legend = F)+  
  scale_color_manual(name = "Protein expression_R406W_2m" , labels= c("Upregulated (n=24)", "Downregulated (n=20)", "Not significant (n=1293)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#for 4m time point
t_test_protein_R406W_4m <- row_t_paired(NPX [,colnames(NPX) %in% mut4m_R406W],
                                        NPX [,colnames(NPX) %in% wt4m_R406W])

t_test_protein_R406W_4m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_R406W_4m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_R406W_4m$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_R406W_4m <- merge(uniprot, t_test_protein_R406W_4m, by = "PanelAssay", all.y = T)

t_test_protein_R406W_4m$diffexpressed <- "Not-significant"  
t_test_protein_R406W_4m$diffexpressed[t_test_protein_R406W_4m$mean.diff>0.0 & t_test_protein_R406W_4m$pvalue<0.05] <- "UP"
t_test_protein_R406W_4m$diffexpressed[t_test_protein_R406W_4m$mean.diff<0.0 & t_test_protein_R406W_4m$pvalue<0.05] <- "DOWN"
t_test_protein_R406W_4m$protein2 <- gsub(".*_", "", t_test_protein_R406W_4m$PanelAssay)
t_test_protein_R406W_4m$delabel <- NA
t_test_protein_R406W_4m$delabel[t_test_protein_R406W_4m$diffexpressed != "Not-significant"] <- t_test_protein_R406W_4m$protein2[t_test_protein_R406W_4m$diffexpressed != "Not-significant"]
t_test_protein_R406W_4m$diffexpressed <- factor(t_test_protein_R406W_4m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 

#enhanced vol plot
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_R406W_4m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "R406W_4m",
                     labels = c("Up (9)", "Down (13)", "Not-sig (1315)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)






ggplot(data = t_test_protein_R406W_4m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression_R406W_4m" , labels= c("Upregulated (n=9)", "Downregulated (n=13)", "Not significant (n=1315)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#for V337M 
mut1_V337M_2m <- OlinkData@NPX [ , c("43C")]
wt1_V337M_2m <- OlinkData@NPX [ , c("41A")]
mut1_V337M_6m <- OlinkData@NPX [ , c("44D")]
wt1_V337M_6m <- OlinkData@NPX [ , c("42B")]


mut2_V337M_2m <- OlinkData@NPX [ , c("45E")]
wt2_V337M_2m <- OlinkData@NPX [ , c("48H")]
mut2_V337M_4m <- OlinkData@NPX [ , c("46F")]
wt2_V337M_4m <- OlinkData@NPX [ , c("49A")]
mut2_V337M_6m <- OlinkData@NPX [ , c("47G")]
wt2_V337M_6m <- OlinkData@NPX [ , c("50B")]


mut3_V337M_2m <- OlinkData@NPX [ , c("54F")]
wt3_V337M_2m <- OlinkData@NPX [ , c("51C")]
mut3_V337M_4m <- OlinkData@NPX [ , c("55G")]
wt3_V337M_4m <- OlinkData@NPX [ , c("53E")]


#to omit NA
mut1_V337M_2m <- na.omit(mut1_V337M_2m)
wt1_V337M_2m <- na.omit(wt1_V337M_2m)
mut1_V337M_6m <- na.omit(mut1_V337M_6m)
wt1_V337M_6m <- na.omit(wt1_V337M_6m)

mut2_V337M_2m <- na.omit(mut2_V337M_2m)
wt2_V337M_2m <- na.omit(wt2_V337M_2m)
mut2_V337M_4m <- na.omit(mut2_V337M_4m)
wt2_V337M_4m <- na.omit(wt2_V337M_4m)
mut2_V337M_6m <- na.omit(mut2_V337M_6m)
wt2_V337M_6m <- na.omit(wt2_V337M_6m)

mut3_V337M_2m <- na.omit(mut3_V337M_2m)
wt3_V337M_2m <- na.omit(wt3_V337M_2m)
mut3_V337M_4m <- na.omit(mut3_V337M_4m)
wt3_V337M_4m <- na.omit(wt3_V337M_4m)



#to calculate the mean 
mut1_V337M_2m <- as.data.frame(mut1_V337M_2m)
mut1_V337M_2m$panelassay <- rownames(mut1_V337M_2m)

wt1_V337M_2m <- as.data.frame(wt1_V337M_2m)
wt1_V337M_2m$panelassay <- rownames(wt1_V337M_2m)

mut1_V337M_6m <- as.data.frame(mut1_V337M_6m)
mut1_V337M_6m$panelassay <- rownames(mut1_V337M_6m)

wt1_V337M_6m <- as.data.frame(wt1_V337M_6m)
wt1_V337M_6m$panelassay <- rownames(wt1_V337M_6m)

mut2_V337M_2m <- as.data.frame(mut2_V337M_2m)
mut2_V337M_2m$panelassay <- rownames(mut2_V337M_2m)

wt2_V337M_2m <- as.data.frame(wt2_V337M_2m)
wt2_V337M_2m$panelassay <- rownames(wt2_V337M_2m)

mut2_V337M_4m <- as.data.frame(mut2_V337M_4m)
mut2_V337M_4m$panelassay <- rownames(mut2_V337M_4m)

wt2_V337M_4m <- as.data.frame(wt2_V337M_4m)
wt2_V337M_4m$panelassay <- rownames(wt2_V337M_4m)

mut2_V337M_6m <- as.data.frame(mut2_V337M_6m)
mut2_V337M_6m$panelassay <- rownames(mut2_V337M_6m)

wt2_V337M_6m <- as.data.frame(wt2_V337M_6m)
wt2_V337M_6m$panelassay <- rownames(wt2_V337M_6m)

mut3_V337M_2m <- as.data.frame(mut3_V337M_2m)
mut3_V337M_2m$panelassay <- rownames(mut3_V337M_2m)

wt3_V337M_2m <- as.data.frame(wt3_V337M_2m)
wt3_V337M_2m$panelassay <- rownames(wt3_V337M_2m)

mut3_V337M_4m <- as.data.frame(mut3_V337M_4m)
mut3_V337M_4m$panelassay <- rownames(mut3_V337M_4m)

wt3_V337M_4m <- as.data.frame(wt3_V337M_4m)
wt3_V337M_4m$panelassay <- rownames(wt3_V337M_4m)


V337M <- as.data.frame(cbind (mut1_V337M_2m$mut1_V337M_2m,
                              wt1_V337M_2m$wt1_V337M_2m,
                              mut1_V337M_6m$mut1_V337M_6m,
                              wt1_V337M_6m$wt1_V337M_6m,
                              mut2_V337M_2m$mut2_V337M_2m,
                              wt2_V337M_2m$wt2_V337M_2m,
                              mut2_V337M_4m$mut2_V337M_4m,
                              wt2_V337M_4m$wt2_V337M_4m,
                              mut2_V337M_6m$mut2_V337M_6m,
                              wt2_V337M_6m$wt2_V337M_6m,
                              mut3_V337M_2m$mut3_V337M_2m,
                              wt3_V337M_2m$wt3_V337M_2m,
                              mut3_V337M_4m$mut3_V337M_4m,
                              wt3_V337M_4m$wt3_V337M_4m))


colnames(V337M) <- c('mut1_V337M_2m', 'wt1_V337M_2m',
                     'mut1_V337M_6m', 'wt1_V337M_6m',
                     'mut2_V337M_2m', 'wt2_V337M_2m',
                     'mut2_V337M_4m', 'wt2_V337M_4m',
                     'mut2_V337M_6m', 'wt2_V337M_6m',
                     'mut3_V337M_2m', 'wt3_V337M_2m',
                     'mut3_V337M_4m', 'wt3_V337M_4m')

rownames(V337M) <- mut1_V337M_2m$panelassay 
V337M$PanelAssay <- mut1_V337M_2m$panelassay 

#t-test for IVS10+16 

#1.adding the UNIPROT accession ID to the proteins 
uniprot <- OlinkData@proteins
uniprot$protein <- paste0(uniprot$panel , "_", uniprot$Assay)


#2.the new dataframe (uniprot) contains all values including those that are below the LOD and those that are NA (missingness)
missingness <- OlinkData@prct.NPX.below.LOD #assigning the prct.NPX.below.LOD to the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m') #separate the wt from the mut


#the thing that I did that here is that I can't combine based on the proteins as some proteins are in duplicates, I need UNIQUE value for each thing so I will combine using the PanelAssay
all_V337M <- merge(uniprot, V337M, by="PanelAssay", all.y=TRUE)
write.csv(all_V337M, "all_V337M.csv")

library ("FactoMineR") #for multivariate exploratory data analysis (MV-EDA), provides set of functions for PCA and clustering
library ("factoextra") #for extracting and visualizing the results of mulivariate data 


pca_data <- t(all_V337M [ ,14:27]) #the means for the 4 samples that I have 
PCA (pca_data, scale.unit = TRUE, ncp = 5, graph = TRUE) #scale unit is important to standardize all data (z-scores normalization)thus ensures that all data will contribute to the PCA making PCA more robust (decreases results skewness)
res.pca <- PCA(pca_data, graph = FALSE) #the results of PCA in a list

eig.val <- get_eigenvalue(res.pca) #eigenvalues represent the total amount of variance that can be explained by PCA. 
eig.val

fviz_eig (res.pca, addlabels = TRUE, ylim = c (0,50)) #to plot the eigenvalues/variance against the no. dimensions
var <- get_pca_var(res.pca)
var

#to identify 3 distinct groups (clusters) based on the similarity betweeen data points 
set.seed(123)
res.km <- kmeans (var$coord, centers = 3, nstart = 25) #to cluster the observations (var$coord) into 3 clusters based upon the nearest mean (centroid)
grp <- as.factor(res.km$cluster)

# coloring the clusters 
fviz_pca_var(res.pca, col.var = grp, palette = c ("#F2AFB4", "#9DC893", "#36454F"), legend.title = "Cluster")

# detailed description of the first 2 PCs  
res.desc <- dimdesc(res.pca, axes = c (1,2), proba = 0.05) #list for each dim contains the p-value and the correlation 

#to get info about the pca but for each dim (we have 2 dims and this will assign each sample with varaibles to specific location)
ind <- get_pca_ind(res.pca)
ind
grp <- as.factor (c('mut1_V337M_2m', 'wt1_V337M_2m',
                    'mut1_V337M_6m', 'wt1_V337M_6m',
                    'mut2_V337M_2m', 'wt2_V337M_2m',
                    'mut2_V337M_4m', 'wt2_V337M_4m',
                    'mut2_V337M_6m', 'wt2_V337M_6m',
                    'mut3_V337M_2m', 'wt3_V337M_2m',
                    'mut3_V337M_4m', 'wt3_V337M_4m'))

fviz_pca_ind(res.pca, fill.ind = grp, pointshape = 21, palette = c("#F2AFB4", "#9DC893", "#36454F", "#005b96", "#265828", "#7293ec", "#790049"
                                                                   , "#0B6623", "#1A2421", "#317873", "#A9BA9D", "#36454F", "#ffc0cb","#59000f"), repel = TRUE)

#t_test 
library(matrixTests)
library(plyr)
library(ggrepel)


all_V337M <- read.csv ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/all_V337M.csv", row.names = 1)
#create NPX df that contains the means for the protein in each sample
NPX <-all_V337M [, (14:27)]
NPX$PanelAssay<-all_V337M$PanelAssay
#define the con cloumns and mut columns separtely as it is paired t-tes 
#this arrangement is critical because we need to define the mean.diff (mean.diff= mean_mut - mean_con)

#to take in ur considerations the background 

mut2m_V337M <- c ('mut1_V337M_2m', 'mut2_V337M_2m', 'mut3_V337M_2m')
wt2m_V337M <- c ('wt1_V337M_2m', 'wt2_V337M_2m', 'wt3_V337M_2m')
mut4m_V337M <- c ('mut2_V337M_4m', 'mut3_V337M_4m')
wt4m_V337M <- c ('wt2_V337M_4m', 'wt3_V337M_4m')
mut6m_V337M <- c ('mut1_V337M_6m', 'mut2_V337M_6m')
wt6m_V337M <- c ('wt1_V337M_6m', 'wt2_V337M_6m')


#I will separate t-test according to the time point
#for 2m time point
t_test_protein_V337M_2m <- row_t_paired(NPX [,colnames(NPX) %in% mut2m_V337M],
                                        NPX [,colnames(NPX) %in% wt2m_V337M])

t_test_protein_V337M_2m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_V337M_2m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_V337M_2m$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_V337M_2m <- merge(uniprot, t_test_protein_V337M_2m, by = "PanelAssay", all.y = T)



t_test_protein_V337M_2m$diffexpressed <- "Not-significant"  
t_test_protein_V337M_2m$diffexpressed[t_test_protein_V337M_2m$mean.diff>0.0 & t_test_protein_V337M_2m$pvalue<0.05] <- "UP"
t_test_protein_V337M_2m$diffexpressed[t_test_protein_V337M_2m$mean.diff<0.0 & t_test_protein_V337M_2m$pvalue<0.05] <- "DOWN"
t_test_protein_V337M_2m$protein2 <- gsub(".*_", "", t_test_protein_V337M_2m$PanelAssay)
t_test_protein_V337M_2m$delabel <- NA
t_test_protein_V337M_2m$delabel[t_test_protein_V337M_2m$diffexpressed != "Not-significant"] <- t_test_protein_V337M_2m$protein2[t_test_protein_V337M_2m$diffexpressed != "Not-significant"]
t_test_protein_V337M_2m$diffexpressed <- factor(t_test_protein_V337M_2m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 


Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_V337M_2m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "V337M_2m",
                     labels = c("Up (10)", "Down (28)", "Not-sig (1299)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)

ggplot(data = t_test_protein_V337M_2m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression_V337M_2m" , labels= c("Upregulated (n=10)", "Downregulated (n=28)", "Not significant (n=1299)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#for 4m timepoint 
t_test_protein_V337M_4m <- row_t_paired(NPX [,colnames(NPX) %in% mut4m_V337M],
                                        NPX [,colnames(NPX) %in% wt4m_V337M])

t_test_protein_V337M_4m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_V337M_4m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_V337M_4m$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_V337M_4m <- merge(uniprot, t_test_protein_V337M_4m, by = "PanelAssay", all.y = T)



t_test_protein_V337M_4m$diffexpressed <- "Not-significant"  
t_test_protein_V337M_4m$diffexpressed[t_test_protein_V337M_4m$mean.diff>0.0 & t_test_protein_V337M_4m$pvalue<0.05] <- "UP"
t_test_protein_V337M_4m$diffexpressed[t_test_protein_V337M_4m$mean.diff<0.0 & t_test_protein_V337M_4m$pvalue<0.05] <- "DOWN"
t_test_protein_V337M_4m$protein2 <- gsub(".*_", "", t_test_protein_V337M_4m$PanelAssay)
t_test_protein_V337M_4m$delabel <- NA
t_test_protein_V337M_4m$delabel[t_test_protein_V337M_4m$diffexpressed != "Not-significant"] <- t_test_protein_V337M_4m$protein2[t_test_protein_V337M_4m$diffexpressed != "Not-significant"]
t_test_protein_V337M_4m$diffexpressed <- factor(t_test_protein_V337M_4m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 


Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_V337M_4m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "V337M_4m",
                     labels = c("Up (35)", "Down (21)", "Not-sig (1281)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)



ggplot(data = t_test_protein_V337M_4m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 20, show.legend = F)+  
  scale_color_manual(name = "Protein expression_V337M_4m" , labels= c("Upregulated (n=35)", "Downregulated (n=21)", "Not significant (n=1281)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()

#for 6m timepoint 
t_test_protein_V337M_6m <- row_t_paired(NPX [,colnames(NPX) %in% mut6m_V337M],
                                        NPX [,colnames(NPX) %in% wt6m_V337M])

t_test_protein_V337M_6m$PanelAssay <- NPX$PanelAssay
p_value <- as.vector(t_test_protein_V337M_6m$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein_V337M_6m$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut_2m',
                           'prct.below.LOD.mut_4m',
                           'prct.below.LOD.mut_6m',
                           'prct.below.LOD.wt_2m',
                           'prct.below.LOD.wt_4m',
                           'prct.below.LOD.wt_6m')
t_test_uniprot_V337M_6m <- merge(uniprot, t_test_protein_V337M_6m, by = "PanelAssay", all.y = T)



t_test_protein_V337M_6m$diffexpressed <- "Not-significant"  
t_test_protein_V337M_6m$diffexpressed[t_test_protein_V337M_6m$mean.diff>0.0 & t_test_protein_V337M_6m$pvalue<0.05] <- "UP"
t_test_protein_V337M_6m$diffexpressed[t_test_protein_V337M_6m$mean.diff<0.0 & t_test_protein_V337M_6m$pvalue<0.05] <- "DOWN"
t_test_protein_V337M_6m$protein2 <- gsub(".*_", "", t_test_protein_V337M_6m$PanelAssay)
t_test_protein_V337M_6m$delabel <- NA
t_test_protein_V337M_6m$delabel[t_test_protein_V337M_6m$diffexpressed != "Not-significant"] <- t_test_protein_V337M_6m$protein2[t_test_protein_V337M_6m$diffexpressed != "Not-significant"]
t_test_protein_V337M_6m$diffexpressed <- factor(t_test_protein_V337M_6m$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 


Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein_V337M_6m, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "V337M_6m",
                     labels = c("Up (34)", "Down (15)", "Not-sig (1288)"),
                     values = c("#880808", "#6495ED", "#A9A9A9")) +
  labs(x = "Log2 Fold Change", y = "-log10(q-value)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  ) +
  geom_hline(yintercept = Significance_level, color = "black", linetype = "dashed") +
  geom_vline(xintercept = logFC_threshold, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -logFC_threshold, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(-log10(pvalue) > Significance_level & abs(mean.diff) > logFC_threshold, delabel, '')),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50'  # Color of the lines connecting points to labels
  )

# Print the plot
print(volcano_plot)




ggplot(data = t_test_protein_V337M_6m, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression_V337M_6m" , labels= c("Upregulated (n=34)", "Downregulated (n=15)", "Not significant (n=1288)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 18))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()




#biological analysis for each mutation 

HS = "org.Hs.eg.db"
#install.packages("BiocManager")
#BiocManager::install(HS, character.only = TRUE)
library(HS, character.only = TRUE)
library(org.Hs.eg.db)
library(clusterProfiler)
library (AnnotationDbi)
#BiocManager::install("clusterProfiler") #to be able to use gseGO function (it took toooooo long!)
#BiocManager::install("AnnotationDbi")
#BiocManager::install("biomaRt")
library("biomaRt")
#install.packages("ggridges")
library(ggridges)
#install.packages("dplyr")
library("dplyr")
#for visualiztion we need 
#install.packages("pheatmap")
#install.packages ("DOSE")
#install.packages("enrichplot")
#install.packages("ggupset")
#install.packages("UpSetR")


#BiocManager::install("pathview")
#install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)
library(ComplexHeatmap)
library(UpSetR)



#for IVS10 2m 

write.csv(t_test_uniprot_IVS10_2m, "GO_test_IVS10_2m.csv")
GO_test_IVS10_2m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_IVS10_2m.csv")
GO_test_IVS10_2m <- na.omit (GO_test_IVS10_2m)
row_protein_list <- GO_test_IVS10_2m$mean.diff #we need to use mean.diff 
protein_list_IVS10_2m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_IVS10_2m) <- GO_test_IVS10_2m$UNIPROT 
GO_test_IVS10_2m_unique <- GO_test_IVS10_2m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_IVS10_2m <- setNames(GO_test_IVS10_2m_unique$mean.diff, as.character(GO_test_IVS10_2m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_IVS10_2m_unique <- GO_test_IVS10_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_IVS10_2m <- setNames(GO_test_IVS10_2m_unique$mean.diff, as.character(GO_test_IVS10_2m_unique$ENTREZID))
protein_list_IVS10_2m <- sort (protein_list_IVS10_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_IVS10_2m <- subset(GO_test_IVS10_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)


#over representation analysis 
GO_test_IVS10_2m_unique <- GO_test_IVS10_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_IVS10_2m <- setNames(GO_test_IVS10_2m_unique$mean.diff, as.character(GO_test_IVS10_2m_unique$ENTREZID))
protein_list_IVS10_2m <- sort (protein_list_IVS10_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_IVS10_2m <- subset(GO_test_IVS10_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_IVS10_2m_unique <- sig_protein_list_IVS10_2m %>% distinct(UNIPROT, .keep_all = TRUE)
sig_protein_list_IVS10_2m_unique <- na.omit(sig_protein_list_IVS10_2m_unique)


go_enrich_IVS10_2m <- enrichGO(gene = as.character(sig_protein_list_IVS10_2m_unique$ENTREZID),
                               universe = as.character(GO_test_IVS10_2m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "MF", #for this one there was no gene can be maped for the PB and the MF 
                               pvalueCutoff = 0.9, #this is not significant  
                               qvalueCutoff = 0.9) 



dotplot(go_enrich_IVS10_2m)

IVS10_2m_go_enrich_MF <- go_enrich_IVS10_2m@result

IVS10_2m_go_enrich_MF <- IVS10_2m_go_enrich_MF[IVS10_2m_go_enrich_MF$pvalue < 0.05, ]
IVS10_2m_go_enrich_MF <- IVS10_2m_go_enrich_MF[order(IVS10_2m_go_enrich_MF$pvalue), ]
level_order <- rev(IVS10_2m_go_enrich_MF$Description)
IVS10_2m_go_enrich_MF$Description <- factor(IVS10_2m_go_enrich_MF$Description, levels = level_order)
IVS10_2m_dotplot <- ggplot(data = IVS10_2m_go_enrich_MF) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("IVS10_2m")
IVS10_2m_dotplot

library(writexl)

write.xlsx(IVS10_2m_go_enrich_MF, "IVS10_2m_go_enrich_MF.xlsx")

#to create other plots for the results of the GO 
wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_IVS10_2m_unique$SYMBOL)


#for IVS10 4m 

write.csv(t_test_uniprot_IVS10_4m, "GO_test_IVS10_4m.csv")
GO_test_IVS10_4m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_IVS10_4m.csv")
GO_test_IVS10_2m <- na.omit (GO_test_IVS10_2m)
row_protein_list_4m <- GO_test_IVS10_4m$mean.diff #we need to use mean.diff 
protein_list_IVS10_4m <- na.omit(row_protein_list_4m) #9proteins have been removed with final list 1328
names(protein_list_IVS10_4m) <- GO_test_IVS10_4m$UNIPROT 
GO_test_IVS10_4m_unique <- GO_test_IVS10_4m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_IVS10_4m <- setNames(GO_test_IVS10_4m_unique$mean.diff, as.character(GO_test_IVS10_4m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_IVS10_4m_unique <- GO_test_IVS10_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_IVS10_4m <- setNames(GO_test_IVS10_4m_unique$mean.diff, as.character(GO_test_IVS10_4m_unique$ENTREZID))
protein_list_IVS10_4m <- sort (protein_list_IVS10_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_IVS10_4m <- subset(GO_test_IVS10_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_IVS10_4m <- na.omit(sig_protein_list_IVS10_4m)


#over representation analysis 
GO_test_IVS10_4m_unique <- GO_test_IVS10_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_IVS10_4m <- setNames(GO_test_IVS10_4m_unique$mean.diff, as.character(GO_test_IVS10_4m_unique$ENTREZID))
protein_list_IVS10_4m <- sort (protein_list_IVS10_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_IVS10_4m <- subset(GO_test_IVS10_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_IVS10_4m_unique <- sig_protein_list_IVS10_4m %>% distinct(UNIPROT, .keep_all = TRUE)
sig_protein_list_IVS10_4m_unique <- na.omit (sig_protein_list_IVS10_4m_unique)


go_enrich_IVS10_4m <- enrichGO(gene = as.character(sig_protein_list_IVS10_4m_unique$ENTREZID),
                               universe = as.character(GO_test_IVS10_4m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "MF", #for this one there was no gene can be maaped for the PB and the MF 
                               pvalueCutoff = 0.7, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_IVS10_4m)



IVS10_4m_go_enrich_MF <- go_enrich_IVS10_4m@result

IVS10_4m_go_enrich_MF <- IVS10_4m_go_enrich_MF[IVS10_4m_go_enrich_MF$pvalue < 0.05, ]
IVS10_4m_go_enrich_MF <- IVS10_4m_go_enrich_MF[order(IVS10_4m_go_enrich_MF$pvalue), ]
level_order <- rev(IVS10_4m_go_enrich_MF$Description)
IVS10_4m_go_enrich_MF$Description <- factor(IVS10_4m_go_enrich_MF$Description, levels = level_order)
IVS10_4m_dotplot <- ggplot(data = IVS10_4m_go_enrich_MF) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("IVS10_4m")
IVS10_4m_dotplot



#to create other plots for the results of the GO 
wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_IVS10_4m_unique$SYMBOL)


#for S305I 2m

write.csv(t_test_uniprot_S305I_2m, "GO_test_S305I_2m.csv")
GO_test_S305I_2m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_S305I_2m.csv")
GO_test_S305I_2m <- na.omit (GO_test_S305I_2m)
row_protein_list <- GO_test_S305I_2m$mean.diff #we need to use mean.diff 
protein_list_S305I_2m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_S305I_2m) <- GO_test_S305I_2m$UNIPROT 
GO_test_S305I_2m_unique <- GO_test_S305I_2m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_S305I_2m <- setNames(GO_test_S305I_2m_unique$mean.diff, as.character(GO_test_S305I_2m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_S305I_2m_unique <- GO_test_S305I_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_S305I_2m <- setNames(GO_test_S305I_2m_unique$mean.diff, as.character(GO_test_S305I_2m_unique$ENTREZID))
protein_list_S305I_2m <- sort (protein_list_S305I_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_S305I_2m <- subset(GO_test_S305I_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
GO_test_S305I_2m_unique <- GO_test_S305I_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_S305I_2m <- setNames(GO_test_S305I_2m_unique$mean.diff, as.character(GO_test_S305I_2m_unique$ENTREZID))
protein_list_S305I_2m <- sort (protein_list_S305I_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_S305I_2m <- subset(GO_test_S305I_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_S305I_2m_unique <- sig_protein_list_S305I_2m %>% distinct(UNIPROT, .keep_all = TRUE)

go_enrich_S305I_2m <- enrichGO(gene = as.character(sig_protein_list_S305I_2m_unique$ENTREZID),
                               universe = as.character(GO_test_S305I_2m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "BP", 
                               pvalueCutoff = 0.7,   
                               qvalueCutoff = 0.9) 

dotplot (go_enrich_S305I_2m)

S305I_2m_go_enrich_BP <- go_enrich_S305I_2m@result

S305I_2m_go_enrich_BP <- S305I_2m_go_enrich_BP[S305I_2m_go_enrich_BP$pvalue < 0.05, ]
S305I_2m_go_enrich_BP <- S305I_2m_go_enrich_BP[order(S305I_2m_go_enrich_BP$pvalue), ]
level_order <- rev(S305I_2m_go_enrich_BP$Description)
S305I_2m_go_enrich_BP$Description <- factor(S305I_2m_go_enrich_BP$Description, levels = level_order)
S305I_2m_dotplot <- ggplot(data = S305I_2m_go_enrich_BP) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("S305I_2m")
S305I_2m_dotplot



wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_S305I_2m_unique$SYMBOL)


#for S305I 4m

write.csv(t_test_uniprot_S305I_4m, "GO_test_S305I_4m.csv")
GO_test_S305I_4m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_S305I_4m.csv")
GO_test_S305I_4m <- na.omit (GO_test_S305I_4m)
row_protein_list <- GO_test_S305I_4m$mean.diff #we need to use mean.diff 
protein_list_S305I_4m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_S305I_4m) <- GO_test_S305I_4m$UNIPROT 
GO_test_S305I_4m_unique <- GO_test_S305I_4m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_S305I_4m <- setNames(GO_test_S305I_4m_unique$mean.diff, as.character(GO_test_S305I_4m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_S305I_4m_unique <- GO_test_S305I_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_S305I_4m <- setNames(GO_test_S305I_4m_unique$mean.diff, as.character(GO_test_S305I_4m_unique$ENTREZID))
protein_list_S305I_4m <- sort (protein_list_S305I_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_S305I_4m <- subset(GO_test_S305I_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
GO_test_S305I_4m_unique <- GO_test_S305I_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_S305I_4m <- setNames(GO_test_S305I_4m_unique$mean.diff, as.character(GO_test_S305I_4m_unique$ENTREZID))
protein_list_S305I_4m <- sort (protein_list_S305I_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_S305I_4m <- subset(GO_test_S305I_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_S305I_4m_unique <- sig_protein_list_S305I_4m %>% distinct(UNIPROT, .keep_all = TRUE)

go_enrich_S305I_4m <- enrichGO(gene = as.character(sig_protein_list_S305I_4m_unique$ENTREZID),
                               universe = as.character(GO_test_S305I_4m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "MF", #for this one there was no gene can be maped for the PB and the MF 
                               pvalueCutoff = 0.9, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_S305I_4m)


S305I_4m_go_enrich_MF <- go_enrich_S305I_4m@result

S305I_4m_go_enrich_MF <- S305I_4m_go_enrich_MF[S305I_4m_go_enrich_MF$pvalue < 0.05, ]
S305I_4m_go_enrich_MF <- S305I_4m_go_enrich_MF[order(S305I_4m_go_enrich_MF$pvalue), ]
level_order <- rev(S305I_4m_go_enrich_MF$Description)
S305I_4m_go_enrich_MF$Description <- factor(S305I_4m_go_enrich_MF$Description, levels = level_order)
S305I_4m_dotplot <- ggplot(data = S305I_4m_go_enrich_MF) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("S305I_4m")
S305I_4m_dotplot



wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_S305I_4m_unique$SYMBOL)


#for R406W 2m 
write.csv(t_test_uniprot_R406W_2m, "GO_test_R406W_2m.csv")
GO_test_R406W_2m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_R406W_2m.csv")
GO_test_R406W_2m <- na.omit (GO_test_R406W_2m)
row_protein_list <- GO_test_R406W_2m$mean.diff #we need to use mean.diff 
protein_list_R406W_2m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_R406W_2m) <- GO_test_R406W_2m$UNIPROT 
GO_test_R406W_2m_unique <- GO_test_R406W_2m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_R406W_2m <- setNames(GO_test_R406W_2m_unique$mean.diff, as.character(GO_test_R406W_2m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_R406W_2m_unique <- GO_test_R406W_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_R406W_2m <- setNames(GO_test_R406W_2m_unique$mean.diff, as.character(GO_test_R406W_2m_unique$ENTREZID))
protein_list_R406W_2m <- sort (protein_list_R406W_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_R406W_2m <- subset(GO_test_R406W_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
GO_test_R406W_2m_unique <- GO_test_R406W_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_R406W_2m <- setNames(GO_test_R406W_2m_unique$mean.diff, as.character(GO_test_R406W_2m_unique$ENTREZID))
protein_list_R406W_2m <- sort (protein_list_R406W_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_R406W_2m <- subset(GO_test_R406W_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_R406W_2m_unique <- sig_protein_list_R406W_2m %>% distinct(UNIPROT, .keep_all = TRUE)

go_enrich_R406W_2m <- enrichGO(gene = as.character(sig_protein_list_R406W_2m_unique$ENTREZID),
                               universe = as.character(GO_test_R406W_2m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "BP", 
                               pvalueCutoff = 0.7, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_R406W_2m)

R406W_2m_go_enrich_BP <- go_enrich_R406W_2m@result

R406W_2m_go_enrich_BP <- R406W_2m_go_enrich_BP[R406W_2m_go_enrich_BP$pvalue < 0.05, ]
R406W_2m_go_enrich_BP <- R406W_2m_go_enrich_BP[order(R406W_2m_go_enrich_BP$pvalue), ]
level_order <- rev(R406W_2m_go_enrich_BP$Description)
R406W_2m_go_enrich_BP$Description <- factor(R406W_2m_go_enrich_BP$Description, levels = level_order)
R406W_2m_dotplot <- ggplot(data = R406W_2m_go_enrich_BP) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("R406W_2m")
R406W_2m_dotplot


wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_R406W_2m_unique$SYMBOL)



#for R406W 4m 
write.csv(t_test_uniprot_R406W_4m, "GO_test_R406W_4m.csv")
GO_test_R406W_4m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_R406W_4m.csv")
GO_test_R406W_4m <- na.omit (GO_test_R406W_4m)
row_protein_list <- GO_test_R406W_4m$mean.diff #we need to use mean.diff 
protein_list_R406W_4m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_R406W_4m) <- GO_test_R406W_4m$UNIPROT 
GO_test_R406W_4m_unique <- GO_test_R406W_4m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_R406W_4m <- setNames(GO_test_R406W_4m_unique$mean.diff, as.character(GO_test_R406W_4m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_R406W_4m_unique <- GO_test_R406W_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_R406W_4m <- setNames(GO_test_R406W_4m_unique$mean.diff, as.character(GO_test_R406W_4m_unique$ENTREZID))
protein_list_R406W_4m <- sort (protein_list_R406W_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_R406W_4m <- subset(GO_test_R406W_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
GO_test_R406W_4m_unique <- GO_test_R406W_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_R406W_4m <- setNames(GO_test_R406W_4m_unique$mean.diff, as.character(GO_test_R406W_4m_unique$ENTREZID))
protein_list_R406W_4m <- sort (protein_list_R406W_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_R406W_4m <- subset(GO_test_R406W_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_R406W_4m_unique <- sig_protein_list_R406W_4m %>% distinct(UNIPROT, .keep_all = TRUE)

go_enrich_R406W_4m <- enrichGO(gene = as.character(sig_protein_list_R406W_4m_unique$ENTREZID),
                               universe = as.character(GO_test_R406W_4m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "BP", 
                               pvalueCutoff = 0.9, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_R406W_4m)

R406W_4m_go_enrich_BP <- go_enrich_R406W_4m@result

R406W_4m_go_enrich_BP <- R406W_4m_go_enrich_BP[R406W_4m_go_enrich_BP$pvalue < 0.05, ]
R406W_4m_go_enrich_BP <- R406W_4m_go_enrich_BP[order(R406W_4m_go_enrich_BP$pvalue), ]
level_order <- rev(R406W_4m_go_enrich_BP$Description)
R406W_4m_go_enrich_BP$Description <- factor(R406W_4m_go_enrich_BP$Description, levels = level_order)
R406W_4m_dotplot <- ggplot(data = R406W_4m_go_enrich_BP) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("R406W_4m")
R406W_4m_dotplot





wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_R406W_4m_unique$SYMBOL)



#for V337M 2m 
write.csv(t_test_uniprot_V337M_2m, "GO_test_V337M_2m.csv")
GO_test_V337M_2m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_V337M_2m.csv")
GO_test_V337M_2m <- na.omit (GO_test_V337M_2m)
row_protein_list <- GO_test_V337M_2m$mean.diff #we need to use mean.diff 
protein_list_V337M_2m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_V337M_2m) <- GO_test_V337M_2m$UNIPROT 
GO_test_V337M_2m_unique <- GO_test_V337M_2m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_V337M_2m <- setNames(GO_test_V337M_2m_unique$mean.diff, as.character(GO_test_V337M_2m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_V337M_2m_unique <- GO_test_V337M_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_V337M_2m <- setNames(GO_test_V337M_2m_unique$mean.diff, as.character(GO_test_V337M_2m_unique$ENTREZID))
protein_list_V337M_2m <- sort (protein_list_V337M_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_V337M_2m <- subset(GO_test_V337M_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
GO_test_V337M_2m_unique <- GO_test_V337M_2m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_V337M_2m <- setNames(GO_test_V337M_2m_unique$mean.diff, as.character(GO_test_V337M_2m_unique$ENTREZID))
protein_list_V337M_2m <- sort (protein_list_V337M_2m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_V337M_2m <- subset(GO_test_V337M_2m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_V337M_2m_unique <- sig_protein_list_V337M_2m %>% distinct(UNIPROT, .keep_all = TRUE)

go_enrich_V337M_2m <- enrichGO(gene = as.character(sig_protein_list_V337M_2m_unique$ENTREZID),
                               universe = as.character(GO_test_V337M_2m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "BP", 
                               pvalueCutoff = 0.8, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_V337M_2m)

V337M_2m_go_enrich_BP <- go_enrich_V337M_2m@result

V337M_2m_go_enrich_BP <- V337M_2m_go_enrich_BP[V337M_2m_go_enrich_BP$pvalue < 0.05, ]
V337M_2m_go_enrich_BP <- V337M_2m_go_enrich_BP[order(V337M_2m_go_enrich_BP$pvalue), ]
level_order <- rev(V337M_2m_go_enrich_BP$Description)
V337M_2m_go_enrich_BP$Description <- factor(V337M_2m_go_enrich_BP$Description, levels = level_order)
V337M_2m_dotplot <- ggplot(data = V337M_2m_go_enrich_BP) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("V337M_2m")
V337M_2m_dotplot



wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_V337M_2m_unique$SYMBOL)




#for V337M 4m 
write.csv(t_test_uniprot_V337M_4m, "GO_test_V337M_4m.csv")
GO_test_V337M_4m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_V337M_4m.csv")
GO_test_V337M_4m <- na.omit (GO_test_V337M_4m)
row_protein_list <- GO_test_V337M_4m$mean.diff #we need to use mean.diff 
protein_list_V337M_4m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_V337M_4m) <- GO_test_V337M_4m$UNIPROT 
GO_test_V337M_4m_unique <- GO_test_V337M_4m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_V337M_4m <- setNames(GO_test_V337M_4m_unique$mean.diff, as.character(GO_test_V337M_4m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_V337M_4m_unique <- GO_test_V337M_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_V337M_4m <- setNames(GO_test_V337M_4m_unique$mean.diff, as.character(GO_test_V337M_4m_unique$ENTREZID))
protein_list_V337M_4m <- sort (protein_list_V337M_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_V337M_4m <- subset(GO_test_V337M_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
#for this timepoint there is no genes can be mapped at all!
GO_test_V337M_4m_unique <- GO_test_V337M_4m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_V337M_4m <- setNames(GO_test_V337M_4m_unique$mean.diff, as.character(GO_test_V337M_4m_unique$ENTREZID))
protein_list_V337M_4m <- sort (protein_list_V337M_4m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_V337M_4m <- subset(GO_test_V337M_4m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_V337M_4m_unique <- sig_protein_list_V337M_4m %>% distinct(UNIPROT, .keep_all = TRUE)


go_enrich_V337M_4m <- enrichGO(gene = as.character(sig_protein_list_V337M_4m_unique$ENTREZID),
                               universe = as.character(GO_test_V337M_4m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "MF", 
                               pvalueCutoff = 0.9, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_V337M_4m) 


V337M_4m_go_enrich_MF <- go_enrich_V337M_4m@result

V337M_4m_go_enrich_MF <- V337M_4m_go_enrich_MF[V337M_4m_go_enrich_MF$pvalue < 0.05, ]
V337M_4m_go_enrich_MF <- V337M_4m_go_enrich_MF[order(V337M_4m_go_enrich_MF$pvalue), ]
level_order <- rev(V337M_4m_go_enrich_MF$Description)
V337M_4m_go_enrich_MF$Description <- factor(V337M_4m_go_enrich_MF$Description, levels = level_order)
V337M_4m_dotplot <- ggplot(data = V337M_4m_go_enrich_MF) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("V337M_4m")
V337M_4m_dotplot 




wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_V337M_4m_unique$SYMBOL)




#for V337M 6m 
write.csv(t_test_uniprot_V337M_6m, "GO_test_V337M_6m.csv")
GO_test_V337M_6m <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_test_V337M_6m.csv")
GO_test_V337M_6m <- na.omit (GO_test_V337M_6m)
row_protein_list <- GO_test_V337M_6m$mean.diff #we need to use mean.diff 
protein_list_V337M_6m <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list_V337M_6m) <- GO_test_V337M_6m$UNIPROT 
GO_test_V337M_6m_unique <- GO_test_V337M_6m %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_V337M_6m <- setNames(GO_test_V337M_6m_unique$mean.diff, as.character(GO_test_V337M_6m_unique$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


GO_test_V337M_6m_unique <- GO_test_V337M_6m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_V337M_6m <- setNames(GO_test_V337M_6m_unique$mean.diff, as.character(GO_test_V337M_6m_unique$ENTREZID))
protein_list_V337M_6m <- sort (protein_list_V337M_6m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_V337M_6m <- subset(GO_test_V337M_6m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)



#over representation analysis 
GO_test_V337M_6m_unique <- GO_test_V337M_6m %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list_V337M_6m <- setNames(GO_test_V337M_6m_unique$mean.diff, as.character(GO_test_V337M_6m_unique$ENTREZID))
protein_list_V337M_6m <- sort (protein_list_V337M_6m, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list_V337M_6m <- subset(GO_test_V337M_6m_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list_V337M_6m_unique <- sig_protein_list_V337M_6m %>% distinct(UNIPROT, .keep_all = TRUE)


go_enrich_V337M_6m <- enrichGO(gene = as.character(sig_protein_list_V337M_6m_unique$ENTREZID),
                               universe = as.character(GO_test_V337M_6m_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "MF", 
                               pvalueCutoff = 0.9, #this is not significant  
                               qvalueCutoff = 0.9) 


dotplot(go_enrich_V337M_6m) #similar to what I got from the IVS10 2m 


V337M_6m_go_enrich_MF <- go_enrich_V337M_6m@result

V337M_6m_go_enrich_MF <- V337M_6m_go_enrich_MF[V337M_6m_go_enrich_MF$pvalue < 0.05, ]
V337M_6m_go_enrich_MF <- V337M_6m_go_enrich_MF[order(V337M_6m_go_enrich_MF$pvalue), ]
level_order <- rev(V337M_6m_go_enrich_MF$Description)
V337M_6m_go_enrich_MF$Description <- factor(V337M_6m_go_enrich_MF$Description, levels = level_order)
V337M_6m_dotplot <- ggplot(data = V337M_6m_go_enrich_MF) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("V337M_6m")
V337M_6m_dotplot


wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

goplot(go_enrich, showCategory = 10)

cnetplot(go_enrich, categorySize="pvalue", foldChange=GO_test_V337M_6m_unique$SYMBOL)




#heatmap 
install.packages("pheatmap")
library(pheatmap)
## try to select the right data - all data from the different mutations
h_IVS10_2m <- t_test_protein_IVS10_2m[ , c('pvalue', 'PanelAssay')]
names(h_IVS10_2m) <- c('IVS10_2m' , 'PanelAssay')


h_IVS10_4m <- t_test_protein_IVS10_4m[ , c('pvalue', 'PanelAssay')]
names(h_IVS10_4m) <- c('IVS10_4m' , 'PanelAssay')


h_S305I_2m <- t_test_protein_S305I_2m[ , c('pvalue', 'PanelAssay')]
names(h_S305I_2m) <- c('S305I_2m' , 'PanelAssay')


h_S305I_4m <- t_test_protein_S305I_4m[ , c('pvalue', 'PanelAssay')]
names(h_S305I_4m) <- c('S305I_4m', 'PanelAssay')


h_R406W_2m <- t_test_protein_R406W_2m[ , c('pvalue', 'PanelAssay')]
names(h_R406W_2m) <- c('R406W_2m', 'PanelAssay')


h_R406W_4m <- t_test_protein_R406W_4m[ , c('pvalue', 'PanelAssay')]
names(h_R406W_4m) <- c('R406W_4m', 'PanelAssay')


h_V337M_2m <- t_test_protein_V337M_2m[ , c('pvalue', 'PanelAssay')]
names(h_V337M_2m) <- c('V337M_2m', 'PanelAssay')


h_V337M_4m <- t_test_protein_V337M_4m[ , c('pvalue', 'PanelAssay')]
names(h_V337M_4m) <- c('V337M_4m', 'PanelAssay')


h_V337M_6m <- t_test_protein_V337M_6m[ , c('pvalue', 'PanelAssay')]
names(h_V337M_6m) <- c('V337M_6m', 'PanelAssay')


#library(dplyr)
#library(purrr)
#library(pheatmap)


mut_list <- list(h_IVS10_2m, h_IVS10_4m, h_S305I_2m, h_S305I_4m, h_R406W_2m, h_R406W_4m, h_V337M_2m, h_V337M_4m, h_V337M_6m)
combined_mut <- reduce(mut_list, full_join, by = "PanelAssay")
rownames(combined_mut) <- combined_mut$PanelAssay
combined_mut <- combined_mut[ , -2]
combined_mut <- as.matrix(combined_mut)
pheatmap(combined_mut)

# to enhance the vizulation of the heatmap and make it readable (didn't enhance it)
pheatmap(combined_mut,
         color = colorRampPalette(c("blue", "white", "red"))(100), #as R couldn't define the colors as a vector of colors (col transition frpm blut to white to red)
         scale = "row", #to scale the row to have a mean = 0 and var =1 
         clustering_distance_rows = "euclidean", #based on the values of the rows/cols, Uses Euclidean distance to cluster the rows
         clustering_distance_cols = "euclidean", 
         cluster_rows = TRUE, #clustering of rows 
         cluster_cols = TRUE, #clustering of cols 
         show_rownames = TRUE, #display row names 
         show_colnames = TRUE, #display col names 
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45, #rotation degree of the col names 
         main = "Heatmap of the organoids mutations",
         cellwidth = 10, # to make the width of the cell a little bit larger that makes better visualization 
         cellhight = 10)
#this heatmap is based on the p-value (mutations are not arranged and the proteins aren't readable)


library(ggplot2)

d_IVS10_2m <- t_test_protein_IVS10_2m[ , c('mean.x', 'PanelAssay')]
names(d_IVS10_2m) <- c('IVS10_2m' , 'PanelAssay')


d_IVS10_4m <- t_test_protein_IVS10_4m[ , c('mean.x', 'PanelAssay')]
names(d_IVS10_4m) <- c('IVS10_4m' , 'PanelAssay')


d_S305I_2m <- t_test_protein_S305I_2m[ , c('mean.x', 'PanelAssay')]
names(d_S305I_2m) <- c('S305I_2m' , 'PanelAssay')


d_S305I_4m <- t_test_protein_S305I_4m[ , c('mean.x', 'PanelAssay')]
names(d_S305I_4m) <- c('S305I_4m', 'PanelAssay')


d_R406W_2m <- t_test_protein_R406W_2m[ , c('mean.x', 'PanelAssay')]
names(d_R406W_2m) <- c('R406W_2m', 'PanelAssay')


d_R406W_4m <- t_test_protein_R406W_4m[ , c('mean.x', 'PanelAssay')]
names(d_R406W_4m) <- c('R406W_4m', 'PanelAssay')


d_V337M_2m <- t_test_protein_V337M_2m[ , c('mean.x', 'PanelAssay')]
names(d_V337M_2m) <- c('V337M_2m', 'PanelAssay')


d_V337M_4m <- t_test_protein_V337M_4m[ , c('mean.x', 'PanelAssay')]
names(d_V337M_4m) <- c('V337M_4m', 'PanelAssay')


d_V337M_6m <- t_test_protein_V337M_6m[ , c('mean.x', 'PanelAssay')]
names(d_V337M_6m) <- c('V337M_6m', 'PanelAssay')


density_list <- list(d_IVS10_2m, d_IVS10_4m, d_S305I_2m, d_S305I_4m, d_R406W_2m, d_R406W_4m, d_V337M_2m, d_V337M_4m, d_V337M_6m)
combined_density_list <- reduce(density_list, full_join, by = "PanelAssay")
rownames(combined_density_list) <- combined_density_list$PanelAssay
combined_density_list <- combined_density_list[ , -2]
combined_density_list$Protein <- rownames(combined_density_list)


library(tidyr)
library(ggplot2)
library(dplyr)

#to covert he data frame (wide) to a long data frame to have a column for the NPX, a column for time point values, and a column for the protein names
combined_density_long <- combined_density_list %>%
  pivot_longer(
    cols = -Protein,  # Exclude the Protein column from the pivot
    names_to = "TimePoint",
    values_to = "NPX"
  )

#across the different time points 
ggplot(combined_density_long, aes(x = NPX, fill = TimePoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of NPX Across Time Points",
       x = "NPX",
       y = "Density") +
  theme_minimal()

#there is a general trend toward the upregulation (shift to the right, with a general normal distribution.

ggplot(combined_density_long, aes(x = NPX, fill = TimePoint)) +
  geom_density(alpha = 0.5, adjust = 0.5) + #to make it smaller (to make sure if there's any trends or changes compared to the snooth one)
  labs(title = "Density of NPX Across Time Points",
       x = "NPX",
       y = "Density") +
  theme_minimal()


ggplot(combined_density_long, aes(x = NPX, fill = TimePoint)) +
  geom_density(alpha = 0.5, adjust = 2) +
  labs(title = "Density of NPX Across Time Points",
       x = "NPX",
       y = "Density") +
  theme_minimal()


#the density plot based on the mean.diff

library(ggplot2)

x_IVS10_2m <- t_test_protein_IVS10_2m[ , c('mean.diff', 'PanelAssay')]
names(x_IVS10_2m) <- c('IVS10_2m' , 'PanelAssay')


x_IVS10_4m <- t_test_protein_IVS10_4m[ , c('mean.diff', 'PanelAssay')]
names(x_IVS10_4m) <- c('IVS10_4m' , 'PanelAssay')


x_S305I_2m <- t_test_protein_S305I_2m[ , c('mean.diff', 'PanelAssay')]
names(x_S305I_2m) <- c('S305I_2m' , 'PanelAssay')


x_S305I_4m <- t_test_protein_S305I_4m[ , c('mean.diff', 'PanelAssay')]
names(x_S305I_4m) <- c('S305I_4m', 'PanelAssay')


x_R406W_2m <- t_test_protein_R406W_2m[ , c('mean.diff', 'PanelAssay')]
names(x_R406W_2m) <- c('R406W_2m', 'PanelAssay')


x_R406W_4m <- t_test_protein_R406W_4m[ , c('mean.diff', 'PanelAssay')]
names(x_R406W_4m) <- c('R406W_4m', 'PanelAssay')


x_V337M_2m <- t_test_protein_V337M_2m[ , c('mean.diff', 'PanelAssay')]
names(x_V337M_2m) <- c('V337M_2m', 'PanelAssay')


x_V337M_4m <- t_test_protein_V337M_4m[ , c('mean.diff', 'PanelAssay')]
names(x_V337M_4m) <- c('V337M_4m', 'PanelAssay')


x_V337M_6m <- t_test_protein_V337M_6m[ , c('mean.diff', 'PanelAssay')]
names(x_V337M_6m) <- c('V337M_6m', 'PanelAssay')


x_density_list <- list(x_IVS10_2m, x_IVS10_4m, x_S305I_2m, x_S305I_4m, x_R406W_2m, x_R406W_4m, x_V337M_2m, x_V337M_4m, x_V337M_6m)
combined2_density_list <- reduce(x_density_list, full_join, by = "PanelAssay")
rownames(combined2_density_list) <- combined2_density_list$PanelAssay
combined2_density_list <- combined2_density_list[ , -2]
combined2_density_list$Protein <- rownames(combined2_density_list)


library(tidyr)
library(ggplot2)
library(dplyr)

#to covert he data frame (wide) to a long data frame to have a column for the NPX, a column for time point values, and a column for the protein names
combined2_density_list <- combined2_density_list %>%
  pivot_longer(
    cols = -Protein,  # Exclude the Protein column from the pivot
    names_to = "TimePoint",
    values_to = "NPX"
  )

#across the different time points 
ggplot(combined2_density_list, aes(x = NPX, fill = TimePoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of NPX Across Time Points",
       x = "NPX",
       y = "Density") +
  theme_minimal()

#there is a general trend toward the upregulation (shift to the right, with a general normal distribution.

ggplot(combined2_density_list, aes(x = NPX, fill = TimePoint)) +
  geom_density(alpha = 0.5, adjust = 0.5) + #to make it smaller (to make sure if there's any trends or changes compared to the snooth one)
  labs(title = "Density of NPX Across Time Points",
       x = "NPX",
       y = "Density") +
  theme_minimal()


ggplot(combined2_density_list, aes(x = NPX, fill = TimePoint)) +
  geom_density(alpha = 0.5, adjust = 2) +
  labs(title = "Density of NPX Across Time Points",
       x = "NPX",
       y = "Density") +
  theme_minimal()


#for the result that I got and according to the t_test that I did, I need to combine the NPX values and take the mean
#for 2m 
I_IVS10_2m <- IVS10[, c('mut1_IVS10_2m', 'wt1_IVS10_2m', 'PanelAssay')]
names(I_IVS10_2m) <- c('mut1_IVS10_2m', 'wt1_IVS10_2m', 'PanelAssay')


II_IVS10_2m <- IVS10[, c('mut2_IVS10_2m', 'wt2_IVS10_2m', 'PanelAssay')]
names(II_IVS10_2m) <- c('mut2_IVS10_2m', 'wt2_IVS10_2m', 'PanelAssay')

significant_pro <- t_test_protein_IVS10_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro) <- c('pvalue', 'PanelAssay')

mut_IVS10 <- IVS10[, c('mut1_IVS10_2m', 'mut2_IVS10_2m', 'PanelAssay')]
names(mut_IVS10) <- c('mut1_IVS10_2m', 'mut2_IVS10_2m', 'PanelAssay')
mut_IVS10 <- IVS10[, c('mut1_IVS10_2m', 'mut2_IVS10_2m')]
mut_IVS10$mut_IVS10 <- rowMeans(mut_IVS10, na.rm = TRUE)
mut_IVS10 <- cbind(mut_IVS10, IVS10$PanelAssay)
names(mut_IVS10)[ncol(mut_IVS10)] <- "PanelAssay"


wt_IVS10 <- IVS10[, c('wt1_IVS10_2m', 'wt2_IVS10_2m', 'PanelAssay')]
names(wt_IVS10) <- c('wt1_IVS10_2m', 'wt2_IVS10_2m', 'PanelAssay')
wt_IVS10 <- IVS10[, c('wt1_IVS10_2m', 'wt2_IVS10_2m')]
wt_IVS10$wt_IVS10 <- rowMeans(wt_IVS10, na.rm = TRUE)
wt_IVS10 <- cbind(wt_IVS10, IVS10$PanelAssay)
names(wt_IVS10)[ncol(wt_IVS10)] <- "PanelAssay"



mut_IVS10 <- mut_IVS10[, c('mut_IVS10', 'PanelAssay')]
names(mut_IVS10) <- c('mut_IVS10', 'PanelAssay')


wt_IVS10 <- wt_IVS10[, c('wt_IVS10', 'PanelAssay')]
names(wt_IVS10) <- c('wt_IVS10', 'PanelAssay')


significant_pro_IVS10_2m <- t_test_protein_IVS10_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_IVS10_2m) <- c('pvalue', 'PanelAssay')


combined_IVS10_2m <- list(mut_IVS10, wt_IVS10, significant_pro_IVS10_2m) 
combined_IVS10_2m <- reduce(combined_IVS10_2m, full_join, by = "PanelAssay")
rownames(combined_IVS10_2m) <- combined_IVS10_2m$PanelAssay
combined_IVS10_2m <- combined_IVS10_2m[ , -2]
combined_IVS10_2m <- as.matrix(combined_IVS10_2m)
sig_matrix_2 <- combined_IVS10_2m[combined_IVS10_2m[, "pvalue"] < 0.05,]
sig_matrix_2 <- sig_matrix_2[, -3] #to remove the p-value column 


pheatmap(sig_matrix_2,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in IVS10 2m",
         cellwidth = 10,
         cellheight = 10)




#for 4m
I_IVS10_4m <- IVS10[, c('mut1_IVS10_4m', 'wt1_IVS10_4m', 'PanelAssay')]
names(I_IVS10_4m) <- c('mut1_IVS10_4m', 'wt1_IVS10_4m', 'PanelAssay')

II_IVS10_4m <- IVS10[, c('mut2_IVS10_4m', 'wt2_IVS10_4m', 'PanelAssay')]
names(II_IVS10_4m) <- c('mut2_IVS10_4m', 'wt2_IVS10_4m', 'PanelAssay')

significant_pro_IVS10_4m <- t_test_protein_IVS10_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_IVS10_4m) <- c('pvalue', 'PanelAssay')

mut_IVS10_4m <- IVS10[, c('mut1_IVS10_4m', 'mut2_IVS10_4m', 'PanelAssay')]
names(mut_IVS10_4m) <- c('mut1_IVS10_4m', 'mut2_IVS10_4m', 'PanelAssay')
mut_IVS10_4m <- IVS10[, c('mut1_IVS10_4m', 'mut2_IVS10_4m')]
mut_IVS10_4m$mut_IVS10_4m <- rowMeans(mut_IVS10_4m, na.rm = TRUE)
mut_IVS10_4m <- cbind(mut_IVS10_4m, IVS10$PanelAssay)
names(mut_IVS10_4m)[ncol(mut_IVS10_4m)] <- "PanelAssay"


wt_IVS10_4m <- IVS10[, c('wt1_IVS10_4m', 'wt2_IVS10_4m', 'PanelAssay')]
names(wt_IVS10_4m) <- c('wt1_IVS10_4m', 'wt2_IVS10_4m', 'PanelAssay')
wt_IVS10_4m <- IVS10[, c('wt1_IVS10_4m', 'wt2_IVS10_4m')]
wt_IVS10_4m$wt_IVS10_4m <- rowMeans(wt_IVS10_4m, na.rm = TRUE)
wt_IVS10_4m <- cbind(wt_IVS10_4m, IVS10$PanelAssay)
names(wt_IVS10_4m)[ncol(wt_IVS10_4m)] <- "PanelAssay"


mut_IVS10_4m <- mut_IVS10_4m[, c('mut_IVS10_4m', 'PanelAssay')]
names(mut_IVS10_4m) <- c('mut_IVS10_4m', 'PanelAssay')


wt_IVS10_4m <- wt_IVS10_4m[, c('wt_IVS10_4m', 'PanelAssay')]
names(wt_IVS10_4m) <- c('wt_IVS10_4m', 'PanelAssay')


significant_pro_IVS10_4m <- t_test_protein_IVS10_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_IVS10_4m) <- c('pvalue', 'PanelAssay')


combined_IVS10_4m <- list(mut_IVS10_4m, wt_IVS10_4m, significant_pro_IVS10_4m) 
combined_IVS10_4m <- reduce(combined_IVS10_4m, full_join, by = "PanelAssay")
rownames(combined_IVS10_4m) <- combined_IVS10_4m$PanelAssay
combined_IVS10_4m <- combined_IVS10_4m[ , -2]
combined_IVS10_4m <- as.matrix(combined_IVS10_4m)
sig_matrix_4 <- combined_IVS10_4m[combined_IVS10_4m[, "pvalue"] < 0.05,]
sig_matrix_4 <- sig_matrix_4[, -3] #to remove the p-value column 


pheatmap(sig_matrix_4,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in IVS10 4m",
         cellwidth = 10,
         cellheight = 10)


#for S305I_2m 

I_S305I_2m <- S305I[, c('mut1_S305I_2m', 'wt1_S305I_2m', 'PanelAssay')]
names(I_S305I_2m) <- c('mut1_S305I_2m', 'wt1_S305I_2m', 'PanelAssay')

II_S305I_2m <- S305I[, c('mut2_S305I_2m', 'wt2_S305I_2m', 'PanelAssay')]
names(II_S305I_2m) <- c('mut2_S305I_2m', 'wt2_S305I_2m', 'PanelAssay')

significant_pro_S305I_2m <- t_test_protein_S305I_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_S305I_2m) <- c('pvalue', 'PanelAssay')

mut_S305I_2m <- S305I[, c('mut1_S305I_2m', 'mut2_S305I_2m', 'PanelAssay')]
names(mut_S305I_2m) <- c('mut1_S305I_2m', 'mut2_S305I_2m', 'PanelAssay')
mut_S305I_2m <- S305I[, c('mut1_S305I_2m', 'mut2_S305I_2m')]
mut_S305I_2m$mut_S305I_2m <- rowMeans(mut_S305I_2m, na.rm = TRUE)
mut_S305I_2m <- cbind(mut_S305I_2m, S305I$PanelAssay)
names(mut_S305I_2m)[ncol(mut_S305I_2m)] <- "PanelAssay"


wt_S305I_2m <- S305I[, c('wt1_S305I_2m', 'wt2_S305I_2m', 'PanelAssay')]
names(wt_S305I_2m) <- c('wt1_S305I_2m', 'wt2_S305I_2m', 'PanelAssay')
wt_S305I_2m <- S305I[, c('wt1_S305I_2m', 'wt2_S305I_2m')]
wt_S305I_2m$wt_S305I_2m <- rowMeans(wt_S305I_2m, na.rm = TRUE)
wt_S305I_2m <- cbind(wt_S305I_2m, S305I$PanelAssay)
names(wt_S305I_2m)[ncol(wt_S305I_2m)] <- "PanelAssay"


mut_S305I_2m <- mut_S305I_2m[, c('mut_S305I_2m', 'PanelAssay')]
names(mut_S305I_2m) <- c('mut_S305I_2m', 'PanelAssay')


wt_S305I_2m <- wt_S305I_2m[, c('wt_S305I_2m', 'PanelAssay')]
names(wt_S305I_2m) <- c('wt_S305I_2m', 'PanelAssay')


significant_pro_S305I_2m <- t_test_protein_S305I_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_S305I_2m) <- c('pvalue', 'PanelAssay')


combined_S305I_2m <- list(mut_S305I_2m, wt_S305I_2m, significant_pro_S305I_2m) 
combined_S305I_2m <- reduce(combined_S305I_2m, full_join, by = "PanelAssay")
rownames(combined_S305I_2m) <- combined_S305I_2m$PanelAssay
combined_S305I_2m <- combined_S305I_2m[ , -2]
combined_S305I_2m <- as.matrix(combined_S305I_2m)
sig_matrix_S305I_2m <- combined_S305I_2m[combined_S305I_2m[, "pvalue"] < 0.05,]
sig_matrix_S305I_2m <- sig_matrix_S305I_2m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_S305I_2m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in S305I 2m",
         cellwidth = 10,
         cellheight = 10)



#for S305I_4m 

I_S305I_4m <- S305I[, c('mut1_S305I_4m', 'wt1_S305I_4m', 'PanelAssay')]
names(I_S305I_4m) <- c('mut1_S305I_4m', 'wt1_S305I_4m', 'PanelAssay')

II_S305I_4m <- S305I[, c('mut2_S305I_4m', 'wt2_S305I_4m', 'PanelAssay')]
names(II_S305I_4m) <- c('mut2_S305I_4m', 'wt2_S305I_4m', 'PanelAssay')

significant_pro_S305I_4m <- t_test_protein_S305I_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_S305I_4m) <- c('pvalue', 'PanelAssay')

mut_S305I_4m <- S305I[, c('mut1_S305I_4m', 'mut2_S305I_4m', 'PanelAssay')]
names(mut_S305I_4m) <- c('mut1_S305I_4m', 'mut2_S305I_4m', 'PanelAssay')
mut_S305I_4m <- S305I[, c('mut1_S305I_4m', 'mut2_S305I_4m')]
mut_S305I_4m$mut_S305I_4m <- rowMeans(mut_S305I_4m, na.rm = TRUE)
mut_S305I_4m <- cbind(mut_S305I_4m, S305I$PanelAssay)
names(mut_S305I_4m)[ncol(mut_S305I_4m)] <- "PanelAssay"


wt_S305I_4m <- S305I[, c('wt1_S305I_4m', 'wt2_S305I_4m', 'PanelAssay')]
names(wt_S305I_4m) <- c('wt1_S305I_4m', 'wt2_S305I_4m', 'PanelAssay')
wt_S305I_4m <- S305I[, c('wt1_S305I_4m', 'wt2_S305I_4m')]
wt_S305I_4m$wt_S305I_4m <- rowMeans(wt_S305I_4m, na.rm = TRUE)
wt_S305I_4m <- cbind(wt_S305I_4m, S305I$PanelAssay)
names(wt_S305I_4m)[ncol(wt_S305I_4m)] <- "PanelAssay"


mut_S305I_4m <- mut_S305I_4m[, c('mut_S305I_4m', 'PanelAssay')]
names(mut_S305I_4m) <- c('mut_S305I_4m', 'PanelAssay')


wt_S305I_4m <- wt_S305I_4m[, c('wt_S305I_4m', 'PanelAssay')]
names(wt_S305I_4m) <- c('wt_S305I_4m', 'PanelAssay')


significant_pro_S305I_4m <- t_test_protein_S305I_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_S305I_4m) <- c('pvalue', 'PanelAssay')


combined_S305I_4m <- list(mut_S305I_4m, wt_S305I_4m, significant_pro_S305I_4m) 
combined_S305I_4m <- reduce(combined_S305I_4m, full_join, by = "PanelAssay")
rownames(combined_S305I_4m) <- combined_S305I_4m$PanelAssay
combined_S305I_4m <- combined_S305I_4m[ , -2]
combined_S305I_4m <- as.matrix(combined_S305I_4m)
sig_matrix_S305I_4m <- combined_S305I_4m[combined_S305I_4m[, "pvalue"] < 0.05,]
sig_matrix_S305I_4m <- sig_matrix_S305I_4m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_S305I_4m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in S305I 4m",
         cellwidth = 10,
         cellheight = 10)



#for R406W_2m 
I_R406W_2m <- R406W[, c('mut1_R406W_2m', 'wt1_R406W_2m', 'PanelAssay')]
names(I_R406W_2m) <- c('mut1_R406W_2m', 'wt1_R406W_2m', 'PanelAssay')

II_R406W_2m <- R406W[, c('mut2_R406W_2m', 'wt2_R406W_2m', 'PanelAssay')]
names(II_R406W_2m) <- c('mut2_R406W_2m', 'wt2_R406W_2m', 'PanelAssay')

significant_pro_R406W_2m <- t_test_protein_R406W_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_R406W_2m) <- c('pvalue', 'PanelAssay')

mut_R406W_2m <- R406W[, c('mut1_R406W_2m', 'mut2_R406W_2m', 'PanelAssay')]
names(mut_R406W_2m) <- c('mut1_R406W_2m', 'mut2_R406W_2m', 'PanelAssay')
mut_R406W_2m <- R406W[, c('mut1_R406W_2m', 'mut2_R406W_2m')]
mut_R406W_2m$mut_R406W_2m <- rowMeans(mut_R406W_2m, na.rm = TRUE)
mut_R406W_2m <- cbind(mut_R406W_2m, R406W$PanelAssay)
names(mut_R406W_2m)[ncol(mut_R406W_2m)] <- "PanelAssay"


wt_R406W_2m <- R406W[, c('wt1_R406W_2m', 'wt2_R406W_2m', 'PanelAssay')]
names(wt_R406W_2m) <- c('wt1_R406W_2m', 'wt2_R406W_2m', 'PanelAssay')
wt_R406W_2m <- R406W[, c('wt1_R406W_2m', 'wt2_R406W_2m')]
wt_R406W_2m$wt_R406W_2m <- rowMeans(wt_R406W_2m, na.rm = TRUE)
wt_R406W_2m <- cbind(wt_R406W_2m, R406W$PanelAssay)
names(wt_R406W_2m)[ncol(wt_R406W_2m)] <- "PanelAssay"


mut_R406W_2m <- mut_R406W_2m[, c('mut_R406W_2m', 'PanelAssay')]
names(mut_R406W_2m) <- c('mut_R406W_2m', 'PanelAssay')


wt_R406W_2m <- wt_R406W_2m[, c('wt_R406W_2m', 'PanelAssay')]
names(wt_R406W_2m) <- c('wt_R406W_2m', 'PanelAssay')


significant_pro_R406W_2m <- t_test_protein_R406W_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_R406W_2m) <- c('pvalue', 'PanelAssay')


combined_R406W_2m <- list(mut_R406W_2m, wt_R406W_2m, significant_pro_R406W_2m) 
combined_R406W_2m <- reduce(combined_R406W_2m, full_join, by = "PanelAssay")
rownames(combined_R406W_2m) <- combined_R406W_2m$PanelAssay
combined_R406W_2m <- combined_R406W_2m[ , -2]
combined_R406W_2m <- as.matrix(combined_R406W_2m)
sig_matrix_R406W_2m <- combined_R406W_2m[combined_R406W_2m[, "pvalue"] < 0.05,]
sig_matrix_R406W_2m <- sig_matrix_R406W_2m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_R406W_2m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in R406W 2m",
         cellwidth = 10,
         cellheight = 10)


#for R406W_4m 
I_R406W_4m <- R406W[, c('mut1_R406W_4m', 'wt1_R406W_4m', 'PanelAssay')]
names(I_R406W_4m) <- c('mut1_R406W_4m', 'wt1_R406W_4m', 'PanelAssay')

II_R406W_4m <- R406W[, c('mut2_R406W_4m', 'wt2_R406W_4m', 'PanelAssay')]
names(II_R406W_4m) <- c('mut2_R406W_4m', 'wt2_R406W_4m', 'PanelAssay')

significant_pro_R406W_4m <- t_test_protein_R406W_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_R406W_4m) <- c('pvalue', 'PanelAssay')

mut_R406W_4m <- R406W[, c('mut1_R406W_4m', 'mut2_R406W_4m', 'PanelAssay')]
names(mut_R406W_4m) <- c('mut1_R406W_4m', 'mut2_R406W_4m', 'PanelAssay')
mut_R406W_4m <- R406W[, c('mut1_R406W_4m', 'mut2_R406W_4m')]
mut_R406W_4m$mut_R406W_4m <- rowMeans(mut_R406W_4m, na.rm = TRUE)
mut_R406W_4m <- cbind(mut_R406W_4m, R406W$PanelAssay)
names(mut_R406W_4m)[ncol(mut_R406W_4m)] <- "PanelAssay"


wt_R406W_4m <- R406W[, c('wt1_R406W_4m', 'wt2_R406W_4m', 'PanelAssay')]
names(wt_R406W_4m) <- c('wt1_R406W_4m', 'wt2_R406W_4m', 'PanelAssay')
wt_R406W_4m <- R406W[, c('wt1_R406W_4m', 'wt2_R406W_4m')]
wt_R406W_4m$wt_R406W_4m <- rowMeans(wt_R406W_4m, na.rm = TRUE)
wt_R406W_4m <- cbind(wt_R406W_4m, R406W$PanelAssay)
names(wt_R406W_4m)[ncol(wt_R406W_4m)] <- "PanelAssay"


mut_R406W_4m <- mut_R406W_4m[, c('mut_R406W_4m', 'PanelAssay')]
names(mut_R406W_4m) <- c('mut_R406W_4m', 'PanelAssay')


wt_R406W_4m <- wt_R406W_4m[, c('wt_R406W_4m', 'PanelAssay')]
names(wt_R406W_4m) <- c('wt_R406W_4m', 'PanelAssay')


significant_pro_R406W_4m <- t_test_protein_R406W_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_R406W_4m) <- c('pvalue', 'PanelAssay')


combined_R406W_4m <- list(mut_R406W_4m, wt_R406W_4m, significant_pro_R406W_4m) 
combined_R406W_4m <- reduce(combined_R406W_4m, full_join, by = "PanelAssay")
rownames(combined_R406W_4m) <- combined_R406W_4m$PanelAssay
combined_R406W_4m <- combined_R406W_4m[ , -2]
combined_R406W_4m <- as.matrix(combined_R406W_4m)
sig_matrix_R406W_4m <- combined_R406W_4m[combined_R406W_4m[, "pvalue"] < 0.05,]
sig_matrix_R406W_4m <- sig_matrix_R406W_4m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_R406W_4m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in R406W 4m",
         cellwidth = 10,
         cellheight = 10)



#for V337M_2m 
I_V337M_2m <- V337M[, c('mut1_V337M_2m', 'wt1_V337M_2m', 'PanelAssay')]
names(I_V337M_2m) <- c('mut1_V337M_2m', 'wt1_V337M_2m', 'PanelAssay')

II_V337M_2m <- V337M[, c('mut2_V337M_2m', 'wt2_V337M_2m', 'PanelAssay')]
names(II_V337M_2m) <- c('mut2_V337M_2m', 'wt2_V337M_2m', 'PanelAssay')

significant_pro_V337M_2m <- t_test_protein_V337M_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_V337M_2m) <- c('pvalue', 'PanelAssay')

mut_V337M_2m <- V337M[, c('mut1_V337M_2m', 'mut2_V337M_2m', 'PanelAssay')]
names(mut_V337M_2m) <- c('mut1_V337M_2m', 'mut2_V337M_2m', 'PanelAssay')
mut_V337M_2m <- V337M[, c('mut1_V337M_2m', 'mut2_V337M_2m')]
mut_V337M_2m$mut_V337M_2m <- rowMeans(mut_V337M_2m, na.rm = TRUE)
mut_V337M_2m <- cbind(mut_V337M_2m, V337M$PanelAssay)
names(mut_V337M_2m)[ncol(mut_V337M_2m)] <- "PanelAssay"


wt_V337M_2m <- V337M[, c('wt1_V337M_2m', 'wt2_V337M_2m', 'PanelAssay')]
names(wt_V337M_2m) <- c('wt1_V337M_2m', 'wt2_V337M_2m', 'PanelAssay')
wt_V337M_2m <- V337M[, c('wt1_V337M_2m', 'wt2_V337M_2m')]
wt_V337M_2m$wt_V337M_2m <- rowMeans(wt_V337M_2m, na.rm = TRUE)
wt_V337M_2m <- cbind(wt_V337M_2m, V337M$PanelAssay)
names(wt_V337M_2m)[ncol(wt_V337M_2m)] <- "PanelAssay"


mut_V337M_2m <- mut_V337M_2m[, c('mut_V337M_2m', 'PanelAssay')]
names(mut_V337M_2m) <- c('mut_V337M_2m', 'PanelAssay')


wt_V337M_2m <- wt_V337M_2m[, c('wt_V337M_2m', 'PanelAssay')]
names(wt_V337M_2m) <- c('wt_V337M_2m', 'PanelAssay')


significant_pro_V337M_2m <- t_test_protein_V337M_2m[, c('pvalue', 'PanelAssay')]
names(significant_pro_V337M_2m) <- c('pvalue', 'PanelAssay')


combined_V337M_2m <- list(mut_V337M_2m, wt_V337M_2m, significant_pro_V337M_2m) 
combined_V337M_2m <- reduce(combined_V337M_2m, full_join, by = "PanelAssay")
rownames(combined_V337M_2m) <- combined_V337M_2m$PanelAssay
combined_V337M_2m <- combined_V337M_2m[ , -2]
combined_V337M_2m <- as.matrix(combined_V337M_2m)
sig_matrix_V337M_2m <- combined_V337M_2m[combined_V337M_2m[, "pvalue"] < 0.05,]
sig_matrix_V337M_2m <- sig_matrix_V337M_2m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_V337M_2m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in V337M 2m",
         cellwidth = 10,
         cellheight = 10)


#for V337M_4m 
II_V337M_4m <- V337M[, c('mut2_V337M_4m', 'wt2_V337M_4m', 'PanelAssay')]
names(II_V337M_4m) <- c('mut2_V337M_4m', 'wt2_V337M_4m', 'PanelAssay')

III_V337M_4m <- V337M[, c('mut3_V337M_4m', 'wt3_V337M_4m', 'PanelAssay')]
names(III_V337M_4m) <- c('mut3_V337M_4m', 'wt3_V337M_4m', 'PanelAssay')

significant_pro_V337M_4m <- t_test_protein_V337M_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_V337M_4m) <- c('pvalue', 'PanelAssay')

mut_V337M_4m <- V337M[, c('mut2_V337M_4m', 'mut3_V337M_4m', 'PanelAssay')]
names(mut_V337M_4m) <- c('mut2_V337M_4m', 'mut3_V337M_4m', 'PanelAssay')
mut_V337M_4m <- V337M[, c('mut2_V337M_4m', 'mut3_V337M_4m')]
mut_V337M_4m$mut_V337M_4m <- rowMeans(mut_V337M_4m, na.rm = TRUE)
mut_V337M_4m <- cbind(mut_V337M_4m, V337M$PanelAssay)
names(mut_V337M_4m)[ncol(mut_V337M_4m)] <- "PanelAssay"


wt_V337M_4m <- V337M[, c('wt2_V337M_4m', 'wt3_V337M_4m', 'PanelAssay')]
names(wt_V337M_4m) <- c('wt2_V337M_4m', 'wt3_V337M_4m', 'PanelAssay')
wt_V337M_4m <- V337M[, c('wt2_V337M_4m', 'wt3_V337M_4m')]
wt_V337M_4m$wt_V337M_4m <- rowMeans(wt_V337M_4m, na.rm = TRUE)
wt_V337M_4m <- cbind(wt_V337M_4m, V337M$PanelAssay)
names(wt_V337M_4m)[ncol(wt_V337M_4m)] <- "PanelAssay"


mut_V337M_4m <- mut_V337M_4m[, c('mut_V337M_4m', 'PanelAssay')]
names(mut_V337M_4m) <- c('mut_V337M_4m', 'PanelAssay')


wt_V337M_4m <- wt_V337M_4m[, c('wt_V337M_4m', 'PanelAssay')]
names(wt_V337M_4m) <- c('wt_V337M_4m', 'PanelAssay')


significant_pro_V337M_4m <- t_test_protein_V337M_4m[, c('pvalue', 'PanelAssay')]
names(significant_pro_V337M_4m) <- c('pvalue', 'PanelAssay')


combined_V337M_4m <- list(mut_V337M_4m, wt_V337M_4m, significant_pro_V337M_4m) 
combined_V337M_4m <- reduce(combined_V337M_4m, full_join, by = "PanelAssay")
rownames(combined_V337M_4m) <- combined_V337M_4m$PanelAssay
combined_V337M_4m <- combined_V337M_4m[ , -2]
combined_V337M_4m <- as.matrix(combined_V337M_4m)
sig_matrix_V337M_4m <- combined_V337M_4m[combined_V337M_4m[, "pvalue"] < 0.05,]
sig_matrix_V337M_4m <- sig_matrix_V337M_4m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_V337M_4m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in V337M 4m",
         cellwidth = 10,
         cellheight = 10)

#for V337M_6m 
I_V337M_6m <- V337M[, c('mut1_V337M_6m', 'wt1_V337M_6m', 'PanelAssay')]
names(I_V337M_6m) <- c('mut1_V337M_6m', 'wt1_V337M_6m', 'PanelAssay')

II_V337M_6m <- V337M[, c('mut2_V337M_6m', 'wt2_V337M_6m', 'PanelAssay')]
names(II_V337M_6m) <- c('mut2_V337M_6m', 'wt2_V337M_6m', 'PanelAssay')

significant_pro_V337M_6m <- t_test_protein_V337M_6m[, c('pvalue', 'PanelAssay')]
names(significant_pro_V337M_6m) <- c('pvalue', 'PanelAssay')

mut_V337M_6m <- V337M[, c('mut1_V337M_6m', 'mut2_V337M_6m', 'PanelAssay')]
names(mut_V337M_6m) <- c('mut1_V337M_6m', 'mut2_V337M_6m', 'PanelAssay')
mut_V337M_6m <- V337M[, c('mut1_V337M_6m', 'mut2_V337M_6m')]
mut_V337M_6m$mut_V337M_6m <- rowMeans(mut_V337M_6m, na.rm = TRUE)
mut_V337M_6m <- cbind(mut_V337M_6m, V337M$PanelAssay)
names(mut_V337M_6m)[ncol(mut_V337M_6m)] <- "PanelAssay"


wt_V337M_6m <- V337M[, c('wt1_V337M_6m', 'wt2_V337M_6m', 'PanelAssay')]
names(wt_V337M_6m) <- c('wt1_V337M_6m', 'wt2_V337M_6m', 'PanelAssay')
wt_V337M_6m <- V337M[, c('wt1_V337M_6m', 'wt2_V337M_6m')]
wt_V337M_6m$wt_V337M_6m <- rowMeans(wt_V337M_6m, na.rm = TRUE)
wt_V337M_6m <- cbind(wt_V337M_6m, V337M$PanelAssay)
names(wt_V337M_6m)[ncol(wt_V337M_6m)] <- "PanelAssay"


mut_V337M_6m <- mut_V337M_6m[, c('mut_V337M_6m', 'PanelAssay')]
names(mut_V337M_6m) <- c('mut_V337M_6m', 'PanelAssay')


wt_V337M_6m <- wt_V337M_6m[, c('wt_V337M_6m', 'PanelAssay')]
names(wt_V337M_6m) <- c('wt_V337M_6m', 'PanelAssay')


significant_pro_V337M_6m <- t_test_protein_V337M_6m[, c('pvalue', 'PanelAssay')]
names(significant_pro_V337M_6m) <- c('pvalue', 'PanelAssay')


combined_V337M_6m <- list(mut_V337M_6m, wt_V337M_6m, significant_pro_V337M_6m) 
combined_V337M_6m <- reduce(combined_V337M_6m, full_join, by = "PanelAssay")
rownames(combined_V337M_6m) <- combined_V337M_6m$PanelAssay
combined_V337M_6m <- combined_V337M_6m[ , -2]
combined_V337M_6m <- as.matrix(combined_V337M_6m)
sig_matrix_V337M_6m <- combined_V337M_6m[combined_V337M_6m[, "pvalue"] < 0.05,]
sig_matrix_V337M_6m <- sig_matrix_V337M_6m[, -3] #to remove the p-value column 


pheatmap(sig_matrix_V337M_6m,
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color gradient
         scale = "row", # Scale by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         main = "Heatmap of Significant Proteins in V337M 6m",
         cellwidth = 10,
         cellheight = 10)



#to extract the significant protein list for each time point 
#for IVS10_2
sig_protein_list_IVS10_2m_unique$UniProt <- row.names(sig_protein_list_IVS10_2m_unique)
pro_list_IVS10_2m <- sig_protein_list_IVS10_2m_unique [c("UNIPROT")]
#write_delim(pro_list_IVS10_2m, "pro_list_IVS10_2m.txt", delim = "\t")

#for IVS10_4m 
sig_protein_list_IVS10_4m$UniProt <- row.names(sig_protein_list_IVS10_4m)
pro_list_IVS10_4m <- sig_protein_list_IVS10_4m [c("UNIPROT")]
#write_delim(pro_list_IVS10_4m, "pro_list_IVS10_4m.txt", delim = "\t")


#for S305I_2m 
sig_protein_list_S305I_2m$UniProt <- row.names(sig_protein_list_S305I_2m)
pro_list_S305I_2m <- sig_protein_list_S305I_2m [c("UNIPROT")]
#write_delim(pro_list_S305I_2m, "pro_list_S305I_2m.txt", delim = "\t")


#for S305I_4m
sig_protein_list_S305I_4m$UniProt <- row.names(sig_protein_list_S305I_4m)
pro_list_S305I_4m <- sig_protein_list_S305I_4m [c("UNIPROT")]
#write_delim(pro_list_S305I_4m, "pro_list_S305I_4m.txt", delim = "\t")

#for R406W_2m 
sig_protein_list_R406W_2m$UniProt <- row.names(sig_protein_list_R406W_2m)
pro_list_R406W_2m <- sig_protein_list_R406W_2m [c("UNIPROT")]
#write_delim(pro_list_R406W_2m, "pro_list_R406W_2m.txt", delim = "\t")


#for R406W_4m
sig_protein_list_R406W_4m$UniProt <- row.names(sig_protein_list_R406W_4m)
pro_list_R406W_4m <- sig_protein_list_R406W_4m [c("UNIPROT")]
#write_delim(pro_list_R406W_4m, "pro_list_R406W_4m.txt", delim = "\t")


#for V337M_2m 
sig_protein_list_V337M_2m$UniProt <- row.names(sig_protein_list_V337M_2m)
pro_list_V337M_2m <- sig_protein_list_V337M_2m [c("UNIPROT")]
#write_delim(pro_list_V337M_2m, "pro_list_V337M_2m.txt", delim = "\t")


#for V337M_4m
sig_protein_list_V337M_4m$UniProt <- row.names(sig_protein_list_V337M_4m)
pro_list_V337M_4m <- sig_protein_list_V337M_4m [c("UNIPROT")]
#write_delim(pro_list_V337M_4m, "pro_list_V337M_4m.txt", delim = "\t")


#for V337M_6m
sig_protein_list_V337M_6m$UniProt <- row.names(sig_protein_list_V337M_6m)
pro_list_V337M_6m <- sig_protein_list_V337M_6m [c("UNIPROT")]
#write_delim(pro_list_V337M_6m, "pro_list_V337M_6m.txt", delim = "\t")


# to generate a data frame that combines all significant proteins
#1. add a column that specifies the mutation and time point
sig_protein_list_IVS10_2m <- na.omit(sig_protein_list_IVS10_2m)
sig_protein_list_S305I_4m <- na.omit(sig_protein_list_S305I_4m)
sig_protein_list_V337M_4m <- na.omit(sig_protein_list_V337M_4m)


sig_protein_list_IVS10_2m$Mutation_TimePoint <- "IVS10_2m"
sig_protein_list_IVS10_4m$Mutation_TimePoint <- "IVS10_4m"
sig_protein_list_S305I_2m$Mutation_TimePoint <- "S305I_2m"
sig_protein_list_S305I_4m$Mutation_TimePoint <- "S305I_4m"
GO_test_R406W_2m_unique$Mutation_TimePoint <- "R406W_2m"
sig_protein_list_R406W_4m$Mutation_TimePoint <- "R406W_4m"
sig_protein_list_V337M_2m$Mutation_TimePoint <- "V337M_2m"
sig_protein_list_V337M_4m$Mutation_TimePoint <- "V337M_4m"
sig_protein_list_V337M_6m$Mutation_TimePoint <- "V337M_6m"

#2.combine all data frames 
combined_sig_proteins <- bind_rows(
  sig_protein_list_IVS10_2m,
  sig_protein_list_IVS10_4m,
  sig_protein_list_S305I_2m,
  sig_protein_list_S305I_4m,
  GO_test_R406W_2m_unique,
  sig_protein_list_R406W_4m,
  sig_protein_list_V337M_2m,
  sig_protein_list_V337M_4m,
  sig_protein_list_V337M_6m)



#3. remove the duplicates from the protein (I had an error for that)
combined_sig_proteins <- combined_sig_proteins[!duplicated(combined_sig_proteins$Protein), ] #314 protein remains out of 405 

rownames(combined_sig_proteins) <- combined_sig_proteins$Protein


combined_sig_proteins <- combined_sig_proteins[, -which(names(combined_sig_proteins) == "delabel")] #to remove the delabel column 

combined_sig_proteins <- combined_sig_proteins[, -which(names(combined_sig_proteins) == "protein2")]# to remove the protein2 column (I already have the names)

combined_sig_proteins <- combined_sig_proteins[, -which(names(combined_sig_proteins) == "diffexpressed")]# to remove the diffexpressed column 


selected_columns <- c ("Mutation_TimePoint", "pvalue", "PanelAssay") #the columns that I need in my new data frame 

subset_combined <- combined_sig_proteins %>% select(all_of(selected_columns))



pheatmap(
  combined_sig_proteins,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = TRUE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Significant Proteins"
)



#I need to extract the NPX values for all mutations and time points 

IVS10_2m_NPX <- subset (all_IVS10, select = c('PanelAssay', 'mut1_IVS10_2m', 'wt1_IVS10_2m', 'mut2_IVS10_2m', 'wt2_IVS10_2m'))

IVS10_4m_NPX <- subset (all_IVS10, select = c('PanelAssay', 'mut1_IVS10_4m', 'wt1_IVS10_4m', 'mut2_IVS10_4m', 'wt2_IVS10_4m'))


S305I_2m_NPX <- subset(all_S305I, select = c('PanelAssay', 'mut1_S305I_2m', 'wt1_S305I_2m', 'mut2_S305I_2m', 'wt2_S305I_2m'))

S305I_4m_NPX <- subset(all_S305I, select = c('PanelAssay', 'mut1_S305I_4m', 'wt1_S305I_4m', 'mut2_S305I_4m', 'wt2_S305I_4m'))


R406W_2m_NPX <- subset(all_R406W, select = c('PanelAssay', 'mut1_R406W_2m', 'wt1_R406W_2m',  'mut2_R406W_2m', 'wt2_R406W_2m'))

R406W_4m_NPX <- subset(all_R406W, select = c('PanelAssay', 'mut1_R406W_4m', 'wt1_R406W_4m', 'mut2_R406W_4m', 'wt2_R406W_4m')) 


#V337M_2m_NPX <- subset(all_V337M, select = c('PanelAssay', 'wt1_V337M_2m', 'wt2_V337M_2m', 'wt3_V337M_2m',
                                                'wt2_V337M_4m', 'wt3_V337M_4m',
                                                'wt1_V337M_6m','wt2_V337M_6m',
                                                'mut1_V337M_2m','mut2_V337M_2m','mut3_V337M_2m',
                                                'mut2_V337M_4m', 'mut3_V337M_4m',
                                                'mut1_V337M_6m','mut2_V337M_6m'))


V337M_2m_NPX <- subset(all_V337M, select = c('PanelAssay', 'mut1_V337M_2m', 'wt1_V337M_2m', 'mut2_V337M_2m', 'wt2_V337M_2m', 'mut3_V337M_2m', 'wt3_V337M_2m'))

V337M_4m_NPX <- subset(all_V337M, select = c('PanelAssay','mut2_V337M_4m', 'wt2_V337M_4m', 'mut3_V337M_4m', 'wt3_V337M_4m'))

V337M_6m_NPX <- subset(all_V337M, select = c('PanelAssay', 'mut1_V337M_6m', 'wt1_V337M_6m', 'mut2_V337M_6m', 'wt2_V337M_6m'))



IVS10_2m_hmap <- cbind (IVS10_2m_NPX, sig_protein_list_IVS10_2m)
#IVS10_2m_hmap <- IVS10_2m_hmap[IVS10_2m_hmap[, "pvalue"] < 0.05, ]
IVS10_2m_hmap <- IVS10_2m_hmap[, !duplicated(colnames(IVS10_2m_hmap))]



IVS10_4m_hmap <- cbind (IVS10_4m_NPX, significant_pro_IVS10_4m)
#IVS10_4m_hmap <- IVS10_4m_hmap[IVS10_4m_hmap[, "pvalue"] < 0.05, ]
IVS10_4m_hmap <- IVS10_4m_hmap[, !duplicated(colnames(IVS10_4m_hmap))]


S305I_2m_hmap <- cbind (S305I_2m_NPX, significant_pro_S305I_2m)
#S305I_2m_hmap <- S305I_2m_hmap[S305I_2m_hmap[, "pvalue"] < 0.05, ]
S305I_2m_hmap <- S305I_2m_hmap[, !duplicated(colnames(S305I_2m_hmap))]


S305I_4m_hmap <- cbind (S305I_4m_NPX, significant_pro_S305I_4m)
#S305I_4m_hmap <- S305I_4m_hmap[S305I_4m_hmap[, "pvalue"] < 0.05, ]
S305I_4m_hmap <- S305I_4m_hmap[, !duplicated(colnames(S305I_4m_hmap))]


R406W_2m_hmap <- cbind (R406W_2m_NPX, significant_pro_R406W_2m)
#R406W_2m_hmap <- R406W_2m_hmap[R406W_2m_hmap[, "pvalue"] < 0.05, ]
R406W_2m_hmap <- R406W_2m_hmap[, !duplicated(colnames(R406W_2m_hmap))]


R406W_4m_hmap <- cbind (R406W_4m_NPX, significant_pro_R406W_4m)
#R406W_4m_hmap <- R406W_4m_hmap[R406W_4m_hmap[, "pvalue"] < 0.05, ]
R406W_4m_hmap <- R406W_4m_hmap[, !duplicated(colnames(R406W_4m_hmap))]


V337M_2m_hmap <- cbind (V337M_2m_NPX, significant_pro_V337M_2m)
#V337M_2m_hmap <- V337M_2m_hmap[V337M_2m_hmap[, "pvalue"] < 0.05, ]
V337M_2m_hmap <- V337M_2m_hmap[, !duplicated(colnames(V337M_2m_hmap))]


V337M_4m_hmap <- cbind (V337M_4m_NPX, significant_pro_V337M_4m)
#V337M_4m_hmap <- V337M_4m_hmap[V337M_4m_hmap[, "pvalue"] < 0.05, ]
V337M_4m_hmap <- V337M_4m_hmap[, !duplicated(colnames(V337M_4m_hmap))]


V337M_6m_hmap <- cbind (V337M_6m_NPX, significant_pro_V337M_6m)
#V337M_6m_hmap <- V337M_6m_hmap[V337M_6m_hmap[, "pvalue"] < 0.05, ]
V337M_6m_hmap <- V337M_6m_hmap[, !duplicated(colnames(V337M_6m_hmap))]



NPX_list <- list(IVS10_2m_hmap, IVS10_4m_hmap, S305I_2m_hmap, S305I_4m_hmap, R406W_2m_hmap, R406W_4m_hmap, V337M_2m_hmap, V337M_4m_hmap, V337M_6m_hmap)
combined_NPX <- reduce(NPX_list, full_join, by = "PanelAssay")
rownames(combined_NPX) <- combined_NPX$PanelAssay
combined_NPX$PanelAssay <- rownames(combined_NPX)


#NPX_list <- list(IVS10_2m_NPX, IVS10_4m_NPX, S305I_2m_NPX, S305I_4m_NPX, R406W_2m_NPX, R406W_4m_NPX, V337M_2m_NPX, V337M_4m_NPX, V337M_6m_NPX)
NPX_list <- list( V337M_4m_NPX, V337M_6m_NPX)

combined_NPX <- V337M_2m_NPX

combined_NPX <- reduce(NPX_list, full_join, by = "PanelAssay")



## make a list of the right protein names V337M_6m
V337M_6m_hmap_pvalue <- V337M_6m_hmap[V337M_6m_hmap[, "pvalue"] < 0.05, ]
V337M_6m_p_value_protein_names <- c(V337M_6m_hmap_pvalue$PanelAssay)

## make a list of the right protein names V337M_4m
V337M_4m_hmap_pvalue <- V337M_4m_hmap[V337M_4m_hmap[, "pvalue"] < 0.05, ]
V337M_4m_p_value_protein_names <- c(V337M_4m_hmap_pvalue$PanelAssay)

## make a list of the right protein names V337M_2m 
V337M_2m_hmap_pvalue <- V337M_2m_hmap[V337M_2m_hmap[, "pvalue"] < 0.05, ]
V337M_2m_p_value_protein_names <- c(V337M_2m_hmap_pvalue$PanelAssay)


p_value_hmap_protein_names_V337M <-  c(V337M_2m_p_value_protein_names,
                                       V337M_4m_p_value_protein_names,
                                       V337M_6m_p_value_protein_names)
dups <- p_value_hmap_protein_names_V337M[duplicated(p_value_hmap_protein_names_V337M)]
p_value_hmap_protein_names_V337M <-  c(V337M_6m_p_value_protein_names)


rownames(combined_NPX) <- combined_NPX$PanelAssay
#combined_NPX$PanelAssay <- rownames(combined_NPX)

combined_NPX <- combined_NPX[combined_NPX$PanelAssay %in% p_value_hmap_protein_names_V337M, ] 
combined_NPX <- combined_NPX[combined_NPX$PanelAssay %in% dups, ] 

combined_NPX <- combined_NPX[, -which(names(combined_NPX) == "PanelAssay")] 

diff_V337M_4m_2 <- combined_NPX$mut2_V337M_4m - combined_NPX$wt2_V337M_4m
diff_V337M_4m_3 <- combined_NPX$mut3_V337M_4m - combined_NPX$wt3_V337M_4m
diff_V337M_6m_1 <- combined_NPX$mut1_V337M_6m - combined_NPX$wt1_V337M_6m
diff_V337M_6m_2 <- combined_NPX$mut2_V337M_6m - combined_NPX$wt2_V337M_6m

diff_V337M <- cbind(diff_V337M_4m_2,diff_V337M_4m_3, diff_V337M_6m_1, diff_V337M_6m_2)

combined_NPX_test <- as.matrix(combined_NPX)

pheatmap(
 combined_NPX,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of V337M"
)


 #what if I combined all the data sets and then selected the sig proteins onlys based on p-vlue
#I have multiple columns for p_value based on the time point and the selected p-value will be the one for V337M_6m as it is the last one 
combined_NPX <- reduce(NPX_list, full_join, by = "PanelAssay")
rownames(combined_NPX) <- combined_NPX$PanelAssay

combined_NPX <- combined_NPX[combined_NPX [, "pvalue"] < 0.05, ]
#combined_NPX <- combined_NPX[ , -2]
#to covert he data frame (wide) to a long data frame to have a column for the NPX, a column for time point values, and a column for the protein names
combined_NPX <- combined_NPX %>%
  pivot_longer(
    cols = -PanelAssay,  # Exclude the Protein column from the pivot
    names_to = "TimePoint",
    values_to = "NPX"
  )






#for IVS10
#to calculate the mean.diff
IVS10_NPX_list <- list (IVS10_2m_NPX, IVS10_4m_NPX)
combined_NPX_IVS10 <- reduce(IVS10_NPX_list, full_join, by = "PanelAssay")


IVS10_2m_1 <- combined_NPX_IVS10$mut1_IVS10_2m - combined_NPX_IVS10$wt1_IVS10_2m
IVS10_2m_2 <- combined_NPX_IVS10$mut2_IVS10_2m - combined_NPX_IVS10$wt2_IVS10_2m
IVS10_4m_1 <- combined_NPX_IVS10$mut1_IVS10_4m - combined_NPX_IVS10$wt1_IVS10_4m
IVS10_4m_2 <- combined_NPX_IVS10$mut2_IVS10_4m - combined_NPX_IVS10$wt1_IVS10_4m

diff_IVS10 <- cbind(IVS10_2m_1, IVS10_2m_2, IVS10_4m_1, IVS10_4m_2)
rownames(diff_IVS10) <- combined_NPX$PanelAssay

IVS10_hmap <- cbind(diff_IVS10, t_test_protein_IVS10_2m, by = "PanelAssay")
#IVS10_hmap <- IVS10_hmap[, c("IVS10_2m_1", "IVS10_2m_2", "IVS10_4m_1", "IVS10_4m_2", "PanelAssay", "pvalue", "pvalue")]

IVS10_hmap <- IVS10_hmap[, c("IVS10_2m_1", "IVS10_4m_1", "IVS10_2m_2", "IVS10_4m_2", "PanelAssay", "pvalue", "pvalue")]


#IVS10_hmap <- IVS10_hmap %>%
  rename(pvalue_2m = pvalue,
         pvalue_4m = pvalue.1)

IVS10_hmap <- IVS10_hmap[, -which(names(IVS10_hmap) == "PanelAssay")] 

#for IVS10_2m 

IVS10_hmap <- IVS10_hmap[IVS10_hmap [, "pvalue"] < 0.05, ]
IVS10_hmap <- IVS10_hmap[, -which(names(IVS10_hmap) == "pvalue")] 
IVS10_hmap_2m <- IVS10_hmap[, -which(names(IVS10_hmap) == "pvalue.1")] 


IVS10_hmap_2m <- as.matrix(IVS10_hmap_2m)

pheatmap(
  IVS10_hmap_2m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid IVS10_2m"
)





#for IVS10_4m 
diff_IVS10 <- cbind(IVS10_2m_1, IVS10_2m_2, IVS10_4m_1, IVS10_4m_2)
rownames(diff_IVS10) <- combined_NPX$PanelAssay

IVS10_hmap <- cbind(diff_IVS10, t_test_protein_IVS10_4m, by = "PanelAssay")
IVS10_hmap <- IVS10_hmap[, c("IVS10_2m_1", "IVS10_2m_2", "IVS10_4m_1", "IVS10_4m_2", "PanelAssay", "pvalue", "pvalue")]

#IVS10_hmap <- IVS10_hmap[, c("IVS10_2m_1", "IVS10_4m_1", "IVS10_2m_2", "IVS10_4m_2", "PanelAssay", "pvalue", "pvalue")]


IVS10_hmap <- IVS10_hmap[, -which(names(IVS10_hmap) == "PanelAssay")] 
IVS10_hmap <- IVS10_hmap[, -which(names(IVS10_hmap) == "pvalue.1")] 

IVS10_hmap <- IVS10_hmap[IVS10_hmap [, "pvalue"] < 0.05, ]

IVS10_hmap_4m <- IVS10_hmap[, -which(names(IVS10_hmap) == "pvalue")] 


IVS10_hmap_4m <- as.matrix(IVS10_hmap_4m)


pheatmap(
  IVS10_hmap_4m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid IVS10_4m"
)



#for S305I mutation 
#for S305I_2m 

S305I_NPX_list <- list (S305I_2m_NPX, S305I_4m_NPX)
combined_NPX_S305I <- reduce(S305I_NPX_list, full_join, by = "PanelAssay")


S305I_2m_1 <- combined_NPX_S305I$mut1_S305I_2m - combined_NPX_S305I$wt1_S305I_2m
S305I_2m_2 <- combined_NPX_S305I$mut2_S305I_2m - combined_NPX_S305I$wt2_S305I_2m
S305I_4m_1 <- combined_NPX_S305I$mut1_S305I_4m - combined_NPX_S305I$wt1_S305I_4m
S305I_4m_2 <- combined_NPX_S305I$mut2_S305I_4m - combined_NPX_S305I$wt1_S305I_4m

diff_S305I <- cbind(S305I_2m_1, S305I_2m_2, S305I_4m_1, S305I_4m_2)
rownames(diff_S305I) <- combined_NPX$PanelAssay

S305I_hmap <- cbind(diff_S305I, t_test_protein_S305I_2m, by = "PanelAssay")
S305I_hmap <- S305I_hmap[, c("S305I_2m_1", "S305I_2m_2", "S305I_4m_1", "S305I_4m_2", "PanelAssay", "pvalue", "pvalue")]

#S305I_hmap <- S305I_hmap[, c("S305I_2m_1", "S305I_4m_1", "S305I_2m_2", "S305I_4m_2", "PanelAssay", "pvalue", "pvalue")]



S305I_hmap <- S305I_hmap[, -which(names(S305I_hmap) == "PanelAssay")] 


S305I_hmap <- S305I_hmap[S305I_hmap [, "pvalue"] < 0.05, ]
S305I_hmap <- S305I_hmap[, -which(names(S305I_hmap) == "pvalue")] 
S305I_hmap_2m <- S305I_hmap[, -which(names(S305I_hmap) == "pvalue.1")] 


S305I_hmap_2m <- as.matrix(S305I_hmap_2m)

pheatmap(
  S305I_hmap_2m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid S305I_2m"
)



#for S305I_4m 
S305I_NPX_list <- list (S305I_2m_NPX, S305I_4m_NPX)
combined_NPX_S305I <- reduce(S305I_NPX_list, full_join, by = "PanelAssay")


S305I_2m_1 <- combined_NPX_S305I$mut1_S305I_2m - combined_NPX_S305I$wt1_S305I_2m
S305I_2m_2 <- combined_NPX_S305I$mut2_S305I_2m - combined_NPX_S305I$wt2_S305I_2m
S305I_4m_1 <- combined_NPX_S305I$mut1_S305I_4m - combined_NPX_S305I$wt1_S305I_4m
S305I_4m_2 <- combined_NPX_S305I$mut2_S305I_4m - combined_NPX_S305I$wt1_S305I_4m

diff_S305I <- cbind(S305I_2m_1, S305I_2m_2, S305I_4m_1, S305I_4m_2)
rownames(diff_S305I) <- combined_NPX$PanelAssay

S305I_hmap <- cbind(diff_S305I, t_test_protein_S305I_4m, by = "PanelAssay")
#S305I_hmap <- S305I_hmap[, c("S305I_2m_1", "S305I_2m_2", "S305I_4m_1", "S305I_4m_2", "PanelAssay", "pvalue", "pvalue")]

S305I_hmap <- S305I_hmap[, c("S305I_2m_1", "S305I_4m_1", "S305I_2m_2", "S305I_4m_2", "PanelAssay", "pvalue", "pvalue")]



S305I_hmap <- S305I_hmap[, -which(names(S305I_hmap) == "PanelAssay")] 


S305I_hmap <- S305I_hmap[S305I_hmap [, "pvalue"] < 0.05, ]
S305I_hmap <- S305I_hmap[, -which(names(S305I_hmap) == "pvalue")] 
S305I_hmap_4m <- S305I_hmap[, -which(names(S305I_hmap) == "pvalue.1")] 


S305I_hmap_4m <- as.matrix(S305I_hmap_4m)

pheatmap(
  S305I_hmap_4m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid S305I_4m"
)



#for R406W_2m 

R406W_NPX_list <- list (R406W_2m_NPX, R406W_4m_NPX)
combined_NPX_R406W <- reduce(R406W_NPX_list, full_join, by = "PanelAssay")


R406W_2m_1 <- combined_NPX_R406W$mut1_R406W_2m - combined_NPX_R406W$wt1_R406W_2m
R406W_2m_2 <- combined_NPX_R406W$mut2_R406W_2m - combined_NPX_R406W$wt2_R406W_2m
R406W_4m_1 <- combined_NPX_R406W$mut1_R406W_4m - combined_NPX_R406W$wt1_R406W_4m
R406W_4m_2 <- combined_NPX_R406W$mut2_R406W_4m - combined_NPX_R406W$wt1_R406W_4m

diff_R406W <- cbind(R406W_2m_1, R406W_2m_2, R406W_4m_1, R406W_4m_2)
rownames(diff_R406W) <- combined_NPX$PanelAssay

R406W_hmap <- cbind(diff_R406W, t_test_protein_R406W_2m, by = "PanelAssay")
#R406W_hmap <- R406W_hmap[, c("R406W_2m_1", "R406W_2m_2", "R406W_4m_1", "R406W_4m_2", "PanelAssay", "pvalue", "pvalue")]

R406W_hmap <- R406W_hmap[, c("R406W_2m_1", "R406W_4m_1", "R406W_2m_2", "R406W_4m_2", "PanelAssay", "pvalue", "pvalue")]



R406W_hmap <- R406W_hmap[, -which(names(R406W_hmap) == "PanelAssay")] 


R406W_hmap <- R406W_hmap[R406W_hmap [, "pvalue"] < 0.05, ]
R406W_hmap <- R406W_hmap[, -which(names(R406W_hmap) == "pvalue")] 
R406W_hmap_2m <- R406W_hmap[, -which(names(R406W_hmap) == "pvalue.1")] 


R406W_hmap_2m <- as.matrix(R406W_hmap_2m)

pheatmap(
  R406W_hmap_2m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid R406W_2m"
)


#R406W_4m 

R406W_NPX_list <- list (R406W_2m_NPX, R406W_4m_NPX)
combined_NPX_R406W <- reduce(R406W_NPX_list, full_join, by = "PanelAssay")


R406W_2m_1 <- combined_NPX_R406W$mut1_R406W_2m - combined_NPX_R406W$wt1_R406W_2m
R406W_2m_2 <- combined_NPX_R406W$mut2_R406W_2m - combined_NPX_R406W$wt2_R406W_2m
R406W_4m_1 <- combined_NPX_R406W$mut1_R406W_4m - combined_NPX_R406W$wt1_R406W_4m
R406W_4m_2 <- combined_NPX_R406W$mut2_R406W_4m - combined_NPX_R406W$wt1_R406W_4m

diff_R406W <- cbind(R406W_2m_1, R406W_2m_2, R406W_4m_1, R406W_4m_2)
rownames(diff_R406W) <- combined_NPX$PanelAssay

R406W_hmap <- cbind(diff_R406W, t_test_protein_R406W_4m, by = "PanelAssay")
R406W_hmap <- R406W_hmap[, c("R406W_2m_1", "R406W_2m_2", "R406W_4m_1", "R406W_4m_2", "PanelAssay", "pvalue", "pvalue")]

#R406W_hmap <- R406W_hmap[, c("R406W_2m_1", "R406W_4m_1", "R406W_2m_2", "R406W_4m_2", "PanelAssay", "pvalue", "pvalue")]



R406W_hmap <- R406W_hmap[, -which(names(R406W_hmap) == "PanelAssay")] 


R406W_hmap <- R406W_hmap[R406W_hmap [, "pvalue"] < 0.05, ]
R406W_hmap <- R406W_hmap[, -which(names(R406W_hmap) == "pvalue")] 
R406W_hmap_4m <- R406W_hmap[, -which(names(R406W_hmap) == "pvalue.1")] 


R406W_hmap_4m <- as.matrix(R406W_hmap_4m)

pheatmap(
  R406W_hmap_4m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid R406W_4m"
)



#for V337M_2m 
V337M_NPX_list <- list (V337M_2m_NPX, V337M_4m_NPX, V337M_6m_NPX)
combined_NPX_V337M <- reduce(V337M_NPX_list, full_join, by = "PanelAssay")


V337M_2m_1 <- combined_NPX_V337M$mut1_V337M_2m - combined_NPX_V337M$wt1_V337M_2m
V337M_2m_2 <- combined_NPX_V337M$mut2_V337M_2m - combined_NPX_V337M$wt2_V337M_2m
V337M_2m_3 <- combined_NPX_V337M$mut3_V337M_2m - combined_NPX_V337M$wt3_V337M_2m
V337M_4m_2 <- combined_NPX_V337M$mut2_V337M_4m - combined_NPX_V337M$wt2_V337M_4m
V337M_4m_3 <- combined_NPX_V337M$mut3_V337M_4m - combined_NPX_V337M$wt3_V337M_4m
V337M_6m_1 <- combined_NPX_V337M$mut1_V337M_6m - combined_NPX_V337M$wt1_V337M_6m
V337M_6m_2 <- combined_NPX_V337M$mut2_V337M_6m - combined_NPX_V337M$wt2_V337M_6m


diff_V337M <- cbind(V337M_2m_1, V337M_2m_2, V337M_2m_3, V337M_4m_2, V337M_4m_3, V337M_6m_1, V337M_6m_2)
rownames(diff_V337M) <- combined_NPX$PanelAssay

V337M_hmap <- cbind(diff_V337M, t_test_protein_V337M_2m, by = "PanelAssay")

#V337M_hmap <- V337M_hmap[, c("V337M_2m_1", "V337M_6m_1", 
                             "V337M_2m_2", "V337M_4m_2","V337M_6m_2",
                             "V337M_2m_3", "V337M_4m_3", 
                              "PanelAssay", "pvalue", "pvalue")]


#V337M_hmap <- V337M_hmap[, c("V337M_2m_1", "V337M_2m_2", "V337M_2m_3", 
                              "V337M_4m_2","V337M_4m_3", 
                             "V337M_6m_1", "V337M_6m_2",
                             "PanelAssay", "pvalue", "pvalue")]





V337M_hmap <- V337M_hmap[, -which(names(V337M_hmap) == "PanelAssay")] 

#for IVS10_2m 

V337M_hmap <- V337M_hmap[V337M_hmap [, "pvalue"] < 0.05, ]
V337M_hmap <- V337M_hmap[, -which(names(V337M_hmap) == "pvalue")] 
V337M_hmap_2m <- V337M_hmap[, -which(names(V337M_hmap) == "pvalue.1")] 


V337M_hmap_2m <- as.matrix(V337M_hmap_2m)

pheatmap(
  V337M_hmap_2m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid V337M_2m"
)



#for V337M_4m 
V337M_NPX_list <- list (V337M_2m_NPX, V337M_4m_NPX, V337M_6m_NPX)
combined_NPX_V337M <- reduce(V337M_NPX_list, full_join, by = "PanelAssay")


V337M_2m_1 <- combined_NPX_V337M$mut1_V337M_2m - combined_NPX_V337M$wt1_V337M_2m
V337M_2m_2 <- combined_NPX_V337M$mut2_V337M_2m - combined_NPX_V337M$wt2_V337M_2m
V337M_2m_3 <- combined_NPX_V337M$mut3_V337M_2m - combined_NPX_V337M$wt3_V337M_2m
V337M_4m_2 <- combined_NPX_V337M$mut2_V337M_4m - combined_NPX_V337M$wt2_V337M_4m
V337M_4m_3 <- combined_NPX_V337M$mut3_V337M_4m - combined_NPX_V337M$wt3_V337M_4m
V337M_6m_1 <- combined_NPX_V337M$mut1_V337M_6m - combined_NPX_V337M$wt1_V337M_6m
V337M_6m_2 <- combined_NPX_V337M$mut2_V337M_6m - combined_NPX_V337M$wt2_V337M_6m


diff_V337M <- cbind(V337M_2m_1, V337M_2m_2, V337M_2m_3, V337M_4m_2, V337M_4m_3, V337M_6m_1, V337M_6m_2)
rownames(diff_V337M) <- combined_NPX$PanelAssay

V337M_hmap <- cbind(diff_V337M, t_test_protein_V337M_4m, by = "PanelAssay")

#V337M_hmap <- V337M_hmap[, c("V337M_2m_1", "V337M_6m_1", 
"V337M_2m_2", "V337M_4m_2","V337M_6m_2",
"V337M_2m_3", "V337M_4m_3", 
"PanelAssay", "pvalue", "pvalue")]


V337M_hmap <- V337M_hmap[, c("V337M_2m_1", "V337M_2m_2", "V337M_2m_3", 
"V337M_4m_2","V337M_4m_3", 
"V337M_6m_1", "V337M_6m_2",
"PanelAssay", "pvalue", "pvalue")]






V337M_hmap <- V337M_hmap[, -which(names(V337M_hmap) == "PanelAssay")] 

#for IVS10_2m 

V337M_hmap <- V337M_hmap[V337M_hmap [, "pvalue"] < 0.05, ]
V337M_hmap <- V337M_hmap[, -which(names(V337M_hmap) == "pvalue")] 
V337M_hmap_4m <- V337M_hmap[, -which(names(V337M_hmap) == "pvalue.1")] 


V337M_hmap_4m <- as.matrix(V337M_hmap_4m)

pheatmap(
  V337M_hmap_4m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid V337M_4m"
)


#for V337M_6m 
V337M_NPX_list <- list (V337M_2m_NPX, V337M_4m_NPX, V337M_6m_NPX)
combined_NPX_V337M <- reduce(V337M_NPX_list, full_join, by = "PanelAssay")


V337M_2m_1 <- combined_NPX_V337M$mut1_V337M_2m - combined_NPX_V337M$wt1_V337M_2m
V337M_2m_2 <- combined_NPX_V337M$mut2_V337M_2m - combined_NPX_V337M$wt2_V337M_2m
V337M_2m_3 <- combined_NPX_V337M$mut3_V337M_2m - combined_NPX_V337M$wt3_V337M_2m
V337M_4m_2 <- combined_NPX_V337M$mut2_V337M_4m - combined_NPX_V337M$wt2_V337M_4m
V337M_4m_3 <- combined_NPX_V337M$mut3_V337M_4m - combined_NPX_V337M$wt3_V337M_4m
V337M_6m_1 <- combined_NPX_V337M$mut1_V337M_6m - combined_NPX_V337M$wt1_V337M_6m
V337M_6m_2 <- combined_NPX_V337M$mut2_V337M_6m - combined_NPX_V337M$wt2_V337M_6m


diff_V337M <- cbind(V337M_2m_1, V337M_2m_2, V337M_2m_3, V337M_4m_2, V337M_4m_3, V337M_6m_1, V337M_6m_2)
rownames(diff_V337M) <- combined_NPX$PanelAssay

V337M_hmap <- cbind(diff_V337M, t_test_protein_V337M_6m, by = "PanelAssay")

V337M_hmap <- V337M_hmap[, c("V337M_2m_1", "V337M_6m_1", 
"V337M_2m_2", "V337M_4m_2","V337M_6m_2",
"V337M_2m_3", "V337M_4m_3", 
"PanelAssay", "pvalue", "pvalue")]


#V337M_hmap <- V337M_hmap[, c("V337M_2m_1", "V337M_2m_2", "V337M_2m_3", 
                             "V337M_4m_2","V337M_4m_3", 
                             "V337M_6m_1", "V337M_6m_2",
                             "PanelAssay", "pvalue", "pvalue")]






V337M_hmap <- V337M_hmap[, -which(names(V337M_hmap) == "PanelAssay")] 



V337M_hmap <- V337M_hmap[V337M_hmap [, "pvalue"] < 0.05, ]
V337M_hmap <- V337M_hmap[, -which(names(V337M_hmap) == "pvalue")] 
V337M_hmap_6m <- V337M_hmap[, -which(names(V337M_hmap) == "pvalue.1")] 


V337M_hmap_6m <- as.matrix(V337M_hmap_6m)

pheatmap(
  V337M_hmap_6m,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colors from blue to red
  scale = "row", # Scale each row
  clustering_distance_rows = "euclidean", # Distance metric for rows
  clustering_distance_cols = "euclidean", # Distance metric for columns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = FALSE, # Cluster columns
  show_rownames = TRUE, # Show row names
  show_colnames = TRUE, # Show column names
  fontsize_row = 6, # Font size for row names
  fontsize_col = 6, # Font size for column names
  angle_col = 45, # Angle for column names
  main = "Heatmap of Organoid V337M_6m"
)



#install.packages(c("shiny", "ggplot2", "dplyr"))
library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  titlePanel("Interactive Volcano Plot for IVS10_2m"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("pvalueThreshold", "P-value Threshold:", min = 0, max = 0.1, value = 0.05, step = 0.01),
      sliderInput("foldChangeThreshold", "Fold Change Threshold:", min = 0, max = 2, value = 0.5, step = 0.1),
      actionButton("updatePlot", "Update Plot")
    ),
    mainPanel(
      plotOutput("volcanoPlot")
    )
  )
)

erver <- function(input, output, session) {
  # Load data or prepare your dataset
  # Let's assume t_test_protein_IVS10_2m is already loaded
  
  # Reactive expression for filtering based on user input
  reactiveData <- reactive({
    df <- t_test_protein_IVS10_2m
    df$diffexpressed <- ifelse(df$mean.diff > input$foldChangeThreshold & df$pvalue < input$pvalueThreshold, "UP",
                               ifelse(df$mean.diff < -input$foldChangeThreshold & df$pvalue < input$pvalueThreshold, "DOWN", "Not-significant"))
    df
  })
  
  # Plot output
  output$volcanoPlot <- renderPlot({
    df <- reactiveData()
    ggplot(data = df, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed, label = delabel)) +
      geom_point(size = 3, alpha = 0.6) +
      geom_text_repel(max.overlaps = 15, show.legend = FALSE) +
      scale_color_manual(name = "Protein expression_IVS10_2m",
                         labels = c("Upregulated", "Downregulated", "Not significant"), 
                         values = c("red", "blue", "grey")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(x = expression(Log[2]~FC), y = expression(-Log[10]~p-value)) +
      theme_bw() +
      theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14))
  })
  
  # Update plot when button is clicked
  observeEvent(input$updatePlot, {
    output$volcanoPlot <- renderPlot({
      df <- reactiveData()
      # Same ggplot code as above
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)





















































