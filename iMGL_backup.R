





#install.packages("readxl")
library(readxl)
#install.packages("writexl")
library(writexl)
#install.packages("xlsx")

## TO DEFINE WHERE THE OUTPUT WILL BE ##
source("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/olinkdata_objectTools_v1.1.R")

## FOR CLEANING TOOL ##
source("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/olinkdata_cleaningTools_v1.1.R")

outputDir = "M:/NL/6_Scientific Results/2_Students/Noha Salem/3 Output/20240212"

project_name = 'MAPT'

## WE NEED TO CREATE THE OLINK OBJECT ##

## TO SUBSET THE MICROGLIA DATA ONLY FROM THE FILES AND SAVE FOR CREATING THE OBJECT ##
sample_info <- read_excel("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/sample_info.xlsx")
microglia <-sample_info[sample_info$group == "iMGL" , ]

olink_npx_data <- read_excel("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/olink_data_pilot_Mapt_contr_removed.xlsx")
olink_npx_data <- olink_npx_data[which(olink_npx_data$SampleID %in% microglia$SampleID), ]

excel_file <-("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/microglia.xlsx")
write_xlsx(x=microglia, path=excel_file)
excel_file <-("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/olink_npx_data.xlsx")
write_xlsx(x=olink_npx_data, path=excel_file)

## DEFINE THE FILES ##
olinkDataFiles <- c("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/olink_npx_data.xlsx")
olinkDataSetNames <- c("data_MAPT")
OlinkData <- createOlinkDataObject(olinkDataFiles=olinkDataFiles,olinkDataSetNames=olinkDataSetNames,olinkLayoutFiles="")

# I GOT EROR HERE (FEB 8TH) #
#OlinkData@samples[['projectName']] = "iMGL"


## WE NEED TO ANNOTATE IN ORDER TO MERGE THE SAMPLE INFOR FILE WITH THE OLINKDATA FILE ##
microglia <- c("M:/NL/6_Scientific Results/2_Students/Noha Salem/1 Data/for R path/microglia.xlsx")
SampleIDcolumnsName <- c("SampleID")
OlinkData = annotateOlinkData(OlinkData,annotationFiles=microglia,OlinkIDcolumnsName=SampleIDcolumnsName)


### 

OlinkData@design[['project']] = OlinkData@samples[['project']]
OlinkData@design[['design']] = OlinkData@samples[['Mutation']]

## VISUALIZE THE DESIGN THAT WILL BE USED FOR ALL CLEANING STEPS ##
#mut_color <- "red"
#con_color <- "blue"
#sample_colors <- ifelse(OlinkDat$design == "Mutation", mut_color, con_color)



#WHAT DOES IT MEAN BY THE ERROR HERE?##

#OlinkData@color.vector <- c(OlinkData@color.vector, setNames("#F2AFB4", "WT"))
#OlinkData@color.vector <- c(OlinkData@color.vector, setNames("#9DC893", "Mutation"))
# CUZ WE HAD AN ERROR. WE HAD TO DEFINE THE PROJECT (iMGL) WITH ANOTHER COLOR#
#OlinkData@color.vector <- c(OlinkData@color.vector, setNames("#36454F", "iMGL"))
OlinkData@color.vector <- c(OlinkData@color.vector, setNames("#36454F", "iMGL"))


project_name = 'MAPT'
OlinkData@samples[['projectName']] = "iMGL"


OlinkData = set.design(OlinkData,
                       projectName= "iMGL",
                       samples.columnName.group = "Mutation", 
                       samples.columnName.time ="time",
                       samples.columnName.time.label="timepointLabel",
                       samples.columnName.sampleID="SampleID",
                       samples.columnName.subjectID="Donor")

designLayoutPlot = plot.design.layout(OlinkData,prct=F, colorbyProject=T)
plot(designLayoutPlot)

## select proteins above LOD. 
install.packages("GGally")
PCA_plot_raw = plot.PCA(OlinkData,slot="NPX.raw",nComponents=3,samples.col.name.color="project")
PCA_plot_raw



#1. we need to detect the NPX below the LOD (measurement bias), sturctural measurments below the LOD will be replaced by NA 
install.packages("diptest")
OlinkData = detectLowNPX(OlinkData,LOD.factor=1,numberCores=1)

#2. in order to effect of the cutoff score (for LOD), create a cutoff plot, the default value for the cutoff is 100% beavuse we only have 4 samples 
#install.packages("fitdistrplus")
cutoffPlot_LOD = makeCutoffPlot(OlinkData, functionName="detectLowNPX",designCol ="project")
cutoffPlot_LOD

#3. Exclude the proteins that do not have an NPX values above the LOD at least 85% 
#due to error (connection), we changed the number cores to be 1 instead of 12
OlinkData = excludeLowNPX(OlinkData, prctPass=100, FilterPerProject=T,numberCores=1)

#explore the cleaning step via PCA 
PCA_plot_LOD = plot.PCA(OlinkData,slot="NPX",nComponents=3,samples.col.name.color="project")
PCA_plot_LOD

# to calculate the CVs

cvcal_mut_1 <- OlinkData@NPX[ , c('1A', '71G', '75C')]
cvcal_mut_2 <- OlinkData@NPX[ , c('3C', '72H', '73A')] 
cvcal_con_1 <- OlinkData@NPX[ , c('2B', '69E', '76D')]
cvcal_con_2 <- OlinkData@NPX[ , c('4D', '70F', '74B')]


#to remove the NA
cvcal_mut_1 <- na.omit(cvcal_mut_1)
cvcal_mut_2 <- na.omit(cvcal_mut_2)
cvcal_con_1 <- na.omit(cvcal_con_1)
cvcal_con_2 <- na.omit(cvcal_con_2)

## calculate the CVs for each sample separetly 
#mut_1
cvcal_mut_1 <- as.data.frame(cvcal_mut_1)                    
cvcal_mut_1$mut_1 <- rowMeans(cvcal_mut_1) 
cvcal_mut_1$Median = apply(cvcal_mut_1, 1, median, na.rm=T)
cvcal_mut_1$SD = apply(cvcal_mut_1, 1, sd, na.rm= T)
cvcal_mut_1$CV_step1 = log(2)*cvcal_mut_1$SD
cvcal_mut_1$CV_mut_1 = sqrt(exp((cvcal_mut_1$CV_step1)^2) - 1)*100
cvcal_mut_1$panelassay <- rownames (cvcal_mut_1)

#con_1
cvcal_con_1 <- as.data.frame(cvcal_con_1)                    
cvcal_con_1$con_1 <- rowMeans(cvcal_con_1) 
cvcal_con_1$Median = apply(cvcal_con_1, 1, median, na.rm=T)
cvcal_con_1$SD = apply(cvcal_con_1, 1, sd, na.rm= T)
cvcal_con_1$CV_step1 = log(2)*cvcal_con_1$SD
cvcal_con_1$CV_con_1 = sqrt(exp((cvcal_con_1$CV_step1)^2) - 1)*100
cvcal_con_1$panelassay <- rownames (cvcal_con_1)

#mut_2
cvcal_mut_2 <- as.data.frame(cvcal_mut_2)         
cvcal_mut_2$mut_2 <- rowMeans(cvcal_mut_2) 
cvcal_mut_2$Median = apply(cvcal_mut_2, 1, median, na.rm=T)
cvcal_mut_2$SD = apply(cvcal_mut_2, 1, sd, na.rm= T)
cvcal_mut_2$CV_step1 = log(2)*cvcal_mut_2$SD
cvcal_mut_2$CV_mut_2 = sqrt(exp((cvcal_mut_2$CV_step1)^2) - 1)*100
cvcal_mut_2$panelassay <- rownames (cvcal_mut_2)

#con_2
cvcal_con_2 <- as.data.frame(cvcal_con_2)                    
cvcal_con_2$con_2 <- rowMeans(cvcal_con_2) 
cvcal_con_2$Median = apply(cvcal_con_2, 1, median, na.rm=T)
cvcal_con_2$SD = apply(cvcal_con_2, 1, sd, na.rm= T)
cvcal_con_2$CV_step1 = log(2)*cvcal_con_2$SD
cvcal_con_2$CV_con_2 = sqrt(exp((cvcal_con_2$CV_step1)^2) - 1)*100
cvcal_con_2$panelassay <- rownames (cvcal_con_2)



cvs_con <- cbind(cvcal_con_1$CV_con_1,cvcal_con_2$CV_con_2)
cvs_con <- as.data.frame (cbind(cvcal_con_1$CV_con_1,cvcal_con_2$CV_con_2))
cvs_con$cv_mean_con <- rowMeans(cvs_con)
cvs_con <- cbind (cvs_con, cvcal_con_1$panelassay)
colnames(cvs_con) [colnames(cvs_con) == "cvcal_con_1$panelassay"] <- "Proteins"


cvs_mut <- cbind(cvcal_mut_1$CV_mut_1, cvcal_mut_2$CV_mut_2)
cvs_mut <- as.data.frame( cbind(cvcal_mut_1$CV_mut_1, cvcal_mut_2$CV_mut_2))
cvs_mut$cv_mean_mut <- rowMeans(cvs_mut)
cvs_mut <- cbind (cvs_mut, cvcal_mut_1$panelassay)
colnames(cvs_mut) [colnames(cvs_mut) == "cvcal_mut_1$panelassay"] <- "Proteins"


# TO DEFINE THE CUTT OFF SCORES FOR THE CVS TO BE LESS THAN 20% #

selected_proteins_con  <- cvs_con[cvs_con$cv_mean_con <20, ]
removed_proteins_con <- setdiff(cvs_con, selected_proteins_con)

selected_proteins_mut  <- cvs_mut[cvs_mut$cv_mean_mut <20, ]
removed_proteins_mut <- setdiff(cvs_mut, selected_proteins_mut)




#to combine the cv and mean values for both mut and control 
cvs_both <- c(selected_proteins_mut$Proteins, selected_proteins_con$Proteins)

#to create a new data frame that contains the means and cvs for all samples and name the rows and columns 
mapt_clean <- as.data.frame (cbind( cvcal_mut_1$mut_1,
                                    cvcal_mut_2$mut_2,
                                    cvcal_con_1$con_1,
                                    cvcal_con_2$con_2,
                                    cvcal_mut_1$CV_mut_1,
                                    cvcal_mut_2$CV_mut_2,
                                    cvcal_con_1$CV_con_1,
                                    cvcal_con_2$CV_con_2,
                                    cvcal_mut_2$cv_mean_mut,
                                    cvcal_mut_1$cv_mean_con))


colnames(mapt_clean) <- c('mut_1','mut_2',
                          'con_1', 'con_2',
                          'cv_mut_1','cv_mut_2',
                          'cv_con_1', 'cv_con_2')

rownames(mapt_clean) <- cvcal_mut_1$panelassay 
mapt_clean$protein <- cvcal_mut_1$panelassay 



#create histogram for the cvs (IT KEPT RUNNING FOR 20 MIN WIHOUT A RESULT)
library (ggplot2)

#ggplot(mapt_clean, aes(x = cv_con_1, fill = "con_1")) +
geom_histogram(binwidth = 1, position = "dodge") +
  geom_histogram(data = mapt_clean, aes(x = cv_con_2, fill = "con_2"), binwidth = 1, position = "dodge") +
  geom_histogram(data = mapt_clean, aes(x = cv_mut_1, fill = "mut_1"), binwidth = 1, position = "dodge") +
  geom_histogram(data = mapt_clean, aes(x = cv_mut_2, fill = "mut_2"), binwidth = 1, position = "dodge") +
  labs(title = "Histogram of CVs for each protein",
       x = "CVs", y = "Frequency") +
  scale_fill_manual(values = c("con_1" = "blue", "con_2" = "green", 
                               "mut_1" = "red", "mut_2" = "orange"),
                    guide = "none") +
  facet_wrap(~ protein, scales = "free")



#mapt_long <- pivot_longer(mapt_clean, cols = starts_with("cv"), names_to = "Sample", values_to = "CV")
#ggplot(mapt_long, aes(x = CV, fill = Sample)) +
geom_histogram(binwidth = 1, position = "dodge") +
  labs(title = "Histogram of CVs for each protein",
       x = "CVs", y = "Protein") +
  facet_wrap(~ protein, scales = "free") +
  scale_fill_manual(values = c("con_1" = "blue", "con_2" = "green", 
                               "mut_1" = "red", "mut_2" = "orange"))



#TO ADD THE UNIPROT PRIMARY ACCESSION NUMBER FOR THE PROTEINS FOR THE T-TEST 
#1. assigning the proteins in the olink to the uniprot variables (to assign each protein to its counter part uniprot accession number) 
uniprot <- OlinkData@proteins #to create a dataframe called uniprot
uniprot$protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #"_" to separate the panel from the assay in the newly created uniprot dataframe 

#2.the new dataframe (uniprot) contains all values incluing those that are below the LOD and those that are NA (missingness)
missingness <- OlinkData@prct.NPX.below.LOD #assigning the prct.NPX.below.LOD to the missingness
colnames(missingness) <- c('prct.below.LOD.mut', 'prct.below.LOD.con') #separate the con from the mut
all_data <- merge(uniprot, mapt_clean, by="protein", all.y=T) #merge the two dataframes uniprot and the mapt_clean by the common column proteins


#I don't know what this code is doing and it gives me an error as LOD isn't in the all_data (actually the length of all_data and OlinkData is not similar) 
#all_data$LOD <- OlinkData@LOD[,1]
#all_data$LOD <- OlinkData@LOD
#all_data <- merge(all_data, OlinkData@LOD, by = "sampleID")

#SELECT THE PROTEINS THAT ARE ONLY FOUND IN BOTH CON AND MUT (RIGHT ONES)
mapt_clean <- mapt_clean [rownames(mapt_clean) %in% cvs_both,]
pca_data <- t(mapt_clean[ ,1:4])
# to save the mapt_clean df that I got
#write.csv(mapt_clean, "mapt_clean.csv")

# HERE I USED LISA'S CODE BUT I WILL KEEP MINE ALSO, BOTH OF THEM WORKED #
#MAPT_clean <- as.data.frame (cbind (cvcal_con_1$Mean, cvcal_con_2$Mean, cvcal_mut_1$Mean, cvcal_mut_2$Mean))
#MAPT_clean <- cbind(MAPT_clean, cvcal_con_1$panelassay)
#colnames (MAPT_clean) [colnames (MAPT_clean) == "cvcal_con_1$panelassay"] <- "Proteins"
#colnames (MAPT_clean) [colnames (MAPT_clean) == "V1"] <- "con_1_mean"
#colnames (MAPT_clean) [colnames (MAPT_clean) == "V1"] <- "con_1_mean"
#colnames (MAPT_clean) [colnames (MAPT_clean) == "V2"] <- "con_2_mean"
#colnames (MAPT_clean) [colnames (MAPT_clean) == "V3"] <- "mut_1_mean"
#colnames (MAPT_clean) [colnames (MAPT_clean) == "V4"] <- "mut_2_mean"


library ("FactoMineR") #for multivariate exploratory data analysis (MV-EDA), provides set of functions for PCA and clustering
library ("factoextra") #for extracting and visualizing the results of mulivariate data 

PCA (pca_data, scale.unit = TRUE, ncp = 5, graph = TRUE) #scale unit is important to standardize all data (z-scores normalization)thus ensures that all data will contribute to the PCA making PCA more robust (decreases results skewness)
res.pca <- PCA(pca_data, graph = FALSE) #the results of PCA in a list

eig.val <- get_eigenvalue(res.pca) #eigenvalues represent the total amount of variance that can be explained by PCA. 
eig.val
# results explanation: 
#1. eigenvalue each one represents the amount of variance in the original data that is captured by each component 
#2. variance.percent represents the total variance of the dataset that explained by each principle component (it is calculated by: (eigenvalue for each PC/ total eigenvalues)* 100)
#3. cumulative.variance.percent represents 63% for the 1st PC, 86% for both PC1 and 2, and 100% for the n PC that captures all variablity in the dataset.

fviz_eig (res.pca, addlabels = TRUE, ylim = c (0,50)) #to plot the eigenvalues/variance against the no. dimensions

var <- get_pca_var(res.pca)
var
#$coord (coordinates for the variables):coordinates are values assigned to variables (proteins) within a dataset. These coordinates help to locate and understand the behavior of the proteins within the system
#$cor (correleation between varaiables): correlation matrix or correlation coefficients between the variables and dimensions in a dataset.
#$cos2 (squared cos for the variables): the squared cos of th angel between the variable and the PC, measures the quality of the representation of each varaible by the PC (provides insight into how much each original variable contributes to the principal components),
# high cos2 means the variable is well-represented by the principal component and vice versa.


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
grp <- as.factor (c ('con_1', 'con_2', 'mut_1', 'mut_2'))
fviz_pca_ind(res.pca, fill.ind = grp, pointshape = 21, palette = c ("#F2AFB4", "#9DC893", "#36454F", "#005b96"), repel = TRUE)



#T-TEST 
library(matrixTests)
library(plyr)
library(ggrepel)



mapt_clean <- read.csv ("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/mapt_clean.csv", row.names = 1)
#create NPX df that contains the means for the protein in each sample
NPX <-mapt_clean [, (1:4)]
#define the con cloumns and mut columns separtely as it is paired t-tes 
mut <- c ('mut_1', 'mut_2') #this arrangement is critical because we need to define the mean.diff (mean.diff= mean_mut - mean_con)
con <- c ('con_1', 'con_2')

#using paied t-test because the samples are taken from the same donor (i.g: con_1 and mut_1 from Donor 1)

t_test_protein <- row_t_paired(NPX[,colnames(NPX) %in% mut], 
                               NPX[,colnames(NPX) %in% con])

t_test_protein$Protein <- rownames (t_test_protein)
p_value <- as.vector(t_test_protein$pvalue) #selecting the p-value in a separate vector 
adjusted_p_value <- p.adjust(p_value, method = "fdr") #to remove type I error (alpha/ false positive)
t_test_protein$adjusted_p_value <- adjusted_p_value


uniprot <- OlinkData@proteins
uniprot$Protein <- paste0(uniprot$Panel, "_", uniprot$Assay) #to assign the protein to the uniprot accession number
missingness <- OlinkData@prct.NPX.below.LOD #to include the missingness
colnames(missingness) <- c('prct.below.LOD.mut', 'prct.below.LOD.wt')
t_test_uniprot <- merge(uniprot, t_test_protein, by = "Protein", all.y = T)


#write.csv(t_test_uniprot, "all_information of iMGL.csv")
#x= mut, y= con
# check this with Lisa if the mean.diff is the same as log2fold change
#pseudocount <- 1 #use constant 
#log2_FC <- log2((t_test_uniprot$mean.y + pseudocount) / (t_test_uniprot$mean.x + pseudocount))
#t_test_protein$log2_FC <- log2_FC

#to create volcano plots 
#1. to create a column in the t_test_protein df called diffexpressed that identifies each protein either into (not, up or down) based on the p_value
t_test_protein$diffexpressed <- "Not-significant"  
t_test_protein$diffexpressed[t_test_protein$mean.diff>0.0 & t_test_protein$pvalue<0.05] <- "UP"
t_test_protein$diffexpressed[t_test_protein$mean.diff<0.0 & t_test_protein$pvalue<0.05] <- "DOWN"
t_test_protein$protein2 <- gsub(".*_", "", t_test_protein$Protein)
t_test_protein$delabel <- NA
t_test_protein$delabel[t_test_protein$diffexpressed != "Not-significant"] <- t_test_protein$protein2[t_test_protein$diffexpressed != "Not-significant"]
t_test_protein$diffexpressed <- factor(t_test_protein$diffexpressed, levels = c("UP", "DOWN", "Not-significant")) 




# -ve mean.diff means the mean of con is higher and +ve mean.diff means the mean mut is higher 


#to create enhnaced volcano plot: 
Significance_level <- -log10(0.05)  # Adjust p-value threshold as necessary for the horizontal line
logFC_threshold <- 0                # Adjust log fold change threshold for the vertical line and labeling

# Create the plot
volcano_plot <- ggplot(data = t_test_protein, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  #geom_text_repel(max.overlaps = 15, show.legend = F)+ 
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "iMGL IVS10+16",
                     labels = c("Upr (73)", "Down (16)", "Not-sign (1271)"),
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



##########################



#2. to create the volcano plot 
ggplot(data = t_test_protein, aes(x=mean.diff, y = -log10(pvalue), color=diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression" , labels= c("Upregulated (n=73)", "Downregulated (n=16)", "Not significant (n=1271)"), values = c("#f03b20","blue","grey")) +
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



library(ggplot2)

# Create the plot
volcano_plot <- ggplot(data = t_test_protein, aes(x = mean.diff, y = -log10(pvalue), color = diffexpressed)) +
  geom_point(size = 3, alpha = 1) +  # Uniform size and no transparency
  scale_color_manual(name = "Protein expression",
                     labels = c("Up (73)", "Down (16)", "Non-sig (1271)"),
                     values = c("#880808", "#6495ED", "#818589")) +  # Red for up, blue for down, grey for not significant
  labs(x = "Log2 Fold Change", y = "-log10(p-value)") +
  theme_bw() +  # Use a base theme that starts with no panel background
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional for a clean look)
  )

# Print the plot
print(volcano_plot)





#to select the significant proteins only 

significant_pro_df <- subset(t_test_protein, pvalue < 0.05)

#3. to identify the top 25 proteins 
t_test_protein$top25 <- "NA"
t_test_protein$top25[t_test_protein$pvalue<0.03] <- "yes" 
t_test_protein$delabel2[t_test_protein$top25 != "NA"] <- t_test_protein$protein2[t_test_protein$top25 != "NA"]
t_test_protein <- t_test_protein[order(-t_test_protein$mean.diff),]
#put new cutoff score to select the top 25 proteins based upon the mean.diff and the p-value
sum(t_test_protein$mean.diff>0.4 & t_test_protein$pvalue<0.03) #All proteins with FC > 0.78 are in the top 10 
t_test_protein$fctop25 <- "NA"
t_test_protein$fctop25[t_test_protein$mean.diff>0.4 & t_test_protein$pvalue<0.03] <- "yes" 
t_test_protein$delabel3[t_test_protein$top25 == "yes" | t_test_protein$fctop25 =="yes"] <- t_test_protein$protein2[t_test_protein$top25 == "yes"]


ggplot(t_test_protein, aes(x=mean.diff, y = -log10(pvalue), color= diffexpressed, label=delabel))+
  geom_point(size=3, alpha=0.6)+
  geom_text_repel(max.overlaps = 15, show.legend = F)+  
  scale_color_manual(name = "Protein expression" , labels= c("Upregulated (n = 73)", "Downregulated (n = 16)", "Not significant (n = 1271)"), values = c("#f03b20","blue","grey")) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.0) +
  labs(x=expression(Log[2]~FC))+
  labs(y=expression(-Log[10]~q-value))+
  theme(axis.title = element_text(size = 16))+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=14))+
  guides()


#GO 
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



write.csv(t_test_uniprot, "GO_t_test_uniprot.csv")
GO_t_test_uniprot <- read.csv("M:/NL/6_Scientific Results/2_Students/Noha Salem/2 Scripts/GO_t_test_uniprot.csv")
GO_t_test_uniprot <- na.omit (GO_t_test_uniprot)
row_protein_list <- GO_t_test_uniprot$mean.diff #we need to use mean.diff 
protein_list <- na.omit(row_protein_list) #9proteins have been removed with final list 1328
names(protein_list) <- GO_t_test_uniprot$UNIPROT 
protein_list_unique <- GO_t_test_uniprot %>% distinct(UNIPROT, .keep_all = TRUE) #to remove the duplicates 
protein_list_unique <- setNames(GO_t_test_uniprot$mean.diff, as.character(GO_t_test_uniprot$UNIPROT))
keytypes(org.Hs.eg.db) #to check what kind of protein IDs are there 


protein_list_unique <- GO_t_test_uniprot %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list <- setNames(protein_list_unique$mean.diff, as.character(protein_list_unique$ENTREZID))
protein_list <- sort (protein_list, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list <- subset(protein_list_unique, pvalue < 0.05) # Exctract significant results (padj < 0.05)


#over representation analysis 
protein_list_unique <- GO_t_test_uniprot %>% distinct(ENTREZID, .keep_all = TRUE)
protein_list <- setNames(GO_t_test_uniprot$mean.diff, as.character(GO_t_test_uniprot$ENTREZID))
protein_list <- sort (protein_list, decreasing = TRUE)# this is the protein list ranked for the GSEA test (based on the differential gene/protein expression between the mut and con)
sig_protein_list <- subset(GO_t_test_uniprot, pvalue < 0.05) # Exctract significant results (padj < 0.05)
sig_protein_list <- sig_protein_list %>% distinct(UNIPROT, .keep_all = TRUE)
sig_protein_list_unique <- na.omit(sig_protein_list)


go_enrich_IVS10 <- enrichGO(gene = as.character(sig_protein_list_unique$ENTREZID),
                               universe = as.character(protein_list_unique$ENTREZID),
                               OrgDb = HS, 
                               keyType = 'ENTREZID',
                               readable = T,
                               ont = "BP", #for this one there was no gene can be maped for the PB and the MF 
                               pvalueCutoff = 0.9, #this is not significant  
                               qvalueCutoff = 0.9) 



dotplot(go_enrich_IVS10)

IVS10_go_enrich_BP <- go_enrich_IVS10@result

IVS10_go_enrich_BP <- IVS10_go_enrich_BP[IVS10_go_enrich_BP$pvalue < 0.05, ]
IVS10_go_enrich_BP <- IVS10_go_enrich_BP[order(IVS10_go_enrich_BP$pvalue), ]
level_order <- rev(IVS10_go_enrich_BP$Description)
IVS10_go_enrich_BP$Description <- factor(IVS10_go_enrich_BP$Description, levels = level_order)
IVS10_dotplot <- ggplot(data = IVS10_go_enrich_BP) +
  geom_point(aes(x = Description, size = Count, y = pvalue, color = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") + # Use a color gradient from blue (low pvalues) to red (high pvalues)
  ggtitle ("iMGL_IVS10")
IVS10_dotplot


#heatmap:

IVS10_NPX <- subset(all_data, select = c ('PanelAssay', 'mut_1', 'con_1', 'mut_2', 'con_2'))

IVS10_NPX <- list(IVS10_NPX)

IVS10_1 <- IVS10_NPX$mut_1 - IVS10_NPX$con_1
IVS10_2 <- IVS10_NPX$mut_2 - IVS10_NPX$con_2

diff_IVS10 <- cbind(IVS10_1, IVS10_2)
rownames(diff_IVS10) <- IVS10_NPX$PanelAssay

IVS10_hmap <- cbind(diff_IVS10, t_test_uniprot, by = "PanelAssay")

IVS10_hmap <- merge(diff_IVS10, t_test_uniprot, by = "PanelAssay")



#density plot


























d_IVS10_2m <- t_test_protein[ , c('mean.x', 'Protein')]
names(d_IVS10_2m) <- c('IVS10_2m' , 'PanelAssay')

ggplot(d_IVS10_2m, aes(x = IVS10_2m)) +
  geom_density(fill = "blue", alpha = 0.5) +  # Fill color and transparency
  labs(title = "Density Plot of IVS10_2m",
       x = "Mean Values",
       y = "Density") +
  theme_minimal()  






























