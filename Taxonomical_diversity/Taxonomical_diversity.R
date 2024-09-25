####READ ME####
#FILE NAME: Taxonomical_diversity.R     
#DATE: 02/03/2023
#TITLE: Addressing marine restoration success: Evidence of species and functional diversity recovery in a ten-year restored macroalgal forest
#AUTHORS: Cristina Galobart, Enric Ballesteros, Raül Golo, Emma Cebrian
#SCRIPT: Raül Golo (raul.gonzalez@udg.edu) & C. Galobart (cgalobart@ceab.csic.es)
#JOURNAL: Frontiers in Marine Science

# DISCLAMER: This script has been developed by an ecologist, not a programmer, 
# so please take into account that the code may have room to be optimised. 
# Positive feedback will always be more than welcome.

# Script Content
# 1. Load data and data arrangements
# 2. Figure 3A - Biomass boxplot
# 3. Figure 3B - Mean species richness boxplot
# 4. Figure 3C - Shannon-Weiner diversity boxplot
# 5. Test differences between Biomass, Mean species richness and Shannon-Weiner diversity
# 6. Figure 4 -  3D non-metric multidimensional scaling (nMDS) ordination analysis
# 7. PERMANOVA to test for differences on the associated species composition

#Load libraries
library(car)
library(ggplot2)
library(dplyr)
library(sf)
library(agricolae)
library(vegan)
library(rgl)


####1. Load data and data arrangements####
#Set working directory 
setwd("DataPath")
setwd("C:/Users/ecologia/Desktop/Data base & R code to Zenodo/Taxonomical_diversity")

raw_data<-read.csv2("Relative_abundance.csv", sep = ";", dec = ".", stringsAsFactors = F)
#61 sp in total

raw_data$Samp <- as.factor(raw_data$Samp)
raw_data$Locality <- as.factor(raw_data$Locality)

#Split the dataframe 
#Data abundance - only the abundance data without the other information - "Samp", "Locality" and "Biomass"
#Data ambient - Sample, Locality and Biomass information
data_abundance <- select(raw_data, -c(Samp, Locality, Biomass))
data_ambient <- select(raw_data, c(Samp, Locality, Biomass))

#Squareroot transformation of the abundances
data_abundance_sqrt <- sqrt(data_abundance) 
range(data_abundance_sqrt)

#Join both dataframes again
data <- cbind(data_ambient, data_abundance_sqrt)

#Order the levels of Locality factor
data$Locality <- factor(data$Locality, levels = c("Non-restored", "Restored", "Expansion_area", "Reference_1", "Reference_2"))

names <- as.vector(c("Non-restored", "Restored", "Expansion_area", "Reference_1", "Reference_2"))
color <- as.vector(c("#9A2E33", "#CC462A", "#E19644", "#EDDE75", "#7A9F84"))

####2. Figure 3A - Biomass boxplot####
Biomass <- ggplot(data, aes(x=Locality, y=Biomass, fill = Locality)) + ggtitle("Biomass") +
  labs(y= "Biomass (dry weight/m2)") + geom_boxplot(width=0.5)+
  theme_bw() + scale_x_discrete(limits = names)+
  scale_fill_manual(values = color) +
  theme(text = element_text(size = 20)) 

print(Biomass)

#Save the plot as a pdf
ggsave("Results/Biomass.pdf", width=20, height=20, units= "cm")


####3. Figure 3B - Mean species richness boxplot####
#Calculate the number of species per sample
S <- specnumber(data_abundance_sqrt) 
S

Richness <- ggplot(data, aes(x=Locality, y=S, fill=Locality)) + ggtitle("Mean species richness") +
  labs(y= "S") + geom_boxplot(width=0.5)+
  theme_bw() + scale_x_discrete(limits = names)+
  scale_fill_manual(values = color) +
  theme(text = element_text(size = 20)) 

print(Richness)

#Save the plot as a pdf
ggsave("Results/Richness.pdf", width=20, height=20, units= "cm")


####4. Figure 3C - Shannon-Weiner diversity boxplot####
#Calculate the diversity index per sample
#We use 2 as the base of the logarithm, but it could be used exp = 1 (natural logarithms) as well.
H <- diversity(data_abundance_sqrt, index = "shannon", base = 2) 
H

Shannon <- ggplot(data, aes(x=Locality, y=H, fill=Locality)) + ggtitle("Shannon-Weiner diversity") +
  labs(y= "H") + geom_boxplot(width=0.5)+ theme_bw() + 
  scale_fill_manual(values = color) + 
  scale_x_discrete(limits = names) +
  theme(text = element_text(size = 20)) 

print(Shannon)

ggsave("Results/Shannon.pdf", width=20, height=20, units= "cm")

#Prepare dataframe with the taxonomic diversity values
taxo_diversity<- data.frame(data$Locality, data$Biomass, S, H)
#Rename the columns
names(taxo_diversity)[names(taxo_diversity) == "data.Locality"] <- "Locality"
names(taxo_diversity)[names(taxo_diversity) == "data.Biomass"] <- "Biomass"
#Export dataframe
write.table(taxo_diversity, file = "Results/taxo_diversity.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


####5. Test differences between Biomass, Mean species richness and Shannon-Weiner diversity####
#One-way ANOVA with locality (five levels: Non-Restored, Restored, Expansion Area, Reference 1 and Reference 2) as a fixed factor
#Biomass
result.aov <- aov(Biomass ~ Locality, data = taxo_diversity)
summary(result.aov)
#Check ANOVA assumptions
#Homogeneity of variances - Levene's test
leveneTest(Biomass ~ Locality, data = taxo_diversity)
#Normality of the residuals - Shapiro Wilk test
aov_residuals <- residuals(object = result.aov)
shapiro.test(x = aov_residuals)

#Pairwise comparsions - Tukey's honestly significant difference (HSD) test
HSD.test(result.aov, "Locality", console=TRUE)


#Mean species richness
result.aov <- aov(S ~ Locality, data = taxo_diversity)
summary(result.aov)
#Check ANOVA assumptions
#Homogeneity of variances - Levene's test
leveneTest(S ~ Locality, data = taxo_diversity)
#Normality of the residuals - Shapiro Wilk test
aov_residuals <- residuals(object = result.aov)
shapiro.test(x = aov_residuals)

#Pairwise comparsions - Tukey's honestly significant difference (HSD) test
HSD.test(result.aov, "Locality", console=TRUE)

#Shannon-Weiner diversity
result.aov <- aov(H ~ Locality, data = taxo_diversity)
summary(result.aov)
#Check ANOVA assumptions
#Homogeneity of variances - Levene's test
leveneTest(H ~ Locality, data = taxo_diversity)
#Normality of the residuals - Shapiro Wilk test
aov_residuals <- residuals(object = result.aov)
shapiro.test(x = aov_residuals)

#Pairwise comparsions - Tukey's honestly significant difference (HSD) test
HSD.test(result.aov, "Locality", console=TRUE)


####6. Figure 4 -  3D non-metric multidimensional scaling (nMDS) ordination analysis####

#Upload the relative abundance data, but without Gongolaria barbata species as we want to visualise the associated community, without considering the structural species
raw_data <- read.csv2("Relative_abundance_without_Gongolaria.csv", sep = ";", dec = ".", stringsAsFactors = F)
#60 sp in total

raw_data$Samp <- as.factor(raw_data$Samp)
raw_data$Locality <- as.factor(raw_data$Locality)

#Only the abundance data without the other information - "Samp", "Locality"
data_abundance <- select(raw_data, -c(Samp, Locality))
data_ambient <- select(raw_data, c(Samp, Locality))

#Squareroot transformation of the abundances
data_abundance_sqrt <- sqrt(data_abundance) 
range(data_abundance_sqrt)

#Join both dataframes again
data <- cbind(data_ambient, data_abundance_sqrt)

#Do the NMDS analyses with three dimensions
Comm_NMDS<-metaMDS(data_abundance_sqrt, k=3, trymax=300, autotransform = F)
# Trymax for the function to do more iterations. 
# k = Number of dimensions 
# Use a bigger number if it does not converge
# Bray-Curtis dissimilarity is used as default

#See NMDS results
Comm_NMDS
Comm_NMDS$diss

#Safe the "stress" value
stress <- Comm_NMDS$stress
stressplot(Comm_NMDS)

#X, Y and Z coordenates to draw the plot, of each sampling site
Comm_NMDS$points

#Create the dataframe with the data for the 3D NMDS
Com_abun<-data.frame(x=Comm_NMDS$point[,1], y=Comm_NMDS$point[,2], z=Comm_NMDS$point[,3], Locality=as.factor(data_ambient$Locality))

Com_abun <- Com_abun %>%
  mutate(color = case_when(Locality == "Non-restored" ~ "#9A2E33",
                          Locality == "Restored" ~ "#CC462A",
                          Locality == "Expansion_area" ~ "#E19644",
                          Locality == "Reference_1" ~ "#EDDE75",
                          Locality == "Reference_2" ~ "#7A9F84"))
  
Com_abun <- as.data.frame(Com_abun)
Com_abun$color <- as.factor(Com_abun$color)
Com_abun

p <- plot3d(x = Com_abun$x, y = Com_abun$y, z = Com_abun$z,
       col = Com_abun$color, type = "s", size = 4,
       xlab = "NMDS axis 1", ylab = "NMDS axis 2", zlab = "NMDS axis 3")

#Rotate the 3D NMDS so that the X, Y and Z axes are orderd and take a screenshot
rgl.snapshot("Results/NMDS_3D.png")

#To create a GIF considering plot rotation
movie3d(spin3d(axis = c(0,1,0), rpm = 4), duration = 15, dir ="./Results")
#Change the c(0,0,0) to set the rotation axis

####7. PERMANOVA to test for differences on the associated species composition####
permanova <- adonis2(data_abundance_sqrt~Locality, data=data, permutations=999, method ="euclidean")
permanova


