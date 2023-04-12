####READ ME####
#FILE NAME: Functional_diversity.R
#DATE: 02/03/2023
#TITLE: Addressing marine restoration success: Evidence of species and functional diversity recovery in a ten-year restored macroalgal forest
#AUTHORS: Cristina Galobart, Enric Ballesteros, Raül Golo, Emma Cebrian
#JOURNAL: Frontiers in Marine Science - under review

#Script obtained from Teixidó et al. 2018 (Nature Communications, DOI: 10.1038/s41467-018-07592-1)
#Modified by: Raül Golo (raul.gonzalez@udg.edu) (automation) & C. Galobart (cgalobart@ceab.csic.es)

# DISCLAMER: This script has been modified by ecologists, not programmers,
# so please take into account that the code may have room to be optimised.
# Positive feedback will always be more than welcome.

# Script content
# 1. Load data and data arrangements
# 2. Overall functional space
# 3. Functional richness
# 4. Figure 5. Species and functional diversity changes among the five localities
# 5. Intersection of functional space between pair of localities
# 6. Load data and data arrangements
# 7. Figure 6 - Distribution of FE abundance across the functional space
# 8. Data arrangements
# 9. Figure 7 - Relative abundances of each trait category in each locality

#Set your working directory here
setwd("C:/Functional_diversity")

#Load required packages
library(FD)
library(tripack)
library(matrixStats)
library(dplyr)
library(RColorBrewer)
library(rcdd)
library(geometry)
library(ggiraphExtra)
library(geosphere)
library(tibble)
library(ggplot2)
library(patchwork)


####1. Load and data arrangements####

#Load abundance data
ab <- read.csv2("Your_data/Data_Ab_example.csv", sep=";", dec=".", row.names = 1)

#To delete columns without abundance, if any
ab <- ab[, colSums(ab != 0) > 0]

#Load Species and their FEs
#Code to automate script - make ab & spe_fes species names match
spe_fes <- read.csv2("Data_Species_FEs.csv", sep=";", dec=".", row.names=)

presence.names<- names(ab)

spe_fes <- lapply(presence.names, function(x) {

  spe_fes <- spe_fes[spe_fes$Species == x,]

})

spe_fes <- do.call(rbind, spe_fes)

Sp <- as.character(spe_fes$Species)

rownames(spe_fes) <- Sp

spe_fes <- as.data.frame( spe_fes[,-1, drop=F])


# Load FEs data
fes <- read.csv2("Data_FEs.csv", sep=";", dec=".", row.names=, stringsAsFactors = T)

str(fes)
fes$Morphological.form <- as.factor(fes$Morphological.form)
fes$Solitary.Colonial <- as.factor(fes$Solitary.Colonial)
fes$Max.Longevity <- as.factor(fes$Max.Longevity)
fes$Height <- as.factor(fes$Height)
fes$Width <- as.factor(fes$Width)
fes$Epibiosis <- as.factor(fes$Epibiosis)
fes$Major.photosynthetic.pigments <- as.factor(fes$Major.photosynthetic.pigments)
fes$Age.reproductive.maturity <- as.factor(fes$Age.reproductive.maturity)
fes$Asexual.Reproduction <- as.factor(fes$Asexual.Reproduction)
fes$Growth.rates <- as.factor(fes$Growth.rates)
fes$Calcification <- as.factor(fes$Calcification)
fes$Chemical.defenses <- as.factor(fes$Chemical.defenses)
fes$Biological.cycle <- as.factor(fes$Biological.cycle)
fes$Vegetative.reproduction <- as.factor(fes$Vegetative.reproduction)
str(fes)

x <- "FE"
presence.names<-as.character(spe_fes[,x])

fes <- lapply(presence.names, function(x) {

  fes <- fes[fes$FE == x,]

})

fes <- do.call(rbind, fes)

fes <- distinct(fes,FE, .keep_all= TRUE)

FE <- as.character(fes[,x])

rownames(fes) <- FE

fes <- fes[,-1]

#Load sites
sites <- read.csv2("Your_data/Data_Site_example.csv", sep=";", dec=".", row.names=1)

# Defining treatments
condition <-  unique( as.character(sites$Treat))

####2. Overall functional space####
#computing  multidimensional functional spaces (2 to 14 D)
#load additional functions

source("Necs_scripts/quality_funct_space.R")

qfs <- quality_funct_space(fes, traits_weights=NULL, nbdim=14, metric="Gower", dendro=FALSE, plot="quality_funct_space")

# quality of spaces (low meanSD = high quality)
round( qfs$meanSD , 4)

# keeping coordiantes of the 5 dimensions, meanSD = 0.0082
fd.coord <- qfs$details_funct_space$mat_coord[,1:5]

write.csv(fd.coord, file="Results/FE_5D_coord.csv") #to use it for further analyses

#see variance explained by the PCoA axes (64%)
gower<-qfs$details_funct_space$mat_dissim

fit <- cmdscale(gower,eig=TRUE, k=5) # PCoA

cumsum(fit$eig[fit$eig>=0]) / sum(fit$eig[fit$eig>0])

####3. Functional richness####

# Data arrangements
ab.conditions <- lapply(condition, function(x) {

  quad <- rownames(sites[sites$Treat == x,])

  colSums(ab[rownames(ab) %in% quad,])

})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

#Calculate convex hull
Fric <- lapply(condition, function (x) {

  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]

  ch <- convhulln(m, options = "FA")

  chg <- convhulln(fd.coord, options = "FA")

  c(length(species), length(species)/length(Sp)*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol) #/chg$vol*100

})#eo lapply

names(Fric) = condition

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP andNbFEsP), and the volume among the 5 localities (Vol5D)
Fric <- do.call(rbind, Fric)
colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol5D")

cols <- as.vector(c("#9A2E33", "#CC462A", "#E19644", "#EBCC82", "#7A9F84"))
colstr <- paste(cols, "70", sep="")
names(cols) <- condition
names(colstr) <- condition

####4. Figure 5 - Species and functional diversity changes among the five localities####
n_axes = 5
labels_fig_cv1 <- c("Sp", "FE")
labels_fig_cv2 <- c(paste0("Vol. ",n_axes,"D"))

pdf("Results/Figure_5.pdf", pointsize=21, height = 10, width = 18)

par(mfrow = c(3,5), mar = c(1, 4, 4, 1)) # change the last number, that matches the number of treatments

Fric1 <- as.matrix(Fric[,c(1,3)])
Fric2 <- as.matrix(Fric[,c(5)])

Null <- c(rep(0, 5))

Fric2 <- cbind(Fric2, Null)

colnames(Fric1) = labels_fig_cv1
colnames(Fric2) = c(labels_fig_cv2, "Null")

for (i in condition) {

  midpoints <-  barplot(Fric1[i,], ylab="Richness", col=cols[i], ylim=c(0,105), main=paste0(i), col.main=cols[i], cex.main=1)

  if(Fric1[1,1] > 1) {

    lab = round(Fric1[i,],0)

  } else {

    lab = round(Fric1[i,],2)

  }#eo ifelse

  text(midpoints, Fric1[i,]+8, labels=lab, col=cols[i], cex=1, font=2)

}#eo for barplot 1

for (i in condition) {

  midpoints <-  barplot(Fric2[i,], ylab="Richness", col=cols[i], ylim=c(0,0.01), cex.axis = 0.5)

  text(midpoints, Fric2[i,]+0.0008, labels=as.character(round(Fric2[i,],4)), col=cols[i], cex=1, font=2)

}#eo for barplot 2

for (i in condition) {

  species <- colnames(ab.conditions)[which(ab.conditions[i,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond, ]

  tr <-tri.mesh(m[,1],m[,2])
  ch <- convex.hull(tr)

  plot(fd.coord[,1], fd.coord[,2], xlab = "PCoA 1", ylab = "PCoA 2", type="n")

  polygon(ch, col=colstr[i], border=cols[i])
  points(m[,1:2], pch = 16, col=cols[i])

}#eo for convex

dev.off()
#This figure has been modified in InkScape to improve visibility and clarity.

####5. Intersection of functional space between pair of localities####

#Load intersect function to compute convex hull (vertices + volume) of two set of points and their intersection
source("Necs_scripts/intersect.R")

mat_int <- lapply(condition, function (x) {

  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]

  return(m)

})#eo lapply

names(mat_int) = condition

#Localities
#1 - No Rest
#2 - Rest
#3 - Exp Area
#4 - Ref 1
#5 - Ref 2

#Intersections represented in Figure 8 are calculated
#Intersect No Rest within Rest
NoRest_int_Rest <- CHVintersect(mat_int[[names(cols[1])]],mat_int[[names(cols[2])]]) 
#Percentage No Rest within Rest
NoRest_int_Rest$vol[3]/NoRest_int_Rest$vol[1]

#Intersect Rest within Ref 1
Rest_int_Ref1 <- CHVintersect(mat_int[[names(cols[2])]],mat_int[[names(cols[4])]]) 
#Percentage Rest within Ref 1
Rest_int_Ref1$vol[3]/Rest_int_Ref1$vol[1]

#Intersect Rest within Ref 2
Rest_int_Ref2 <- CHVintersect(mat_int[[names(cols[2])]],mat_int[[names(cols[5])]]) 
#Percentage Rest within Ref 2
Rest_int_Ref2$vol[3]/Rest_int_Ref2$vol[1]

#Intersect Rest within Exp Area
Rest_int_Exp <- CHVintersect(mat_int[[names(cols[2])]],mat_int[[names(cols[3])]])
#Percentage Rest within Exp Area
Rest_int_Exp$vol[3]/Rest_int_Exp$vol[1]


####6. Load data and data arrangements####

#Load Species and FE data again
spe_fes <- read.csv2("Data_Species_FEs.csv", sep=";", dec=".", row.names=, stringsAsFactors = T)

presence.names<- names(ab)

spe_fes <- lapply(presence.names, function(x) {

  spe_fes <- spe_fes[spe_fes$Species == x,]

})

spe_fes <- do.call(rbind, spe_fes)

Sp <- as.character(spe_fes$Species)

ab.conditions <- lapply(condition, function(x) {

  quad <- rownames(sites[sites$Treat == x,])

  colSums(ab[rownames(ab) %in% quad,])

})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

#Compute abundance of FEs for the five localities

fes_names <- levels(spe_fes$FE)

ab.fe.conditions <- lapply(condition, function (z) {
  
  abund.fes <-  sapply(fes_names, function (x) {
    
    spec <- as.character(spe_fes[which(spe_fes$FE == x),]$Species)
    
    sum(ab.conditions[z,spec])
    
  })#eo sapply
  
  abund.fes
  
})#eo lapply

names(ab.fe.conditions) = condition

ab.fe.conditions <- do.call(rbind, ab.fe.conditions)

####7. Figure 6 - Distribution of FE abundance across the functional space####

n_axes = 5

labels_fig_cv <- c("NbSp", "NbFE", paste0("Vol. ",n_axes,"D"))

pdf("Results/Figure_6.pdf", height=12, width=25, pointsize=24)

par(mfrow=c(1,5)) # The last number matches the number of treatments

ab.fe.conditions <- ab.fe.conditions[, colSums(ab.fe.conditions != 0) > 0]

#Non-Restored
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main=names(cols[1]), col.main=cols[1])
points(fd.coord[,1], fd.coord[,2], pch=21, col=cols[1] , bg=colstr[1], cex=sqrt(ab.fe.conditions[names(cols[1]),]))#/100*30)*3)

#Restored
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main=names(cols[2]), col.main=cols[2])
points(fd.coord[,1], fd.coord[,2], pch=21, col=cols[2] , bg=colstr[2], cex=sqrt(ab.fe.conditions[names(cols[2]),]))#/100*30)*3)

#Expansion Area
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main=names(cols[3]), col.main=cols[3])
points(fd.coord[,1], fd.coord[,2], pch=21, col=cols[3] , bg=colstr[3], cex=sqrt(ab.fe.conditions[names(cols[3]),]))#/100*30)*3)

#Reference 1
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main=names(cols[4]), col.main=cols[4])
points(fd.coord[,1], fd.coord[,2], pch=21, col=cols[4] , bg=colstr[4], cex=sqrt(ab.fe.conditions[names(cols[4]),]))#/100*30)*3)

#Reference 2
plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main=names(cols[5]), col.main=cols[5])
points(fd.coord[,1], fd.coord[,2], pch=21, col=cols[5] , bg=colstr[5], cex=sqrt(ab.fe.conditions[names(cols[5]),]))

dev.off()


####8. Data arrangements####

#List to store results
ab.trait_condition<-list()

#Ordering FEs in trait table as in abundance table
fes_traits<-fes[ colnames(ab.fe.conditions), ]

#For each trait, computing total relative abundance among localities
for (t in colnames(fes_traits) )
{
  # levels of trait t
  levels_t<-as.character( sort( unique(fes_traits[,t]) ) )
  
  # empty vectors to store results
  trait_t<-c()
  condition_t<-c()
  rel_ab_t<-c()
  
  # computing relative abundance of each trait modality in each locality
  for (i in condition)
    for (j in levels_t)
    {
      trait_t<-c(trait_t, j)
      condition_t<-c(condition_t, i)
      rel_ab_t<-c(rel_ab_t, sum( ab.fe.conditions[ i , which(fes_traits[,t]==j) ] )  )
    }# end of i,j
  
  # setting correct order for levels of localities
  condition_t <- factor(condition_t, levels = condition )
  
  # storing results in a dataframe
  ab.trait_condition[[t]]<-data.frame(  trait_val= trait_t,  condition= condition_t, rel_ab_FE=rel_ab_t )
  
}# end of t

####9. Figure 7 - Relative abundances of each trait category in each locality####
#ggplot does not like list of dataframes; thus we code 14 individual plots and later we grouped them together

#Set the colors as needed
cc8<-c("#DF775C","#EFBAA6","#FFE9C3","#EFBA67", "#C4D27B", "#225F6B", "#077962","#70B89A") 
cc7<-c("#DF775C","#EFBAA6","#FFE9C3", "#C4D27B", "#225F6B", "#077962","#70B89A") 
cc6<-c("#DF775C","#EFBAA6","#FFE9C3", "#225F6B", "#077962","#70B89A") 
cc5<-c("#DF775C","#EFBAA6","#FFE9C3", "#077962","#70B89A") 
cc4<-c("#DF775C","#EFBAA6","#FFE9C3", "#70B89A") 
cc3<-c("#DF775C","#EFBAA6", "#70B89A") 
cc2<-c("#DF775C", "#70B89A") 

#plot1: Morphological.form, 8 trait values, cc8
Ftrait<-as.data.frame(ab.trait_condition$Morphological.form)

Ftrait$trait_val<-as.character(Ftrait$trait_val)

Ftrait$trait_val[Ftrait$trait_val == "a"] <- "A"
Ftrait$trait_val[Ftrait$trait_val == "b"] <- "B"
Ftrait$trait_val[Ftrait$trait_val == "e"] <- "E"
Ftrait$trait_val[Ftrait$trait_val == "f"] <- "F"
Ftrait$trait_val[Ftrait$trait_val == "g"] <- "G"
Ftrait$trait_val[Ftrait$trait_val == "j"] <- "J"
Ftrait$trait_val[Ftrait$trait_val == "l"] <- "L"
Ftrait$trait_val[Ftrait$trait_val == "m"] <- "M"

Ftrait$trait_val<-as.factor(Ftrait$trait_val)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_1M_form.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot1 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Morphological form")+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values= cc8)

plot1

#plot2: Solitary.Colonial, 3 trait values, cc3
Ftrait<-as.data.frame(ab.trait_condition$Solitary.Colonial)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_2Colonial.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot2 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Solitary-Colonial")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values= cc3)

plot2

#plot 3 Max.Longevity, 7 trait values, cc7
Ftrait<-as.data.frame(ab.trait_condition$Max.Longevity)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_3Max_long.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot3 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Max Longevity")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values= cc7)

plot3

#plot 4 Height, 5 trait values, cc5
Ftrait<-as.data.frame(ab.trait_condition$Height)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_4Height.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot4 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Height")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values= cc5)

plot4

#plot 5 Width, 5 trait values, cc5
Ftrait<-as.data.frame(ab.trait_condition$Width)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_5Width.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot5 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Width")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values= cc5)

plot5

#plot 6 Epibiosis, 3 trait values, cc3
Ftrait<-as.data.frame(ab.trait_condition$Epibiosis)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_6Epibiosis.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot6 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Epibiosis")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values= cc3)

plot6

#plot 7 Major.photosynthetic.pigments, 5 trait values, cc5
Ftrait<-as.data.frame(ab.trait_condition$Major.photosynthetic.pigments)

Ftrait$trait_val<-as.character(Ftrait$trait_val)

Ftrait$trait_val[Ftrait$trait_val == "a"] <- "A"
Ftrait$trait_val[Ftrait$trait_val == "b"] <- "B"
Ftrait$trait_val[Ftrait$trait_val == "c"] <- "C"

Ftrait$trait_val<-as.factor(Ftrait$trait_val)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_7Pigments.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot7 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(y="Relative Abundance of Functional Categories (%)", fill="Traits")+
  ggtitle("Major photosynthetic pigments")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc3)

plot7

#plot 8 Age.reproductive.maturity, 7 trait values, cc7
Ftrait<-as.data.frame(ab.trait_condition$Age.reproductive.maturity)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_8Maturity.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot8 <-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Age reproductive maturity")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc5)

plot8

#plot9: Asexual.Reproduction, 2 trait values, cc2
Ftrait<-as.data.frame(ab.trait_condition$Asexual.Reproduction)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_9Asexual.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot9<-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(x="",  fill="Traits")+
  ggtitle("Asexual Reproduction")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc2)

plot9

#plot10: Growth.rates, 5 trait values, cc5
Ftrait<-as.data.frame(ab.trait_condition$Growth.rates)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_10Growth_Rate.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot10<-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Growth rates")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc5)

plot10

#plot11: Calcification, 3 trait values, cc3
Ftrait<-as.data.frame(ab.trait_condition$Calcification)

Ftrait$trait_val<-as.character(Ftrait$trait_val)

Ftrait$trait_val[Ftrait$trait_val == "a"] <- "A"
Ftrait$trait_val[Ftrait$trait_val == "c"] <- "C"
Ftrait$trait_val[Ftrait$trait_val == "d"] <- "D"

Ftrait$trait_val<-as.factor(Ftrait$trait_val)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_11Calcificat.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot11<-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Calcification")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc3)

plot11

#plot12: Chemical.defenses, 2 trait values, cc2
Ftrait<-as.data.frame(ab.trait_condition$Chemical.defenses)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_12Defences.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot12<-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Chemical defenses")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc2)

plot12

#plot13: Biological.cycle, 4 trait values, cc4
Ftrait<-as.data.frame(ab.trait_condition$Biological.cycle)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_13Cycle.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot13<-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Biological cycle")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc4)

plot13

#plot14: Vegetative.reproduction, 4 trait values, cc4
Ftrait<-as.data.frame(ab.trait_condition$Vegetative.reproduction)

Ftrait_Non_Rest<-Ftrait[Ftrait$condition == condition[1],]
FT_Non_Rest_percent <- Ftrait_Non_Rest$rel_ab_FE / sum(Ftrait_Non_Rest$rel_ab_FE) * 100
Ftrait_Non_Rest<-Ftrait_Non_Rest[,-3]
Ftrait_Non_Rest$rel_ab_FE <- NA
Ftrait_Non_Rest$rel_ab_FE <- FT_Non_Rest_percent

Ftrait_Rest<-Ftrait[Ftrait$condition == condition[2],]
FT_Rest_percent <- Ftrait_Rest$rel_ab_FE / sum(Ftrait_Rest$rel_ab_FE) * 100
Ftrait_Rest<-Ftrait_Rest[,-3]
Ftrait_Rest$rel_ab_FE <- NA
Ftrait_Rest$rel_ab_FE <- FT_Rest_percent

Ftrait_Exp<-Ftrait[Ftrait$condition == condition[3],]
FT_Exp_percent <- Ftrait_Exp$rel_ab_FE / sum(Ftrait_Exp$rel_ab_FE) * 100
Ftrait_Exp<-Ftrait_Exp[,-3]
Ftrait_Exp$rel_ab_FE <- NA
Ftrait_Exp$rel_ab_FE <- FT_Exp_percent

Ftrait_Ref1<-Ftrait[Ftrait$condition == condition[4],]
FT_Ref1_percent <- Ftrait_Ref1$rel_ab_FE / sum(Ftrait_Ref1$rel_ab_FE) * 100
Ftrait_Ref1<-Ftrait_Ref1[,-3]
Ftrait_Ref1$rel_ab_FE <- NA
Ftrait_Ref1$rel_ab_FE <- FT_Ref1_percent

Ftrait_Ref2<-Ftrait[Ftrait$condition == condition[5],]
FT_Ref2_percent <- Ftrait_Ref2$rel_ab_FE / sum(Ftrait_Ref2$rel_ab_FE) * 100
Ftrait_Ref2<-Ftrait_Ref2[,-3]
Ftrait_Ref2$rel_ab_FE <- NA
Ftrait_Ref2$rel_ab_FE <- FT_Ref2_percent

Ftrait2<-rbind(Ftrait_Non_Rest, Ftrait_Rest, Ftrait_Exp, Ftrait_Ref1, Ftrait_Ref2)

# write.table(Ftrait2, file = "Results/Ftrait_14Vegetative.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, dec = ".")

plot14<-ggplot(data=Ftrait2, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Vegetative reproduction")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = cc4)

plot14

#Group all the plots
all.plots.vert<-(plot1|plot2|plot3) / (plot4|plot5|plot6) /
  (plot7|plot8|plot9) / (plot10|plot11|plot12) /  (plot13|plot14|plot14)
all.plots.vert

ggsave("Results/Figure_7.pdf", plot= all.plots.vert, device="pdf", height=40, width=35, units="cm", dpi=300)

