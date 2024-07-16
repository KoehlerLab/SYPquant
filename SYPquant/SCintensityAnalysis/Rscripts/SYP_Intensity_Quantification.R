library(ggplot2)
library(viridis)
library(sp)
library(autoimage)
library(mgcv) #required for gam function
library(tidyverse)
library(readxl)
library(tidyquant)
library(ggbio)
library(grid)
library(gridExtra)
library(scales)
library(cowplot)
library(zoo)
library(dplyr)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
library(rcompanion)
library(ggh4x)

setwd('./SC_along_Axis_includesHexanediol')

metadata_file <- './SYP4mutants_SC_quantification_w-o_Hex.xlsx'
metadata <- read_xlsx(metadata_file, col_names = TRUE)
metadata <- metadata[metadata$USE & metadata$Hexanediol==0,]
metadata <- data.frame(lapply(metadata, function(x) {gsub("/g/", '/Volumes/',x)}))

###Choose filters
NucleusVolume.um3.UpperCut <- 60
NucleusVolume.um3.LowerrCut <- 10
SphericityLowerCut <- 0.4

XY_pixel_size <- as.numeric(29.5385) #for 60x with SoRA 3.2x amplification

##Color scheme
WT_Color<-"grey"
D114_Color <-"#38a3a5"
AN90_Color <-"#9a031e" 
Deletion_Color <-"#1520a6"
Phospho_Color <- "#fb8b24"

##### Starts here 
poligon_info <- NULL
data <- NULL

for (file in metadata$File){
  print(file)
  tmp_poligon <- read.csv(metadata$Polygon_csv_path[metadata$File==file], header=TRUE)
  poligon_per_file <- tmp_poligon[grepl(paste(metadata$Gonad[metadata$File==file],"_",sep=""),tmp_poligon$File),]
  print(poligon_per_file)
  poligon_per_file$full_path <- file
  poligon_info <- rbind(poligon_info, poligon_per_file)
  
  #print(poligon_info)
  
  tmp_SC <- read.csv(metadata$HTP3_output_csv[metadata$File==file], header=TRUE)
  SC_per_file <- tmp_SC[grepl(paste(metadata$Gonad[metadata$File==file],"_",sep=""), tmp_SC$Tile),]
  SC_per_file$Folder <- metadata$Folder[metadata$File==file]
  SC_per_file$full_path <- file
  SC_per_file$Slide <- metadata$Sample[metadata$File==file]
  SC_per_file$GT <- metadata$Genotype[metadata$File==file]
  data <- rbind(data, SC_per_file)
}

poligon_info$pol.x.pxl <- poligon_info$pol.x*XY_pixel_size
poligon_info$pol.y.pxl <- poligon_info$pol.y*XY_pixel_size

  
data$Cx_rotated=0
data$Cy_rotated=0
data$poligon=0

#Poligon starts at the top of the gonad begining of transition zone, orientation of the gonad from left to right!

for (file in unique(data$full_path)) {
  data$poligon[data$full_path==file]=point.in.polygon(data$Cx[data$full_path==file], data$Cy[data$full_path==file], poligon_info$pol.x.pxl[poligon_info$full_path==file], poligon_info$pol.y.pxl[poligon_info$full_path==file], mode.checked=FALSE)
  center_point=c((poligon_info$pol.x.pxl[poligon_info$full_path==file][length(poligon_info$pol.x.pxl[poligon_info$full_path==file])-1]+ poligon_info$pol.x.pxl[poligon_info$full_path==file][2])/2,
                (poligon_info$pol.y.pxl[poligon_info$full_path==file][length(poligon_info$pol.y.pxl[poligon_info$full_path==file])-1] + poligon_info$pol.y.pxl[poligon_info$full_path==file][2])/2) #finds the center point of the line corresponding to end of TZ
  poligon_info$pol.x.pxl[poligon_info$full_path==file]=poligon_info$pol.x.pxl[poligon_info$full_path==file]-center_point[1]
  poligon_info$pol.y.pxl[poligon_info$full_path==file]=poligon_info$pol.y.pxl[poligon_info$full_path==file]-center_point[2]
  data$Cx_rotated[data$full_path==file]=data$Cx[data$full_path==file]-center_point[1]
  data$Cy_rotated[data$full_path==file]=data$Cy[data$full_path==file]-center_point[2]
  
  theta = atan(-1/((poligon_info$pol.y.pxl[poligon_info$full_path==file][length(poligon_info$pol.y.pxl[poligon_info$full_path==file])-1] - poligon_info$pol.y.pxl[poligon_info$full_path==file][2])/
                (poligon_info$pol.x.pxl[poligon_info$full_path==file][length(poligon_info$pol.x.pxl[poligon_info$full_path==file])-1] - poligon_info$pol.x.pxl[poligon_info$full_path==file][2])))
  beta = atan2((poligon_info$pol.y.pxl[poligon_info$full_path==file][length(poligon_info$pol.y.pxl[poligon_info$full_path==file])-1] - poligon_info$pol.y.pxl[poligon_info$full_path==file][2]),
               (poligon_info$pol.x.pxl[poligon_info$full_path==file][length(poligon_info$pol.x.pxl[poligon_info$full_path==file])-1] - poligon_info$pol.x.pxl[poligon_info$full_path==file][2]))
  #print(c(file, theta*(180/pi),beta*(180/pi)))

  theta[theta<0]=theta[theta<0] + 2*pi #radians
  theta[beta>0]=pi+theta[beta>0] #radians pi+theta[beta>0]
  rotated_coord=autoimage::rotate(cbind(data$Cx_rotated[data$full_path==file], data$Cy_rotated[data$full_path==file]), -theta, pivot = c(0,0))
  data$Cx_rotated[data$full_path==file]=rotated_coord[,1]*(-1) ##I have to multiply this because the gonads are oriented from right to left
  data$Cy_rotated[data$full_path==file]=rotated_coord[,2]
  
  #print(c(file, theta*(180/pi)))
}

df<-subset(data, data$poligon>=1)
df <- df %>% 
  group_by(full_path) %>% 
  arrange(Cx_rotated) 

#All data
p0 <- ggplot(data, aes(x = Cx, y = Cy, color=MeanI_SC_SYP.Bkg)) + geom_point() +facet_wrap(~ full_path, ncol=4) +
  scale_colour_viridis(option="D") + coord_fixed()
p0

#Only data inside the poligon
InPoligon <- ggplot(df, aes(x = Cx, y = Cy, color=MeanI_SC_SYP.Bkg)) + geom_point() +facet_wrap(~ full_path, ncol=4) +
  scale_colour_viridis(option="D") + coord_fixed()
InPoligon

#Only data inside the poligon after reorientation (gonad starts from left to right)
Reorientated <- ggplot(df, aes(x = Cx_rotated, y = Cy_rotated, color=MeanI_SC_SYP.Bkg)) + geom_point(alpha=0.3) + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed() 
Reorientated

#Check nuclei volumes for each gonad
Vol1 <- ggplot(df) + geom_boxplot(aes(x=VolumeNucleus.um3., y=full_path)) + scale_y_discrete(labels=c(unique(df$Tile)))
Vol1

Vol2 <- ggplot(df[df$full_path %in% unique(df$full_path)[1:8] & df$VolumeNucleus.um3.<10,]) + geom_point(aes(x = Cx_rotated, y = Cy_rotated, color=VolumeNucleus.um3.)) + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed()
Vol2

Vol3 <- ggplot(df[df$full_path %in% unique(df$full_path)[1:8] & df$VolumeNucleus.um3.>10 & df$VolumeNucleus.um3.<60,]) + geom_point(aes(x = Cx_rotated, y = Cy_rotated, color=VolumeNucleus.um3.)) + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed()
Vol3

table(df$VolumeNucleus.um3.>NucleusVolume.um3.UpperCut)
table(df$VolumeNucleus.um3.<NucleusVolume.um3.LowerrCut)

#Filter nuclei with volumes bigger than 60 um3 (a pachytene nucleus should be around that volume)
df <- df[df$VolumeNucleus.um3.<NucleusVolume.um3.UpperCut,]
df <- df[df$VolumeNucleus.um3.>NucleusVolume.um3.LowerrCut,]

#Check sphericity for each gonad
Spheri <- ggplot(df) + geom_boxplot(aes(x=Sphericity, y=full_path)) + scale_y_discrete(labels=c(unique(df$Tile)))
Spheri

before <- ggplot(df[df$full_path %in% unique(df$full_path)[1:8],], aes(x = Cx_rotated, y = Cy_rotated)) + geom_point(alpha=0.3, color='darkred') +
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed()

after <- ggplot(df[df$full_path %in% unique(df$full_path)[1:8] & df$Sphericity>0.4,], aes(x = Cx_rotated, y = Cy_rotated)) + geom_point(alpha=0.3, color='darkred') + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed()

grid.arrange(before, after, ncol=2)

table(df$Sphericity<SphericityLowerCut)

df <- df[df$Sphericity>SphericityLowerCut,]

###After filtering
Filtered <- ggplot(df[df$full_path %in% unique(df$full_path)[1:8],], aes(x = Cx_rotated, y = Cy_rotated, color=MeanI_SC_SYP.Bkg)) + geom_point() +facet_wrap(~ full_path, ncol=4) +
  scale_colour_viridis(option="D") + coord_fixed()
Filtered


###Calculate smoothing line
df$Y_fitted=0
df$GonadLengthpxl=0

for (tile in unique(df$full_path)) {
  sub<-subset(df, df$full_path==tile)
  #calculate loess fit
  loess_fit=loess(Cy_rotated ~ Cx_rotated, data=sub, span=0.5)
  df$Y_fitted[df$full_path==tile]=loess_fit$fitted
  
  #calculate eucledian distances
  m<-matrix(c(sub$Cx_rotated, loess_fit$fitted), ncol = 2, byrow = F)
  
  x1<-m[-length(m[,1]),]
  x2<-m[-1,]
  
  dist<-sqrt(rowSums((x1 - x2)^2))
  
  d<-c(0)
  for (value in 1:length(dist)){
    sum=d[value]+dist[value]
    d<-c(d, sum)
  }
  df$GonadLengthpxl[df$full_path==tile]=d
}


###Checking splines for all gonads
Fitted <- ggplot(df[df$full_path %in% unique(df$full_path)[1:8],], aes(x = Cx_rotated, y = Cy_rotated, color=MeanI_SC_SYP.Bkg)) + geom_point(alpha=0.3) + geom_line(aes(x=Cx_rotated, y=Y_fitted)) + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed() 
Fitted

Fitted <- ggplot(df[df$full_path %in% unique(df$full_path)[9:16],], aes(x = Cx_rotated, y = Cy_rotated, color=MeanI_SC_SYP.Bkg)) + geom_point(alpha=0.3) + geom_line(aes(x=Cx_rotated, y=Y_fitted)) + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed() 
Fitted

Fitted <- ggplot(df[df$full_path %in% unique(df$full_path)[17:length( unique(df$full_path))],], aes(x = Cx_rotated, y = Cy_rotated, color=MeanI_SC_SYP.Bkg)) + geom_point(alpha=0.3) + geom_line(aes(x=Cx_rotated, y=Y_fitted)) + 
  facet_wrap(~ full_path, ncol=4) + scale_colour_viridis(option="D") + coord_fixed() 
Fitted

####################################
#Normalise Gonad length (TZ+Pachytene)

df$GonadLengthpxlNorm <- 0
for (file in unique(df$full_path)){
  df$GonadLengthpxlNorm[df$full_path==file]<-(df$GonadLengthpxl[df$full_path==file]-min(df$GonadLengthpxl[df$full_path==file]))/(max(df$GonadLengthpxl[df$full_path==file])-min(df$GonadLengthpxl[df$full_path==file])) #set 0 to begining of TZ
} 

NormalisedLength <- ggplot(data=df) + geom_point(aes(x=GonadLengthpxlNorm, y=MeanI_SC_SYP.Bkg), alpha=0.5) + geom_smooth(aes(x=GonadLengthpxlNorm, y=MeanI_SC_SYP.Bkg), method="loess") +
  facet_wrap(~ full_path, ncol=4) + xlab('Normalised Pachytene Length')
NormalisedLength


#Normalise SC signal

df$MeanI_SC_SYP.Bkg.DAPInorm <- df$MeanI_SC_SYP.Bkg/df$I_Total_DAPI.Bkg  #per nucleus the mean is not really informative 


NormalisedLoading <- ggplot(data=df) + geom_point(aes(x=GonadLengthpxlNorm, y=MeanI_SC_SYP.Bkg.DAPInorm, color=GT), alpha=0.1) + geom_smooth(aes(x=GonadLengthpxlNorm, y=MeanI_SC_SYP.Bkg.DAPInorm, color=GT, group=full_path), method="loess")
NormalisedLoading


#Normalise normalised SC signal per slide
WT <- "ie29-10"

df$I_Total_SYPSC.Bkg.per.slide <- 0
for (folder in unique(df$Folder)){
  for (sample in unique(df$Slide[df$Folder==folder])) {
    print(folder)
    print(sample)
    df$I_Total_SYPSC.Bkg.per.slide[df$Folder==folder & df$Slide==sample] <- df$I_Total_SYPSC.Bkg[df$Folder==folder & df$Slide==sample]/mean(df$I_Total_SYPSC.Bkg[df$Folder==folder & df$Slide==sample & df$GT==WT])
  }
}


#ratio between sum intensities of SYP on SC over sum intensities of SYP on the nucleus
df$ratioI_Total_SYPSC.Bkg.I_Total_SYP.Bkg <- df$I_Total_SYPSC.Bkg/df$I_Total_SYP.Bkg  

#ratio between sum intensities of SYP on SC over sum intensities of SYP on the nucleoplasm
df$ratioI_Total_SYPSC.Bkg.I_Total_SYPNuc.Bkg <- df$I_Total_SYPSC.Bkg/df$I_Total_SYPNuc.Bkg


#Binning 
df$GonadLengthpxlNorm_bin <- floor(df$GonadLengthpxlNorm*11)/10
df[df$GonadLengthpxlNorm_bin>1,'GonadLengthpxlNorm_bin'] <- 1


df.plot <- data.frame(df) %>% group_by(full_path,GonadLengthpxlNorm_bin,GT) %>% summarise(mean_I_Total_SYPSC.Bkg.per.slide = mean(I_Total_SYPSC.Bkg.per.slide),
                                                                                          mean_ratioI_Total_SYPSC.Bkg.I_Total_SYP.Bkg=mean(ratioI_Total_SYPSC.Bkg.I_Total_SYP.Bkg),
                                                                                          mean_ratioI_Total_SYPSC.Bkg.I_Total_SYPNuc.Bkg=mean(ratioI_Total_SYPSC.Bkg.I_Total_SYPNuc.Bkg))

df.plot$GT <- factor(df.plot$GT, levels=c("ie29-10", "D114", "AN90"))

ggplot(df.plot) + geom_jitter(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),height=0,width=.02, alpha=0.4, size=1) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='line', linewidth=0.5, alpha=0.8) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='errorbar', width=0.05) +
  labs(x="Normalised gonad length", y="Normalised sum intensity of HA::SYP-1 on SC to mean(sum intensity wt)") + theme_cowplot() +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) 

ggplot(df.plot) + geom_jitter(aes(x=GonadLengthpxlNorm_bin, y=mean_ratioI_Total_SYPSC.Bkg.I_Total_SYP.Bkg , color=GT),height=0,width=.02, alpha=0.4, size=1) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_ratioI_Total_SYPSC.Bkg.I_Total_SYP.Bkg , color=GT),geom='line', linewidth=0.5, alpha=0.8) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_ratioI_Total_SYPSC.Bkg.I_Total_SYP.Bkg , color=GT),geom='errorbar', width=0.05) +
  labs(x="Normalised gonad length", y="Ratio of SYP-4::HA sum intensity on SC vs Nucleus") + theme_cowplot() +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) 

ggplot(df.plot) + geom_jitter(aes(x=GonadLengthpxlNorm_bin, y=mean_ratioI_Total_SYPSC.Bkg.I_Total_SYPNuc.Bkg , color=GT),height=0,width=.02, alpha=0.4, size=1) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_ratioI_Total_SYPSC.Bkg.I_Total_SYPNuc.Bkg , color=GT),geom='line', linewidth=0.5, alpha=0.8) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_ratioI_Total_SYPSC.Bkg.I_Total_SYPNuc.Bkg , color=GT),geom='errorbar', width=0.05) +
  labs(x="Normalised gonad length", y="Ratio of SYP-4::HA sum intensity on SC vs Nucleoplasm") + theme_cowplot() +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) 



final <- ggplot(df.plot) + geom_jitter(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),height=0,width=.02, alpha=0.4, size=1) +
          stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='line', linewidth=0.5, alpha=0.8) +
          stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='errorbar', width=0.05) +
          labs(x="Normalised gonad length", y="Normalised sum intensity of SYP-4::HA on SC") + theme_cowplot() +
          theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) + ylim(0,8) +
          scale_color_manual(values=c(WT_Color, D114_Color, AN90_Color)) +
          force_panelsizes(rows = unit(1.6, "in"), cols = unit(1.6, "in"))
  
  
pdf(file=paste(Sys.Date(),"_D114_9FA_SCloading_MEANSD_errorbar_v4.pdf",sep = ""),useDingbats = TRUE) #in inches
final
dev.off()


final <- ggplot(df.plot[df.plot$GT%in%c('ie29-10', 'D114'),]) + geom_jitter(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),height=0,width=.02, alpha=0.4, size=1) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='line', linewidth=0.5, alpha=0.8) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='errorbar', width=0.05) +
  labs(x="Normalised gonad length", y="Normalised sum intensity of SYP-4::HA on SC") + theme_cowplot() +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) + ylim(0,8) +
  scale_color_manual(values=c(WT_Color, D114_Color)) +
  force_panelsizes(rows = unit(1.6, "in"), cols = unit(1.6, "in"))


pdf(file=paste(Sys.Date(),"_D114_SCloading_MEANSD_errorbar_v4.pdf",sep = ""),useDingbats = TRUE) #in inches
final
dev.off()


final <- ggplot(df.plot[df.plot$GT%in%c('ie29-10', 'AN90'),]) + geom_jitter(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),height=0,width=.02, alpha=0.4, size=1) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='line', linewidth=0.5, alpha=0.8) +
  stat_summary(aes(x=GonadLengthpxlNorm_bin, y=mean_I_Total_SYPSC.Bkg.per.slide , color=GT),geom='errorbar', width=0.05) +
  labs(x="Normalised gonad length", y="Normalised sum intensity of SYP-4::HA on SC") + theme_cowplot() +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) + ylim(0,8) +
  scale_color_manual(values=c(WT_Color, AN90_Color)) +
  force_panelsizes(rows = unit(1.6, "in"), cols = unit(1.6, "in"))


pdf(file=paste(Sys.Date(),"_9FA_SCloading_MEANSD_errorbar_v4.pdf",sep = ""),useDingbats = TRUE) #in inches
final
dev.off()
