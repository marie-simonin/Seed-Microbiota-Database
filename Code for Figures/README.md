---
title: "Meta-Analysis Seed Microbiome - Sccript to generate FIGURES"
author: "Marie Simonin"
date: "18/01/2021"
output: html_document
---

# Figure 1 Map colored by number of samples

## load map info
```{r}
map <- read.table("Fig1_map.txt", header=TRUE, check.names = FALSE, sep="\t")

head(map)
dim(map)
```

## Make final map - Figure 1A
```{r}
library(ggplot2)
library(dplyr)
require(maps)
require(viridis)
world_map <- map_data("world")
worldSubset <- inner_join(world_map, map, by = "region")
head(worldSubset)

### Define countries outlines as a geom_polygon layer
country.layer <- geom_polygon(aes(x = long, y = lat, group = group),
                             data = world_map, fill = NA, color = "black", size=0.3)


## Let's ditch many of the unnecessary elements
plain <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  plot.title = element_text(hjust = 0.5)
)

Fig1_map <-ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = Sample_number))+ labs(fill = "Number of\nSeed Samples")+ scale_fill_gradient2(trans = "log",  breaks = c(1,10,100,500,2000), midpoint=3, low="white", mid="#ffeda0", high="#f03b20")	+ theme(legend.text = element_text(color="black", size=9, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold")) + scale_y_continuous(limits=c(-55,80))+ scale_x_continuous(limits=c(-170,195)) +
  plain+country.layer


Fig1_map
```


# Figure 2 ASV richness across all plants - based on Subset 2

## 16S V4 + V5/V6 datasets
```{r}
Div16S <- read.table("Metadata_16S_V4_V5V6_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep="\t")

head(Div16S)
dim(Div16S)
```
## gyrB dataset
```{r}
DivgyrB <- read.table("Metadata_gyrB_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep="\t")

head(DivgyrB)
dim(DivgyrB)
```

## ITS dataset
```{r}
DivITS <- read.table("Metadata_ITS1_ITS2_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep="\t")

head(DivITS)
dim(DivITS)
```

## 16S -Average SV richness all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(Div16S, measurevar="observed_otus", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
Div16S <- Div16S %>% group_by(Plant) %>%mutate(median_observed_otus_plant = median(observed_otus))

```


## Figure 2A- 16S Ordered based on Median plant richness
```{r warning=FALSE}
library(ggplot2)
Div16S$Plant<-ordered(Div16S$Plant, levels=c("Cauliflower","Tomato","Broccoli","Tobacco","Cabbage","Eyebright","Turnip","Willow Gentian","Lolium arundinacea","Phelipanche ramosa", "Garden rocket","Great Masterwort","Pincushion Flower","Festuca rubra","Setaria pumila","Wheat","Setaria viridis","Rapeseed","Capsella bursa-pastoris","Grass of Parnassus", "Lolium perenne","Cardamine hirsuta","Melon","Brassica nigra","Erophila verna ", "Alliaria petiolata","Arabidopsis thaliana","Pea", "Heliosperma alpestre","Rorippa sylvestris","Siberian wildrye", "Carrot","Barbarea vulgaris ","Rice","Rhinanthus glacialis","Medicago truncatula","Oat","Dahurian wildrye", "Sinapis arvensis","Radish","Bean","Hairy vetch", "Sunflower","Chiltern Gentian","Berteroa incana", "Alfalfa", "Oak"))
color= c("#5892ae","#ff2600","#aea9ca","#ff9b7c","#811770","#00bb9f","#c5d6d6","#15cd7e","#3c884c","#08ffda","#67bde8","#f0810F","#e09aa5","#2E4600","#c5cd77","#fdf351","#edf050","#008bc4","#6b66a3","#d0425d","#72ae4f","#1b2b7d","#FA6775","#a96699","#5e87c5","#4a1777","#835e9f","#68104d","#ffceda","#3e8a89","#9ff76e","#ee8332","#baa5c8","#b2d24f","#8eddad","#b554a6","#66ff00","#2fbd03","#9abdbb","#4cb5f5","#e03581","#e035ac","#fec767","#00a7b5","#ceaac4","#e035d7","#000000")
p4=ggplot(data=Div16S, aes(x=observed_otus, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Observed ASV richness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(3, 1,2, 16))+geom_vline(aes(xintercept = median_observed_otus_plant, group = Plant, colour = Plant))+ labs(shape = "Gene Region", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_x_log10()+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.9, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"))+ggtitle("16S rRNA gene - Bacteria & Archaea") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+ theme(panel.spacing.y = unit(0.2, "lines"))+ guides(col = guide_legend(ncol = 4))
p4
```


## gyrB Average SV richness all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(DivgyrB, measurevar="observed_otus", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
DivgyrB <- DivgyrB %>% group_by(Plant) %>%mutate(median_observed_otus_plant = median(observed_otus))
```


## Figure 2B gyrB - Ordered based on median plant richness
```{r warning=FALSE}
library(ggplot2)
DivgyrB$Plant<-ordered(DivgyrB$Plant, levels=c("Cauliflower","Melon","Broccoli","Tomato","Turnip","Medicago truncatula","Rapeseed","Carrot","Cabbage","Garden rocket","Brassica nigra","Radish","Cardamine hirsuta","Erophila verna ","Capsella bursa-pastoris","Alliaria petiolata","Barbarea vulgaris ","Berteroa incana","Arabidopsis thaliana","Bean","Sinapis arvensis"))
color= c("#5892ae","#FA6775","#aea9ca","#ff2600","#c5d6d6","#b554a6","#008bc4","#ee8332","#811770","#67bde8","#a96699","#4cb5f5","#1b2b7d","#5e87c5","#6b66a3","#4a1777","#baa5c8","#ceaac4","#835e9f","#e03581","#9abdbb")
p5=ggplot(data=DivgyrB, aes(x=observed_otus, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Observed ASV richness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold"))+scale_shape_manual(values=c(1, 16))+geom_vline(aes(xintercept = median_observed_otus_plant, group = Plant, colour = Plant))+ labs(shape = "Gene Region", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_x_log10()+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"),strip.text = element_text(colour = "black", face = "bold"))+ theme(legend.position = "none")+ggtitle("gyrB gene - Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))
p5
```


## ITS Average SV richness all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(DivITS, measurevar="observed_otus", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
DivITS <- DivITS %>% group_by(Plant) %>%mutate(median_observed_otus_plant = median(observed_otus))
#write.table(DivITS, "DivITS.txt")
```


## Figure 2C - ITS Ordered based on median plant richness
```{r warning=FALSE}
library(ggplot2)
DivITS$Plant<-ordered(DivITS$Plant, levels=c("Medicago truncatula","Great Masterwort","Garden rocket","Broccoli","Cabbage","Turnip","Creeping Bentgrass", "Eyebright","Pincushion Flower","Tomato","Rapeseed","Bean","Grass of Parnassus","Phelipanche ramosa", "Willow Gentian","Cauliflower","Heliosperma alpestre","Radish","Carrot","Pea", "White Clover", "Chiltern Gentian","Rice","Rorripa sylvestris", "Barbarea vulgaris ","Sinapis arvensis","Wheat","Alliaria petiolata","Cardamine hirsuta","Brassica nigra","Dahurian wildrye","Alfalfa", "Hairy vetch", "Red Clover","Oak",  "Rhinanthus glacialis","Erophila verna ","Siberian wildrye", "Oat", "Arabidopsis thaliana","Berteroa incana","Capsella bursa-pastoris","Sunflower"))
color= c("#b554a6","#f0810F","#67bde8","#aea9ca","#811770","#c5d6d6","#bd9103","#00bb9f","#e09aa5","#ff2600","#008bc4","#e03581","#d0425d","#08ffda","#15cd7e","#5892ae","#ffceda","#4cb5f5","#ee8332","#68104d","#8f868c","#00a7b5","#b2d24f","#c473fa","#baa5c8","#9abdbb","#fdf351","#4a1777","#1b2b7d","#a96699","#2fbd03","#e035d7","#e035ac","#dad7d9","#000000","#8eddad","#5e87c5","#9ff76e","#66ff00","#835e9f","#ceaac4","#6b66a3","#fec767")
p6=ggplot(data=DivITS, aes(x=observed_otus, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Observed ASV richness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(1, 16))+geom_vline(aes(xintercept = median_observed_otus_plant, group = Plant, colour = Plant))+ labs(shape = "Gene Region", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_x_log10()+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.9, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"))+ theme(legend.position = "none")+ggtitle("ITS region - Fungi") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+ theme(panel.spacing.y = unit(0.2, "lines"))
p6
```

## Make combined Figure 2 plot with the 3 genes
```{r warning=FALSE}
#http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
#plot: the plot to place (ggplot2 or a gtable)
#x: The x location of the lower left corner of the plot.
#y: The y location of the lower left corner of the plot.
#width, height: the width and the height of the plot
library(cowplot)
legend <- get_legend(p4)
p4=p4+theme(legend.position = "none")
```

```{r warning=FALSE}
library(gridExtra)
figure2_Median=ggdraw() +
  draw_plot(p4, 0, 0, 0.35, 1) +
  draw_plot(p5, 0.35, 0.25, .31, .75) +
  draw_plot(p6, .65, 0.1, .31, .9) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.4, 0.7), c(1, 1, 1), size = 15)
figure2_Median
```


###############################################################################

# Supplementary Figure S4 -Shannon index all genes

## 16S-Average Shannon all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(Div16S, measurevar="shannon", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
Div16S <- Div16S %>% group_by(Plant) %>%mutate(median_shannon_plant = median(shannon))
#write.table(Div16S, "Div16S.txt")
```


## Figure S4A - 16S-Ordered based on Median plant richness
```{r warning=FALSE}
library(ggplot2)
Div16S$Plant<-ordered(Div16S$Plant, levels=c("Cauliflower","Tomato","Broccoli","Phelipanche ramosa","Tobacco","Lolium arundinacea","Eyebright","Cabbage","Pea","Wheat","Festuca rubra","Siberian wildrye","Pincushion Flower","Turnip","Dahurian wildrye","Setaria pumila", "Willow Gentian","Oat","Setaria viridis","Lolium perenne","Brassica nigra","Hairy vetch","Great Masterwort", "Rice","Grass of Parnassus","Alfalfa","Melon","Cardamine hirsuta","Alliaria petiolata","Carrot","Radish","Sunflower","Garden rocket","Heliosperma alpestre","Rapeseed","Capsella bursa-pastoris","Erophila verna ","Sinapis arvensis","Arabidopsis thaliana", "Bean",   "Rorippa sylvestris","Rhinanthus glacialis","Barbarea vulgaris ","Chiltern Gentian","Berteroa incana","Oak","Medicago truncatula"))
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
p4=ggplot(data=Div16S, aes(x=shannon, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Shannon Index")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=9, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold"))+scale_shape_manual(values=c(3, 1,2, 16))+geom_vline(aes(xintercept = median_shannon_plant, group = Plant, colour = Plant))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=9, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"),strip.text = element_text(colour = "black", face = "bold"))+ggtitle("16S rRNA gene - Bacteria & Archaea") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ theme(panel.spacing.y = unit(0.2, "lines"))+scale_x_continuous(limits = c(0,10))+ guides(col = guide_legend(ncol = 4))
p4


```


## GyrB- Average shannon all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(DivgyrB, measurevar="shannon", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
DivgyrB <- DivgyrB %>% group_by(Plant) %>%mutate(median_shannon_plant = median(shannon))
#write.table(DivgyrB, "DivgyrB.txt")
```



## Figure S4B- gyrB Ordered based on median plant richness
```{r warning=FALSE}
library(ggplot2)
DivgyrB$Plant<-ordered(DivgyrB$Plant, levels=c("Cauliflower","Broccoli","Melon","Turnip","Tomato","Carrot","Rapeseed","Cabbage","Cardamine hirsuta","Brassica nigra","Medicago truncatula","Erophila verna ","Sinapis arvensis","Capsella bursa-pastoris","Arabidopsis thaliana","Radish","Alliaria petiolata","Berteroa incana","Barbarea vulgaris ","Garden rocket","Bean"))
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
p5=ggplot(data=DivgyrB, aes(x=shannon, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Shannon Index")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold"))+scale_shape_manual(values=c(1, 16))+geom_vline(aes(xintercept = median_shannon_plant, group = Plant, colour = Plant))+ labs(shape = "Gene Region", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"),strip.text = element_text(colour = "black", face = "bold"))+ theme(legend.position = "none")+ggtitle("gyrB gene - Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+scale_x_continuous(limits = c(0,10))
p5
```


## ITS Average Shannon all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(DivITS, measurevar="shannon", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
DivITS <- DivITS %>% group_by(Plant) %>%mutate(median_shannon_plant = median(shannon))
```



## Figure S4C - ITS Ordered based on median plant richness
```{r warning=FALSE}
library(ggplot2)
DivITS$Plant<-ordered(DivITS$Plant, levels=c('Great Masterwort','Broccoli','Garden rocket','Medicago truncatula','Pincushion Flower','Tomato','Turnip','Cabbage','Pea','Phelipanche ramosa','Creeping Bentgrass','Grass of Parnassus','Heliosperma alpestre','Rapeseed','Willow Gentian','Eyebright','Radish','Bean','Cardamine hirsuta','Rice','White Clover','Rorripa sylvestris','Alliaria petiolata','Chiltern Gentian','Brassica nigra','Carrot','Sinapis arvensis','Capsella bursa-pastoris','Arabidopsis thaliana','Barbarea vulgaris ','Red Clover','Siberian wildrye','Dahurian wildrye','Wheat','Hairy vetch','Erophila verna ','Oak','Alfalfa','Oat','Cauliflower','Rhinanthus glacialis','Sunflower','Berteroa incana'))
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
p6=ggplot(data=DivITS, aes(x=shannon, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Shannon Index")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold"))+scale_shape_manual(values=c(1, 16))+geom_vline(aes(xintercept = median_shannon_plant, group = Plant, colour = Plant))+ labs(shape = "Gene Region", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"),strip.text = element_text(colour = "black", face = "bold"))+ theme(legend.position = "none")+ggtitle("ITS region - Fungi") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ theme(panel.spacing.y = unit(0.2, "lines"))+scale_x_continuous(limits = c(0,10))
p6
```



## Make combined Figure S4 plot with the 3 genes
```{r warning=FALSE}
#http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
#plot: the plot to place (ggplot2 or a gtable)
#x: The x location of the lower left corner of the plot.
#y: The y location of the lower left corner of the plot.
#width, height: the width and the height of the plot
library(cowplot)
legend <- get_legend(p4)
p4=p4+theme(legend.position = "none")
```

```{r warning=FALSE}
library(gridExtra)
library(gridExtra)
figureS4_MedianShannon=ggdraw() +
  draw_plot(p4, 0, 0, 0.35, 1) +
  draw_plot(p5, 0.35, 0.25, .31, .75) +
  draw_plot(p6, .65, 0.1, .31, .9) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.4, 0.7), c(1, 1, 1), size = 15)
figureS4_MedianShannon
```


###############################################################

###############################################################################

# Supplementary Figure S5 -Evenness all genes

## 16S-Average pielou_e all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(Div16S, measurevar="pielou_e", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
Div16S <- Div16S %>% group_by(Plant) %>%mutate(median_pielou_e_plant = median(pielou_e))
#write.table(Div16S, "Div16S.txt")
```


## Figure S5A - 16S-Ordered based on Median plant richness
```{r warning=FALSE}
library(ggplot2)
Div16S$Plant<-ordered(Div16S$Plant, levels=c('Phelipanche ramosa',
'Cauliflower','Pea','Siberian wildrye','Tomato','Dahurian wildrye','Lolium arundinacea','Broccoli','Wheat','Hairy vetch','Alfalfa','Oat','Rice','Tobacco','Cabbage','Festuca rubra','Setaria pumila','Pincushion Flower','Eyebright','Lolium perenne','Cardamine hirsuta','Brassica nigra','Willow Gentian','Setaria viridis','Sunflower','Turnip','Melon','Radish','Grass of Parnassus','Bean','Great Masterwort','Carrot','Alliaria petiolata','Sinapis arvensis','Heliosperma alpestre','Rapeseed','Erophila verna ','Arabidopsis thaliana','Rorippa sylvestris','Capsella bursa-pastoris','Garden rocket','Oak','Chiltern Gentian','Rhinanthus glacialis','Barbarea vulgaris ','Berteroa incana','Medicago truncatula'))
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 

p4=ggplot(data=Div16S, aes(x=pielou_e, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Pielou's Evenness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(3, 1,2, 16))+geom_vline(aes(xintercept = median_pielou_e_plant, group = Plant, colour = Plant))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.9, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"))+
  theme(legend.position = "none")+ggtitle("16S rRNA gene - Bacteria & Archaea") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+ theme(panel.spacing.y = unit(0.2, "lines"))
p4
```


## gyrB Average pielou_e all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(DivgyrB, measurevar="pielou_e", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
DivgyrB <- DivgyrB %>% group_by(Plant) %>%mutate(median_pielou_e_plant = median(pielou_e))
#write.table(DivgyrB, "DivgyrB.txt")
```


## Figure S5B - gyrB -Ordered based on median plant richness
```{r warning=FALSE}
library(ggplot2)
DivgyrB$Plant<-ordered(DivgyrB$Plant, levels=c('Cardamine hirsuta','Cauliflower','Broccoli','Carrot','Alliaria petiolata','Brassica nigra','Tomato','Cabbage','Arabidopsis thaliana','Erophila verna ','Melon','Turnip','Capsella bursa-pastoris','Sinapis arvensis','Rapeseed','Berteroa incana','Barbarea vulgaris ','Medicago truncatula','Radish','Bean','Garden rocket'))
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 

p5=ggplot(data=DivgyrB, aes(x=pielou_e, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Pielou's Evenness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold"))+scale_shape_manual(values=c(1, 16))+geom_vline(aes(xintercept = median_pielou_e_plant, group = Plant, colour = Plant))+ labs(shape = "Gene Region", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"),strip.text = element_text(colour = "black", face = "bold"))+ theme(legend.position = "none")+ggtitle("gyrB gene - Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))
p5
```

## ITS-Average pielou_e all studies
```{r warning=FALSE}
library(Rmisc)
Diversity_stat <- summarySE(DivITS, measurevar="pielou_e", groupvars=c("Plant"), na.rm = TRUE)
library(dplyr) # Calculate median by plant
DivITS <- DivITS %>% group_by(Plant) %>%mutate(median_pielou_e_plant = median(pielou_e))
#write.table(DivITS, "DivITS.txt")
```

## Figure S5C - ITS Ordered based on median plant richness
```{r warning=FALSE}
library(ggplot2)
DivITS$Plant<-ordered(DivITS$Plant, levels=c('Cardamine hirsuta','Capsella bursa-pastoris','Brassica nigra','Tomato','Pea','Great Masterwort','Alliaria petiolata','Erophila verna ','Arabidopsis thaliana','Rorripa sylvestris','Phelipanche ramosa','Pincushion Flower','Grass of Parnassus','White Clover','Rice','Broccoli','Radish','Heliosperma alpestre','Rapeseed','Creeping Bentgrass','Willow Gentian','Garden rocket','Medicago truncatula','Turnip','Eyebright','Cabbage','Bean','Barbarea vulgaris ','Siberian wildrye','Sinapis arvensis','Chiltern Gentian','Red Clover','Oak','Hairy vetch','Carrot','Dahurian wildrye','Alfalfa','Oat','Wheat','Rhinanthus glacialis','Sunflower','Cauliflower','Berteroa incana'))
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 

p6=ggplot(data=DivITS, aes(x=pielou_e, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Pielou's Evenness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold"))+scale_shape_manual(values=c(1, 16))+geom_vline(aes(xintercept = median_pielou_e_plant, group = Plant, colour = Plant))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"),strip.text = element_text(colour = "black", face = "bold"))+ theme(legend.position = "none")+ggtitle("ITS region - Fungi") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+ theme(panel.spacing.y = unit(0.2, "lines"))
p6
```



## Make combined Figure S5 plot with the 3 genes
```{r warning=FALSE}
#http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
#plot: the plot to place (ggplot2 or a gtable)
#x: The x location of the lower left corner of the plot.
#y: The y location of the lower left corner of the plot.
#width, height: the width and the height of the plot
library(cowplot)
legend <- get_legend(p4)
p4=p4+theme(legend.position = "none")
```

```{r warning=FALSE}
library(gridExtra)
figureS5_MedianEvenness=ggdraw() +
  draw_plot(p4, 0, 0, 0.35, 1) +
  draw_plot(p5, 0.35, 0.25, .31, .75) +
  draw_plot(p6, .65, 0.1, .31, .9) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.4, 0.7), c(1, 1, 1), size = 15)
figureS5_MedianEvenness
```



############################################################################################################
#########################################################################################################

# Figure S6 Effect of seed fraction on ASV richness
### Calculate median by Seed fraction
```{r warning=FALSE}
library(dplyr) # Calculate median by plant
Div16S <- Div16S %>% group_by(Seed_fraction) %>%mutate(median_richness_seed_fraction = median(observed_otus))
#write.table(DivITS, "DivITS.txt")
```

### Seed Fraction -16S Richness
```{r warning=FALSE}
F1=ggplot(data=Div16S) + geom_jitter(alpha=0.8,  aes(x=observed_otus, y=Seed_fraction, color=Plant, shape=Seed_fraction)) +xlab("Observed ASVs Richness")+ylab("Seed Fraction")+ theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(3, 1,2, 16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+
  theme(legend.position = "none")+ggtitle("16S rRNA gene - Bacteria & Archaea") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+
  geom_vline(aes(xintercept = median_richness_seed_fraction, group = Seed_fraction, colour = "black"))+scale_x_log10()+geom_boxplot(data=Div16S, aes(x=observed_otus, y=Seed_fraction), alpha=0.3)
F1
```

### Seed Fraction -gyrB Richness
```{r warning=FALSE}
F2=ggplot(data=DivgyrB) + geom_jitter(alpha=0.8,  aes(x=observed_otus, y=Seed_fraction, color=Plant, shape=Seed_fraction)) +xlab("Observed ASVs Richness")+ylab("Seed Fraction")+ theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(1,16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+
  theme(legend.position = "none")+ggtitle("gyrB gene - Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+scale_x_log10()+geom_boxplot(data=DivgyrB, aes(x=observed_otus, y=Seed_fraction), alpha=0.3)
F2
```


### Seed Fraction -ITS - Richness
```{r warning=FALSE}
F3=ggplot(data=DivITS) + geom_jitter(alpha=0.8,  aes(x=observed_otus, y=Seed_fraction, color=Plant, shape=Seed_fraction)) +xlab("Observed ASVs Richness")+ylab("Seed Fraction")+ theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(1,16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+
  theme(legend.position = "none")+ggtitle("ITS Region - Fungi") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+scale_x_log10()+geom_boxplot(data=DivITS, aes(x=observed_otus, y=Seed_fraction), alpha=0.3)
F3
```
## Figure S6 with 3 plots combined
```{r}
library(ggpubr)
figureS6=ggarrange(F1, F2, F3,labels = c("A", "B", "C"), ncol=2, nrow = 2)
figureS6
```




############################################################################################################
#########################################################################################################
# Figure S7 - Effect of Seed Preparation (grinding - soaking) on ASV richness
### Seed Prep -16S Richness
```{r warning=FALSE}
P1=ggplot(data=Div16S) + geom_jitter(alpha=0.8,  aes(x=observed_otus, y=Microbe_collection, color=Plant, shape=Seed_fraction)) +xlab("Observed ASVs Richness")+ylab("Seed Preparation")+ theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(3, 1,2, 16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+
  theme(legend.position = "none")+ggtitle("16S rRNA gene - Bacteria & Archaea") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+
  geom_vline(aes(xintercept = median_richness_seed_fraction, group = Seed_fraction, colour = "black"))+scale_x_log10()+geom_boxplot(data=Div16S, aes(x=observed_otus, y=Microbe_collection), alpha=0.3)
P1
```

### Seed Prep -gyrB Richness
```{r warning=FALSE}
P2=ggplot(data=DivgyrB) + geom_jitter(alpha=0.8,  aes(x=observed_otus, y=Microbe_collection, color=Plant, shape=Seed_fraction)) +xlab("Observed ASVs Richness")+ylab("Seed Fraction")+ theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(1,16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+
  theme(legend.position = "none")+ggtitle("gyrB gene - Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+scale_x_log10()+geom_boxplot(data=DivgyrB, aes(x=observed_otus, y=Microbe_collection), alpha=0.3)
P2
```


### Seed Prep -ITS - Richness
```{r warning=FALSE}
P3=ggplot(data=DivITS) + geom_jitter(alpha=0.8,  aes(x=observed_otus, y=Microbe_collection, color=Plant, shape=Seed_fraction)) +xlab("Observed ASVs Richness")+ylab("Seed Fraction")+ theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(1,16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+
  theme(legend.position = "none")+ggtitle("ITS Region - Fungi") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+scale_x_log10()+geom_boxplot(data=DivITS, aes(x=observed_otus, y=Microbe_collection), alpha=0.3)
P3
```
## Make Figure S7
```{r}
library(ggpubr)
figureS7=ggarrange(P1, P2, P3,labels = c("A", "B", "C"), ncol=2, nrow = 2)
figureS7
```



############################################################################################################
#########################################################################################################
# Figure S1 -16S Basic graphs to describe the meta-analysis dataset (Subset 2)
### 16S Seed Fraction
```{r warning=FALSE}
A1=Div16S %>% 
	group_by(Seed_fraction) %>% 
	summarise(count = n()) %>% 
	ggplot(aes(x = reorder(Seed_fraction,(-count)), y = count)) + 
		geom_bar(stat = 'identity', fill="#bcd4e6") +xlab("Seed Fraction")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
A1
```
### Microbe Collection
```{r warning=FALSE}
A2=Div16S %>% 
	group_by(Microbe_collection) %>% 
	summarise(count = n()) %>% 
	ggplot(aes(x = reorder(Microbe_collection,(-count)), y = count)) + 
		geom_bar(stat = 'identity', fill="#bcd4e6") +xlab("Seed Preparation")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))
A2
```
### Plant Species
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600","Asteraceae"="grey", "Cucurbitaceae"="darkred") 
library(forcats)
A3=ggplot(Div16S,aes(x = fct_infreq(Plant), fill=Family)) + 
	geom_bar(stat = 'count')+scale_fill_manual(values =color)+xlab("Plant Species")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+guides(fill=guide_legend(nrow=5,byrow=TRUE))+ theme(legend.position = c(0.65, 0.55))	
A3
```


### Plant Family
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600","Asteraceae"="grey", "Cucurbitaceae"="darkred") 
library(forcats)
A4=ggplot(Div16S,aes(x = fct_infreq(Family), fill=Family)) + 
	geom_bar(stat = 'count')+scale_fill_manual(values =color)+xlab("Plant Family")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A4
```

### Countries
```{r warning=FALSE}
library(forcats)
A5=ggplot(Div16S,aes(x = fct_infreq(Country))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Country of Origin of the Seed Samples")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A5
```

### Research Institutes
```{r warning=FALSE}
library(forcats)
A6=ggplot(Div16S,aes(x = fct_infreq(Study_Origin))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Research Institutes - Universities")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 65, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A6
```

### Site
```{r warning=FALSE}
A7=ggplot(Div16S,aes(x = fct_infreq(Site))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Production Condition")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A7
```

### Seed Number
```{r warning=FALSE}
A8=ggplot(Div16S, aes(Seed_number))+geom_histogram(binwidth=5, fill="#bcd4e6") +xlab("Number of Seeds in Sample")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
A8
```

# Figure S2 - gyrB
### Seed Fraction
```{r warning=FALSE}
A1=DivgyrB %>% 
	group_by(Seed_fraction) %>% 
	summarise(count = n()) %>% 
	ggplot(aes(x = reorder(Seed_fraction,(-count)), y = count)) + 
		geom_bar(stat = 'identity', fill="#bcd4e6") +xlab("Seed Fraction")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))
A1
```
### Microbe Collection
```{r warning=FALSE}
A2=DivgyrB %>% 
	group_by(Microbe_collection) %>% 
	summarise(count = n()) %>% 
	ggplot(aes(x = reorder(Microbe_collection,(-count)), y = count)) + 
		geom_bar(stat = 'identity', fill="#bcd4e6") +xlab("Seed Preparation")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))
A2
```
### Plant Species
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600","Asteraceae"="grey", "Cucurbitaceae"="darkred") 
library(forcats)
A3=ggplot(DivgyrB,aes(x = fct_infreq(Plant), fill=Family)) + 
	geom_bar(stat = 'count')+scale_fill_manual(values =color)+xlab("Plant Species")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+guides(fill=guide_legend(nrow=5,byrow=TRUE))+ theme(legend.position = c(0.65, 0.55))	
A3
```


### Plant Family
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600","Asteraceae"="grey", "Cucurbitaceae"="darkred") 
library(forcats)
A4=ggplot(DivgyrB,aes(x = fct_infreq(Family), fill=Family)) + 
	geom_bar(stat = 'count')+scale_fill_manual(values =color)+xlab("Plant Family")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A4
```

### Countries
```{r warning=FALSE}
library(forcats)
A5=ggplot(DivgyrB,aes(x = fct_infreq(Country))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Country of Origin of the Seed Samples")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A5
```

### Research Institutes
```{r warning=FALSE}
library(forcats)
A6=ggplot(DivgyrB,aes(x = fct_infreq(Study_Origin))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Research Institutes - Universities")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 65, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A6
```

### Site
```{r warning=FALSE}
A7=ggplot(DivgyrB,aes(x = fct_infreq(Site))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Production Condition")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A7
```


### Seed Number
```{r warning=FALSE}
A8=ggplot(DivgyrB, aes(Seed_number))+geom_histogram(binwidth=5, fill="#bcd4e6") +xlab("Number of Seeds in Sample")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
A8
```

# Figure S3 - ITS
### Seed Fraction
```{r warning=FALSE}
A1=DivITS %>% 
	group_by(Seed_fraction) %>% 
	summarise(count = n()) %>% 
	ggplot(aes(x = reorder(Seed_fraction,(-count)), y = count)) + 
		geom_bar(stat = 'identity', fill="#bcd4e6") +xlab("Seed Fraction")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))
A1
```
### Microbe Collection
```{r warning=FALSE}
A2=DivITS %>% 
	group_by(Microbe_collection) %>% 
	summarise(count = n()) %>% 
	ggplot(aes(x = reorder(Microbe_collection,(-count)), y = count)) + 
		geom_bar(stat = 'identity', fill="#bcd4e6") +xlab("Seed Preparation")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))
A2
```
### Plant Species
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600","Asteraceae"="grey", "Cucurbitaceae"="darkred") 
library(forcats)
A3=ggplot(DivITS,aes(x = fct_infreq(Plant), fill=Family)) + 
	geom_bar(stat = 'count')+scale_fill_manual(values =color)+xlab("Plant Species")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+guides(fill=guide_legend(nrow=5,byrow=TRUE))+ theme(legend.position = c(0.65, 0.55))	
A3
```


### Plant Family
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600","Asteraceae"="grey", "Cucurbitaceae"="darkred") 
library(forcats)
A4=ggplot(DivITS,aes(x = fct_infreq(Family), fill=Family)) + 
	geom_bar(stat = 'count')+scale_fill_manual(values =color)+xlab("Plant Family")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A4
```

### Countries
```{r warning=FALSE}
library(forcats)
A5=ggplot(DivITS,aes(x = fct_infreq(Country))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Country of Origin of the Seed Samples")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A5
```

### Research Institutes
```{r warning=FALSE}
library(forcats)
A6=ggplot(DivITS,aes(x = fct_infreq(Study_Origin))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Research Institutes - Universities")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 65, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A6
```

### Site
```{r warning=FALSE}
A7=ggplot(DivITS,aes(x = fct_infreq(Site))) + 
	geom_bar(stat = 'count', fill="#bcd4e6")+xlab("Production Condition")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black", size=9, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+theme(legend.position = "none")
A7
```

### Seed Number
```{r warning=FALSE}
A8=ggplot(DivITS, aes(Seed_number))+geom_histogram(binwidth=5, fill="#bcd4e6") +xlab("Number of Seeds in Sample")+ylab("Number of Seed Samples")+ theme_classic()+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
A8
```


############################################################################################################
#########################################################################################################

# Figure 3 - Relationship between 16S and ITS ASV richness
```{r}
ITS_16S<-read.table("16S_and_ITS_merged_Subset2_metadata.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(ITS_16S)
dim(ITS_16S)
```

## Figure 3A - Relationship ASV Richness between 16S and ITS datasets
```{r warning=FALSE}
color=c("Apiaceae"="#000000", "Brassicaceae"="#81edff", "Caprifoliaceae"="#8338ec", "Caryophyllaceae"="#ff7c43", "Celastraceae"="#a05195", "Fabaceae"="#f69cbd", "Fagaceae"=	"#90665f", "Gentianaceae"="#6f9d4b", "Orobanchaceae"="#affc41", "Poaceae"=	"#ffdd00", "Solanaceae"="#ff2600") 

c1=ggplot(data=ITS_16S) + geom_point(aes(x=observed_otus_16S, y=observed_otus_ITS, color=Family, shape=Seed_fraction), alpha=0.7, size=2) + theme_classic()+xlab("Bacterial & Archaeal ASV Richness")+ylab("Fungal ASV Richness")+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.title = element_text(colour="black", size=12, face="bold"))+ geom_abline(intercept = 0, slope=1, aes(linetype="dashed", color="red"))+scale_x_log10()+scale_y_log10()+scale_shape_manual(values=c(1,16))+ labs(shape = "Seed Fraction", color = "Plant Family")+ guides(col = guide_legend(ncol = 4))+scale_color_manual(values=color)+guides(color=guide_legend(ncol=3,byrow=TRUE))
c1
```

## Figure 3B - Plot ratio 16S-ITS richness by Plant Species
```{r warning=FALSE}
# Calculate median by plant
ITS_16S <- ITS_16S %>% group_by(Plant) %>%mutate(median_ratio_observed_otus_plant = median(Ratio_observed_otus))
c5=ggplot(data=ITS_16S) + geom_jitter(aes(x=reorder(Plant, median_ratio_observed_otus_plant), y=Ratio_observed_otus, color=Family, shape=Seed_fraction), alpha=0.4, size=1.5) + theme_classic()+xlab("Plant Species")+ylab("Ratio Prokaryotic / Fungal\nObserved ASVs Richness")+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=10, face="bold"))+scale_shape_manual(values=c(1,16))+ labs(shape = "Seed Fraction", color = "Plant Family")+scale_color_manual(values=color)+scale_y_log10()+geom_hline(yintercept = 1)+ geom_point(aes(x=Plant, y=median_ratio_observed_otus_plant, color=Family), fill="black", size=3)+ theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black", size=9, face="bold"))+ theme(legend.position = "none")
c5
```

# Make Figure 3
```{r}
library(ggpubr)
figure3=ggarrange(c1, c5, nrow = 2, labels = c("A", "B"), heights = c(1,1.3)) 
figure3
```


############################################################################################################
#########################################################################################################

# Figure 4: Abundance Occupancy curves (and Figure S8 and S9)
##Calculations and graphs for 16S - V4 gene On Subset 3
```{r}
SV<-read.table("Subset3-16S-V4-table-FINAL-rarefied.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names=1)
head(SV)
dim(SV)
```



## Calculate prevalence and relative abundance for each ASV
```{r}
#Code Shade lab: https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/Core_prioritizing_script.R
library(tidyverse)
library(reshape2)
library(vegan)
#presence-absence data
SV_PA <- 1*((SV>0)==1)                                              
# occupancy calculation
SV_Prevalence <- rowSums(SV_PA)/ncol(SV_PA) 
# relative abundance  
SV_relative_abundance <- apply(decostand(SV, method="total", MARGIN=2),1, mean)     

# combining occupancy and relative abundance of each SV in a table
SVprev_rel <- add_rownames(as.data.frame(cbind(SV_Prevalence, SV_relative_abundance)),'SV') 
head(SVprev_rel)

```

## Merge prevalence/rel abund data with taxonomic info for each SV
```{r}
taxo<-read.table("Subset1-2_All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
SVprev_rel_taxo<-merge(SVprev_rel,taxo,by="SV")
dim(SVprev_rel_taxo)

head(SVprev_rel_taxo)

```


## Open Transposed 16S-V4 ASV table and metadata - Subset 3
```{r}
SV1 <- read.table("Subset3-16S-V4-table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
meta1<-read.table("Metadata_16S_V4_V5V6_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(SV1)
head(meta1)
dim(meta1)
dim(SV1)
SV_use1<-merge(meta1,SV1,by="SampleID")
dim(SV_use1)

head(SV_use1)
```

## Script for counting the number of Plants where each SV is observed
```{r}
#Keep only the Plant column as metadata for collapsing table by plant
SV_use2 <- SV_use1[ -c(1:9,11:38) ]
head (SV_use2)
```
```{r}
SV_use1_Plant=SV_use2 %>%
    group_by(Plant) %>% 
    summarise_each(funs(sum))
head(SV_use1_Plant)
dim(SV_use1_Plant)
###Make table as presence absence of SV for each plant species
SV_use1_Plant_PA <- 1*((SV_use1_Plant>0)==1) 
dim(SV_use1_Plant_PA)
```
## the number of plant species where each SV was observed - do colsums and transpose results
```{r}
SVobs_byPlant=data.frame(colSums(SV_use1_Plant_PA))
SVobs_byPlant=na.omit(SVobs_byPlant)
SVobs_byPlant <- tibble::rownames_to_column(SVobs_byPlant, "SV")
names(SVobs_byPlant)[names(SVobs_byPlant) == 'colSums.SV_use1_Plant_PA.'] <- 'NbPlantObs'
head(SVobs_byPlant)

##Merging SV info on prevalence, rel abund and taxo with the info on the number of plant species where each SV was observed
SVprev_rel_taxo_Plantobs<-merge(SVprev_rel_taxo,SVobs_byPlant,by="SV")
#write.table(SVprev_rel_taxo_Plantobs, file = "Subset3-16S-V4-SV_prevalence_abundance_taxo.txt", sep = "\t")
```



# Figure S8 - Plot Abundance-occupancy graph - 16S V4 - All phyla (for SI)
```{r}
g0=ggplot(data=SVprev_rel_taxo_Plantobs, aes(x=SV_relative_abundance, y=SV_Prevalence)) +
    geom_point(aes(color=NbPlantObs), alpha=0.7, size=0.7) + xlab("Log10(Mean ASV relative abundance)") + ylab("ASV Prevalence")+scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+ scale_y_log10(labels = scales::percent_format(accuracy = 1))+facet_wrap(~Phylum, ncol=7)+ theme_classic()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold")) + guides(color=guide_legend(title="Number of Plant\nSpecies Detected"))+scale_color_gradient2(midpoint=13, low="#c1e7ff", mid="#f2c057",high="red", breaks = c(1,5,10,15,20,25))+ theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour = "black", face = "bold"))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid"), panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed',colour = "black"))
g0
```

## Figure 4A - Plot Abundance-occupancy graph - 16S V4 - Selected phyla 
```{r}
SVprev_rel_taxo_Plantobs <- read.table("/Users/msimonin/OneDrive/INRA/Projets/Meta-Analyse Seed Microbiome/Figures_Stats/Final Figures 2021/16S/Subset3-16S-V4-SV_prevalence_abundance_taxo.txt", header=TRUE, check.names = FALSE, sep = "\t")
SVprev_rel_taxo_Plantobs_subset=subset(SVprev_rel_taxo_Plantobs, Most_abundant_phyla=="yes")
dim(SVprev_rel_taxo_Plantobs_subset)
g1=ggplot(data=SVprev_rel_taxo_Plantobs_subset, aes(x=SV_relative_abundance, y=SV_Prevalence)) +
    geom_point(aes(color=NbPlantObs), alpha=0.9, size=0.7) + xlab("Log10(Mean ASV relative abundance)") + ylab("ASV Prevalence")+scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+ scale_y_log10(labels = scales::percent_format(accuracy = 1))+facet_wrap(~Phylum, ncol=3)+ theme_classic()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold"))+ theme(legend.position = c(0.93, 0.2)) + guides(color=guide_legend(title="Number of Plant\nSpecies Detected"))+scale_color_gradient2(midpoint=14, low="#c1e7ff", mid="#f2c057",high="red", breaks = c(1,5,10,15,20,30))+ theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour = "black", face = "bold"))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid"), panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed',colour = "black"))+theme(legend.position = "none")+ggtitle("16S rRNA gene - Bacteria & Archaea")+theme(plot.title = element_text(hjust = 0.5, face="bold"))
g1
```


## Calculations and graphs for ITS1 dataset
```{r}
SV_ITS<-read.table("Subset3-ITS1_table-FINAL-rarefied.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names=1)
head(SV_ITS)
dim(SV_ITS)
```

## Calculate prevalence and relative abundance for each ASV
```{r}
library(dplyr)
#Code Shade lab: https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/Core_prioritizing_script.R
#presence-absence data
SV_ITS_PA <- 1*((SV_ITS>0)==1)                                              
# occupancy calculation
SV_ITS_Prevalence <- rowSums(SV_ITS_PA)/ncol(SV_ITS_PA) 
# relative abundance  
library(vegan)
SV_ITS_relative_abundance <- apply(decostand(SV_ITS, method="total", MARGIN=2),1, mean)     

# combining occupancy and relative abundance of each SV_ITS in a table
SV_ITSprev_rel <- add_rownames(as.data.frame(cbind(SV_ITS_Prevalence, SV_ITS_relative_abundance)),'SV') 
head(SV_ITSprev_rel)
dim(SV_ITSprev_rel)

```

## Merge prevalence/rel abund data with taxonomic info for each SV
```{r}
taxo_ITS<-read.table("Subset1-2_All_studies_merged_ITS1_taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
SV_ITSprev_rel_taxo<-merge(SV_ITSprev_rel,taxo_ITS,by="SV")
dim(SV_ITSprev_rel_taxo)

head(SV_ITSprev_rel_taxo)
#write.table(SV_ITSprev_rel_taxo, file = "Subset3-ITS1-SV_prevalence_abundance_taxo.txt", sep = "\t")
```


## Open Transposed ITS1 ASV table and metadata - Subset 3 
```{r}
meta2 <- read.table("Metadata_ITS1_ITS2_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
SV_ITS1<-read.table("Subset3-ITS1_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(SV_ITS1)
head(meta2)
dim(meta2)
dim(SV_ITS1)
SV_ITS_use1<-merge(meta2,SV_ITS1,by="SampleID")
dim(SV_ITS_use1)

head(SV_ITS_use1)
```

## Script for counting the number of Plants where each SV is observed
```{r}
#Keep only the Plant column as metadata for collapsing table by plant
SV_ITS_use2 <- SV_ITS_use1[ -c(1:9,11:33) ]
head (SV_ITS_use2)
```
```{r}
SV_ITS_use1_Plant=SV_ITS_use2 %>%
    group_by(Plant) %>% 
    summarise_each(funs(sum))
head(SV_ITS_use1_Plant)
dim(SV_ITS_use1_Plant)
###Make table as presence absence of SV for each plant species
SV_ITS_use1_Plant_PA <- 1*((SV_ITS_use1_Plant>0)==1) 
dim(SV_ITS_use1_Plant_PA)
```
## the number of plant species where each SV was observed - do colsums and transpose results
```{r}
SV_ITSobs_byPlant=data.frame(colSums(SV_ITS_use1_Plant_PA))
SV_ITSobs_byPlant=na.omit(SV_ITSobs_byPlant)
SV_ITSobs_byPlant <- tibble::rownames_to_column(SV_ITSobs_byPlant, "SV")
names(SV_ITSobs_byPlant)[names(SV_ITSobs_byPlant) == 'colSums.SV_ITS_use1_Plant_PA.'] <- 'NbPlantObs'
head(SV_ITSobs_byPlant)

##Merging SV info on prevalence, rel abund and taxo with the info on the number of plant species where each SV was observed
SV_ITSprev_rel_taxo_Plantobs<-merge(SV_ITSprev_rel_taxo,SV_ITSobs_byPlant,by="SV")
head(SV_ITSprev_rel_taxo_Plantobs)
#write.table(SV_ITSprev_rel_taxo_Plantobs, file = "Subset1-ITS1-SV_prevalence_abundance_taxo.txt", sep = "\t")
```


## Figure 4B -Plot Abundance-occupancy graph - ITS1 - All phyla 
```{r}
library(ggplot2)
g2=ggplot(data=SV_ITSprev_rel_taxo_Plantobs, aes(x=SV_ITS_relative_abundance, y=SV_ITS_Prevalence)) +
    geom_point(aes(color=NbPlantObs), alpha=0.9, size=0.7) + xlab("Log10(Mean ASV relative abundance)") + ylab("ASV Prevalence")+scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+ scale_y_log10(labels = scales::percent_format(accuracy = 1))+facet_wrap(~Phylum, ncol=4)+ theme_classic()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold")) +scale_color_gradient2(midpoint=16, low="#c1e7ff", mid="#f2c057",high="red", breaks = c(1,5,10,15,20,31))+ theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour = "black", face = "bold"))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid"), panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed',colour = "black"))+ theme(legend.position = c(0.9, 0.15)) + guides(colour = guide_legend(override.aes = list(size=2), title="Number of Plant\nSpecies Detected"))+ggtitle("ITS1 Region - Fungi")+theme(plot.title = element_text(hjust = 0.5, face="bold"))
g2
```



## Figure S9 Calculations and graphs for gyrB dataset - Subset 3
```{r}
SV_gyrB<-read.table("Subset3-gyrB-MiSeq_table-FINAL-rarefied.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names=1)
head(SV_gyrB)
dim(SV_gyrB)
```

## Calculate prevalence and relative abundance for each ASV
```{r}
#Code Shade lab: https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/Core_prioritizing_script.R
#presence-absence data
SV_gyrB_PA <- 1*((SV_gyrB>0)==1)                                              
# occupancy calculation
SV_gyrB_Prevalence <- rowSums(SV_gyrB_PA)/ncol(SV_gyrB_PA) 
# relative abundance  
library(vegan)
library(dplyr)
SV_gyrB_relative_abundance <- apply(decostand(SV_gyrB, method="total", MARGIN=2),1, mean)     

# combining occupancy and relative abundance of each SV_gyrB in a table
SV_gyrBprev_rel <- add_rownames(as.data.frame(cbind(SV_gyrB_Prevalence, SV_gyrB_relative_abundance)),'SV') 
head(SV_gyrBprev_rel)
dim(SV_gyrBprev_rel)

```

## Merge prevalence/rel abund data with taxonomic info for each SV
```{r}
taxo_gyrB<-read.table("Subset1-2_All_studies_merged_gyrB-rep-seqs-FINAL-filtered-taxonomy-final.tsv", header=TRUE, check.names = FALSE, sep = "\t")
SV_gyrBprev_rel_taxo<-merge(SV_gyrBprev_rel,taxo_gyrB,by="SV")
dim(SV_gyrBprev_rel_taxo)

head(SV_gyrBprev_rel_taxo)
#write.table(SV_gyrBprev_rel_taxo, file = "Subset1-gyrB1-SV_prevalence_abundance_taxo.txt", sep = "\t")
```


## Open Transposed gyrB1 ASV table and metadata - Subset 3 
```{r}
meta2 <- read.table("Metadata_gyrB_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
SV_gyrB1<-read.table("Subset3-gyrB-MiSeq_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(SV_gyrB1)
head(meta2)
dim(meta2)
dim(SV_gyrB1)
SV_gyrB_use1<-merge(meta2,SV_gyrB1,by="SampleID")
dim(SV_gyrB_use1)

head(SV_gyrB_use1)
```

## Script for counting the number of Plants where each SV is observed
```{r}
#Keep only the Plant column as metadata for collapsing table by plant
SV_gyrB_use2 <- SV_gyrB_use1[ -c(1:8,10:39) ]
head (SV_gyrB_use2)
```
```{r}
SV_gyrB_use1_Plant=SV_gyrB_use2 %>%
    group_by(Plant) %>% 
    summarise_each(funs(sum))
head(SV_gyrB_use1_Plant)
dim(SV_gyrB_use1_Plant)
###Make table as presence absence of SV for each plant species
SV_gyrB_use1_Plant_PA <- 1*((SV_gyrB_use1_Plant>0)==1) 
dim(SV_gyrB_use1_Plant_PA)
```
## the number of plant species where each SV was observed - do colsums and transpose results
```{r}
SV_gyrBobs_byPlant=data.frame(colSums(SV_gyrB_use1_Plant_PA))
SV_gyrBobs_byPlant=na.omit(SV_gyrBobs_byPlant)
SV_gyrBobs_byPlant <- tibble::rownames_to_column(SV_gyrBobs_byPlant, "SV")
names(SV_gyrBobs_byPlant)[names(SV_gyrBobs_byPlant) == 'colSums.SV_gyrB_use1_Plant_PA.'] <- 'NbPlantObs'
head(SV_gyrBobs_byPlant)

##Merging SV info on prevalence, rel abund and taxo with the info on the number of plant species where each SV was observed
SV_gyrBprev_rel_taxo_Plantobs<-merge(SV_gyrBprev_rel_taxo,SV_gyrBobs_byPlant,by="SV")
head(SV_gyrBprev_rel_taxo_Plantobs)
#write.table(SV_gyrBprev_rel_taxo_Plantobs, file = "Subset3-gyrB1-SV_prevalence_abundance_taxo_Plant.txt", sep = "\t")
```


## Figure S9 - Plot Abundance-occupancy graph - gyrB - All phyla
```{r}
g3=ggplot(data=SV_gyrBprev_rel_taxo_Plantobs, aes(x=SV_gyrB_relative_abundance, y=SV_gyrB_Prevalence)) +
    geom_point(aes(color=NbPlantObs), alpha=0.9, size=0.7) + xlab("Log10(Mean ASV relative abundance)") + ylab("ASV Prevalence")+scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+ scale_y_log10(labels = scales::percent_format(accuracy = 1))+facet_wrap(~Phylum, ncol=3)+ theme_classic()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=8, face="bold")) +scale_color_gradient2(midpoint=13, low="#c1e7ff", mid="#f2c057",high="red", breaks = c(1,5,10,15,20,25))+ theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour = "black", face = "bold"))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid"), panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed',colour = "black")) + guides(colour = guide_legend(override.aes = list(size=2), title="Number of Plant\nSpecies Detected"))+ggtitle("gyrB gene - Bacteria")+theme(plot.title = element_text(hjust = 0.5, face="bold"))
g3
```


## Figure 4C - Make small table of number of SVs observed in single or multiple species
### For 16S
```{r}
# set up cut-off values 
breaks <- c(1,2,5,10,20,27)
# specify interval/bin labels
tags <- c("1 species","2-4 species", "5-9 species", "10-19 species", "20-27 species")
# bucketing values into bins
Bac_bins <- cut(SVprev_rel_taxo_Plantobs$NbPlantObs, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(Bac_bins)
```

### For ITS
```{r}
# set up cut-off values 
breaks <- c(1,2,5,10,20,31)
# specify interval/bin labels
tags <- c("1 species","2-4 species", "5-9 species", "10-19 species", "20-31 species")
# bucketing values into bins
ITS_bins <- cut(SV_ITSprev_rel_taxo_Plantobs$NbPlantObs, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(ITS_bins)

```


## Make final Figure 4 with 16S and ITS results
```{r warning=FALSE}
library(gridExtra)
library(cowplot)


#plot: the plot to place (ggplot2 or a gtable)
#x: The x location of the lower left corner of the plot.
#y: The y location of the lower left corner of the plot.
#width, height: the width and the height of the plot
figure4=ggdraw() +
  draw_plot(g1, 0, 0, 0.45, 1) +
  draw_plot(g2, 0.45, 0.2, .55, .8) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), size = 15)
figure4
```

######################################################################################################
#########################################################################################

# Figure 5 
## For 16S V4- Import SV table with only the most prevalent bacterial taxa (>20 species)
```{r}
meta2 <- read.table("Metadata_16S_V4_V5V6_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
Top_Bac_subset3<-read.table("Subset3-16S-V4-table-FINAL-rarefied-13mostprevalent.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Top_Bac_subset3)
head(meta2)
dim(meta2)
dim(Top_Bac_subset3)

```

### Transform matrix into long table format
```{r}
library(reshape2)
Top_Bac_subset3_long=setNames(melt(Top_Bac_subset3), c('SampleID', 'Taxon', 'Relative_Abundance'))
head(Top_Bac_subset3_long)
dim(Top_Bac_subset3_long)
```
```{r}
##Merging long format table with metadata
Top_Bac_subset3_use<-merge(meta2,Top_Bac_subset3_long,by="SampleID")
dim(Top_Bac_subset3_use)

head(Top_Bac_subset3_use)
```

## Merge long format table with taxonomic info for each SV
```{r}
taxo<-read.table("Subset1-2_All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
Top_Bac_subset3_use_taxo<-merge(Top_Bac_subset3_use,taxo,by="Taxon")
dim(Top_Bac_subset3_use_taxo)

head(Top_Bac_subset3_use_taxo)

```
## Figure 5A - Plotting 16S most prevalent taxa (n=13) with ggridge
```{r}
Top_Bac_subset3_use_taxo$Plant<-ordered(Top_Bac_subset3_use_taxo$Plant, levels=c("Carrot", "Great Masterwort", "Sunflower", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed", "Rorippa sylvestris", "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus", "Melon","Alfalfa", "Bean","Hairy vetch", "Medicago truncatula","Pea", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright", "Phelipanche ramosa", "Rhinanthus glacialis", "Dahurian wildrye", "Festuca rubra", "Lolium arundinacea", "Lolium perenne", "Oat", "Rice", "Setaria pumila", "Setaria viridis", "Siberian wildrye", "Wheat", "Tobacco", "Tomato"))

Top_Bac_subset3_use_taxo$Taxon<-ordered(Top_Bac_subset3_use_taxo$Taxon, levels=c("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium-2", "Pseudomonas-4", "Pseudomonas-5", "Methylobacterium", "Unclassified Microbacteriaceae", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium-1", "Unclassified Enterobacteriaceae", "Pseudomonas-3", "Paenibacillus", "Pseudomonas-2", "Sphingomonas", "Pseudomonas-1", "Pantoea"))


tax_colors_16S <-  c('Acidobacteria'='#ffbb94','Actinobacteria'='#cf7773' ,'Bacteroidetes'='#82e4de', 'Cyanobacteria'='#547dae','Firmicutes'='#d69ec9' ,'Proteobacteria'='#d4e79c', 'Spirochaetes'='#4f6457','Unclassified'='#b4b4b4' ,'Tenericutes'='#c05805','Other'='black', 'Verrucomicrobia'= "#f4cc95", "WPS-2"="#cfa7da", "Epsilonbacteraeota"="#acd0c0", "Chloroflexi"="#b76c3b", "Deinococcus-Thermus"="#b29e54", "Thaumarchaeota"="#90639f" )

library(ggplot2)
library(ggridges)
plot_ridge_16S <- ggplot(Top_Bac_subset3_use_taxo, aes(x=SampleID, y=Taxon, group=Taxon, height=Relative_Abundance, color=Phylum, fill=Phylum)) +
    geom_density_ridges2(stat = "identity", scale=0.9,size=0.7) +scale_color_manual(values = tax_colors_16S)+scale_fill_manual(values = tax_colors_16S)+facet_grid(.~Plant, scales="free", space = "free",labeller=label_wrap_gen(multi_line = TRUE))+
theme(legend.position="bottom") + theme_classic()+ theme(axis.title.y = element_text(color="black", size=11, face="bold"))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+ggtitle("16S rRNA gene - Bacteria & Archaea")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+theme(strip.text.x = element_text(size=7, angle=90, face = "bold")) +xlab("Samples")+ylab("ASVs present in >20 Plant Species")+theme(legend.position = "bottom") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.05, linetype="solid"))+ theme(panel.spacing = unit(0.1, "lines"))
plot_ridge_16S
```

### Make color band for plant species that goes on top of Fig5 for 16S
```{r}
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
Top_Bac_subset3_use_taxo$Plant<-ordered(Top_Bac_subset3_use_taxo$Plant, levels=c("Carrot", "Great Masterwort", "Sunflower", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed", "Rorippa sylvestris", "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus", "Melon","Alfalfa", "Bean","Hairy vetch", "Medicago truncatula","Pea", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright", "Phelipanche ramosa", "Rhinanthus glacialis", "Dahurian wildrye", "Festuca rubra", "Lolium arundinacea", "Lolium perenne", "Oat", "Rice", "Setaria pumila", "Setaria viridis", "Siberian wildrye", "Wheat", "Tobacco", "Tomato"))

plant_band=ggplot(Top_Bac_subset3_use_taxo, aes(x=SampleID, y=Group, color=Plant))+geom_point(shape=15)+facet_grid(.~Plant, scales="free", space = "free",labeller=label_wrap_gen(multi_line = TRUE))+scale_color_manual(values = color)+ theme_classic()+ theme(axis.title.y = element_text(color="black", size=11, face="bold"))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+theme(strip.text.x = element_text(size=7, angle=90, face = "bold")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.05, linetype="solid"))+ theme(panel.spacing = unit(0.1, "lines"))+theme(legend.position = "none") 
plant_band
```


## For ITS
### Import SV table with only the most prevalent bacterial taxa (>20 species)
```{r}
meta2_ITS <- read.table("Metadata_ITS1_ITS2_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
Top_Fun_subset3<-read.table("Subset3-ITS1_table-FINAL-rarefied-mostprevalent16.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Top_Fun_subset3)
head(meta2_ITS)
dim(meta2_ITS)
dim(Top_Fun_subset3)

```

### Transform matrix into long table format
```{r}
library(reshape2)
Top_Fun_subset3_long=setNames(melt(Top_Fun_subset3), c('SampleID', 'Taxon', 'Relative_Abundance'))
head(Top_Fun_subset3_long)
dim(Top_Fun_subset3_long)
```
```{r}
##Merging long format table with metadata
Top_Fun_subset3_use<-merge(meta2_ITS,Top_Fun_subset3_long,by="SampleID")
dim(Top_Fun_subset3_use)

head(Top_Fun_subset3_use)
```

## Merge long format table with taxonomic info for each SV
```{r}
taxo<-read.table("Subset1-2_All_studies_merged_ITS1_taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
Top_Fun_subset3_use_taxo<-merge(Top_Fun_subset3_use,taxo,by="Taxon")
dim(Top_Fun_subset3_use_taxo)

head(Top_Fun_subset3_use_taxo)
```

## Figure 5C - Plotting ITS most prevalent taxa (n=14) with ggridge
```{r}
Top_Fun_subset3_use_taxo$Plant<-ordered(Top_Fun_subset3_use_taxo$Plant, levels=c("Carrot", "Great Masterwort", "Sunflower", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed","Rorripa sylvestris",  "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus","Alfalfa", "Bean","Hairy vetch", "Medicago truncatula","Pea", "Red Clover", "White Clover", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright","Phelipanche ramosa", "Rhinanthus glacialis","Creeping Bentgrass","Dahurian wildrye", "Oat", "Rice","Siberian wildrye","Wheat", "Tomato"))

Top_Fun_subset3_use_taxo$Taxon<-ordered(Top_Fun_subset3_use_taxo$Taxon, levels=c('Alternaria metachromatica-3','Alternaria metachromatica-2','Gibberella acuminata','Vishniacozyma tephrensis','Unclassified Sclerotiniaceae','Aureobasidium pullulans', 'Cladosporium perangustum-2', 'Vishniacozyma victoriae-2','Epicoccum nigrum','Sporobolomyces roseus', 'Vishniacozyma victoriae-1', 'Filobasidium sp', "Alternaria sp", "Alternaria metachromatica-1", "Unclassified Capnodiales", "Cladosporium perangustum-1"))

tax_colors_ITS <-  c('Agaricomycetes'='#9089b7','Archaeorhizomycetes'='#d35c37' ,'Dothideomycetes'='#fcc875', 'Eurotiomycetes'='#ecdbd2','Leotiomycetes'='#e99787','Malasseziomycetes'='#7e8aab','Microbotryomycetes'='#805778','Saccharomycetes'='#ab4f5f','Sordariomycetes'='#e58f66','Other'='black','Tremellomycetes'='#97b8c2','Wallemiomycetes'='#d4c8d1')

library(ggplot2)
library(ggridges)
plot_ridge_ITS <- ggplot(Top_Fun_subset3_use_taxo, aes(x=SampleID, y=Taxon, group=Taxon, height=Relative_Abundance, color=Class, fill=Class)) +
    geom_density_ridges2(stat = "identity", scale=0.9,size=0.7) +scale_color_manual(values = tax_colors_ITS)+scale_fill_manual(values = tax_colors_ITS)+facet_grid(.~Plant, scales="free", space = "free",labeller=label_wrap_gen(multi_line = TRUE))+
theme(legend.position="bottom") + theme_classic()+ theme(axis.title.y = element_text(color="black", size=11, face="bold"))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+ggtitle("ITS1 Region - Fungi")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+theme(strip.text.x = element_text(size=7, angle=90, face = "bold")) +xlab("Samples")+ylab("ASVs present in >20 Plant Species")+theme(legend.position = "bottom") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.05, linetype="solid"))+ theme(panel.spacing = unit(0.1, "lines"))
plot_ridge_ITS

```

### Make color band for plant species that goes on top of Fig5 for ITS
```{r}
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103")
Top_Fun_subset3_use_taxo$Plant<-ordered(Top_Fun_subset3_use_taxo$Plant, levels=c("Carrot", "Great Masterwort", "Sunflower", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed","Rorripa sylvestris",  "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus","Alfalfa", "Bean","Hairy vetch", "Medicago truncatula","Pea", "Red Clover", "White Clover", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright","Phelipanche ramosa", "Rhinanthus glacialis","Creeping Bentgrass","Dahurian wildrye", "Oat", "Rice","Siberian wildrye","Wheat", "Tomato"))

plant_band=ggplot(Top_Fun_subset3_use_taxo, aes(x=SampleID, y=Group, color=Plant))+geom_point(shape=15)+facet_grid(.~Plant, scales="free", space = "free",labeller=label_wrap_gen(multi_line = TRUE))+scale_color_manual(values = color)+ theme_classic()+ theme(axis.title.y = element_text(color="black", size=11, face="bold"))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+theme(strip.text.x = element_text(size=7, angle=90, face = "bold")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.05, linetype="solid"))+ theme(panel.spacing = unit(0.1, "lines"))+theme(legend.position = "none") 
plant_band
```


## For gyrB
### Import SV table with only the most prevalent bacterial taxa (>12 species)
```{r}
meta2_gyrB <- read.table("Metadata_gyrB_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
Top_Bac2_subset3<-read.table("Subset3-gyrB-MiSeq_table-FINAL-rarefied-18mostprevalent.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Top_Bac2_subset3)
head(meta2_gyrB)
dim(meta2_gyrB)
dim(Top_Bac2_subset3)

```

### Transform matrix into long table format
```{r}
library(reshape2)
Top_Bac2_subset3_long=setNames(melt(Top_Bac2_subset3), c('SampleID', 'Taxon', 'Relative_Abundance'))
head(Top_Bac2_subset3_long)
dim(Top_Bac2_subset3_long)
```
```{r}
##Merging long format table with metadata
Top_Bac2_subset3_use<-merge(meta2_gyrB,Top_Bac2_subset3_long,by="SampleID")
dim(Top_Bac2_subset3_use)

head(Top_Bac2_subset3_use)
```

## Merge long format table with taxonomic info for each SV
```{r}
taxo<-read.table("Subset1-2_All_studies_merged_gyrB-rep-seqs-FINAL-filtered-taxonomy-final.tsv", header=TRUE, check.names = FALSE, sep = "\t")
Top_Bac2_subset3_use_taxo<-merge(Top_Bac2_subset3_use,taxo,by="Taxon")
dim(Top_Bac2_subset3_use_taxo)

head(Top_Bac2_subset3_use_taxo)
```



## Figure 5B - Plotting gyrB most prevalent taxa (n=18) with ggridge
```{r}
Top_Bac2_subset3_use_taxo$Plant<-ordered(Top_Bac2_subset3_use_taxo$Plant, levels=c("Carrot", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed", "Sinapis arvensis", "Turnip", "Melon", "Bean", "Medicago truncatula", "Tomato"))

Top_Bac2_subset3_use_taxo$Taxon<-ordered(Top_Bac2_subset3_use_taxo$Taxon, levels=c('Pseudomonas sp-2','Methylobacterium sp','Bacillus altitudinis complex', 'Pseudomonas sp-1', 'Unclassified Enterobacteriaceae','Pseudomonas antarctica','Pseudomonas poae', 'Paenibacillus sp-2', "Rhizobium sp", "Paenibacillus sp-1", "Pseudomonas viridiflava-3", "Stenotrophomonas sp", "Pseudomonas fluorescens-3", "Pseudomonas fluorescens-2", "Erwinia persicina", "Pseudomonas fluorescens-1", "Pantoea agglomerans-3", "Pseudomonas viridiflava-2", "Pseudomonas viridiflava-1","Cutibacterium acnes", "Unclassified Clostridiales", "Pantoea agglomerans-2", "Pantoea agglomerans-1"))



tax_colors_gyrB <- c('Acidobacteria'='#ffbb94','Actinobacteria'='#cf7773' ,'Bacteroidetes'='#82e4de', 'Cyanobacteria'='#547dae','Firmicutes'='#d69ec9' ,'Proteobacteria'='#d4e79c', 'Spirochaetes'='#4f6457','Unclassified'='#b4b4b4' ,'Tenericutes'='#c05805','Other'='black', 'Verrucomicrobia'= "#f4cc95", "WPS-2"="#cfa7da", "Epsilonbacteraeota"="#acd0c0", "Chloroflexi"="#b76c3b", "Deinococcus-Thermus"="#b29e54", "Thaumarchaeota"="#90639f" )


library(ggplot2)
library(ggridges)
plot_ridge_gyrB <- ggplot(Top_Bac2_subset3_use_taxo, aes(x=SampleID, y=Taxon, group=Taxon, height=Relative_Abundance, color=Phylum, fill=Phylum)) +
    geom_density_ridges2(stat = "identity", scale=0.9,size=0.7) +scale_color_manual(values = tax_colors_gyrB)+scale_fill_manual(values = tax_colors_gyrB)+facet_grid(.~Plant, scales="free", space = "free",labeller=label_wrap_gen(multi_line = TRUE))+
theme(legend.position="bottom") + theme_classic()+ theme(axis.title.y = element_text(color="black", size=11, face="bold"))+ theme(axis.text.y = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+ggtitle("gyrB gene - Bacteria")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+theme(strip.text.x = element_text(size=7, angle=90, face = "bold")) +xlab("Samples")+ylab("ASVs present in >12 Plant Species")+theme(legend.position = "bottom") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.05, linetype="solid"))+ theme(panel.spacing = unit(0.1, "lines"))
plot_ridge_gyrB

```

### Make color band for plant species that goes on top of Fig5 for gyrB
```{r}
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103")
Top_Bac2_subset3_use_taxo$Plant<-ordered(Top_Bac2_subset3_use_taxo$Plant, levels=c("Carrot", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed", "Sinapis arvensis", "Turnip", "Melon", "Bean", "Medicago truncatula", "Tomato"))

plant_band=ggplot(Top_Bac2_subset3_use_taxo, aes(x=SampleID, y=Group, color=Plant))+geom_point(shape=15)+facet_grid(.~Plant, scales="free", space = "free",labeller=label_wrap_gen(multi_line = TRUE))+scale_color_manual(values = color)+ theme_classic()+ theme(axis.title.y = element_text(color="black", size=11, face="bold"))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+theme(strip.text.x = element_text(size=7, angle=90, face = "bold")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.05, linetype="solid"))+ theme(panel.spacing = unit(0.1, "lines"))+theme(legend.position = "none") 
plant_band
```





##########################################################################################################
###########################################################################################################


# Figure 6 - Plant-specific patterns in community composition/structure 

## Figure 6 Panels C and D - Bargarphs taxonomy by plant
## 16S V4: Import Subset 2 dataset (rarefied by study)
```{r}
meta2 <- read.table("Metadata_16S_V4_V5V6_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
Subset2_16S<-read.table("Subset2-All_studies_merged_16S_table-FINAL-V4-MiSeq-filtered-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset2_16S)
head(meta2)
dim(meta2)
dim(Subset2_16S)
SV_16S_use1<-merge(meta2,Subset2_16S,by="SampleID")
dim(SV_16S_use1)

head(SV_16S_use1)
```


### Make matrix of just SVs without metadata 
```{r}
dim(SV_16S_use1)
matrix_16S<-SV_16S_use1[c(39:31046)]
dim(matrix_16S)
head(matrix_16S)
```

### Prepare phyloseq object
```{r}
library(microbiome)
##Convert data as phyloseq object
matrix_16St=t(matrix_16S)
taxo_16S2 <- read.table(file="Subset1-2_All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-taxonomy.tsv", sep='\t', header=TRUE,check.names=FALSE,row.names=1) 
taxo_16S2=as.matrix(taxo_16S2)
TAXO_16S = tax_table(taxo_16S2)
OTU_16S = otu_table(matrix_16St, taxa_are_rows = TRUE)
meta_16S=SV_16S_use1[c(1:35)]
META_16S=sample_data(meta_16S)
physeq_16S = phyloseq(OTU_16S, TAXO_16S,META_16S)
physeq_16S
```

```{r}
library(phyloseq)
#For each plant species, Merge all samples from the same study together 
physeq_16S_merged=merge_samples(physeq_16S, "PlantbyStudy")
physeq_16S_merged
```

```{r}
#Turn all OTUs into class (or phylum or order level) counts
glom16S <- tax_glom(physeq_16S_merged, taxrank = 'Phylum')
glom16S # should list # taxa as # phyla
glom16S2 = transform_sample_counts(glom16S, function(x) x / sum(x) )
glom16S2
```
```{r}
data_glom16S<- psmelt(glom16S2) # create dataframe from phyloseq object
data_glom16S$Phylum <- as.character(data_glom16S$Phylum) #convert to character

### Recreate the Plant and Study_ID columns lost during sample merging
library(stringr)
col16S=str_split_fixed(data_glom16S$Sample, "--", 2)
col16S=data.frame(col16S)
names(col16S)[1] <- "Plant2"
names(col16S)[2] <- "Study_ID2"
data_glom16S=cbind(data_glom16S, col16S)
```

```{r}
#simple way to rename phyla with < 1% abundance
data_glom_16S3=data_glom16S
data_glom_16S3$Phylum[data_glom_16S3$Abundance < 0.003] <- "Other"

#Count # phyla to set color palette
Count = length(unique(data_glom_16S3$Phylum))
Count
#write.table(data_glom16S, file = "Subset2-16S-V4-table_glom_phylum.txt", sep = "\t")
```

# Figure 6C: Bar graph 16S taxonomy
```{r}
data_glom_16S3$Plant2<-ordered(data_glom_16S3$Plant2, levels=c("Carrot", "Great Masterwort", "Sunflower", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed", "Rorippa sylvestris", "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus", "Melon","Alfalfa", "Bean","Hairy vetch", "Medicago truncatula","Pea", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright", "Phelipanche ramosa", "Rhinanthus glacialis", "Dahurian wildrye", "Festuca rubra", "Lolium arundinacea", "Lolium perenne", "Oat", "Rice", "Setaria pumila", "Setaria viridis", "Siberian wildrye", "Wheat", "Tobacco", "Tomato"))

data_glom_16S3$Phylum<-ordered(data_glom_16S3$Phylum, levels=c('Other','Unclassified','Acidobacteria','Actinobacteria','Bacteroidetes',"Chloroflexi",'Cyanobacteria',"Deinococcus-Thermus","Epsilonbacteraeota",'Firmicutes','Proteobacteria','Spirochaetes', 'Tenericutes', 'Microbotryomycetes','Tremellomycetes',"Verrucomicrobia","WPS-2",'Wallemiomycetes', "Thaumarchaeota"))

tax_colors_16S <-  c('Acidobacteria'='#ffbb94','Actinobacteria'='#cf7773' ,'Bacteroidetes'='#82e4de', 'Cyanobacteria'='#547dae','Firmicutes'='#d69ec9' ,'Proteobacteria'='#d4e79c', 'Spirochaetes'='#4f6457','Unclassified'='#b4b4b4' ,'Tenericutes'='#c05805','Other'='black', 'Verrucomicrobia'= "#f4cc95", "WPS-2"="#cfa7da", "Epsilonbacteraeota"="#acd0c0", "Chloroflexi"="#b76c3b", "Deinococcus-Thermus"="#b29e54", "Thaumarchaeota"="#90639f" )

#plot with condensed phyla into "unknown" category
f5_16S <- ggplot(data=data_glom_16S3, aes(x=Abundance, y=Study_ID2, fill=Phylum))+facet_grid(Plant2~., scales="free", space = "free")+ geom_bar(aes(), stat="identity", position="stack") +
theme(legend.position="bottom") + theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=8, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid"))+ggtitle("16S rRNA gene - Bacteria & Archaea")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold")) +xlab("Relative Abundance")+ylab("Studies")+ scale_x_continuous(labels = scales::percent_format(accuracy = 1))+theme(legend.position = "bottom") +scale_fill_manual(values=tax_colors_16S)+guides(fill=guide_legend(nrow=3,byrow=TRUE))+ theme(panel.spacing.y = unit(0.01, "lines"))

f5_16S
```


## For ITS Import Subset 2 dataset (rarefied by study)
```{r}
meta2 <- read.table("Metadata_ITS1_ITS2_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
Subset2_ITS1<-read.table("Subset2-ITS1region_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset2_ITS1)
head(meta2)
dim(meta2)
dim(Subset2_ITS1)
SV_ITS_use1<-merge(meta2,Subset2_ITS1,by="SampleID")
dim(SV_ITS_use1)

head(SV_ITS_use1)
```


### Make matrix of just SVs without metadata 
```{r}
dim(SV_ITS_use1)
matrix_ITS<-SV_ITS_use1[c(34:ncol(SV_ITS_use1))]
dim(matrix_ITS)
head(matrix_ITS)
```

### Prepare phyloseq object
```{r}
library(microbiome)
##Convert data as phyloseq object
matrix_ITSt=t(matrix_ITS)
taxo_ITS2 <- read.table(file="Subset1-2_All_studies_merged_ITS1_taxonomy.tsv", sep='\t', header=TRUE,check.names=FALSE,row.names=1) 
taxo_ITS2=as.matrix(taxo_ITS2)
TAXO_ITS = tax_table(taxo_ITS2)
OTU_ITS = otu_table(matrix_ITSt, taxa_are_rows = TRUE)
meta_ITS=SV_ITS_use1[c(1:35)]
META_ITS=sample_data(meta_ITS)
physeq_ITS = phyloseq(OTU_ITS, TAXO_ITS,META_ITS)
physeq_ITS

```


```{r}
library(phyloseq)
#For each plant species, Merge all samples from the same study together 
physeq_ITS_merged=merge_samples(physeq_ITS, "PlantbyStudy")
physeq_ITS_merged
```

```{r}
#Turn all OTUs into class (or phylum or order level) counts
glomITS <- tax_glom(physeq_ITS_merged, taxrank = 'Class')
glomITS # should list # taxa as # phyla
glomITS2 = transform_sample_counts(glomITS, function(x) x / sum(x) )
glomITS2
data_glomITS<- psmelt(glomITS2) # create dataframe from phyloseq object
data_glomITS$Class <- as.character(data_glomITS$Class) #convert to character

#Recreate the Plant and Study_ID columns lost during sample merging
library(stringr)
colITS=str_split_fixed(data_glomITS$Sample, "--", 2)
colITS=data.frame(colITS)
names(colITS)[1] <- "Plant2"
names(colITS)[2] <- "Study_ID2"
data_glomITS=cbind(data_glomITS, colITS)
```

```{r}
#simple way to rename phyla with < 1% abundance
data_glomITS2=data_glomITS
data_glomITS2$Class[data_glomITS2$Abundance < 0.01] <- "Other"

#Count # phyla to set color palette
Count = length(unique(data_glomITS2$Class))
Count

```

# Figure 6D - Bar graph IST1 taxonomy
```{r}
data_glomITS2$Plant2<-ordered(data_glomITS2$Plant2, levels=c("Carrot", "Great Masterwort", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden_rocket", "Radish", "Rapeseed", "Rorripa sylvestris", "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus", "Bean", "Medicago truncatula","White Clover", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright", "Phelipanche ramosa", "Rhinanthus glacialis", "Wheat", "Tomato"))

data_glomITS2$Class<-ordered(data_glomITS2$Class, levels=c( 'Other',"Unclassified",'Archaeorhizomycetes','Dothideomycetes','Eurotiomycetes','Leotiomycetes','Saccharomycetes','Sordariomycetes','Agaricomycetes',"Agaricostilbomycetes", 'Malasseziomycetes', 'Microbotryomycetes', "Mortierellomycetes","Pezizomycetes","Taphrinomycetes",'Tremellomycetes','Wallemiomycetes'))

tax_colors_ITS <-  c('Agaricomycetes'='#9089b7', "Agaricostilbomycetes"="gray", 'Archaeorhizomycetes'='#d35c37' ,'Dothideomycetes'='#fcc875', 'Eurotiomycetes'='#ecdbd2','Leotiomycetes'='#a5a58d','Malasseziomycetes'='#7e8aab','Microbotryomycetes'='#805778',"Mortierellomycetes"="#a8dadc", "Pezizomycetes"="#cb997e", 'Saccharomycetes'='#ab4f5f','Sordariomycetes'='#e58f66','Other'='black',"Taphrinomycetes"="pink", 'Tremellomycetes'='#219ebc',"Unclassified"="black", 'Wallemiomycetes'='#d4c8d1')


#plot with condensed phyla into "unknown" category
f5_ITS <- ggplot(data=data_glomITS2, aes(x=Abundance, y=Study_ID2, fill=Class))+facet_grid(Plant2~., scales="free", space = "free")+ geom_bar(aes(), stat="identity", position="stack") +
theme(legend.position="bottom") + theme_classic()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=8, face="bold")) + theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid"))+ggtitle("ITS Region - Fungi")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold")) +xlab("Relative Abundance")+ylab("Studies")+ scale_x_continuous(labels = scales::percent_format(accuracy = 1))+theme(legend.position = "bottom") + theme(panel.spacing = unit(0.1, "lines"))+scale_fill_manual(values=tax_colors_ITS)+guides(fill=guide_legend(nrow=4,byrow=TRUE))

f5_ITS
```


## Figure 6 panels C & D combined bargraph for 16S and ITS
```{r}
library(ggpubr)
figure6CD=ggarrange(f5_16S, f5_ITS,labels = c("A", "B"), ncol = 2, widths = c(1.2, 1))
figure6CD
```


# Figure 6 Part 1 - Beta div - Performed on Subset 3 because needs to be perfomed on common distance matrix

### Import Subset 3 dataset (rarefied across all studies) - ITS
```{r}
meta2 <- read.table("Metadata_ITS1_ITS2_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
Subset3_ITS1<-read.table("Subset3-ITS1_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset3_ITS1)
head(meta2)
dim(meta2)
dim(Subset3_ITS1)
Subset3_ITS1_use<-merge(meta2,Subset3_ITS1,by="SampleID")
dim(Subset3_ITS1_use)

head(Subset3_ITS1_use)
Subset3_ITS1_use$Plant<-ordered(Subset3_ITS1_use$Plant, levels=c("Carrot", "Great Masterwort", "Sunflower", "Alliaria petiolata", "Arabidopsis thaliana", "Barbarea vulgaris ", "Berteroa incana", "Brassica nigra", "Broccoli", "Cabbage", "Capsella bursa-pastoris", "Cardamine hirsuta", "Cauliflower", "Erophila verna ", "Garden rocket", "Radish", "Rapeseed","Rorripa sylvestris",  "Sinapis arvensis", "Turnip", "Pincushion Flower", "Heliosperma alpestre", "Grass of Parnassus","Alfalfa", "Bean","Hairy vetch", "Medicago truncatula","Pea", "Red Clover", "White Clover", "Oak", "Chiltern Gentian", "Willow Gentian", "Eyebright","Phelipanche ramosa", "Rhinanthus glacialis","Creeping Bentgrass","Dahurian wildrye", "Oat", "Rice","Siberian wildrye","Wheat", "Tomato"))
```

### Make matrix of just SVs without metadata 
```{r}
dim(Subset3_ITS1_use)
matrix3_ITS<-Subset3_ITS1_use[c(34:ncol(Subset3_ITS1_use))]
dim(matrix3_ITS)
head(matrix3_ITS)
```


## Make phyloseq object
```{r}
library(microbiome)
##Convert data as phyloseq object
matrix3_ITSt=t(matrix3_ITS)
OTU_ITS_subset3 = otu_table(matrix3_ITSt, taxa_are_rows = TRUE)
meta_ITS_Subset3=Subset3_ITS1_use[c(1:33)]
META_ITS_Subset3=sample_data(meta_ITS_Subset3)
physeq_ITS_Subset3 = phyloseq(OTU_ITS_subset3,META_ITS_Subset3)
physeq_ITS_Subset3

```

## Figure 6B - ITS Plot PCoA on all samples - Bray-Curtis
```{r}
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
ITS.ord <- ordinate(physeq_ITS_Subset3, "PCoA", "bray")
b3 = plot_ordination(physeq_ITS_Subset3, ITS.ord, type="samples", color="Plant", shape="Seed_fraction", title="ITS1 region - Fungi (Bray-Curtis)")+ scale_color_manual(values=color)+theme_classic(base_size = 12)+xlab("PCoA1 = 18.8%")+ylab("PCoA2 = 11%")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=7, face="bold"))+theme(legend.position = "none")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+scale_shape_manual(values=c(1, 16))
b3
```

### Import Subset 3 dataset (rarefied across all studies) - 16S
```{r}
Subset3_16S <- read.table("Subset3-16S-V4-table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
meta2 <- read.table("Metadata_16S_V4_V5V6_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset3_16S)
head(meta2)
dim(meta2)
dim(Subset3_16S)
Subset3_16S_use<-merge(meta2,Subset3_16S,by="SampleID")
dim(Subset3_16S_use)
```

### Make matrix of just SVs without metadata 
```{r}
dim(Subset3_16S_use)
matrix3_16S<-Subset3_16S_use[c(39:ncol(Subset3_16S_use))]
dim(matrix3_16S)
head(matrix3_16S)
```

## Make phyloseq object
```{r}
library(microbiome)
##Convert data as phyloseq object
matrix3_16St=t(matrix3_16S)
OTU_16S_subset3 = otu_table(matrix3_16St, taxa_are_rows = TRUE)
meta_16S_Subset3=Subset3_16S_use[c(1:33)]
META_16S_Subset3=sample_data(meta_16S_Subset3)
physeq_16S_Subset3 = phyloseq(OTU_16S_subset3,META_16S_Subset3)
physeq_16S_Subset3
```

## Figure 6A - 16S Plot PCoA on all samples - Bray-Curtis
```{r}
bac.ord <- ordinate(physeq_16S_Subset3, "PCoA", "bray")
b1z = plot_ordination(physeq_16S_Subset3, bac.ord, type="samples", color="Plant", shape="Seed_fraction", title="16S rRNA gene region - Bacteria & Archaea (Bray-Curtis)")+ scale_color_manual(values=color)+theme_classic(base_size = 12)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=7, face="bold"))+theme(legend.position = "none")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+scale_shape_manual(values=c(3, 1,2, 16))+xlab("PCoA1 = 15.8%")+ylab("PCoA2 = 11.2%")+scale_y_continuous(limits = c(-0.12, 0.2))
print(b1z)
```


# Make Fig 6 panel A & B ordinations
```{r}
library(ggpubr)
figure6AB=ggarrange(b1z, b3,labels = c("A", "B"), ncol = 2)
figure6AB
```

# Figure S10 Picrust2 on 16S and gyrB datasets
### Import Picrust2 KOs Subset 1 dataset - 16S-V4
```{r}
Subset3_16S_picrust <- read.table("16S_Subset1_KOs_pred_metagenome_unstrat.tsv", header=TRUE, check.names = FALSE, sep = "\t")
meta2 <- read.table("Metadata_16S_V4_V5V6_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset3_16S_picrust)
head(meta2)
dim(meta2)
dim(Subset3_16S_picrust)
Subset3_16S_picrust_use<-merge(meta2,Subset3_16S_picrust,by="SampleID")
dim(Subset3_16S_picrust_use)


```

### Make matrix of just SVs without metadata 
```{r}
dim(Subset3_16S_picrust_use)
matrix3_16S_picrust<-Subset3_16S_picrust_use[c(39:ncol(Subset3_16S_picrust_use))]
dim(matrix3_16S_picrust)
head(matrix3_16S_picrust)
```


## Make phyloseq object
```{r}
library(microbiome)
##Convert data as phyloseq object
matrix3_16S_picrustt=t(matrix3_16S_picrust)
OTU_16S_picrust_subset3 = otu_table(matrix3_16S_picrustt, taxa_are_rows = TRUE)
meta_16S_picrust_Subset3=Subset3_16S_picrust_use[c(1:33)]
META_16S_picrust_Subset3=sample_data(meta_16S_picrust_Subset3)
physeq_16S_picrust_Subset3 = phyloseq(OTU_16S_picrust_subset3,META_16S_picrust_Subset3)
physeq_16S_picrust_Subset3
```

## Figure S10A - Picrust 16S Plot PCoA on all samples - Bray-Curtis
```{r}
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
bac.pic <- ordinate(physeq_16S_picrust_Subset3, "NMDS", "bray")
pic2 = plot_ordination(physeq_16S_picrust_Subset3, bac.pic, type="samples", color="Plant", shape="Seed_fraction", title="16S rRNA gene - Bacteria & Archaea (Picrust2 KOs)")+ scale_color_manual(values=color)+theme_classic(base_size = 12)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=7, face="bold"))+theme(legend.position = "none")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+scale_shape_manual(values=c(3, 1,2, 16))
print(pic2)
```


## 16S KO richness by sample
```{r}
# KO richness (S) 
library(vegan)
S <- specnumber(matrix3_16S_picrust) 
Subset3_16S_picrust_use=cbind(Subset3_16S_picrust_use,S)
```

### Plot Richness (S) 
```{r warning=FALSE}
library(ggplot2)
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 

p4=ggplot(data=Subset3_16S_picrust_use, aes(x=S, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Observed KO richness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(3, 1,2, 16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.9, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"))+
  theme(legend.position = "none")+ggtitle("16S rRNA gene - Bacteria & Archaea") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+ theme(panel.spacing.y = unit(0.2, "lines"))
p4
```



### Import Picrust2 KOs Subset 3 dataset (rarefied across all studies) - gyrB
```{r}
Subset3_gyrB_picrust <- read.table("Subset1_gyrB_KO_pred_metagenome_unstrat.tsv", header=TRUE, check.names = FALSE, sep = "\t")
meta2 <- read.table("Metadata_gyrB_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset3_gyrB_picrust)
head(meta2)
dim(meta2)
dim(Subset3_gyrB_picrust)
Subset3_gyrB_picrust_use<-merge(meta2,Subset3_gyrB_picrust,by="SampleID")
dim(Subset3_gyrB_picrust_use)
```

### Make matrix of just SVs without metadata 
```{r}
dim(Subset3_gyrB_picrust_use)
matrix3_gyrB_picrust<-Subset3_gyrB_picrust_use[c(40:ncol(Subset3_gyrB_picrust_use))]
dim(matrix3_gyrB_picrust)
head(matrix3_gyrB_picrust)
```

## Make phyloseq object
```{r}
library(microbiome)
##Convert data as phyloseq object
matrix3_gyrB_picrustt=t(matrix3_gyrB_picrust)
OTU_gyrB_picrust_subset3 = otu_table(matrix3_gyrB_picrustt, taxa_are_rows = TRUE)
meta_gyrB_picrust_Subset3=Subset3_gyrB_picrust_use[c(1:33)]
META_gyrB_picrust_Subset3=sample_data(meta_gyrB_picrust_Subset3)
physeq_gyrB_picrust_Subset3 = phyloseq(OTU_gyrB_picrust_subset3,META_gyrB_picrust_Subset3)
physeq_gyrB_picrust_Subset3
```

## Subset with 3 main plant species
```{r}
#create subset with only plant species with min 3 independent studies and 50 samples
Species = c("Radish", "Rapeseed", "Bean")
physeq_gyrB_picrust_Subset3_Core=subset_samples(physeq_gyrB_picrust_Subset3, Plant %in% Species)
physeq_gyrB_picrust_Subset3_Core2=prune_taxa(taxa_sums(physeq_gyrB_picrust_Subset3_Core) >= 1, physeq_gyrB_picrust_Subset3_Core)
physeq_gyrB_picrust_Subset3_Core2
```

## Figure S10B - Picrust gyrB Plot PCoA on all samples - Bray-Curtis
```{r}
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 
bac.pic <- ordinate(physeq_gyrB_picrust_Subset3, "NMDS", "bray")
pic2 = plot_ordination(physeq_gyrB_picrust_Subset3, bac.pic, type="samples", color="Plant", shape="Seed_fraction", title="gyrB gene - Bacteria (Picrust2 KOs)")+ scale_color_manual(values=color)+theme_classic(base_size = 12)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=7, face="bold"))+theme(legend.position = "none")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+scale_shape_manual(values=c(1,16))
print(pic2)
```


## gyrB KO richness by sample
```{r}
# KO richness (S) 
library(vegan)
S <- specnumber(matrix3_gyrB_picrust) 
Subset3_gyrB_picrust_use=cbind(Subset3_gyrB_picrust_use,S)

```

### Plot Richness (S)  KO gyrB
```{r warning=FALSE}
library(ggplot2)
color=c("Carrot"="#ee8332", "Great Masterwort"="#f0810F", "Sunflower"="#fec767", "Alliaria petiolata"="#4a1777", "Arabidopsis thaliana"="#835e9f", "Barbarea vulgaris "="#baa5c8", "Berteroa incana"=	"#ceaac4", "Brassica nigra"="#a96699", "Broccoli"="#aea9ca", "Cabbage"=	"#811770", "Capsella bursa-pastoris"="#6b66a3", "Cardamine hirsuta"="#1b2b7d", "Cauliflower"="#5892ae", "Erophila verna "=	"#5e87c5", "Garden rocket"="#67bde8", "Radish"="#4cb5f5", "Rapeseed"="#008bc4", "Rorippa sylvestris"="#3e8a89", "Sinapis arvensis"="#9abdbb", "Turnip"=	"#c5d6d6", "Pincushion Flower"="#e09aa5", "Heliosperma alpestre"="#ffceda", "Grass of Parnassus"="#d0425d", "Melon"=	"#FA6775","Alfalfa"=	"#e035d7", "Bean"=	"#e03581", "Hairy vetch"=	"#e035ac",
"Medicago truncatula"=	"#b554a6", "Pea"=	"#68104d", "Oak"=	"#000000", "Chiltern Gentian"=	"#00a7b5", "Willow Gentian"="#15cd7e", "Eyebright"=	"#00bb9f", "Phelipanche ramosa"=	"#08ffda", "Rhinanthus glacialis"="#8eddad", "Dahurian wildrye"="#2fbd03","Festuca rubra"="#2E4600", "Lolium arundinacea"=	"#3c884c", "Lolium perenne"	="#72ae4f", "Oat"=	"#66ff00", "Rice"=	"#b2d24f", "Setaria pumila"=	"#c5cd77", "Setaria viridis"="#edf050", "Siberian wildrye"="#9ff76e", "Wheat"="#fdf351", "Tobacco"="#ff9b7c", "Tomato"="#ff2600", "Rorripa sylvestris"=	"#c473fa", "White Clover"="#8f868c", "Red Clover"="#dad7d9", "Creeping Bentgrass"="#bd9103") 

p4=ggplot(data=Subset3_gyrB_picrust_use, aes(x=S, y=Study_ID, color=Plant, shape=Seed_fraction)) + geom_jitter( alpha=0.8) +xlab("Observed KO richness")+ylab("Studies")+ theme_gray()+ theme(axis.title = element_text(color="black", size=9, face="bold"))+ theme(axis.text = element_text(color="black", size=7, face="bold"))+facet_grid(Plant~., scales="free", space = "free")+ theme(legend.text = element_text(color="black", size=8, face="bold"))+theme(strip.text.y = element_text(size=8, angle=0, face = "bold",margin = margin( b = 2, t = 2)))+scale_shape_manual(values=c(1, 16))+ labs(shape = "Seed Fraction", color = "Plant Species")+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(panel.background = element_rect(fill = "#eeeeee",colour = "#eeeeee",size = 0.9, linetype="solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "#eeeeee"))+ theme(strip.background = element_rect(fill = "#d9dbdb"))+
  theme(legend.position = "none")+ggtitle("gyrB gene - Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))+ theme(panel.spacing.y = unit(0.2, "lines"))
p4
```




# Figure 7 Core & Flexible taxa by plant species

## Load 16S-V4 Subset 3 and transform to long format for DB Browser

```{r}
Subset3_16S<-read.table("Subset3-16S-V4-table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names=1)
head(Subset3_16S)
dim(Subset3_16S)

```
## Convert Subset 3 table to long format for DB Browser
```{r}
library(reshape2)
Subset3_16S_long=setNames(melt(Subset2_16S), c("SampleID", 'SV', 'Count'))
head(Subset3_16S_long)
dim(Subset3_16S_long)
#write.table(Subset3_16S_long, "Subset3_16S_SVtable_long.txt")
```

## gyrB SubSet3- core taxa analysis
```{r}
Subset3_gyrB<-read.table("Subset3-gyrB-MiSeq_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset3_gyrB)

dim(Subset3_gyrB)
```
## Convert Subset 3 table to long format for DB Browser
```{r}
library(reshape2)
Subset3_gyrB_long=setNames(melt(Subset3_gyrB), c("SampleID", 'SV', 'Count'))
head(Subset3_gyrB_long)
dim(Subset3_gyrB_long)
#write.table(Subset3_gyrB_long, "Subset3_gyrB_SVtable_long.txt")
```


## ITS1 Subset3 - core taxa analysis
```{r}
Subset3_ITS<-read.table("Subset3-ITS1_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Subset3_ITS)
dim(Subset3_ITS)
```
## Convert Subset 2 table to long format for DB Browser
```{r}
library(reshape2)
Subset3_ITS_long=setNames(melt(Subset3_ITS), c("SampleID", 'SV', 'Count'))
head(Subset3_ITS_long)
dim(Subset3_ITS_long)
#write.table(Subset3_ITS_long, "Subset3_ITS_SVtable_long.txt")
```

### In DB Browser SQl, run the following code to get the table for core taxa analysis
```{r}
#For example for the 16S-V4 dataset (16S_SV is the long table we have just prepared, 16S_meta is the metadata table)
select t1.*, t2.sample_number_Plant, t2.count_sum_Plant  from
(
select SV, Plant, sum(Count) as count_sum_SV, count(*) as sample_number_SV,  count(DISTINCT Study_ID) as Nb_studies
from "16S_SV"
inner join "16S_meta"
on "16S_SV".SampleID = "16S_meta".SampleID
where count > 0
group by SV, Plant
) t1
inner join 
(
select Plant, sum(Count) as count_sum_Plant, count(distinct  "16S_SV".SampleID) as sample_number_Plant
from "16S_SV"
inner join "16S_meta"
on "16S_SV".SampleID = "16S_meta".SampleID
where count > 0
group by Plant
) t2
on t1.Plant = t2.Plant
```



# Figure 7 - Part 1: SV rel abund and prev
## Back in R: 16S-Import list core & flexible taxa for each plant species
```{r}
Core_Flex_16S<-read.table("16S-V4_core_flexible_abund_prev.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Core_Flex_16S)
dim(Core_Flex_16S)
```

```{r}
library(ggplot2)
color_core= c(Core="#ffa600", Flexible="#4F4A45")
p16S_core_flex=ggplot(data=Core_Flex_16S, aes(x=SV_Rel_Abund, y=SV_Prev,color=Type))+geom_point(size=3, alpha=0.7)+theme_classic(base_size = 12)+xlab("ASV Relative Abundance (%)")+ylab("ASV Prevalence (%)")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.title = element_text(colour="black", size=10, face="bold"))+ggtitle("16S rRNA gene - Bacteria & Archaea") +theme(plot.title = element_text(hjust = 0.5, face="bold", size = 14))+scale_x_log10(breaks =c(10, 0.1, 0.001), labels = c(10, 0.1, 0.001))+facet_grid(.~Plant)+ theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+theme(strip.text.x = element_text(face = "bold", size = 13))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.2, linetype="solid"), panel.grid.major.x = element_line(size = 0.3, linetype = 'dashed',colour = "grey"))+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = c(0.09, 0.85))+scale_color_manual(values=color_core)	+ labs( color = "ASV Type")
p16S_core_flex
```


## gyrB-Import list core & flexible taxa for each plant species
```{r}
Core_Flex_gyrB<-read.table("gyrB_core_flexible_abund_prev.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Core_Flex_gyrB)
dim(Core_Flex_gyrB)
```

```{r}
library(ggplot2)
color_core= c(Core="#ffa600", Flexible="#4F4A45")
pgyrB_core_flex=ggplot(data=Core_Flex_gyrB, aes(x=SV_Rel_Abund, y=SV_Prev,color=Type))+geom_point(size=3, alpha=0.7)+theme_classic(base_size = 12)+xlab("ASV Relative Abundance (%)")+ylab("ASV Prevalence (%)")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.title = element_text(colour="black", size=12, face="bold"))+ggtitle("gyrB gene - Bacteria") +theme(plot.title = element_text(hjust = 0.5, face="bold", size = 14))+scale_x_log10(breaks =c(10, 0.1, 0.001), labels = c(10, 0.1, 0.001))+facet_grid(.~Plant)+ theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+theme(strip.text.x = element_text(face = "bold", size = 13))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.2, linetype="solid"), panel.grid.major.x = element_line(size = 0.3, linetype = 'dashed',colour = "grey"))+ theme(panel.spacing = unit(0.1, "lines"))+scale_color_manual(values=color_core)	+theme(legend.position = "none")
pgyrB_core_flex
```

## ITS-Import list core & flexible taxa for each plant species
```{r}
Core_Flex_ITS<-read.table("ITS_core_flexible_taxa_abund_prev.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Core_Flex_ITS)
dim(Core_Flex_ITS)
```

```{r}
library(ggplot2)
color_core= c(Core="#ffa600", Flexible="#4F4A45")
pITS_core_flex=ggplot(data=Core_Flex_ITS, aes(x=SV_Rel_Abund, y=SV_Prev,color=Type))+geom_point(size=3, alpha=0.7)+theme_classic(base_size = 12)+xlab("ASV Relative Abundance (%)")+ylab("ASV Prevalence (%)")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.title = element_text(colour="black", size=12, face="bold"))+ggtitle("ITS1 Region - Fungi") +theme(plot.title = element_text(hjust = 0.5, face="bold", size = 14))+scale_x_log10(breaks =c(10, 0.1, 0.001), labels = c(10, 0.1, 0.001))+facet_grid(.~Plant)+ theme(strip.background = element_rect(fill = "#f0eeec"),strip.text = element_text(colour = "black", face = "bold"))+theme(strip.text.x = element_text(face = "bold", size = 13))+ theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.2, linetype="solid"), panel.grid.major.x = element_line(size = 0.3, linetype = 'dashed',colour = "grey"))+ theme(panel.spacing = unit(0.1, "lines"))+scale_color_manual(values=color_core)	+theme(legend.position = "none")
pITS_core_flex
```

# Figure 7 - Part 2: Cumulative Rel Abund of core & flexible
## 16S- Cumulative values for core and flexible
```{r}
Core_Flex_16S_Cum<-read.table("16S-V4_core_flexible_taxa_cumulative.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Core_Flex_16S_Cum)
dim(Core_Flex_16S_Cum)
```

```{r}
library(ggplot2)
color_core= c(Core="#ffa600", Flexible="#4F4A45")
color_num= c(Core="#de425b", Flexible="white")
Core_Flex_16S_Cum$Type<-ordered(Core_Flex_16S_Cum$Type, levels=c("Flexible", "Core"))
p16S_core_flex2=ggplot(data=Core_Flex_16S_Cum, aes(x=Plant, y=Rel_abund_SV,fill=Type))+ geom_bar(aes(), stat="identity", position="stack")+theme_classic(base_size = 14)+ylab("Cumulative ASV\nRelative Abundance (%)")+xlab("Plant Species")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold")) + theme(legend.position = c(0.05, 0.85))+scale_fill_manual(values=color_core)	+ labs( color = "ASV Type")+theme(legend.position = "none")+geom_text(aes(label = nb_SV, fontface = "bold", colour = Type, size=12),position = position_stack(vjust = .5))+scale_color_manual(values=color_num)
p16S_core_flex2
```

## gyrB- Cumulative values for core and flexible
```{r}
Core_Flex_gyrB_Cum<-read.table("gyrB_core_flexible_taxa_cumulative.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Core_Flex_gyrB_Cum)
dim(Core_Flex_gyrB_Cum)
```

```{r}
library(ggplot2)
color_core= c(Core="#ffa600", Flexible="#4F4A45")
color_num= c(Core="#de425b", Flexible="white")
Core_Flex_gyrB_Cum$Type<-ordered(Core_Flex_gyrB_Cum$Type, levels=c("Flexible", "Core"))
pgyrB_core_flex2=ggplot(data=Core_Flex_gyrB_Cum, aes(x=Plant, y=Rel_abund_SV,fill=Type))+ geom_bar(aes(), stat="identity", position="stack")+theme_classic(base_size = 14)+ylab("Cumulative ASVs Relative Abundance (%)")+xlab("Plant Species")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.position = c(0.05, 0.85))+scale_fill_manual(values=color_core)	+ labs( color = "ASV Type")+theme(legend.position = "none")+geom_text(aes(label = nb_SV, fontface = "bold", colour = Type, size=12),position = position_stack(vjust = .5))+scale_color_manual(values=color_num)+ theme(axis.title.y = element_blank())
pgyrB_core_flex2
```

## ITS- Cumulative values for core and flexible
```{r}
Core_Flex_ITS_Cum<-read.table("ITS_core_flexible_taxa_cumulative.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(Core_Flex_ITS_Cum)
dim(Core_Flex_ITS_Cum)
```


```{r}
library(ggplot2)
color_core= c(Core="#ffa600", Flexible="#4F4A45")
color_num= c(Core="#de425b", Flexible="white")
Core_Flex_ITS_Cum$Type<-ordered(Core_Flex_ITS_Cum$Type, levels=c("Flexible", "Core"))
pITS_core_flex2=ggplot(data=Core_Flex_ITS_Cum, aes(x=Plant, y=Rel_abund_SV,fill=Type))+ geom_bar(aes(), stat="identity", position="stack")+theme_classic(base_size = 14)+ylab("Cumulative ASVs Relative Abundance (%)")+xlab("Plant Species")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.position = c(0.05, 0.85))+scale_fill_manual(values=color_core)	+ labs( color = "ASV Type")+theme(legend.position = "none")+geom_text(aes(label = nb_SV, fontface = "bold", colour = Type, size=12),position = position_stack(vjust = .5))+scale_color_manual(values=color_num)+ theme(axis.title.y = element_blank())
pITS_core_flex2
```


# Figure 7 - Part 3: Taxo Bar graph Core & Flexible

## 16S-Merge Core & Flexible data with taxonomic info for each SV
```{r}
taxo<-read.table("Subset1-2_All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
Core_Flex_16S_taxo<-merge(Core_Flex_16S,taxo,by="SV")
dim(Core_Flex_16S_taxo)
head(Core_Flex_16S_taxo)
```



## Figure 16S-V4 for CORE taxa
```{r}
#Subset just core taxa
library(ggplot2)
library(dplyr)
Core_Flex_16S_taxo_core=subset(Core_Flex_16S_taxo, Type=="Core")

#Aggregate/Sum Rel abund by Plant and SV Species 
Core_Flex_16S_taxo_core_ag=Core_Flex_16S_taxo_core %>% 
  group_by(Plant, Genus) %>%                            # multiple group columns
  summarise(SV_Rel_Abund = sum(SV_Rel_Abund))  # multiple summary columns

head(Core_Flex_16S_taxo_core_ag)
dim(Core_Flex_16S_taxo_core_ag)
#Core_Flex_16S_taxo_core_ag$Genus[Core_Flex_16S_taxo_core_ag$SV_Rel_Abund < 1] <- "Other"
```

```{r}
tax_colors_16S_core <-  c('Blastomonas'='#ffbb94','Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium'='yellow' ,'Aureimonas'='#cfa7da', 'Bacillus'='#547dae','Chryseobacterium'='blue' ,'Curtobacterium'='#acd0c0', 'Escherichia-Shigella'='#4f6457','Enhydrobacter'='#b4b4b4' ,'Methylobacterium'='#c05805','Sphingobium'='black', 'Paenibacillus'= "#f4cc95", "Pantoea"="#5fa8d3", "Pseudomonas"="#d4e79c", "Sphingomonas"="orange", "Stenotrophomonas"="#b29e54", "Xanthomonas"="#90639f", "Unclassified Enterobacteriaceae"="lightgreen", "Unclassified Microbacteriaceae"="#cf7773","Streptococcus"="pink","Staphylococcus"="#b76c3b" )

Core_Flex_16S_taxo_core_ag$Genus<-ordered(Core_Flex_16S_taxo_core_ag$Genus, levels=c("Staphylococcus","Streptococcus", "Unclassified Microbacteriaceae","Stenotrophomonas",'Methylobacterium','Blastomonas','Escherichia-Shigella',"Enhydrobacter",'Curtobacterium', "Sphingobium",'Chryseobacterium','Bacillus','Aureimonas','Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium','Paenibacillus',"Xanthomonas", "Unclassified Enterobacteriaceae","Sphingomonas","Pseudomonas","Pantoea"))

p16S_core_flex_tax=ggplot(data=Core_Flex_16S_taxo_core_ag, aes(x=Plant, y=SV_Rel_Abund,fill=Genus))+ geom_bar(aes(), stat="identity", position="stack")+theme_classic(base_size = 14)+ylab("Core ASV\nRelative Abundance (%)")+xlab("Plant Species")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=10, face="bold"))+guides(fill=guide_legend(nrow=5))+theme(legend.position = "bottom", legend.spacing.x = unit(0.05, 'cm'), 
        legend.spacing.y = unit(0.05, 'cm'), 
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.justification="left")+
  scale_fill_manual(labels = c("Staphylococcus","Streptococcus", "Unclassified\nMicrobacteriaceae","Stenotrophomonas",'Methylobacterium','Blastomonas','Escherichia-Shigella',"Enhydrobacter",'Curtobacterium', "Sphingobium",'Chryseobacterium','Bacillus','Aureimonas','Allorhizobium-Neorhizobium\nPararhizobium-Rhizobium','Paenibacillus',"Xanthomonas", "Unclassified\nEnterobacteriaceae","Sphingomonas","Pseudomonas","Pantoea"), values=tax_colors_16S_core)
p16S_core_flex_tax
```

## gyrB-Merge Core & Flexible data with taxonomic info for each SV
```{r}
taxo_gyrB<-read.table("Subset1-2_All_studies_merged_gyrB-rep-seqs-FINAL-filtered-taxonomy-final.tsv", header=TRUE, check.names = FALSE, sep = "\t")
Core_Flex_gyrB_taxo<-merge(Core_Flex_gyrB,taxo_gyrB,by="SV")
dim(Core_Flex_gyrB_taxo)
head(Core_Flex_gyrB_taxo)
```


## Figure gyrB for CORE taxa
```{r}
#Subset just core taxa
library(ggplot2)
library(dplyr)
Core_Flex_gyrB_taxo_core=subset(Core_Flex_gyrB_taxo, Type=="Core")

#Aggregate/Sum Rel abund by Plant and SV Species 
Core_Flex_gyrB_taxo_core_ag=Core_Flex_gyrB_taxo_core %>% 
  group_by(Plant, Taxon) %>%                            # multiple group columns
  summarise(SV_Rel_Abund = sum(SV_Rel_Abund))  # multiple summary columns

head(Core_Flex_gyrB_taxo_core_ag)
dim(Core_Flex_gyrB_taxo_core_ag)
#Core_Flex_gyrB_taxo_core_ag$Genus[Core_Flex_gyrB_taxo_core_ag$SV_Rel_Abund < 1] <- "Other"
```

```{r}
tax_colors_gyrB <-  c('Erwinia persicina'='#e9c46a','Erwinia persicina-2'='#f4a261' ,'Erwinia sp'='#e76f51', 'Paenibacillus sp-1'='#5f0f40','Pantoea agglomerans-1'='#5fa8d3','Pantoea agglomerans-2'='#03045e','Pantoea agglomerans-3'='#90e0ef','Pantoea agglomerans-4'='#3a86ff','Pseudomonas fluorescens-5'='#40916c','Pseudomonas fluorescens-6'='#b7e4c7','Pseudomonas sp-3'='#1b4332','Pseudomonas syringae'='#c71f37', "Pseudomonas viridiflava-1"="#5f0f40", "Pseudomonas viridiflava-2"="#b392ac", "Serratia marcescens"="#5a189a" )


Core_Flex_gyrB_taxo_core_ag$Taxon<-ordered(Core_Flex_gyrB_taxo_core_ag$Taxon, levels=c("Pseudomonas viridiflava-1", "Pseudomonas viridiflava-2",'Pseudomonas fluorescens-5','Pseudomonas fluorescens-6','Pseudomonas sp-3', 'Paenibacillus sp-1','Pseudomonas syringae', "Serratia marcescens",'Erwinia sp','Erwinia persicina-2','Erwinia persicina', 'Pantoea agglomerans-4','Pantoea agglomerans-3','Pantoea agglomerans-2','Pantoea agglomerans-1'))
pgyrB_core_flex_tax=ggplot(data=Core_Flex_gyrB_taxo_core_ag, aes(x=Plant, y=SV_Rel_Abund,fill=Taxon))+ geom_bar(aes(), stat="identity", position="stack")+theme_classic(base_size = 14)+ylab("ASV Relative Abundance (%)")+xlab("Plant Species")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=10, face="bold"))+scale_fill_manual(values=tax_colors_gyrB)+guides(fill=guide_legend(nrow=4))+theme(legend.position = "bottom", legend.spacing.x = unit(0.05, 'cm'), 
        legend.spacing.y = unit(0.05, 'cm'), 
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.justification="left") + theme(axis.title.y = element_blank())
pgyrB_core_flex_tax
```


## ITS-Merge Core & Flexible data with taxonomic info for each SV
```{r}
taxo_ITS<-read.table("Subset1-2_All_studies_merged_ITS1_taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
Core_Flex_ITS_taxo<-merge(Core_Flex_ITS,taxo_ITS,by="SV")
dim(Core_Flex_ITS_taxo)
head(Core_Flex_ITS_taxo)
```



## Figure ITS1 for CORE taxa
```{r}
#Subset just core taxa
library(ggplot2)
library(dplyr)
Core_Flex_ITS_taxo_core=subset(Core_Flex_ITS_taxo, Type=="Core")

#Aggregate/Sum Rel abund by Plant and SV Species 
Core_Flex_ITS_taxo_core_ag=Core_Flex_ITS_taxo_core %>% 
  group_by(Plant, Genus) %>%                            # multiple group columns
  summarise(SV_Rel_Abund = sum(SV_Rel_Abund))  # multiple summary columns

head(Core_Flex_ITS_taxo_core_ag)
dim(Core_Flex_ITS_taxo_core_ag)
#Core_Flex_ITS_taxo_core_ag$Genus[Core_Flex_ITS_taxo_core_ag$SV_Rel_Abund < 1] <- "Other"
```

```{r}
tax_colors_ITS <-  c('Alternaria'='#9089b7','Bensingtonia'='#ecdbd2' ,'Bulleromyces'='#fcc875', 'Cladosporium'='#d35c37','Dioszegia'='#e99787','Filobasidium'='#ffbb94','Fusarium'='#805778','Gibberella'='#ab4f5f','Holtermanniella'='#e58f66','Mortierella'='black','Itersonilia'='#cf7773','Plectosphaerella'='#d4c8d1', "Sporobolomyces"="#b76c3b", "Stemphylium"="#b29e54", "Acremonium"="pink",'Unclassified Basidiomycota'='darkgray','Unclassified Capnodiales'='darkorange' ,'Unclassified Sclerotiniaceae'='#82e4de', 'Unclassified Nectriaceae'='#547dae', "Vishniacozyma"="light green", "Bullera"="gray", "Mycosphaerella"="#97b8c2", "Phaeosphaeria"="brown" )


Core_Flex_ITS_taxo_core_ag$Genus<-ordered(Core_Flex_ITS_taxo_core_ag$Genus, levels=c('Bensingtonia','Bulleromyces','Dioszegia','Fusarium','Gibberella','Holtermanniella','Mortierella','Itersonilia','Plectosphaerella', "Sporobolomyces", "Stemphylium", "Acremonium",'Unclassified Basidiomycota','Unclassified Sclerotiniaceae', 'Unclassified Nectriaceae', "Bullera", "Mycosphaerella", "Phaeosphaeria",'Unclassified Capnodiales', "Vishniacozyma",'Filobasidium', 'Alternaria', 'Cladosporium'))
pITS_core_flex_tax=ggplot(data=Core_Flex_ITS_taxo_core_ag, aes(x=Plant, y=SV_Rel_Abund,fill=Genus))+ geom_bar(aes(), stat="identity", position="stack")+theme_classic(base_size = 14)+ylab("ASV Relative Abundance (%)")+xlab("Plant Species")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 8, face = "bold"))+ theme(legend.title = element_text(colour="black", size=10, face="bold"))+scale_fill_manual(values=tax_colors_ITS)+guides(fill=guide_legend(nrow=5))+theme(legend.position = "bottom", legend.spacing.x = unit(0.05, 'cm'), 
        legend.spacing.y = unit(0.05, 'cm'), 
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.justification="left") + theme(axis.title.y = element_blank())
pITS_core_flex_tax
```

# Final Figure 7 combined
```{r}
library(ggpubr)
figure7=ggarrange(p16S_core_flex,pgyrB_core_flex,pITS_core_flex,p16S_core_flex2,pgyrB_core_flex2,pITS_core_flex2,p16S_core_flex_tax,pgyrB_core_flex_tax,pITS_core_flex_tax,labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), ncol = 3, nrow=3, heights = c(1, 1, 1.5))
figure7
```





