library(dplyr)
library(ggplot2)

setwd("~/Documents/GitHub/Photorespiration-temperate-species/Data/Literature data/Final Literature values")


df <- read.csv("Compiled-literature-ALL-DATA.csv", stringsAsFactors = T, header = T)


df <- subset(df, tleaf == 25)


boxplot(df$phi ~ df$duration)

boxplot(df$phi ~ df$pc)

percent.table <- df %>%
    group_by(pc) %>%
    summarise(pr.percent.avg = mean(phi),
              se = sd(phi)/ sqrt(length(phi)))



ggplot(df, aes(x = pc , y=phi)) + 
    geom_boxplot(aes(x = pc, y=phi), 
                 outlier.colour="darkred", 
                 outlier.shape=21,
                 outlier.size=3, 
                 outlier.alpha = 1, fill = "white")+ theme_bw() +
    geom_point(inherit.aes = F, size = 2, stroke = 0.85, aes(x = pc, y = phi, color = species, shape = species), position = position_dodge(width=0.5)) 

+
    scale_y_continuous(limits = c(0,1))+
    scale_shape_manual(values = c(9,10,19,12,13,17,18))+
    labs(x = "Temperate Tree Species", y = "", color = "Species", shape = "Species")+
    theme(legend.text = element_text(family = "serif", face = "italic", size = 10) )

