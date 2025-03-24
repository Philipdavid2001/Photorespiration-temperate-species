

library(ggplot2)
library(stats)
library(base)
library(dplyr)
library(ggplot2)




df <- read.csv("Compiled-literature-rates.csv", header = T, stringsAsFactors = T, sep = ";")




svg("PR-lit.svg", width = 4.5, height = 4)
ggplot(df, aes(x = cat.1, y=phi, fill = duration)) + 
  geom_boxplot(aes(x = cat.1, y=phi, fill = duration), 
               outlier.colour="darkred", 
               outlier.shape=21,
               outlier.size=3, 
               outlier.alpha = 0.5) + theme_bw() +
  geom_point(size = 2, stroke = 0.85, aes(color = genus, shape = genus)) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,14,15,16)) + ylim(0, 1) +
dev.off()




###Our values


setwd("C:/Users/Phili/Desktop/Github/Photorespiration-temperate-species/output")

df <- read.csv("species-output2.csv", header = T, stringsAsFactors = T, sep = ";")


percent.table <- df %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.percent.avg = mean(pr.percent),
            se = sd(pr.percent)/ sqrt(length(pr.percent)))

svg("PR-our.svg", width = 4.5, height = 4)
ggplot(percent.table, aes(x = 1 , y=pr.percent.avg)) + 
  geom_boxplot(aes(x = 1 , y=pr.percent.avg), 
               outlier.colour="darkred", 
               outlier.shape=21,
               outlier.size=3, 
               outlier.alpha = 0.5, fill = "cornflowerblue")+ theme_bw() +
  geom_point(inherit.aes = F, size = 2, stroke = 0.85, aes(x = 1, y = pr.percent.avg, color = sp, shape = sp), position = position_dodge(width=0.5)) +
  scale_x_discrete(expand = c(0.5,0.5))+
  scale_y_continuous(limits = c(0,1))+
  scale_shape_manual(values = c(9,10,19,12,13,17,18))+
  labs(x = "Temperate Tree species", y = "", color = "Species", shape = "Species")+
  theme(legend.text = element_text(family = "serif", face = "italic", size = 10) )
dev.off()
    


### Beans!

setwd("C:/Users/Phili/Desktop/Github/Photorespiration-temperate-species/Data/Coolbeans")

df <- read.csv("BeansRP.csv", header = T, stringsAsFactors = T, sep = ";")


percent.tableBean <- df %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.percent.avg = mean(pr.percent),
            se = sd(pr.percent)/ sqrt(length(pr.percent)))

svg("PR-our.svg", width = 4.5, height = 4)
ggplot(percent.tableBean, aes(x = 1 , y=pr.percent.avg)) + 
  geom_boxplot(aes(x = 1 , y=pr.percent.avg), 
               outlier.colour="darkred", 
               outlier.shape=21,
               outlier.size=3, 
               outlier.alpha = 0.5, fill = "cornflowerblue")+ theme_bw() +
  geom_point(inherit.aes = F, size = 2, stroke = 0.85, aes(x = 1, y = pr.percent.avg, color = sp, shape = sp), position = position_dodge(width=0.5)) +
  scale_x_discrete(expand = c(0.5,0.5))+
  scale_y_continuous(limits = c(0,1))+
  scale_shape_manual(values = c(9,10,19))+
  labs(x = "Temperate Tree species", y = "", color = "Species", shape = "Species")+
  theme(legend.text = element_text(family = "serif", face = "italic", size = 10) )
dev.off()



