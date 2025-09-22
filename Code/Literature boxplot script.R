

library(ggplot2)
library(stats)
library(base)
library(dplyr)





df <- read.csv("Data/Literature data/Compiled-literature-ALL-DATA.csv", header = T, stringsAsFactors = T, sep = ",")


df <- subset(df,tl2==1)
# excludes our data
df$cat <- factor(df$cat, levels = c("Short",  "Tree", "Current 25", "Current 30", "Current 35"))

# 21 Sept 2025 box figure
svg("PR-lit.svg", width = 8, height = 5)

ggplot(df, aes(x = cat, y = phi, fill = cat)) + 
  geom_boxplot(
    outlier.colour = "darkred", 
    outlier.shape = 21,
    outlier.size = 3, 
    outlier.alpha = 0
  ) +
  geom_point(
    aes(color = species, shape = species),
    position = position_jitter(width = 0.2, height = 0),
    size = 3, 
    stroke = 0.85, 
    alpha = 0.8
  ) +
  scale_shape_manual(values = c(0,0, 1, 2,0,4,5,1, 2, 4, 1, 5, 6, 20, 2,6, 4, 20, 23,  5)) + 
scale_color_manual(values = c("red","blue", "red", "red","darkgreen","red","red","blue", "blue", "blue", "darkgreen", "blue", "blue", "blue", "darkgreen","red", "darkgreen", "red", "blue",  "darkgreen")) + 
  scale_fill_manual(values = c("Short" = "#86EBFD", 
                               "Tree" = "#51FB6E", 
                               "Current 25" = "grey80", 
                               "Current 30" = "grey60", 
                               "Current 35" = "grey40")) +
  xlab("") +
  ylim(0, 1.3) +
  ylab(expression(Phi(R[p] / A[net]))) +
  ggthemes::theme_base() +
  theme(
    legend.text = element_text(family = "Times New Roman", face = "italic"),
    legend.title = element_text(family = "Times New Roman", face = "italic")
  ) + coord_flip()

dev.off()





se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

phi_summary <- df %>%
  group_by(cat) %>%
  summarise(
    mean_phi   = mean(phi, na.rm = TRUE),
    median_phi = median(phi, na.rm = TRUE),
    se_phi     = se(phi)
  ) %>%
  arrange(factor(cat, levels = c("Short", "Tree", "Current 25", "Current 30", "Current 35")))
print(phi_summary)





###Our values


setwd("C:/Users/Phili/Desktop/Github/Photorespiration-temperate-species/Data/Literature data/Final Literature values")

###Species
df <- read.csv("species-output2.csv", header = T, stringsAsFactors = T, sep = ";")

###Ecotypes
df <- read.csv("ecotypes-output.csv", header = T, stringsAsFactors = T, sep = ";")


percent.table <- outs %>%
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
  labs(x = "Temperate Tree Species", y = "", color = "Species", shape = "Species")+
  theme(legend.text = element_text(family = "serif", face = "italic", size = 10) )
dev.off()
    



#### ecotype

svg("PR-our-ecotypes.svg", width = 4.5, height = 4)
ggplot(percent.table, aes(x = 1 , y=pr.percent.avg)) + 
  geom_boxplot(aes(x = 1 , y=pr.percent.avg), 
               outlier.colour="darkred", 
               outlier.shape=21,
               outlier.size=3, 
               outlier.alpha = 0.5, fill = "cornflowerblue")+ theme_bw() +
  geom_point(inherit.aes = F, size = 2, stroke = 0.85, aes(x = 1, y = pr.percent.avg, color = sp, shape = sp), position = position_dodge(width=0.5)) +
  scale_x_discrete(expand = c(0.5,0.5))+
  scale_y_continuous(limits = c(0,1))+
  scale_shape_manual(values = c(20,21,22,23,24,25,-1))+
  labs(x = "Ecotypes", y = "", color = "Species", shape = "Species")+
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



