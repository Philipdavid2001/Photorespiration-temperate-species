library(dplyr)
library(ggplot2)
library(nlme)
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
                 outlier.alpha = 1, fill = "white")+ theme_minimal() +
    geom_point(inherit.aes = F, size = 2, stroke = 0.85, aes(x = pc, y = phi, color = species, shape = species), position = position_dodge(width=0.5)) +
    scale_y_continuous(limits = c(0,1))+
    scale_shape_manual(values = c(9,10,19,12,13,17,18, 20,21,22,23,24,25, 8)) +
    labs(x = "", y = "Phi", color = "Species", shape = "Species")+
    theme(legend.text = element_text(family = "serif", face = "italic", size = 10) ) + coord_flip()

--
    
    
    mod4 <- nlme::lme(pr.real ~  sp * setTleaf, data = outs, 
                      random = ~1|treeid, 
                      method = "REML", 
                      na.action=na.omit) ; anova(mod4)


adf <- read.csv("~/Documents/GitHub/Photorespiration-temperate-species/output/species-output-annotated.csv",
                comment.char = "#" )


adf <- subset(adf, ETR.delta >=-150)

adf <- subset(adf, pr.percent <=1.1)


ggplot(adf, aes(y = ETR.delta, x=pr.percent, color = sp)) +
    geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
    scale_x_continuous(limits = c(0, 1.2), 
                       name = expression(pr.percent)) + 
    ggthemes::theme_base() +
    geom_smooth(se = T, method="lm", col = "grey80") +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          panel.border = element_rect(color = "grey70")) +
    scale_shape_manual(values = c(21,22,23,24,25,1,2))+
    facet_wrap(~tleaf)


ggplot(adf, aes(y = ETR.delta, x=pr, color = sp)) +
    geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
    scale_x_continuous(limits = c(-1, 15), 
                       name = expression(pr)) + 
    ggthemes::theme_base() +
    geom_smooth(se = T, method="lm", col = "grey80") +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          panel.border = element_rect(color = "grey70")) +
    scale_shape_manual(values = c(21,22,23,24,25,1,2))+
    facet_wrap(~tleaf)







pt <- adf %>%
    group_by(sp, tleaf) %>%
    summarise(prp = mean(pr.percent),
              prpse = sd(pr.percent)/sqrt(length(pr.percent)),
              dETR = mean(ETR.delta),
              dETRse = sd(ETR.delta)/ sqrt(length(ETR.delta)))



mod4 <- nlme::lme(pr.percent ~  ETR.delta * tleaf, data = adf, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(pr.percent ~  ETR.percent * tleaf, data = adf, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(ETR.delta ~ pr.percent + tleaf + sp, data = adf, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)


mod4 <- nlme::lme(ETR.delta ~ pr.percent * tleaf * sp, data = adf, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

pt <- subset(pt,  sp != "Scandosorbus intermedia")

svg(filename = "dETR-pr.svg", width = 16, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = dETR, x=prp)) +
    geom_point(size = 3, stroke = 0.85, aes(color = tleaf, shape = sp)) +
    scale_x_continuous(limits = c(0, 1.1), 
                       name = expression(phi == R[p] / A[net])) + 
    scale_y_continuous(limits = c(-65, 15), 
                       name = expression(ETR[0*"%"] - ETR[21*"%"])) + 
    ggthemes::theme_base() +
    geom_smooth(se = F, method="lm", col = "red") +
    geom_errorbar(aes(ymin=dETR-dETRse, ymax=dETR+dETRse), width=.2) +  
    geom_errorbarh(aes(xmin=prp-prpse, xmax=prp+prpse), width=.2) +  
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          legend.position = "none", 
          panel.border = element_rect(color = "grey70")) +
    scale_shape_manual(values = c(21,22,23,24,25,1,2))+
    facet_wrap(~ sp, nrow = 1)
dev.off()


pt$tleaf <- as.factor(pt$tleaf)

svg(filename = "dETR-pr.svg", width = 14, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = dETR, x=prp)) +
    geom_point(size = 3, stroke = 1, aes(color = tleaf), pch = 21) +
    scale_x_continuous(limits = c(0, 1.1), 
                       name = expression(phi == R[p] / A[net])) + 
    scale_y_continuous(limits = c(-65, 15), 
                       name = expression(ETR[0*"%"] - ETR[21*"%"])) + 
    ggthemes::theme_base() +
    geom_smooth(se = F, method="lm", col = "red") +
    geom_errorbar(aes(ymin=dETR-dETRse, ymax=dETR+dETRse), width=.2, col = "grey") +  
    geom_errorbarh(aes(xmin=prp-prpse, xmax=prp+prpse), col = "grey") +  
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          legend.position = "none", 
          panel.border = element_rect(color = "grey70")) +
    scale_color_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20"))+
    geom_point(size = 3, stroke = 1.5, aes(color = tleaf), pch = 21) +
    facet_wrap(~ sp, nrow = 1)

dev.off()
