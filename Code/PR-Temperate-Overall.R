library(dplyr)
library(ggplot2)
library(nlme)
library(patchwork)


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

adf$tlf <- as.factor(adf$tleaf)


# adf <- subset(adf, ETR.delta >=-150)
# adf <- subset(adf, pr.percent <=1.1)
###

  geom_point(size = 3, stroke = 0.85, aes(color = tleaf, shape = sp)) +
    
ggplot(adf, aes(y = pr.percent, x=tlf)) +
  scale_y_continuous(limits = c(0, 1.75), 
                     name = expression(pr)) +
  geom_boxplot(aes()) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  ggthemes::theme_base() +
  geom_smooth(se = T, method="lm", col = "grey80") +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp)

ggplot(adf, aes(y = pr.percent, x = tleaf)) +
  scale_y_continuous(limits = c(0, 1.75), 
                     name = expression(pr)) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), se = TRUE, method = "lm", color = "grey80") + # Smooth line
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp)



ggplot(adf, aes(y = pr.percent, x = tleaf)) +
  scale_y_continuous(limits = c(0, 1.75), 
                     name = expression(pr)) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              se = TRUE, 
              method = "lm", 
              color = "grey0", 
              fill = NA) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp)




adf <- subset(adf,  sp != "Scandosorbus intermedia")


adf$sp <- factor(adf$sp, levels = c("Betula pendula", 
                                  "Fagus sylvatica","Betula pubescens", "Acer platanoides", "Tilia cordata", "Corylus avellana"))

p3 <- ggplot(adf, aes(y = pr.percent, x = tleaf)) +
  scale_y_continuous(limits = c(0, 1.75), 
                     name = expression(phi == R[p] / A[net])) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(pr)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p3



  
  
  

p1 <- ggplot(adf, aes(y = anet.21p, x = tleaf)) +
  scale_y_continuous(limits = c(0, 25), 
                     name = expression(paste(A[net]*" ("~mu~mol~m^{-2}~s^{-1}~")"))) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(paste(T[leaf]~"("~degree~"C)"))) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p1


p1b <- ggplot(adf, aes(y = anet.0p, x = tleaf)) +
  scale_y_continuous(limits = c(0, 30), 
                     name = expression(paste(A[net]*" ("~mu~mol~m^{-2}~s^{-1}~")"))) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(paste(T[leaf]~"("~degree~"C)"))) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p1b




p2 <- ggplot(adf, aes(y = pr, x = tleaf)) +
  scale_y_continuous(limits = c(0, 25), 
                     name = expression(paste(R[p]*" ("~mu~mol~m^{-2}~s^{-1}~")"))) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(tleaf)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p2


###





combined_plot <- p1 / p2 / p3 
svg(filename = "Temperate-PR-boxplot.svg", width = 10, height = 12, bg = "transparent")
combined_plot
dev.off()

combined_plot <- p1 / p1b / p2 
svg(filename = "Temperate-PR-boxplot-anet.svg", width = 10, height = 12, bg = "transparent")
combined_plot
dev.off()




p4 <- ggplot(adf, aes(y = -ETR.delta, x = tleaf)) +
  scale_y_continuous(limits = c(-20, 100), 
                     name = expression("|"~ETR[0*"%"] - ETR[21*"%"]~"|")) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(tleaf)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p4



p5 <- ggplot(adf, aes(y = ETR.21p, x = tleaf)) +
  scale_y_continuous(limits = c(-20, 250), 
                     name =expression(ETR["21%"], ' (', mu * 'mol m'^{-2}*'s'^{-1}*')')) + 
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(tleaf)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p5






p6 <- ggplot(adf, aes(y = ETR.0p, x = tleaf)) +
  scale_y_continuous(limits = c(-20, 250), 
                     name =expression(ETR["0%"], ' (', mu * 'mol m'^{-2}*'s'^{-1}*')')) + 
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(tleaf)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p6



combined_plot <- p5 / p6 / p4 
svg(filename = "Temperate-PR-ETR-boxplot.svg", width = 10, height = 12, bg = "transparent")
combined_plot
dev.off()




p6 <- ggplot(adf, aes(y = gsw.21p, x = tleaf)) +
  scale_y_continuous(limits = c(0, .8), 
                     name =expression(g["s @21%"], ' (', mu * 'mol m'^{-2}*'s'^{-1}*')')) + 
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(tleaf)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p6





p7 <- ggplot(adf, aes(y = gsw.0p, x = tleaf)) +
  scale_y_continuous(limits = c(0, .8), 
                     name =expression(g["s @0%"], ' (', mu * 'mol m'^{-2}*'s'^{-1}*')')) + 
  geom_boxplot(aes(group = cut_width(tleaf, 1))) + # Grouping for boxplot
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(tleaf)) +
  ggthemes::theme_base() +
  geom_smooth(aes(group = 1), 
              method = "lm", 
              se = T, 
              color = "grey20") + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  facet_wrap(~sp, nrow = 1); p7













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


pt$tleaf <- as.factor(pt$tleaf)

svg(filename = "dETR-pr.svg", width = 16, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = -ETR.delta, x=prp)) +
  geom_point(size = 3, stroke = 0.85, aes(color = tleaf), pch = 21) +
  scale_x_continuous(limits = c(0, 1.2), 
                     name = expression(phi == R[p] / A[net]),
                     labels = scales::number_format(accuracy = 0.1)) + 
  scale_y_continuous(limits = c(-10,70), 
                     name = expression("|"~ETR[0*"%"] - ETR[21*"%"]~"|")) + 
  ggthemes::theme_base() +
  geom_smooth(se = F, method="lm", col = "red") +
  geom_errorbar(aes(ymin=-ETR.delta-ETR.delta.se, ymax=-ETR.delta+ETR.delta.se), width=.2) +  
  geom_errorbarh(aes(xmin=prp-prp.se, xmax=prp+prp.se), width=.2) +  
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "none", 
        panel.border = element_rect(color = "grey70")) +
  scale_color_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20"))+
  geom_point(size = 3, stroke = 1.5, aes(color = tleaf), pch = 21) +
  facet_wrap(~ sp, nrow = 1)

dev.off()


pt$tlfac <- as.factor(pt$tleaf)


pt <- read.csv("pr-temperate-2024-means-se.csv", stringsAsFactors = T)




ggplot(pt, aes(y = anet21p, x=tleaf)) + 
  geom_point(size = 3, stroke = 0.85, pch = 21, aes(color = tlfac)) + 
  scale_x_continuous(limits = c(20, 40), 
                     name = expression(paste(italic(T[leaf])~~"("~degree~"C)")), 
                     labels = scales::number_format(accuracy = 5)) + 
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(A[net])*" ("~mu~mol[CO[2]]~m^{-2}~s^{-1}~")"))) + 
  ggthemes::theme_base() + 
  geom_smooth(se = F, method="lm", col = "red") + 
  geom_errorbar(aes(ymin=anet21p-anet21p.se, ymax=anet21p+anet21p.se), width=.2) + 
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        legend.position = "none", 
        panel.border = element_rect(color = "grey70")) + 
  scale_color_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20")) + 
  geom_point(size = 3, stroke = 1.5, aes(color = tlfac), pch = 21) + 
  facet_wrap(~ sp, nrow = 1)

anet21p_plot <- ggplot(pt, aes(y = anet21p, x = tleaf)) + 
  geom_point(size = 3, stroke = 0.85, pch = 21, aes(color = tlfac, fill = tlfac)) + 
  scale_x_continuous(
    limits = c(20, 40), 
    name = expression(paste(italic(T[leaf])~~"("~degree~"C)")), 
    labels = scales::number_format(accuracy = 5)
  ) + 
  scale_y_continuous(
    limits = c(0, 25), 
    name = expression(paste(italic(A[net])*" ("~mu~mol[CO[2]]~m^{-2}~s^{-1}~")"))
  ) + 
  ggthemes::theme_base() + 
  geom_smooth(se = F, method = "lm", col = "red") + 
  geom_errorbar(aes(ymin = anet21p - anet21p.se, ymax = anet21p + anet21p.se), width = 0.2) + 
  theme(
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(size = 15), 
    legend.position = "none", 
    panel.border = element_rect(color = "grey70")
  ) + 
  scale_color_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20")) + 
  scale_fill_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20")) + 
  
  facet_wrap(~ sp, nrow = 1); anet21p_plot


anet0p_plot <- ggplot(pt, aes(y = anet0p, x = tleaf)) + 
  geom_point(size = 3, stroke = 0.85, pch = 21, aes(color = tlfac)) + 
  scale_x_continuous(
    limits = c(20, 40), 
    name = expression(paste(italic(T[leaf])~~"("~degree~"C)")), 
    labels = scales::number_format(accuracy = 5)
  ) + 
  scale_y_continuous(
    limits = c(0, 25), 
    name = expression(paste(italic(A[net])*" ("~mu~mol[CO[2]]~m^{-2}~s^{-1}~")"))
  ) + 
  ggthemes::theme_base() + 
  geom_smooth(se = F, method = "lm", col = "blue") + 
  geom_errorbar(aes(ymin = anet0p - anet0p.se, ymax = anet0p + anet0p.se), width = 0.2) + 
  theme(
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(size = 15), 
    legend.position = "none", 
    panel.border = element_rect(color = "grey70")
  ) + 
  scale_color_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20")) + 
  facet_wrap(~ sp, nrow = 1); anet0p_plot



pt <- subset(pt,  sp != "Scandosorbus intermedia")


pt$sp <- factor(pt$sp, levels = c("Betula pendula", 
  "Fagus sylvatica","Betula pubescens", "Acer platanoides", "Tilia cordata", "Corylus avellana"))





p1 <- ggplot(pt, aes(x = tleaf)) +
  geom_point(aes(y = anet21p, color = "anet21p"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = anet21p - anet21p.se, ymax = anet21p + anet21p.se, color = "anet21p"), width = 0.2) +
  geom_smooth(aes(y = anet21p, color = "anet21p"), se = FALSE, method = "lm") +
  
  geom_point(aes(y = anet0p, color = "anet0p"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = anet0p - anet0p.se, ymax = anet0p + anet0p.se, color = "anet0p"), width = 0.2) +
  geom_smooth(aes(y = anet0p, color = "anet0p"), se = FALSE, method = "lm") +

  geom_point(aes(y = pr, color = "pr"), size = 3, stroke = 0.85, pch = 20) +
  geom_errorbar(aes(ymin = pr - pr.se, ymax = pr + pr.se, color = "pr"), width = 0.2) +
  geom_smooth(aes(y = pr, color = "pr"), se = FALSE, method = "lm") +

    
  scale_x_continuous(
    limits = c(20, 40), 
    name = expression(paste(italic(T[leaf])~~"("~degree~"C)")), 
    labels = scales::number_format(accuracy = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 25), 
    name = expression(paste(A[net]*" ("~mu~mol[CO[2]]~m^{-2}~s^{-1}~")"))
  ) +
  
  scale_color_manual(
    name = "Dataset",
    values = c("anet21p" = "black", "anet0p" = "grey", "pr" = "blue"), # Assign colors to each dataset
    labels = c("anet21p", "anet0p", "pr")
  ) +
  ggthemes::theme_base() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    legend.position = "bottom", 
    panel.border = element_rect(color = "grey70")
  ) +
  facet_wrap(~ sp, nrow = 1); p1


p2 <- ggplot(pt, aes(x = tleaf)) +
  geom_point(aes(y = ETR.21p, color = "ETR.21p"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = ETR.21p - ETR.21p.se, ymax = ETR.21p + ETR.21p.se, color = "ETR.21p"), width = 0.2) +
  geom_smooth(aes(y = ETR.21p, color = "ETR.21p"), se = FALSE, method = "lm") +
  
  geom_point(aes(y = ETR.0p, color = "ETR.0p"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = ETR.0p - ETR.0p.se, ymax = ETR.0p + ETR.0p.se, color = "ETR.0p"), width = 0.2) +
  geom_smooth(aes(y = ETR.0p, color = "ETR.0p"), se = FALSE, method = "lm") +
  
  scale_x_continuous(
    limits = c(20, 40), 
    name = expression(paste(italic(T[leaf])~~"("~degree~"C)")), 
    labels = scales::number_format(accuracy = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 170), 
    name = expression(paste(ETR*" ("~mu~mol~e^{"-"}~m^{-2}~s^{-1}~")"))) +
  
  scale_color_manual(
    name = "Dataset",
    values = c("ETR.21p" = "black", "ETR.0p" = "grey"), # Assign colors to each dataset
    labels = c("ETR.21p", "ETR.0p")
  ) +
  ggthemes::theme_base() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    legend.position = "bottom", 
    panel.border = element_rect(color = "grey70")
  ) +
  facet_wrap(~ sp, nrow = 1); p2



p3 <- ggplot(pt, aes(x = tleaf)) +
  geom_point(aes(y = gsw.21p, color = "gsw.21p"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = gsw.21p - gsw.21p.se, ymax = gsw.21p + gsw.21p.se, color = "gsw.21p"), width = 0.2) +
  geom_smooth(aes(y = gsw.21p, color = "gsw.21p"), se = FALSE, method = "lm") +
  
  geom_point(aes(y = gsw.0p, color = "gsw.0p"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = gsw.0p - gsw.0p.se, ymax = gsw.0p + gsw.0p.se, color = "gsw.0p"), width = 0.2) +
  geom_smooth(aes(y = gsw.0p, color = "gsw.0p"), se = FALSE, method = "lm") +
  
  scale_x_continuous(
    limits = c(20, 40), 
    name = expression(paste(T[leaf]~"("~degree~"C)")), 
    labels = scales::number_format(accuracy = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 0.4), 
    name = expression(italic(g[s])~mol[H[2]*O]~m^{-2}~s^{-1})) +
  
  scale_color_manual(
    name = "Dataset",
    values = c("gsw.21p" = "black", "gsw.0p" = "grey"), # Assign colors to each dataset
    labels = c("gsw.21p", "gsw.0p")
  ) +
  ggthemes::theme_base() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    legend.position = "bottom", 
    panel.border = element_rect(color = "grey70")
  ) +
  facet_wrap(~ sp, nrow = 1); p3




p4 <- ggplot(pt, aes(x = tleaf)) +
  geom_point(aes(y = prp, color = "prp"), size = 3, stroke = 0.85, pch = 21) +
  geom_errorbar(aes(ymin = prp - prp.se, ymax = prp + prp.se, color = "prp"), width = 0.2) +
  geom_smooth(aes(y = prp, color = "prp"), se = FALSE, method = "lm") +
  

  scale_x_continuous(
    limits = c(20, 40), 
    name = expression(paste(T[leaf]~"("~degree~"C)")), 
    labels = scales::number_format(accuracy = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 1.2), 
    name = expression(phi == R[p] / A[net])) +
  
  scale_color_manual(
    name = "Dataset",
    values = c("prp" = "black"), # Assign colors to each dataset
    labels = c("prp")
  ) +
  ggthemes::theme_base() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    legend.position = "bottom", 
    panel.border = element_rect(color = "grey70")
  ) +
  facet_wrap(~ sp, nrow = 1); p4


p5 <- ggplot(pt, aes(y = -ETR.delta, x=prp)) +
  geom_point(size = 3, stroke = 0.85, aes(color = tlfac), pch = 21) +
  scale_x_continuous(limits = c(0, 1.2), 
                     name = expression(phi == R[p] / A[net]),
                     labels = scales::number_format(accuracy = 0.1)) + 
  scale_y_continuous(limits = c(-10,70), 
                     name = expression("|"~ETR[0*"%"] - ETR[21*"%"]~"|")) + 
  ggthemes::theme_base() +
  geom_smooth(se = F, method="lm", col = "red") +
  geom_errorbar(aes(ymin=-ETR.delta-ETR.delta.se, ymax=-ETR.delta+ETR.delta.se), width=.2) +  
  geom_errorbarh(aes(xmin=prp-prp.se, xmax=prp+prp.se), width=.2) +  
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "none", 
        panel.border = element_rect(color = "grey70")) +
  scale_color_manual(values = c("25" = "grey85", "30" = "grey60", "35" = "grey20"))+
  geom_point(size = 3, stroke = 1.5, aes(color = tlfac), pch = 21) +
  facet_wrap(~ sp, nrow = 1)





library(patchwork)
combined_plot <- p1 / p2 / p3 / p4 / p5



svg(filename = "Temperate-combined.svg", width = 14, height = 18, bg = "transparent")
combined_plot
dev.off()





