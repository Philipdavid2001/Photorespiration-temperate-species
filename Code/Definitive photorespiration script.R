
library(stats)
library(base)
library(dplyr)
library(patchwork)
library(lme4)
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(ggpubr)


# Species #######################
df <- read.csv("Data/Photorespiration temperate/Uppsala-2024-Summer-Photorespiration-SpotMes-TreeSpp.csv", header = T, stringsAsFactors = T, sep = ";")


###### Photorespiration rate calculation loop -----------------------------------------------------

## use the following filters as necessary
df              <- subset(df, dq == "yes") # eliminates unreliable values.
df$pairid       <- paste0(df$sprp, "-", df$setTleaf)
df$pairid       <- as.factor(df$pairid)

length(unique(df$pairid))       # individual runs + season
dflist          <- split(df,list(df$pairid))
dflist          <- split(df,unique(list(df$pairid)))
dflist          <- dflist[sapply(dflist, nrow)>0] 

setwd("output")
# tweak the code below for photorespiration rate calculation. 
correct_RD <- function(data, output_path){
  ColumnNames <- c("sp","treeid","setTleaf",
             "anet.21p",    
             "anet.0p",
             "anet.delta",
             "pr.CO2", 
             "pr.real",
             "pr.percent",    
             "gsw.21p",    
             "gsw.0p",     
             "gsw.delta",  
             "gsw.percent", 
             "E.21p",     
             "E.0p",      
             "E.delta",   
             "E.percent", 
             "ETR.21p",     
             "ETR.0p",
             "ETR.delta",   
             "ETR.percent",
             "JT",
             "JO1",
             "JC1",
             "JO2",
             "JC2",
             "Rp.ETR",
             "JO1.percent",
             "JC1.percent",
             "Rp.ETR.percent",
             "JC2.percent",
             "NPQ.21p",
             "NPQ.0p",
             "NPQ.delta",
             "NPQ.percent")
  
  output_data <- data.frame(matrix(nrow = 0, ncol = length(ColumnNames)))
  names(output_data) <- ColumnNames

  for(dataIdx in 1:length(data)){
    
    speciesdata <- data[[dataIdx]]
    if(nrow(speciesdata) != 2){
      print(names(dflist)[dataIdx])
      print(nrow(speciesdata))
      next
    }
    
    sp               <-     speciesdata$sp[1]
    treeid           <-     speciesdata$treeid[1]
    setTleaf         <-     speciesdata$setTleaf[1]
    p21              <-     subset(speciesdata, olev == "21p")    
    p0               <-     subset(speciesdata, olev == "0p")
    anet.21p         <-     p21$A
    anet.0p          <-     p0$A
     
    
    # Correcting for dark respiration - NOT DONE
    # Since we assume mitochondrial respiration under ambient and non-photorespiratory conditions are assumed constant - we will not correct for Rd. hence, this part of the code is ignored. 
    # In case used: The values for slope and intercept calculated using literature values. 
    #Rd               <-      0.05554*anet.21p+ 0.11395 
    #corranet         <-      anet.21p-Rd          
    # "corranet" is Anet corrected for Rd
    
    anet.delta       <-      anet.0p - anet.21p
    ## uses Walker 2017 temperature function for lambda
    ## lambda = 0.389 + 0.00876  Tleaf 
    ## to use default lambda: pr.CO2 <- anet.delta * 0.5
  
    pr.CO2           <-      anet.delta * (0.38926+(0.008765*setTleaf))
    pr.real          <-      anet.delta + pr.CO2
    pr.percent       <-      pr.real/anet.21p
    
    gsw.21p          <-      p21$gsw
    gsw.0p           <-      p0$gsw
    gsw.delta        <-      p0$gsw - p21$gsw
    gsw.percent      <-      gsw.delta/p21$gsw
    
    E.21p            <-      p21$E
    E.0p             <-      p0$E
    E.delta          <-      p0$E - p21$E
    E.percent        <-      E.delta/p21$E
    
    
    ETR.21p          <-      p21$ETR
    ETR.0p           <-      p0$ETR
    ETR.delta        <-      p0$ETR - p21$ETR
    ETR.percent      <-      ETR.delta/p21$ETR
    JT               <-      ETR.21p + ETR.0p
    JO1              <-      pr.percent * JT
    JC1              <-      JT-JO1
    JO2              <-      ((2/3)*(JT-4*(anet.21p)))
    JC2              <-      ((1/3)*(JT+8*(anet.21p)))
    Rp.ETR           <-      ((1/12)*(JT-4*(anet.21p)))
    JO1.percent      <-      JO1/JT
    JC1.percent      <-      JC1/JT
    Rp.ETR.percent   <-      Rp.ETR/JT
    JC2.percent      <-      JC2/JT
    
    NPQ.21p          <-      p21$NPQ
    NPQ.0p           <-      p0$NPQ
    NPQ.delta        <-      NPQ.0p - NPQ.21p
    NPQ.percent      <-      NPQ.delta/NPQ.21p
    
    new_data = data.frame(as.character(sp),
                          as.character(treeid),
                          as.character(setTleaf), 
                          anet.21p,    
                          anet.0p,     
                          anet.delta,
                          pr.CO2,
                          pr.real,
                          pr.percent,
                          gsw.21p,    
                          gsw.0p,     
                          gsw.delta,  
                          gsw.percent,
                          E.21p,     
                          E.0p,      
                          E.delta,   
                          E.percent, 
                          ETR.21p,     
                          ETR.0p,   
                          ETR.delta,   
                          ETR.percent,
                          JT,
                          JO1,
                          JC1,
                          JO2,
                          Rp.ETR,
                          JC2,
                          JO1.percent,
                          JC1.percent,
                          Rp.ETR.percent,
                          JC2.percent,
                          NPQ.21p,
                          NPQ.0p,
                          NPQ.delta, 
                          NPQ.percent)
    
    names(new_data) <- names(output_data)
    output_data <- rbind(output_data, new_data)
    
    
  }
  
  file_path = paste(output_path, "Species-output.csv", sep = "/")
  if(!dir.exists(output_path)){
    dir.create(output_path, recursive = T)
  }
  
  
  index <- 1
  while(file.exists(file_path)){
    file_path <- paste(output_path, paste("Species-output", as.character(index), ".csv", sep =""), sep = "/")
    index = index + 1
  }
 
  

  write.table(output_data, file = file_path, 
              row.names = FALSE, col.names = T, sep = ";")
}
correct_RD(dflist, "./")


####### plotting output -------

outs <- read.csv("Species-output.csv", stringsAsFactors = T, sep = ";")


###Species order
outs$sp <- factor(outs$sp, levels = c(
  "Betula pendula",         
  "Fagus sylvatica",         
  "Betula pubescens",       
  "Acer platanoides", 
  "Tilia cordata",
  "Corylus avellana",       
  "Scandosorbus intermedia"))



###Ecotype order
# outs$sp <- factor(outs$sp, levels = c(
#   "Docksta",         
#   "Dunker",         
#   "Nobbele",       
#   "Ryninngsholm", 
#   "Skelleftea",
#   "Uppsala",       
#   "Ystad"))




# Making the table for the concatenated PR values for each species and each temperature point.





## Main Figure 1 ---------------------------------

## plots updateds to add anova p value for each panel


pval_df <- outs %>%
  group_by(sp) %>%
  filter(n_distinct(setTleaf) >= 2) %>%  # Ensure valid ANOVA
  summarise(
    pval = tryCatch({
      summary(aov(anet.21p ~ factor(setTleaf)))[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA)
  ) %>%
  mutate(label = paste0("P = ", signif(pval, 3)))

annot_df <- outs %>%
  group_by(sp) %>%
  summarise(ypos = max(anet.21p, na.rm = TRUE) + 1) %>%
  left_join(pval_df, by = "sp")

plot1 <- ggplot(outs, aes(x = setTleaf, y = anet.21p)) + 
  geom_boxplot(aes(group = factor(setTleaf)), fill = "#6666ff", alpha = .4, 
               outlier.colour = "blue", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, 
               fill = NA, color =  "grey30", stroke = 1) + 
  geom_smooth(method = "lm", col = "blue", se = FALSE) +  
  scale_y_continuous(limits = c(0, 25), 
                     name = expression(paste(italic(A)[net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  scale_x_continuous(limits = c(22, 37)) +  
  ggthemes::theme_base() +
  facet_wrap(~sp, ncol = 7) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    panel.border = element_rect(color = "grey70"),
    plot.margin = margin(5, 5, 5, 5),
    strip.text = element_blank()
  ) +
  geom_text(data = annot_df, aes(x = 30, y = 24, label = label), inherit.aes = FALSE, size = 5);  plot1


pval_df <- outs %>%
  group_by(sp) %>%
  filter(n_distinct(setTleaf) >= 2) %>%  # Ensure valid ANOVA
  summarise(
    pval = tryCatch({
      summary(aov(pr.real ~ factor(setTleaf)))[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA)
  ) %>%
  mutate(label = paste0("P = ", signif(pval, 3)))

annot_df <- outs %>%
  group_by(sp) %>%
  summarise(ypos = max(pr.real, na.rm = TRUE) + 1) %>%
  left_join(pval_df, by = "sp")


plot2 <- ggplot(outs, aes(x = setTleaf, y = pr.real)) + 
  geom_boxplot(aes(group = factor(setTleaf)), fill = "#FF9900", alpha = 0.7, 
               outlier.colour = "#CC0000", outlier.shape = 16) +
  geom_smooth(method = "lm", col = "#990000", se = FALSE) +  
  scale_y_continuous(limits = c(0, 25), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  stat_summary(fun = mean, geom = "point", shape = 23, 
               size = 2, fill = NA, color =  "grey30", stroke = 1) + 
  scale_x_continuous(limits = c(22, 37)) +  
  ggthemes::theme_base() +
  facet_wrap(~sp, ncol = 7) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    panel.border = element_rect(color = "grey70"),
    strip.text = element_text(family = "Times New Roman", face = "italic", size = 13)
  )  +
  geom_text(data = annot_df, aes(x = 30, y = 24, label = label), inherit.aes = FALSE, size = 5); plot2

pval_df <- outs %>%
  group_by(sp) %>%
  filter(n_distinct(setTleaf) >= 2) %>%  # Ensure valid ANOVA
  summarise(
    pval = tryCatch({
      summary(aov(pr.percent ~ factor(setTleaf)))[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA)
  ) %>%
  mutate(label = paste0("P = ", signif(pval, 3)))

annot_df <- outs %>%
  group_by(sp) %>%
  summarise(ypos = max(pr.percent, na.rm = TRUE) + 1) %>%
  left_join(pval_df, by = "sp")


plot3 <- ggplot(outs, aes(x = setTleaf, y = pr.percent)) + 
  geom_boxplot(aes(group = factor(setTleaf)), fill = "#FFFF99", alpha = 0.5, 
               outlier.colour = "grey30", outlier.shape = 16) +
  geom_smooth(method = "lm", col = "grey30", se = FALSE) +  
  stat_summary(fun = mean, geom = "point", shape = 23, 
               size = 2, fill = NA, color =  "grey30", stroke = 1) + 
  scale_y_continuous(limits = c(0, 2), 
                     breaks = seq(0, 2, 1),  
                     name = expression(italic(R)[p] / italic(A)[net])) +
  scale_x_continuous(limits = c(22, 37), name = "Leaf temperature (°C)") +  
  ggthemes::theme_base() +
  facet_wrap(~sp, ncol = 7) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    panel.border = element_rect(color = "grey70"),
    strip.text = element_blank()
  ) +
  geom_text(data = annot_df, aes(x = 30, y = 2, label = label), inherit.aes = FALSE, size = 5); plot3

svg(filename = "MAIN-Figure1a.svg", width = 10, height = 10, bg = "transparent")

(plot2 / plot1 / plot3) + plot_layout(nrow = 3)

dev.off()


###------- Mixed effects models: Photorespiration ----

mod4 <- nlme::lme(pr.real ~  sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(pr.real ~  setTleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(pr.real ~  sp * setTleaf, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

###------- Mixed effects models: Rp/Anet (Phi) ----


mod6 <- nlme::lme(pr.percent ~  sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)

mod6 <- nlme::lme(pr.percent ~  setTleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


mod6 <- nlme::lme(pr.percent ~  sp * setTleaf, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)

###------- Mixed effects models: Photosynthesis (Anet) ----

mod5 <- nlme::lme(anet.21p ~  sp, data = outs, 
                  random = ~1| treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)

mod5 <- nlme::lme(anet.21p ~  setTleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)


mod6 <- nlme::lme(anet.21p ~  sp * setTleaf, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


## exclude B. pubsecense and SCAINT 
range(outs$anet.21p)
outsx <- subset(outs, sp !="Scandosorbus intermedia")
outsx <- subset(outsx, sp !="Betula pubescens")


mod5 <- nlme::lme(pr.real ~  sp * setTleaf, data = outsx, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)

###------- Mixed effects models: Photosynthesis and Photorespiration ----

mod5 <- nlme::lme(pr.real ~   anet.21p, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)

mod4a <- nlme::lme(pr.real ~   anet.21p * setTleaf, data = outs, 
                   random = ~1|sp, 
                   method = "REML", 
                   na.action=na.omit) ; anova(mod4a)


mod4 <- nlme::lme(pr.real ~   anet.21p * setTleaf * sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

anova(mod4, mod4a)

### plotting averages ----

pt <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr = mean(pr.real),
            prse = sd(pr.real)/sqrt(length(pr.real)),
            ps = mean(anet.21p),
            psse = sd(anet.21p)/ sqrt(length(anet.21p)),
            phi = mean(pr.percent),
            phise = sd(pr.percent)/sqrt(length(pr.percent)),
            ETR.0p = mean(ETR.0p),
            ETR.0pse = sd(ETR.0p)/sqrt(length(ETR.0p)),
            ETR.21p = mean(ETR.21p),
            ETR.21pse = sd(ETR.21p)/ sqrt(length(ETR.21p)),
            ETR.delta = mean((ETR.21p-ETR.0p)),
            ETR.deltase = sd(ETR.delta)/sqrt(length(ETR.percent)),
            ETR.percent = mean(ETR.percent),
            ETR.percentse = sd(ETR.percent)/sqrt(length(ETR.percent))
  )


pt25 <- subset(outs, setTleaf == 25)
pt30 <- subset(outs, setTleaf == 30)
pt35 <- subset(outs, setTleaf == 35)

m25 <- nlme::lme(pr.real ~  anet.21p , data = pt25, 
                 random = ~1|sp, 
                 method = "REML", 
                 na.action=na.omit) ; anova(m25)

m30 <- nlme::lme(pr.real ~  anet.21p , data = pt30, 
                 random = ~1|sp, 
                 method = "REML", 
                 na.action=na.omit) ; anova(m30)

m35 <- nlme::lme(pr.real ~  anet.21p , data = pt35, 
                 random = ~1|sp, 
                 method = "REML", 
                 na.action=na.omit) ; anova(m35)

AIC(m25, m30, m35)  

model <- lmer(pr.real ~ anet.21p + setTleaf + (1 | sp), data = outs)
summary(model)

# pearson correlation coefficient
cor25 <- cor(pt25$anet.21p, pt25$pr.real, method = 'pearson'); round(cor25, 3) 
cor30 <- cor(pt30$anet.21p, pt30$pr.real, method = 'pearson'); round(cor30, 3)
cor35 <- cor(pt35$anet.21p, pt35$pr.real, method = 'pearson'); round(cor35, 3)

# lm r2 
rsq25 <- summary(lm(pr.real ~ anet.21p, pt25))["r.squared"]; round(rsq25$r.squared, 3) 
rsq30 <- summary(lm(pr.real ~ anet.21p, pt30))["r.squared"]; round(rsq30$r.squared, 3) 
rsq35 <- summary(lm(pr.real ~ anet.21p, pt35))["r.squared"]; round(rsq35$r.squared, 3) 

# slope
slp25 <- summary(lm(pr.real ~ anet.21p, pt25))["coefficients"][1]$coefficients[2,1]; round(slp25, 3) 
slp30 <- summary(lm(pr.real ~ anet.21p, pt30))["coefficients"][1]$coefficients[2,1]; round(slp30, 3) 
slp35 <- summary(lm(pr.real ~ anet.21p, pt35))["coefficients"][1]$coefficients[2,1]; round(slp35, 3) 

# slope error
slp25se <- summary(lm(pr.real ~ anet.21p, pt25))["coefficients"][1]$coefficients[2,2]; round(slp25se, 3) 
slp30se <- summary(lm(pr.real ~ anet.21p, pt30))["coefficients"][1]$coefficients[2,2]; round(slp30se, 3) 
slp35se <- summary(lm(pr.real ~ anet.21p, pt35))["coefficients"][1]$coefficients[2,2]; round(slp35se, 2) 


results <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(pr.real ~ setTleaf, data = .)))

results2 <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(pr.percent ~ setTleaf, data = .)))


results3 <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(anet.21p ~ setTleaf, data = .)))


outs_summary <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(anet_mean = mean(anet.21p, na.rm = TRUE),
            anet_se = sd(anet.21p, na.rm = TRUE) / sqrt(n()),
            pr_mean = mean(pr.real, na.rm = TRUE),
            pr_se = sd(pr.real, na.rm = TRUE) / sqrt(n()),
            phi_mean = mean(pr.percent, na.rm = TRUE),
            phi_se = sd(pr.percent, na.rm = TRUE) / sqrt(n()))





















### misc code for scavenging later -----
##---------------------------------------------------------


#### Plots of ETR over temperature #####--------------------------------------------
svg(filename = "ETR.21p-plot.svg", width = 16, height = 4.5, bg = "transparent")


ETR.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(ETR.21p.avg = mean(ETR.21p),
            se = sd(ETR.21p)/ sqrt(length(ETR.21p)))


ggplot(ETR.table, aes(y = ETR.21p.avg, x=setTleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = setTleaf, y = ETR.21p, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=ETR.21p.avg-se, ymax=ETR.21p.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "limegreen", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7)+
  scale_y_continuous(limits = c(0, 200), name = 
                       expression(paste(ETR), ' (', mu * ~'mol'~ " photon"[2]~' m'^{-2}*' s'^{-1}*')'))+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) 
dev.off()

###Plot of O2-free ETR over temperature##-----------------------------------------
svg(filename = "ETR.0p-plot.svg", width = 16, height = 4.5, bg = "transparent")

ETRRp.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(ETR.0p.avg = mean(ETR.0p),
            se = sd(ETR.0p)/ sqrt(length(ETR.0p)))

ggplot(ETRRp.table, aes(y = ETR.0p.avg, x=setTleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = setTleaf, y = gsw.0p, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=ETR.0p.avg-se, ymax=ETR.0p.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "cornflowerblue", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7)+
  scale_y_continuous(limits = c(0, 200), name = 
                       expression(paste(ETR)))+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) 
dev.off()


###Plot of percentage change in ambient vs O2-free gsw over temperature##-----------------------------------------###
svg(filename = "ETR.percent-plot.svg", width = 16, height = 4.5, bg = "transparent")

ETRpercent.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(ETR.percent.avg = mean(ETR.percent),
            se = sd(ETR.percent)/ sqrt(length(ETR.percent)))
ggplot(ETRpercent.table, aes(y = ETR.percent.avg, x=setTleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = setTleaf, y = ETR.percent, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=ETR.percent.avg-se, ymax=ETR.percent.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "red3", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7)+
  scale_y_continuous(limits = c(-0.7, 0.7), name = 
                       expression(ETR0p/ETR21p))+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) 
dev.off()

##---------------------------------------------------------

## Figure 1A ----
svg(filename = "pr-raw-new.svg", width = 16, height = 4.5, bg = "transparent")
absolute.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.real.avg = mean(pr.real),
            se = sd(pr.real)/ sqrt(length(pr.real)))

ggplot(absolute.table, aes(y = pr.real.avg, x=setTleaf, group=sp))+
  geom_errorbar(aes(ymin=pr.real.avg-se, ymax=pr.real.avg+se),col ="grey70", width= 2.5)+ #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  geom_point(col="grey30", size = 3, pch = 21, fill = "cornflowerblue", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7) +
  scale_y_continuous(limits = c(0, 15), name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  scale_x_continuous(limits = c(22,37), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70"))
dev.off()

### Plot showing Anet 21p over temp##-------------------------------###
## Figure 1B ----
svg(filename = "anet21-raw-new.svg", width = 16, height = 4.5, bg = "transparent")
anet.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(anet.21p.avg = mean(anet.21p),
            se = sd(anet.21p) / sqrt(length(anet.21p)))

ggplot(anet.table, aes(y = anet.21p.avg, x=setTleaf, group=sp))+
  geom_smooth(se = F, method="lm", col = "grey80")+
  geom_errorbar(aes(ymin=anet.21p.avg-se, ymax=anet.21p.avg+se),col ="grey70", width= 2.5) +
  geom_point(col="grey30", size = 3, pch = 21, fill = "limegreen", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=anet.21p, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7) +
  scale_y_continuous(limits = c(0, 25), name = expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70"))
dev.off()

###Plot of pr.percent##-----------------------------------------###
svg(filename = "pr-percent-new.svg", width = 16, height = 4.5, bg = "transparent")

percent.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.percent.avg = mean(pr.percent),
            se = sd(pr.percent)/ sqrt(length(pr.percent)))

ggplot(percent.table, aes(y = pr.percent.avg, x=setTleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = setTleaf, y = pr.percent, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=pr.percent.avg-se, ymax=pr.percent.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "red3", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7)+
  scale_y_continuous(limits = c(0, 1.3), name = 
                       expression(paste(italic(R)[p]/italic(A)[Net])))+
  scale_x_continuous(limits = c(22,37), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) 
dev.off()


##---------------------------------------------------------##


###########

outs$tl <- as.factor(outs$setTleaf)

## photorespiration across species boxplot (ignoring temp). First one is for just the relationship, the second one shows the temperature spread.

svg(filename = "Rp-species.svg", width = 12, height = 6, bg = "transparent")

ggplot(outs, aes(y = pr.real, x=sp)) +
  geom_boxplot(aes(x = sp, y = pr.real), outliers =  F)  +
  geom_jitter(aes(size = tl),  color = "grey", alpha = 0) +
  
  xlab(" ") +
  scale_y_continuous(limits = c(0,15), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  coord_flip()
dev.off()



svg(filename = "Rp-species-spotted.svg", width = 12, height = 6, bg = "transparent")

ggplot(outs, aes(y = pr.real, x=sp)) +
  geom_boxplot(aes(x = sp, y = pr.real), outliers =  F)  +
  geom_jitter(aes(size = tl), color = as.factor(outs$setTleaf), alpha = 0.5) +
  xlab(" ") +
  scale_y_continuous(limits = c(0,15), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        panel.border = element_rect(color = "grey70")) +
  coord_flip()
dev.off()




#------------------------------#


## photosynthesis across species boxplot (ignoring temp). First one is for just the relationship, the second one shows the temperature spread.

svg(filename = "Anet-species.svg", width = 12, height = 6, bg = "transparent")

ggplot(outs, aes(y = anet.21p, x=sp)) +
  geom_boxplot(aes(x = sp, y = anet.21p), outliers =  F)  +
  geom_jitter(aes(size = tl),  color = "grey", alpha = 0) +
  
  xlab(" ") +
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(A)[net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  coord_flip()
dev.off()


svg(filename = "Anet-species-spotted.svg", width = 12, height = 6, bg = "transparent")

ggplot(outs, aes(y = anet.21p, x=sp)) +
  geom_boxplot(aes(x = sp, y = anet.21p), outliers =  F)  +
  geom_jitter(aes(size = tl), color = as.factor(outs$setTleaf), alpha = 0.5) +
  xlab(" ") +
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(A)[net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        panel.border = element_rect(color = "grey70")) +
  coord_flip()
dev.off()



ggplot(pt, aes(y = pr, x=setTleaf)) +
  geom_point(aes())  +
  geom_errorbar(aes(ymin=pr-prse, ymax=pr+prse),col ="grey70", width= 2.5)+
  geom_smooth(se = F, method="lm", col = "grey80") +
  facet_wrap(~sp, ncol = 7)+
  scale_y_continuous(limits = c(0, 1.3), name = 
                       expression(paste(italic(R)[p]/italic(A)[Net])))+
  scale_x_continuous(limits = c(22,37), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) 

  
  
  geom_errorbar()
  xlab(" ") +
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(A)[net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        panel.border = element_rect(color = "grey70")) +
  coord_flip()





### Species Rp over Anet graph original

svg(filename = "Rp-over-Anet-alltemp.svg", width = 10, height = 3.3, bg = "transparent")

### use the next one---
ggplot(outs, aes(y = pr.real, x=anet.21p, color = sp)) +
  geom_smooth(se = T, method="lm", col = "grey80") +
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(0, 20), 
                     name = expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  scale_y_continuous(limits = c(0,20), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  ggthemes::theme_base() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") + 
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  facet_wrap(~setTleaf)
  geom_errorbar(aes(ymin=pr.real-prse, ymax=pr.real+prse)) +
  geom_errorbarh(aes(xmin=pr.real-psse, xmax=pr.real+psse)) + 

dev.off()

come here
  
  
  ggplot(outs, aes(y = pr.real, x=setTleaf)) +
    geom_smooth(se = T, method="lm", col = "grey80") +
    geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
    
    scale_y_continuous(limits = c(0,20), 
                       name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
    ggthemes::theme_base() +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          panel.border = element_rect(color = "grey70")) +
    scale_shape_manual(values = c(21,22,23,24,25,1,2))+
    geom_abline(slope = 1, linetype = "dashed", color = "grey50") + 
    geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
    facet_wrap(~sp, ncol = 7) + geom_errorbar(aes(ymin=pr.real-prse, ymax=pr.real+prse)) +
    geom_errorbarh(aes(xmin=pr.real-psse, xmax=pr.real+psse)) + 
    
    
  
  
  
  
###### main figure 2 ------

svg(filename = "Figure-2.svg", width = 9, height = 3, bg = "transparent")
ggplot(outs_summary, aes(x = anet_mean, y = pr_mean, color = sp)) +
  geom_smooth(se = TRUE, method = "lm", col = "grey50") +  
  geom_point(size = 2, stroke = 0.5, aes(color = sp, shape = sp)) +
  geom_point(data = outs_summary, aes(x = anet_mean, y = pr_mean, shape = sp), 
             size = 2, color = "black", fill = "white", stroke = 1) +
  geom_errorbar(data = outs_summary, aes(ymin = pr_mean - pr_se, ymax = pr_mean + pr_se), width = 0, color = "black") +
  geom_errorbarh(data = outs_summary, aes(xmin = anet_mean - anet_se, xmax = anet_mean + anet_se), height = 0, color = "black") +
  
  scale_x_continuous(limits = c(0, 20), 
                     name = expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  scale_y_continuous(limits = c(0, 20), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  ggthemes::theme_base() +
  facet_wrap(~setTleaf) +
  
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50")

dev.off()

#--------------------------------------------------------------------#

###ETR analysis starts here



### Photorespiration ratio plotted over ratio of ETR used in each environment

svg(filename = "Phi-ETRdelta.svg", width = 10, height = 3.3, bg = "transparent")

ggplot(pt, aes(y = ETR.delta, x=pr, color = sp)) +
  geom_boxplot(aes(group = cut_width(setTleaf, 1))) +
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(0, 20), 
                     name = expression(Rp)) +
  scale_y_continuous(limits = c(0,65), 
                     name = expression(paste(ETR.delta))) +
  ggthemes::theme_base() +
  geom_errorbar(aes(ymin=ETR.delta-ETR.deltase, ymax=ETR.delta + ETR.deltase)) +
  geom_errorbarh(aes(xmin=pr-prse, xmax=pr+prse)) + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") + 
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  facet_wrap(~setTleaf)

dev.off()



#--------------------- test ETR over ETR
svg(filename = "ETR.percent-over-Rp.svg", width = 7, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = ETR.percent, x=pr, color = sp)) +
  geom_boxplot(aes(group = cut_width(setTleaf, 1))) +
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(0, 20), 
                     name = expression(Rp)) +
  scale_y_continuous(limits = c(-0.45,0.27), 
                     name = expression(paste(ETR.percent))) +
  ggthemes::theme_base() +
  geom_errorbar(aes(ymin=ETR.percent-ETR.percentse, ymax=ETR.percent + ETR.percentse)) +
  geom_errorbarh(aes(xmin=pr-prse, xmax=pr + prse)) + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") + 
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  facet_wrap(~setTleaf)
dev.off()


##Anet

svg(filename = "ETR.21p-tleaf.svg", width = 7, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = ETR.21p, x=setTleaf, color = sp)) +
  geom_boxplot(aes(group = cut_width(setTleaf, 1))) +
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(22.5,37.5), 
                     name = expression(Anet.21p)) +
  scale_y_continuous(limits = c(0,150), 
                     name = expression(paste(ETR.21p))) +
  ggthemes::theme_base() +
  geom_errorbar(aes(ymin=ETR.21p-ETR.21pse, ymax=ETR.21p + ETR.21pse)) +
  geom_errorbarh(aes(xmin=ps-psse, xmax=ps + psse)) + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") + 
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp))
dev.off()

###Rp

svg(filename = "ETR.0p-tleaf.svg", width = 7, height = 4.5, bg = "transparent")
ggplot(pt, aes(y = ETR.0p, x=setTleaf, color = sp)) +
  geom_boxplot(aes(group = cut_width(setTleaf, 1))) +
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(22.5, 37.5), 
                     name = expression(Rp)) +
  scale_y_continuous(limits = c(0,150), 
                     name = expression(paste(ETR.0p))) +
  ggthemes::theme_base() +
  geom_errorbar(aes(ymin=ETR.0p-ETR.0pse, ymax=ETR.0p + ETR.0pse)) +
  geom_errorbarh(aes(xmin=pr-prse, xmax=pr + prse)) + 
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") + 
  geom_point(size = 3, stroke = 0.85, aes(color = sp, shape = sp)) 

dev.off()






svg(filename = "ETR21p-over-ETR0p.svg", width = 7, height = 4.5, bg = "transparent")
outs$st <- as.factor(outs$setTleaf)
ggplot(outs, aes(y = anet.21p, x=ETR.21p)) +
  geom_point(size = 3, 
             stroke = 0.85, aes(color = sp, shape = sp)) +
  geom_smooth(se = T, method="lm", col = "grey50") +
  scale_x_continuous(limits = c(0, 200), 
                     name = expression(paste(ETR))) +
  scale_y_continuous(limits = c(0,20), 
                     name = expression(paste(Anet)))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") 
dev.off()



svg(filename = "Rp-over-ETR.svg", width = 7, height = 4.5, bg = "transparent")
outs$st <- as.factor(outs$setTleaf)
ggplot(outs, aes(y = pr.real, x=ETR.0p)) +
  geom_point(size = 3, 
             stroke = 0.85, aes(color = sp, shape = sp)) +
  geom_smooth(se = T, method="lm", col = "grey50") +
  scale_x_continuous(limits = c(0, 200), 
                     name = expression(paste(ETR))) +
  scale_y_continuous(limits = c(0,20), 
                     name = expression(paste(Rp)))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") 
dev.off()



svg(filename = "Phi-over-ETRdelta.svg", width = 7, height = 4.5, bg = "transparent")
outs$st <- as.factor(outs$setTleaf)
ggplot(outs, aes(y = pr.percent, x=ETR.21p)) +
  geom_point(size = 3, 
             stroke = 0.85, aes(color = sp, shape = sp)) +
  geom_smooth(se = T, method="lm", col = "grey50") +
  scale_x_continuous(limits = c(0, 200), 
                     name = expression(paste(ETR))) +
  scale_y_continuous(limits = c(0,1), 
                     name = expression(paste(Phi)))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") 
dev.off()



mod6 <- nlme::lme(ETR.percent ~  phi * setTleaf, data = pt, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6) ### significant!

mod6 <- nlme::lme(pr.real ~  ETR.percent, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6) ### Nothing 0.76




# Does higher photorespiration rates = higher photosynthesis rates??? (spoiler, yes)



#### all replicates plotted - donot use. 

svg(filename = "Rp-over-Anet.svg", width = 7, height = 4.5, bg = "transparent")
outs$st <- as.factor(outs$setTleaf)
ggplot(outs, aes(y = pr.real, x=anet.21p)) +
  geom_point(size = 3, 
             stroke = 0.85, aes(color = sp, shape = sp)) +
  geom_smooth(se = T, method="lm", col = "grey50") +
  scale_x_continuous(limits = c(0, 25), 
                     name = expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") 
dev.off()

#---------------------------------------------------------------------#





# svg(filename = "Rp-over-Anet-alltemp.svg", width = 16, height = 4.5, bg = "transparent")

#### all replicates plotted - donot use. 


ggplot(outs, aes(y = pr.real, x=anet.21p, color = sp)) +
  geom_smooth(data = outs, mapping = aes(y=pr.real, x=anet.21p), 
              se = T, method="lm", col = "grey80") +
  geom_point(size = 3, 
             stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(0, 25), 
                     name = expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~setTleaf)
dev.off()


######  Creating the dataframes for each temperature for the models ----

outs25 <- subset(outs, setTleaf == 25)
outs30 <- subset(outs, setTleaf == 30)
outs35 <- subset(outs, setTleaf == 35)


######  Rp against Anet split by temperature ----

mod6 <- nlme::lme(pr.real ~  anet.21p, data = outs, 
                 random = ~1|sp, 
                 method = "REML", 
                 na.action=na.omit) ; anova(mod6)


mod7 <- nlme::lme(pr.real ~  anet.21p, data = outs25, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod7)

mod7 <- nlme::lme(pr.real ~  anet.21p, data = outs30, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod7)

mod7 <- nlme::lme(pr.real ~  anet.21p, data = outs35, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod7)










#-------------------------------#


### I WOULD LIKE TO DISCUSS THIS TO CONFIRM IM DOING THINGS CORRECTLY

# Checking the individual significance of species temperature response.

pt <- outs %>%
  group_by(sp, treeid, setTleaf) %>%
  summarise(pr = mean(pr.real),
            prse = sd(pr.real)/sqrt(length(pr.real)),
            ps = mean(anet.21p),
            psse = sd(anet.21p)/ sqrt(length(anet.21p)),
            phi = mean(pr.percent),
            phise = sd(pr.percent)/sqrt(length(pr.percent)))


###PHOTORESPIRATION RATES / TEMP

Betulapen <- subset(pt, sp == "Betula pendula")

mod6 <- nlme::lme(pr ~  setTleaf, data = Betulapen, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Fagus <- subset(pt, sp == "Fagus sylvatica")

mod6 <- nlme::lme(pr ~  setTleaf, data = Fagus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Acer <- subset(pt, sp == "Acer platanoides")

mod6 <- nlme::lme(pr ~  setTleaf, data = Acer, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Batulapub <- subset(pt, sp == "Betula pubescens")

mod6 <- nlme::lme(pr ~  setTleaf, data = Batulapub, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Tilia <- subset(pt, sp == "Tilia cordata")

mod6 <- nlme::lme(pr ~  setTleaf, data = Tilia, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Corylus <- subset(pt, sp == "Corylus avellana")

mod6 <- nlme::lme(pr ~  setTleaf, data = Corylus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Sorbus <- subset(pt, sp == "Scandosorbus intermedia")


mod6 <- nlme::lme(pr ~  setTleaf, data = Sorbus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)




###PHOTOSYNTHESIS RATES / TEMP

Betulapen <- subset(pt, sp == "Betula pendula")

mod6 <- nlme::lme(ps ~  setTleaf, data = Betulapen, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Fagus <- subset(pt, sp == "Fagus sylvatica")

mod6 <- nlme::lme(ps ~  setTleaf, data = Fagus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Acer <- subset(pt, sp == "Acer platanoides")

mod6 <- nlme::lme(ps ~  setTleaf, data = Acer, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Batulapub <- subset(pt, sp == "Betula pubescens")

mod6 <- nlme::lme(ps ~  setTleaf, data = Batulapub, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Tilia <- subset(pt, sp == "Tilia cordata")

mod6 <- nlme::lme(ps ~  setTleaf, data = Tilia, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Corylus <- subset(pt, sp == "Corylus avellana")

mod6 <- nlme::lme(ps ~  setTleaf, data = Corylus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Sorbus <- subset(pt, sp == "Scandosorbus intermedia")


mod6 <- nlme::lme(ps ~  setTleaf, data = Sorbus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


# Creating a sheet with each temperature point containing one compiled value for each species. This table is used when gathering the highest/lowest mean differences in Anet and Rp across Tleaf. 

pt <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr = mean(pr.real),
            prse = sd(pr.real)/sqrt(length(pr.real)),
            ps = mean(anet.21p),
            psse = sd(anet.21p)/ sqrt(length(anet.21p)))

Sorbus <- subset(pt, sp == "Scandosorbus intermedia")


mod6 <- nlme::lme(pr ~  setTleaf, data = Sorbus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


-----#####Stomatal conductance for O2-free vs ambient across temperatures#---
svg(filename = "gsw-ratio-split-temp.svg", width = 16, height = 4.5, bg = "transparent")

#### all replicates plotted 


ggplot(outs, aes(y = gsw.0p, x=gsw.21p, color = sp)) +
  geom_smooth(data = outs, mapping = aes(y=gsw.0p, x=gsw.21p), 
              se = T, method="lm", col = "grey80") +
  geom_point(size = 3, 
             stroke = 0.85, aes(color = sp, shape = sp)) +
  scale_x_continuous(limits = c(0, 0.4), 
                     name = expression(paste(italic(gsw)["21p"]))) +
  scale_y_continuous(limits = c(0,0.4), 
                     name = expression(paste(italic(gsw)["0p"])))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~setTleaf)
dev.off()













#-------------------------------#



#Box plots of pr.percent across the leaf temperatures

ggplot(outs, aes(y = pr.percent, x=tl, color=tl)) +
  geom_boxplot(aes(x = tl, y = pr.percent, fill = setTleaf)) +
  geom_jitter() +
  scale_x_continuous(limits = c(20, 40), 
                     name = "Leaf Temperature (°C)") +
  scale_y_continuous(limits = c(0,15), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70"))


#-------- boxplot of Rp over temperature ----_#
svg(filename = "Rp-over-temp-box.svg", width = 16, height = 4.5, bg = "transparent")


ggplot(outs, aes(y = pr.real, x=tl)) +
  geom_boxplot(aes(x = tl, y = pr.real)) +
  geom_jitter(aes(color = sp)) +
  xlab("Leaf Temperature (°C)") +
  scale_y_continuous(limits = c(0,15), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70"))
dev.off()







###### STATISTICAL TESTS #########


## example model fits nlme ---- 

# Here we are running a for loop testing     the relationship between species, setTleaf and the combined effect of them (I.e the temperature response of species), against pr.real (Absolute photorespiration rates ("10")), pr.percent (Proportional photorespiration rates ("11")), anet.21p (Ambient photosynthesis rates ("4")), anet.0p (O2-free photosynthesis rates ("5")) and JO2.percent (proportional oxygenation rates calculated using ETR only ("31)).


colindeces <- c(10, 11, 4, 5, 31)
colindeces <- c(10)
method.slope <- c()
value.slope <- c()
species.slope <- c()

names(outs)[colindeces[1]]

for (index in 1:length(colindeces)){
  
  
  independentVariable = outs[,colindeces[index]] 
  
  mod1 <- nlme::lme(independentVariable ~  sp , data = outs, 
                    random = ~1|treeid, 
                    method = "REML", 
                    na.action=na.omit) ; anova(mod1)
  print(paste("This analysis is based on :", names(outs)[colindeces[index]]))
  print(summary(mod1))
  
  mod1 <- nlme::lme(independentVariable ~  setTleaf , data = outs, 
                    random = ~1|sp, 
                    method = "REML", 
                    na.action=na.omit) ; anova(mod1)
  
  print(paste("This analysis is based on :", names(outs)[colindeces[index]]))
  print(summary(mod1))
  
  mod1 <- nlme::lme(independentVariable ~  sp * setTleaf, 
                    data = outs, 
                    random = ~1|treeid, 
                    method = "REML", 
                    na.action=na.omit) ; anova(mod1)
  
  print(paste("This analysis is based on :", names(outs)[colindeces[index]]))
  print(summary(mod1))
  
  
  print("----------------------------------------------------------------")
}


outs25 <- subset(outs, setTleaf == 25)
outs30 <- subset(outs, setTleaf == 30)
outs35 <- subset(outs, setTleaf == 35)



######----


######  Rp against species split by treeid nested for each temp point ----
 mod4 <- nlme::lme(pr.real ~  sp, data = outs25, 
                   random = ~1|treeid, 
                   method = "REML", 
                   na.action=na.omit) ; anova(mod4)
 
 mod4 <- nlme::lme(pr.real ~  sp, data = outs30, 
                   random = ~1|treeid, 
                   method = "REML", 
                   na.action=na.omit) ; anova(mod4)
 
 mod4 <- nlme::lme(pr.real ~  sp, data = outs35, 
                   random = ~1|treeid, 
                   method = "REML", 
                   na.action=na.omit) ; anova(mod4)
 

mod4 <- nlme::lme(pr.real ~  anet.21p * setTleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)




long_data <- outs %>%
  gather(key = "method", value = "y", pr.percent, JO2.percent)

# Fit a model with an interaction term between method and setTleaf
long_model <- lm(y ~ sp * setTleaf * method, data = long_data)

# Check the model summary
summary(long_model)


# Fit the model with and without the interaction term between method and setTleaf
model_interaction <- lm(y ~ sp * setTleaf * method, data = long_data)
model_no_interaction <- lm(y ~ sp * setTleaf + method, data = long_data)

# Compare the models using a likelihood ratio test
anova(model_no_interaction, model_interaction)

# Visualize the comparison

ggplot(long_data, aes(x = setTleaf, y = y, color = method)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~sp)+
  labs(title = "Comparison of Slopes for y1 and y2") +
  theme_minimal()


### post hoc test ------- HSD, honest significant difference 

m7 <- aov(anet.0p ~ sp, data = outs, random = ~1|treeid, method = "REML", na.action=na.omit) ; anova(m7)
par(las=2)
par(mar=c(8,8,2,2)) 
plot(TukeyHSD(m7))


t1 <- TukeyHSD(m7)
t1 <- as.data.frame(t1$sp)
t1$com <- rownames(t1)

rownames(t1) <- NULL
colnames(t1) <- c("prdiff","lwr", "upr","padj", "com")
blah <- subset(t1, padj <= 0.05)





