#install.packages("Hmisc", dependencies = T)
#
# packages_to_load <- c("drc", "bbmle", "labelled", "MALDIquant", 
#                      "magicfor", "ggplot2", "ggpubr", "ggpmisc", "tidyverse", #"broom", 
#                     "Hmisc", "plotly", "PairedData", "DescTools", 
#                     "generalhoslem", "graphics", "nls2", "plyr",
#                     "ggthemes", "ggpubr", "signal", "segmented")
#
# lapply(packages_to_load, library, character.only = TRUE)

## lable ----
ytleaf          <- expression("T"["leaf"]~degree~"C")
photorate       <- expression("" ~mu ~ mol ~ "CO"[2]~m^{-2} ~ "s"^{-1})
labrsq          <- expression(R^{2}~"=")
labvpdl         <- "VPD (kPa)"
lab_photo       <- expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))
lab_photo1      <- expression(paste(italic(A)[Net], ' (', mu * 'mol m'^{-2}*'s'^{-1}*')'))
lab_Tleaf       <- expression(paste("Leaf temperature (", degree,"C)"))
labcond         <- expression("g"["s"]~~"in"~mol~H[2]*O~m^{-2}~s^{-1})
labcond1        <- expression("Conductance to "~H[2]*O~"in"~mol~H[2]*O~m^{-2}~s^{-1})
labtrmmol       <- expression("E"~~"in"~mol~H[2]*O~m^{-2}~s^{-1})
etrrate         <- expression("" ~mu ~ mol ~"photon" ~m^{-2} ~ "s"^{-1})
photorate       <- expression("" ~mu ~ mol ~ "CO"[2]~m^{-2} ~ "s"^{-1})
labgs           <- expression(paste(italic(g[s])*" ("~mol[H[2]*O]~m^{-2}~s^{-1}~")"))


## data processing ---- Load ONE of the following csv's


# Species
df <- read.csv("Data/Uppsala-2024-Summer-Photorespiration-SpotMes-TreeSpp.csv", header = T, stringsAsFactors = T, sep = ";")


# Ecotypes

df <- read.csv("Data/Uppsala-2024-Summer-Photorespiration-SpotMes-Birch-Ecotypes.csv", header = T, stringsAsFactors = T, sep = ";")


#-----------------------------------------------------#


## use the following filters as necessary
df     <- subset(df, dq == "yes") # eliminates unreliable values.

df$pairid <- paste0(df$sprp, "-", df$setTleaf)

df$pairid <- as.factor(df$pairid)


length(unique(df$pairid))       # individual runs + season

dflist    <-   split(df,list(df$pairid))
dflist    <-   split(df,unique(list(df$pairid)))
dflist    <-   dflist[sapply(dflist, nrow)>0] 

setwd("output")


correct_RD <- function(data, output_path){
  
  ColumnNames <- c("sp","treeid","setTleaf",
             "anet.21p",    
             "anet.0p",
             "Rd",
             "corranet",
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
             "JO1.percent",
             "JC1.percent",
             "JO2.percent",
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
     
    
    # Calculate corrected Anet. The values for slope and intercept       where calculated using values taken from the literature. 
    
    Rd               <-      0.05554*anet.21p+ 0.11395 
    corranet         <-      anet.21p-Rd          
    # "corranet" is Anet corrected for Rd
    
    
    anet.delta       <-      anet.0p - anet.21p
    pr.CO2           <-      anet.delta * 0.5
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
    JO1.percent      <-      JO1/JT
    JC1.percent      <-      JC1/JT
    JO2.percent      <-      JO2/JT
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
                          Rd, 
                          corranet,
                          anet.delta,
                          pr.CO2,
                          pr.real,
                          pr.percent,
                          gsw.21p,    
                          gsw.0p,     
                          gsw.delta,  
                          gsw.percent,
                          p21$E,     
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
                          JC2,
                          JO1.percent,
                          JC1.percent,
                          JO2.percent,
                          JC2.percent,
                          NPQ.21p,
                          NPQ.0p,
                          NPQ.delta, 
                          NPQ.percent)
    
    names(new_data) <- names(output_data)
    output_data <- rbind(output_data, new_data)
    
    
  }
  
  file_path = paste(output_path, "species-output.csv", sep = "/")
  if(!dir.exists(output_path)){
    dir.create(output_path, recursive = T)
  }
  
  
  index <- 1
  while(file.exists(file_path)){
    file_path <- paste(output_path, paste("species-output", as.character(index), ".csv", sep =""), sep = "/")
    index = index + 1
  }
 
  

  write.table(output_data, file = file_path, 
              row.names = FALSE, col.names = T, sep = ";")
}

correct_RD(dflist, "./")


outs <- read.csv("species-output.csv", stringsAsFactors = T, sep = ";")


# Making the table for the concatenated PR values for each species and each temperature point.

library(dplyr)
library(stats)
library(base)
library(ggplot2)

s.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.real.avg = mean(pr.real),
            se = sd(pr.real))


#svg(filename = "pr-raw.svg", width = 16, height = 4.5, bg = "transparent")

ggplot(s.table, aes(y = pr.real.avg, x=setTleaf, group=sp))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(aes(ymin=pr.real.avg-se, ymax=pr.real.avg+se), width=.2)+ #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(data = outs, mapping = aes(y=pr.real, x=setTleaf), se = F, method="lm")+
  ggpmisc::stat_poly_eq(
    data = outs,
    formula = y ~ x,
    aes(x = setTleaf, y = pr.real, label = paste(after_stat(p.value.label), sep = "")),
    label.x = 20,
    label.y = 1.5,
    digits = 2
  )+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7) +
  scale_y_continuous(limits = c(0, 15), name = "Net CO2 assimilation rate (0% O2) - Net CO2 assimilation rate (21 % O2)")+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base()
#dev.off()


###Plot of pr.percent

library(dplyr)
s.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.percent.avg = mean(pr.percent),
            se = sd(pr.percent))

svg(filename = "pr-percent.svg", width = 16, height = 4.5, bg = "transparent")
ggplot(s.table, aes(y = pr.percent.avg, x=setTleaf, group=sp))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(aes(ymin=pr.percent.avg-se, ymax=pr.percent.avg+se), width=.2)+ #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(data = outs, mapping = aes(y=pr.percent, x=setTleaf), se = F, method="lm")+
  stat_poly_eq(
    data = outs,
    formula = y ~ x,
    aes(x = setTleaf, y = pr.percent, label = paste(after_stat(p.value.label), sep = "~~~")),
    label.x = 22,
    label.y = 1.5,
    digits = 2
  )+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7)+
  scale_y_continuous(limits = c(0, 1.5), name = "Rp / Net CO2 assimilation rate (21 % O2) (%)")+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  theme_base()
dev.off()



### Plot showing Anet 21p over temp


s.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(anet.21p.avg = mean(anet.21p),
            se = sd(anet.21p))


svg(filename = "anet21-raw.svg", width = 16, height = 4.5, bg = "transparent")

ggplot(s.table, aes(y = anet.21p.avg, x=setTleaf, group=sp))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(aes(ymin=anet.21p.avg-se, ymax=anet.21p.avg+se), width=.2)+ #pmin can be used to cap the ymax. (e.g., pmin(anet.21p.avg+se, 1.0))
  geom_smooth(data = outs, mapping = aes(y=anet.21p, x=setTleaf), se = F, method="lm")+
  ggpmisc::stat_poly_eq(
    data = outs,
    formula = y ~ x,
    aes(x = setTleaf, y = anet.21p, label = paste(after_stat(p.value.label), sep = "")),
    label.x = 20,
    label.y = 1.5
  )+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=anet.21p, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp, ncol = 7) +
  scale_y_continuous(limits = c(0, 25), name = "Net CO2 assimilation rate (21 % O2)")+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature (°C)")+
  ggthemes::theme_base()
dev.off()


### Plot showing JO2 and JC2

library(dplyr)
s.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(JO2.percent.avg = mean(JO2.percent),
            se = sd(JO2.percent))

svg(filename = "JO2.percent.svg", width = 7, height = 4, bg = "transparent")
ggplot(s.table, aes(y = JO2.percent.avg, x=setTleaf, group=sp))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(aes(ymin=JO2.percent.avg-se, ymax=JO2.percent.avg+se), width=.2)+ #pmin can be used to cap the ymax. (e.g., pmin(JO2.percent.avg+se, 1.0))
  geom_smooth(data = outs, mapping = aes(y=JO2.percent, x=setTleaf), se = F, method="lm")+
  stat_poly_eq(
    data = outs,
    formula = y ~ x,
    aes(x = setTleaf, y = JO2.percent, label = paste(after_stat(eq.label),
                                                 after_stat(rr.label),
                                                 after_stat(p.value.label), sep = "~~~")),
    label.x = 22,
    label.y = 1.5,
    digits = 2
  )+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp)+
  scale_y_continuous(limits = c(0, 1), name = "Oxygenation proportion (ETR)")+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature [°C]")+
  theme_minimal()
dev.off()




library(dplyr)
s.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(JC2.percent.avg = mean(JC2.percent),
            se = sd(JC2.percent))

svg(filename = "JC2.percent.svg", width = 7, height = 4, bg = "transparent")
ggplot(s.table, aes(y = JC2.percent.avg, x=setTleaf, group=sp))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(aes(ymin=JC2.percent.avg-se, ymax=JC2.percent.avg+se), width=.2)+ #pmin can be used to cap the ymax. (e.g., pmin(JC2.percent.avg+se, 1.0))
  geom_smooth(data = outs, mapping = aes(y=JC2.percent, x=setTleaf), se = F, method="lm")+
  stat_poly_eq(
    data = outs,
    formula = y ~ x,
    aes(x = setTleaf, y = JC2.percent, label = paste(after_stat(eq.label),
                                                     after_stat(rr.label),
                                                     after_stat(p.value.label), sep = "~~~")),
    label.x = 22,
    label.y = 1.5,
    digits = 2
  )+
  # stat_regline_equation(data = outs, aes(x = setTleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
  # formula = y~x,
  # label.x =22, label.y = 1.5)+
  facet_wrap(~sp)+
  scale_y_continuous(limits = c(0, 1), name = "Carboxylation proportion (ETR)")+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature [°C]")+
  theme_minimal()
dev.off()




# Original pr ggplot script

ggplot(outs, aes(x = setTleaf, y = pr.percent)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() + facet_wrap(~sp) + ylim(0,1) + xlim(20,40) 
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  #geom_smooth(method = "lm", se = F)


# tweaked 
ggplot(outs, aes(x = setTleaf, y = pr.percent)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_boxplot(mapping = aes(group = setTleaf)) + facet_wrap(~sp) + ylim(0,1) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)



## photorespiration proportion ----- 

anova(aov(pr.real ~  sp, data = outs))
anova(aov(pr.real ~  corranet, data = outs))
anova(aov(pr.real ~  setTleaf, data = outs))
anova(aov(pr.real ~  sp * setTleaf, data = outs))

m1 <- nlme::lme(pr.real ~  sp * setTleaf, data = outs, random = ~1|treeid)
anova(aov(pr.real ~  sp * setTleaf, data = outs))

m1



anova(aov(JO2.percent ~  sp, data = outs))             # Very significant
anova(aov(JO2.percent ~  setTleaf, data = outs))       #  Very significant
anova(aov(pr.percent ~  sp, data = outs))              #  Very significant
anova(aov(pr.percent ~  setTleaf, data = outs))        #  Very significant
anova(aov(JO2.percent ~  sp * setTleaf, data = outs))  #  Not significant
anova(aov(pr.percent ~  sp * setTleaf, data = outs))   #  Not significant


###### STATISTICAL TESTS #########


## example model fits nlme ---- 

# Here we are running a for loop testing     the relationship between species, setTleaf and the combined effect of them (I.e the temperature response of species), against pr.real (Absolute photorespiration rates ("10")), pr.percent (Proportional photorespiration rates ("11")), anet.21p (Ambient photosynthesis rates ("4")), anet.0p (O2-free photosynthesis rates ("5")) and JO2.percent (proportional oxygenation rates calculated using ETR only ("31)).


colindeces <- c(10, 11, 4, 5, 31)
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


library(tidyr)

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
library(ggplot2)
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





