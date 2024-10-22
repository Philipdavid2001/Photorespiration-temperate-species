#install.packages("Hmisc", dependencies = T)

packages_to_load <- c("drc", "bbmle", "labelled", "MALDIquant", 
                      "magicfor", "ggplot2", "ggpubr", "ggpmisc", "tidyverse", "broom", 
                      "Hmisc", "plotly", "PairedData", "DescTools", 
                      "generalhoslem", "graphics", "nls2", "plyr",
                      "ggthemes", "ggpubr", "signal", "segmented")

lapply(packages_to_load, library, character.only = TRUE)

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
             "pr.real",    
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
     
    
    # Calculate corrected Anet. The values for slope and intercept where           calculated using values taken from the literature.  
    Rd <- 0.05554*anet.21p+ 0.11395 
    corranet <- anet.21p-Rd          # "corranet" is Anet corrected for Rd
    anet.delta <- anet.0p - corranet
    
    #####TODO correct for photorespiration C02 release######
    
    
    pr.real <- anet.delta/corranet
    
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
                          pr.real,    
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
s.table <- outs %>%
  group_by(sp, setTleaf) %>%
  summarise(pr.real.avg = mean(pr.real),
            se = sd(pr.real))

ggplot(s.table, aes(y = pr.real.avg, x=setTleaf, group=sp))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(aes(ymin=pr.real.avg-se, ymax=pr.real.avg+se), width=.2)+ #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(data = outs, mapping = aes(y=pr.real, x=setTleaf), se = F, method="lm")+
  stat_poly_eq(
    data = outs,
    formula = y ~ x,
    aes(x = setTleaf, y = pr.real, label = paste(after_stat(eq.label),
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
  scale_y_continuous(limits = c(0, 1.8), name = "Photorespiration Rate")+
  scale_x_continuous(limits = c(20,40), name = "Leaf Temperature [°C]")+
  theme_minimal()




# Original pr ggplot script

ggplot(outs, aes(x = setTleaf, y = pr.real)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() + facet_wrap(~sp) + ylim(0,1) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  #geom_smooth(method = "lm", se = F)


# tweaked 
ggplot(outs, aes(x = setTleaf, y = pr.real)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_boxplot(mapping = aes(group = setTleaf)) + facet_wrap(~sp) + ylim(0,1) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)

#--------------------------------------------------------------#
svg(filename = "Figures/anet21p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = corranet)) +
    theme_bw() +
    xlab(lab_Tleaf) +
    ylab("Anet") +
    geom_point() + facet_wrap(~sp) + ylim(0,30) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()


svg(filename = "Figures/anet0p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = anet.0p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("Anet") +
  geom_point() + facet_wrap(~sp) + ylim(0,30) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()



svg(filename = "Figures/ETR21p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = ETR.21p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("ETR") +
  geom_point() + facet_wrap(~sp) + ylim(0,160) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()


svg(filename = "Figures/ETR0p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = ETR.0p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("ETR") +
  geom_point() + facet_wrap(~sp) + ylim(0,160) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()


svg(filename = "Figures/NPQ21p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = NPQ.21p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("NPQ") +
  geom_point() + facet_wrap(~sp) + ylim(-1,-0.997) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()

svg(filename = "Figures/NPQ0p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = NPQ.0p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("NPQ") +
  geom_point() + facet_wrap(~sp) + ylim(-1,-0.997) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()



## photorespiration proportion ----- 

anova(aov(pr.real ~  sp, data = outs))
anova(aov(corranet ~  sp, data = outs))
anova(aov(pr.real ~  corranet, data = outs))
anova(aov(pr.real ~  setTleaf+sp, data = outs))
anova(aov(corranet ~  setTleaf, data = outs))
anova(aov(pr.real ~  sp * setTleaf, data = outs))
anova(aov(corranet ~  sp * setTleaf, data = outs))

m1 <- nlme::lme(pr.real ~  sp * setTleaf, data = outs, random = ~1|treeid)
anova(aov(pr.real ~  sp * setTleaf, data = outs))

m1




# Things I wanna test


# Is higher Anet correlated with higher PR? 
# Is lower ETR correlated to higher PR?
# 

