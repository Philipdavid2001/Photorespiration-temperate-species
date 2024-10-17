#install.packages("Hmisc", dependencies = T)

packages_to_load <- c("drc", "bbmle", "labelled", "MALDIquant", 
                      "magicfor", "ggplot2", "tidyverse", "broom", 
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


## data processing ----

df <- read.csv("Data/Uppsala-2024-Summer-Photorespiration-SpotMes-TreeSpp.csv", header = T, stringsAsFactors = T, sep = ";")



## use the following filters as necessary
df     <- subset(df, dq == "yes") # eliminates unreliable values.

df$pairid <- paste0(df$sprp, "-", df$setTleaf)

df$pairid <- as.factor(df$pairid)


length(unique(df$pairid))       # individual runs + season

dflist    <-   split(df,list(df$pairid))
dflist    <-   split(df,unique(list(df$pairid)))
dflist    <-   dflist[sapply(dflist, nrow)>0] 

setwd("output")

outs = NULL


for(i in 1:length(dflist)) {
    tempdf          <-          (dflist[[i]])
    dlength         <-          nrow(tempdf)
    if( dlength != 2 ) { 
        next 
        }
    sp              <-          tempdf$sp[1]; sp
    treeid          <-          tempdf$treeid[1]; treeid
    setTleaf        <-          tempdf$setTleaf[1]; setTleaf
    p21             <-          subset(tempdf, olev == "21p")    
    p0              <-          subset(tempdf, olev != "21p")

    anet.21p        <-          p21$A
    anet.0p         <-          p0$A
    anet.delta      <-          (anet.0p - anet.21p)
    pr.proxy        <-          anet.delta/anet.21p
    
    gsw.21p        <-          p21$gsw
    gsw.0p         <-          p0$gsw
    gsw.delta      <-          (gsw.0p - gsw.21p)
    gsw.percent    <-          gsw.delta/gsw.21p

    
    E.21p          <-          p21$E
    E.0p           <-          p0$E
    E.delta        <-          (E.0p - E.21p)
    E.percent      <-          E.delta/E.21p
    
    ETR.21p        <-          p21$ETR
    ETR.0p         <-          p0$ETR
    ETR.delta      <-          (ETR.0p - ETR.21p)
    ETR.percent    <-          ETR.delta/ETR.21p
    
  
    Ca.21p          <-          p21$Ca
    Ci.21p          <-          p21$Ci
    Ca.0p           <-          p0$Ca 
    Ci.0p           <-          p21$Ci
    C.21p           <-          Ca.21p - Ci.21p
    C.0p            <-          Ca.0p - Ci.0p
    Ca.delta        <-          (Ca.0p - Ca.21p) 
    Ci.delta        <-          (Ci.0p - Ci.21p)
    C.delta         <-          (C.0p - C.21p)
    Ca.prop         <-          Ca.delta/Ca.21p
    Ci.prop         <-          Ci.delta/Ci.21p
    C.prop          <-          C.delta/C.21p
    
    NPQ.21p         <-          p21$NPQ
    NPQ.0p          <-          p0$NPQ
    
    params <-  rbind(c(as.character(sp),
                       as.character(treeid),
                       as.character(setTleaf), 
                       anet.21p,    
                       anet.0p,     
                       anet.delta,  
                       pr.proxy,    
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
                       Ca.21p,
                       Ci.21p,
                       C.21p,
                       Ca.0p,
                       Ci.0p,
                       C.0p,
                       Ca.delta,
                       Ci.delta,
                       C.delta,
                       Ca.prop,
                       Ci.prop,
                       C.prop,
                       NPQ.21p,
                       NPQ.0p))
    
    colnames(params) <- c("sp","treeid","setTleaf", 
                       "anet.21p",    
                       "anet.0p",     
                       "anet.delta",  
                       "pr.proxy",    
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
                       "Ca.21p",
                       "Ci.21p",
                       "C.21p",
                       "Ca.0p",
                       "Ci.0p",
                       "C.0p",
                       "Ca.delta",
                       "Ci.delta",
                       "C.delta",
                       "Ca.prop",
                       "Ci.prop",
                       "C.prop",
                       "NPQ.21p",
                       "NPQ.0p")
    outs = rbind(outs, params)
    print(i)
    
}



write.table(outs, file = "Species-output.csv", 
            row.names = FALSE, col.names = T, sep = ",")
outs <- read.csv("Species-output.csv", stringsAsFactors = T)

outs <- read.csv("Species-output.csv", stringsAsFactors = T, sep = ";")



outs$se

# Original pr ggplot script

ggplot(outs, aes(x = setTleaf, y = pr.proxy)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() + facet_wrap(~sp) + ylim(0,1) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)

#--------------------------------------------------------------#
svg(filename = "anet21p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = anet.21p)) +
    theme_bw() +
    xlab(lab_Tleaf) +
    ylab("Anet") +
    geom_point() + facet_wrap(~sp) + ylim(0,30) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()


svg(filename = "anet0p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = anet.0p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("Anet") +
  geom_point() + facet_wrap(~sp) + ylim(0,30) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()



svg(filename = "ETR21p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = ETR.21p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("ETR") +
  geom_point() + facet_wrap(~sp) + ylim(0,160) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()


svg(filename = "ETR0p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = ETR.0p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("ETR") +
  geom_point() + facet_wrap(~sp) + ylim(0,160) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()


svg(filename = "NPQ21p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = NPQ.21p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("NPQ") +
  geom_point() + facet_wrap(~sp) + ylim(-1,-0.997) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()

svg(filename = "NPQ0p.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = setTleaf, y = NPQ.0p)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("NPQ") +
  geom_point() + facet_wrap(~sp) + ylim(-1,-0.997) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)
dev.off()



## photorespiration proportion ----- 
anova(aov(pr.proxy ~  sp * setTleaf, data = outs))
anova(aov(pr.proxy ~  sp, data = outs))
anova(aov(pr.proxy ~  setTleaf, data = outs))

m1 <- nlme::lme(pr.proxy ~  sp * setTleaf, data = outs, random = ~1|treeid)
anova(aov(pr.proxy ~  sp * setTleaf, data = outs))




