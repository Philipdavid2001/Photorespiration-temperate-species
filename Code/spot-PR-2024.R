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

setwd("C:/Data/Masterdata/Excel/output")

#setwd("~/GoogleDriveRT/Uppsala-2023/Philip-photorespiration/Uppsala-2024-June-Campaign/output/")

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



write.table(outs, file = "Uppsala-PR-spot-2024-output.csv", 
            row.names = FALSE, col.names = T, sep = ",")
outs <- read.csv("Uppsala-PR-spot-2024-output.csv", stringsAsFactors = T)

outs <- read.csv("Uppsala-PR-spot-2024-output.csv", stringsAsFactors = T, sep = ";")



outs$se

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


ggplot(outs, aes(x = setTleaf, y = Ca.prop)) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("Ca.prop") +
  geom_point() + facet_wrap(~sp) + ylim(-0.008,0.005) + xlim(20,40) +
  #geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + 
  geom_smooth(method = "lm", se = F)



## photorespiration proportion ----- 
anova(aov(pr.proxy ~  sp * setTleaf, data = outs))
anova(aov(pr.proxy ~  sp, data = outs))
anova(aov(pr.proxy ~  setTleaf, data = outs))

m1 <- nlme::lme(pr.proxy ~  sp * setTleaf, data = outs, random = ~1|treeid)
anova(aov(pr.proxy ~  sp * setTleaf, data = outs))

## etr proportion ----- 

ggplot(outs, aes(x = setTleaf, y = E.delta)) +
    theme_bw() +
    xlab(lab_Tleaf) +
    ylab("PR") +
    geom_boxplot(outlier.shape = NA) + facet_wrap(~sp) + ylim(0,1) + xlim(20,40)  


## etr proportion ----- 

outs2 <- read.csv("Uppsala-PR-spot-2024-output.csv", stringsAsFactors = T)


length(unique(outs2$sp))       # individual runs + season

dflist    <-   split(outs2,list(outs2$sp))
dflist    <-   split(outs2,unique(list(outs2$sp)))
dflist    <-   dflist[sapply(dflist, nrow)>0] 

cat <- NULL

for(i in 1:length(dflist)) {
  tempdf          <-          (dflist[[i]])
  sp              <-          tempdf$sp[1]; sp
  mod <- lm(pr.proxy ~ setTleaf, tempdf)
  modsum <- summary(mod)
  
  
  int        <- modsum$coefficients[1,1]; int
  int.er     <- modsum$coefficients[1,2]; int.er  
  int.tval   <- modsum$coefficients[1,3]; int.tval 
  int.pval   <- modsum$coefficients[1,4]; int.pval 
   
  slope      <- modsum$coefficients[2,1]; slope
  slope.er   <- modsum$coefficients[2,2]; slope.er  
  slope.tval <- modsum$coefficients[2,3]; slope.tval 
  slope.pval <- modsum$coefficients[2,4]; slope.pval 
  
  adj.r2     <- modsum$adj.r.squared; adj.r2
  
  params <- rbind(c(as.character(sp),
                    int,       
                    int.er,    
                    int.tval,  
                    int.pval,  
                    slope,     
                    slope.er,  
                    slope.tval,
                    slope.pval,
                    adj.r2))
  
  colnames(params) <- c("sp","int",       
                        "int.er",    
                        "int.tval",  
                        "int.pval",  
                        "slope",     
                        "slope.er",  
                        "slope.tval",
                        "slope.pval",
                        "adj.r2")
  cat = rbind(cat, params)
  print(i)
  
}

cat <- as.data.frame(cat)

# outs <- subset(outs, sp != "scaint")

ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() + geom_smooth(method = "lm", se = F) + 
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + facet_wrap(~sp) 
  
ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() +  
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2)


svg(filename = "p1g.svg", width = 5, height = 4, bg = "transparent")
ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() +  
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2)
dev.off()  


svg(filename = "p2g.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() +  
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2)+ facet_wrap(~sp) 
dev.off()  


svg(filename = "p3g.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() + geom_smooth(method = "lm", se = T) + 
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + facet_wrap(~sp) 
dev.off() 

svg(filename = "p4g.svg", width = 5, height = 4, bg = "transparent")
ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() +  
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) +  
  geom_smooth(method = "lm", se = F)
dev.off() 


svg(filename = "p5g.svg", width = 7, height = 4, bg = "transparent")
ggplot(outs, aes(x = tleaf, y = pr, col = sp) ) +
  theme_bw() +
  xlab(lab_Tleaf) +
  ylab("PR") +
  geom_point() + geom_smooth(method = "lm", se = F) + 
  geom_errorbar(aes(ymin=pr-se, ymax=pr+se), width=.2) + facet_wrap(~sp) 
dev.off() 

