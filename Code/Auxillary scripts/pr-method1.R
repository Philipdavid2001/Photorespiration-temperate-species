
packages_to_load <- c("drc", "bbmle", "labelled", "MALDIquant", 
                      "magicfor", "ggplot2", "tidyverse", "broom", 
                      "Hmisc", "plotly", "PairedData", "DescTools", 
                      "generalhoslem", "graphics", "nls2", "plyr",
                      "ggthemes", "ggpubr", "signal", "segmented")

lapply(packages_to_load, library, character.only = TRUE)

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


setwd("") ### add this --- 

setwd("~/GoogleDriveRT/Uppsala-2023/Philip-photorespiration/PR-TreeSpp-Uppsala-2024-Spot/output/")

# df <- read.csv("Uppsala-2024-Summer-Photorespiration-SpotMes-TreeSpp (1).csv", header = T, stringsAsFactors = T)

df              <- subset(df, dq == "yes") # eliminates unreliable values.
df$pairid       <- paste0(df$sprp, "-", df$setTleaf)
df$pairid       <- as.factor(df$pairid)
length(unique(df$pairid))       # diagnostic
dflist          <-   split(df,list(df$pairid))
dflist          <-   split(df,unique(list(df$pairid)))
dflist          <-   dflist[sapply(dflist, nrow)>0] 


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
  p21             <-          subset(tempdf, olev == "21p"); p21    
  p0              <-          subset(tempdf, olev != "21p"); p0
  
  anet.21p        <-          p21$A; anet.21p   
  anet.0p         <-          p0$A; anet.0p    
  anet.delta      <-          (p0$A - p21$A); anet.delta 
  pr              <-          (0.5 * anet.delta) + anet.delta ; pr
  
  dETR           <-          (p0$ETR - p21$ETR) ; dETR  
  dgs            <-          (p0$gsw - p21$gsw) ; dgs  
  dE             <-          (p0$E - p21$E) ; dETR  
  
  
  params <- rbind(c(as.character(sp),
                    as.character(treeid),
                    as.character(setTleaf), 
                    as.numeric(p21$A),    
                    p0$A,     
                    (p0$A - p21$A),  
                    pr,
                    p21$gsw,
                    p0$gsw,
                    p21$ETR,
                    p0$ETR,
                    p21$E,
                    p0$E)); params
  
  colnames(params) <- c("sp","treeid","setTleaf", 
                        "anet.21p",    
                        "anet.0p",     
                        "anet.delta",  
                        "pr",    
                        "gsw.21p",    
                        "gsw.0p",     
                        "ETR.21p",     
                        "ETR.0p",      
                        "E.21p",     
                        "E.0p")
  outs = rbind(outs, params)
  print(i)
}


write.table(outs, file = "Uppsala-PR-spot-2024-output.csv", 
            row.names = FALSE, col.names = T, sep = ",")

outs <- read.csv("Uppsala-PR-spot-2024-output.csv", stringsAsFactors = T)

outs$setTleaf1 <- as.factor(outs$setTleaf)

outs$sp <- factor(outs$sp, levels = c("Corylus avellana","Tilia cordata", "Acer platanoides",
                                      "Fagus sylvatica","Scandosorbus intermedia",
                                      "Betula pubescens","Betula pendula"))


## example model fits nlme ----

mod1 <- nlme::lme(prr ~  sp , data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod1)

mod1 <- nlme::lme(prr ~  setTleaf , data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod1)


mod1 <- nlme::lme(prr ~  sp * setTleaf, 
                  data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod1)


### post hoc test ------- HSD honest significant difference 

m7 <- aov(pr ~ sp, data = outs35, random = ~1|treeid, method = "REML", na.action=na.omit) ; anova(m7)
par(las=2)
par(mar=c(8,8,2,2)) 
plot(TukeyHSD(m7))


t1 <- TukeyHSD(m7)
t1 <- as.data.frame(t1$sp)
t1$com <- rownames(t1)

rownames(t1) <- NULL
colnames(t1) <- c("prdiff","lwr", "upr","padj", "com")
blah <- subset(t1, padj <= 0.05)



