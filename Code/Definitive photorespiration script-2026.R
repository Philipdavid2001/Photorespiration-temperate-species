
library(stats)
library(base)
library(dplyr)
library(patchwork)
library(tidyr)
library(lme4)
library(broom)
library(ggplot2)
library(ggpubr)
library(purrr)



# Species #######################
setwd("~/Documents/GitHub/Photorespiration-temperate-species/output2026")
df <- read.csv("~/Documents/GitHub/Photorespiration-temperate-species/Data/Photorespiration temperate/Uppsala-2024-Summer-Photorespiration-SpotMes-TreeSpp.csv", header = T, stringsAsFactors = T, sep = ";")

###### Photorespiration rate calculation loop -----------------------------------------------------

## use the following filters as necessary
df              <- subset(df, dq == "yes") # eliminates unreliable values.
df$pairid       <- paste0(df$sprp, "-", df$tleaf)
df$pairid       <- as.factor(df$pairid)

length(unique(df$pairid))       # individual runs + season
dflist          <- split(df,list(df$pairid))
dflist          <- split(df,unique(list(df$pairid)))
dflist          <- dflist[sapply(dflist, nrow)>0] 

setwd("output")

correct_RD <- function(data, output_path){
  ColumnNames <- c("sp","treeid","Tleaf",
                   "anet.21p",    
                   "anet.0p",
                   "Rp.M1",
                   "phi.M1",
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
                   "Rp.ETR",
                   "phi.ETR")
  
  output_data <- data.frame(matrix(nrow = 0, ncol = length(ColumnNames)))
  names(output_data) <- ColumnNames
  
  get1 <- function(x) ifelse(length(x) == 0, NA, x[1])
  
  for(dataIdx in 1:length(data)){
    
    speciesdata <- data[[dataIdx]]
    
    if(nrow(speciesdata) != 2){
      print(names(data)[dataIdx])
      print(nrow(speciesdata))
      next
    }
    
    p21 <- subset(speciesdata, olev == "21p")
    p0  <- subset(speciesdata, olev == "0p")
    
    if(nrow(p21) != 1 | nrow(p0) != 1){
      print(paste("Skipping pair", names(data)[dataIdx],
                  "because 21p rows =", nrow(p21),
                  "and 0p rows =", nrow(p0)))
      next
    }
    
    # use the factor column directly, then convert to character
    sp     <- as.character(speciesdata$sp[1])
    treeid <- as.character(speciesdata$treeid[1])
    Tleaf  <- get1(speciesdata$Tleaf)
    
    anet.21p <- get1(p21$A)
    anet.0p  <- get1(p0$A)
    
    Rp.M1  <- anet.0p - anet.21p
    phi.M1 <- Rp.M1 / anet.21p
    
    gsw.21p     <- get1(p21$gsw)
    gsw.0p      <- get1(p0$gsw)
    gsw.delta   <- gsw.0p - gsw.21p
    gsw.percent <- gsw.delta / gsw.21p
    
    E.21p     <- get1(p21$E)
    E.0p      <- get1(p0$E)
    E.delta   <- E.0p - E.21p
    E.percent <- E.delta / E.21p
    
    ETR.21p     <- get1(p21$ETR)
    ETR.0p      <- get1(p0$ETR)
    ETR.delta   <- ETR.0p - ETR.21p
    ETR.percent <- ETR.delta / ETR.21p
    
    JT         <- ETR.21p + ETR.0p
    
    Rd_assumed <- 1.0
    Rp.ETR     <- (JT - 4 * (anet.21p + Rd_assumed)) / 12
    phi.ETR    <- Rp.ETR / anet.21p
    
    new_data <- data.frame(
      sp         = sp,
      treeid     = treeid,
      Tleaf      = Tleaf,
      anet.21p   = anet.21p,
      anet.0p    = anet.0p,
      Rp.M1      = Rp.M1,
      phi.M1     = phi.M1,
      gsw.21p    = gsw.21p,
      gsw.0p     = gsw.0p,
      gsw.delta  = gsw.delta,
      gsw.percent= gsw.percent,
      E.21p      = E.21p,
      E.0p       = E.0p,
      E.delta    = E.delta,
      E.percent  = E.percent,
      ETR.21p    = ETR.21p,
      ETR.0p     = ETR.0p,
      ETR.delta  = ETR.delta,
      ETR.percent= ETR.percent,
      JT         = JT,
      Rp.ETR     = Rp.ETR,
      phi.ETR    = phi.ETR
    )
    
    names(new_data) <- names(output_data)
    output_data <- rbind(output_data, new_data)
  }
  
  file_path <- paste(output_path, "Species-output.csv", sep = "/")
  if(!dir.exists(output_path)){
    dir.create(output_path, recursive = TRUE)
  }
  
  index <- 1
  while(file.exists(file_path)){
    file_path <- paste(output_path,
                       paste0("Species-output", index, ".csv"),
                       sep = "/")
    index <- index + 1
  }
  
  write.table(output_data, file = file_path, 
              row.names = FALSE, col.names = TRUE, sep = ",")
}

correct_RD(dflist, "./")


####### plotting output -------


outs <- read.csv("Species-output3.csv", stringsAsFactors = T, sep = ",")
# SLA data from GIFT database 
# Denelle, Pierre, Patrick Weigelt, and Holger Kreft. 2023. “GIFT—An R Package to Access the Global Inventory of Floras and Traits.” Methods in Ecology and Evolution 14 (11): 2738–48.
# Weigelt, Patrick, Christian König, and Holger Kreft. 2020. “GIFT – A Global Inventory of Floras and Traits for Macroecology and Biogeography.” Journal of Biogeography 47 (1): 16–43.



### Species order
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



mod4 <- nlme::lme(pr.real ~  sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)


colnames(outs)[3] <- "tleaf"

sdf <- outs %>%
  group_by(sp, tleaf) %>%                     
  summarise(mean_pr_real = mean(pr.percent, na.rm = TRUE)) %>%  
  arrange(sp, tleaf)                          


# sdf <- subset(sdf, tleaf == "35")

lm_phi_tleaf <- lm(pr.percent ~ tleaf, data = outs)
summary(lm_phi_tleaf)

# diagnostic misc plots
ggplot(outs, aes(x = tleaf, y = pr.percent)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(
       x = "Leaf Temperature (tleaf)",
       y = "Phi") +
  theme_minimal()

ggplot(outs, aes(x = tleaf, y = pr.percent, color = sp)) +
  geom_point(alpha = 0.006) +
  geom_smooth(method = "lm", se = FALSE) 


lm_no_interaction <- lm(pr.percent ~ tleaf + sp, data = outs)
lm_interaction <- lm(pr.percent ~ tleaf * sp, data = outs)
summary(lm_no_interaction)
anova(lm_no_interaction, lm_interaction)

prslp <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(pr.percent ~ tleaf, data = .))) %>%
  filter(term == "tleaf") %>%
  dplyr::select(sp, estimate, std.error, p.value) %>%
  rename(slope = estimate, se = std.error) %>%
  arrange(desc(abs(slope)))  

anetslp <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(anet.21p ~ tleaf, data = .))) %>%
  filter(term == "tleaf") %>%
  dplyr::select(sp, estimate, std.error, p.value) %>%
  rename(slope = estimate, se = std.error) %>%
  arrange(desc(abs(slope)))  



# Combine by species
combined_slopes <- left_join(prslp, anetslp, by = "sp")

sla_data <- tribble(
  ~sp, ~sla,
  "Betula pubescens",        144,
  "Scandosorbus intermedia", 157,
  "Fagus sylvatica",         174,
  "Betula pendula",          179,
  "Corylus avellana",        194,
  "Acer platanoides",        197,
  "Tilia cordata",           272
) # data from GIFT

cdb <- combined_slopes %>%
  left_join(sla_data, by = "sp")


cdb$sp <- as.factor(cdb$sp)

mod4 <- nlme::lme(pr_slope ~  sla, data = cdb, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(anet_slope ~  sla, data = cdb, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

plot(anet_slope ~  sla, data = cdb)
plot(pr_slope ~  sla, data = cdb)

plot(pr_slope ~  sla, data = subset(cdb, sla < 250))

cc <- subset(cdb, sp != "Scandosorbus intermedia")

mod4 <- nlme::lme(pr_slope ~  sla, data = subset(cdb, sla < 210), 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)



anova(lm(pr_slope ~  sla, data = subset(cdb, sla < 210)))
cor.test(~ pr_slope + sla, data = subset(cdb, sla < 210))

# Pearsons R = 0.882, P = 0.019 

anova(lm(pr_slope ~  sla, data = subset(cdb, sla < 2100)))
cor.test(~ pr_slope + sla, data = subset(cdb, sla < 2100))
# Pearsons R = 0.138, P = 0.76 
# the relationship does not show up if we include Tilia cordata which has the highest SLA of 272

anova(lm(anet_p ~  sla, data = subset(cdb, sla < 210)))
cor.test(~ anet_p + sla, data = subset(cdb, sla < 210))
anova(lm(anet_p ~  sla, data = subset(cdb, sla < 2100)))
cor.test(~ anet_p + sla, data = subset(cdb, sla < 2100))
# 
# No relationship with slope of Anet

anova(lm(pr.real ~  sla, data = subset(outs, sla < 2100)))
cor.test(~ pr.real + sla, data = subset(outs, sla < 2100))
# Pearsons R = -0.21, P = 0.012 

anova(lm(pr.real ~  sla, data = subset(outs, sla < 210)))
cor.test(~ pr.real + sla, data = subset(outs, sla < 210))
# the relationship weakens if we exclude Tilia cordata which has the highest SLA of 272


anova(lm(anet.21p ~  sla, data = subset(outs, sla < 2100)))
cor.test(~ anet.21p + sla, data = subset(outs, sla < 2100))
# Pearsons R = -0.235, P = 0.0061 

anova(lm(anet.21p ~  sla, data = subset(outs, sla < 210)))
cor.test(~ anet.21p + sla, data = subset(outs, sla < 210))
# the relationship weakens if we exclude Tilia cordata which has the highest SLA of 272

ggplot(outs, aes(x = sla, y = pr.real)) +
  geom_point(alpha = 0.6, color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    x = "SLA (cm²/g)",
    y = "Anet (21p)",
    title = "Scatterplot of Anet vs SLA with Linear Fit"
  ) +
  theme_minimal(base_size = 14)

cor_val <- cor(outs$sla, outs$pr.real, use = "complete.obs")

ggplot(outs, aes(x = sla, y = pr.real, color = tleaf)) +
  geom_point(alpha = 0.6,  size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  annotate("text", x = max(outs$sla) * 0.8, y = max(outs$anet.21p), 
           label = paste0("r = ", round(cor_val, 3)), size = 5, color = "darkred") +
  labs(
    x = "SLA (cm²/g)",
    y = "Anet (21p)"
  ) +
  theme_bw(base_size = 14)

outs <- outs %>%
  mutate(tleaf = factor(tleaf, levels = c(25, 30, 35), labels = c("Low", "Medium", "High")))


ggplot(outs, aes(x = sla, y = pr.real, color = tleaf)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  annotate("text", x = max(outs$sla) * 0.8, y = max(outs$pr.real), 
           label = paste0("r = ", round(cor_val, 3)), size = 5, color = "darkred") +
  scale_color_manual(values = c("blue", "grey20", "orange")) +
  labs(
    x = "SLA (cm²/g)",
    y = "Anet (21p)",
    color = "Leaf Temp"
  ) + facet_wrap(~tleaf) +
  theme_bw(base_size = 14)


# Step 1: Compute r and slope per group
labels_df <- outs %>%
  group_by(tleaf) %>%
  summarise(
    r = cor(sla, pr.real, use = "complete.obs"),
    lm_fit = list(lm(pr.real ~ sla)),
    .groups = "drop"
  ) %>%
  mutate(
    slope = sapply(lm_fit, function(mod) coef(mod)[2]),
    pval = sapply(lm_fit, function(mod) summary(mod)$coefficients[2, 4]),
    label = paste0(
      "r = ", round(r, 2),
      "\nslope = ", round(slope, 3),
      "\np = ", format.pval(pval, digits = 2, eps = 0.001)
    ),
    x = max(outs$sla) * 0.7,
    y = max(outs$pr.real) * 0.9
  )

outs$tleaf <- as.factor(outs$tleaf)
ggplot(outs, aes(x = sla, y = pr.real, color = tleaf)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  geom_text(data = labels_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 5, color = "darkred") +
  scale_color_manual(values = c("blue", "grey20", "orange")) +
  labs(
    x = "SLA (cm²/g)",
    y = "Anet (21p)",
    color = "Leaf Temp"
  ) +
  facet_wrap(~tleaf) +
  theme_bw(base_size = 14)

# interaction model
interaction_model <- lm(pr.real ~ sla * tleaf, data = outs)
summary(interaction_model)






mod_interaction <- lm(pr.percent ~ tleaf * sp, data = outs)

mod_main_effects <- lm(pr.percent ~ tleaf + sp, data = outs)

anova(mod_main_effects, mod_interaction)


## Main Figure 1 ---------------------------------


come back here

## plots updateds to add anova p value for each panel

# next one works - this does not. 
outs$tleaf <- as.numeric(outs$tleaf)
pval_df <- outs %>%
  group_by(sp) %>%
  filter(n_distinct(tleaf) >= 2) %>%  # Ensure valid ANOVA
  summarise(
    pval = tryCatch({
      summary(aov(anet.21p ~ factor(tleaf)))[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA)
  ) %>%
  mutate(label = paste0("P = ", signif(pval, 3)))

outs %>%
  group_by(sp) %>%
  summarise(test = 1)




pval_df <- outs %>%
  group_by(sp) %>%
  filter(n_distinct(tleaf) >= 2) %>%
  summarise(
    pval = tryCatch({
      sm <- summary(aov(anet.21p ~ factor(tleaf)))
      # Only extract p-value if column exists
      if ("Pr(>F)" %in% colnames(sm[[1]])) {
        sm[[1]][["Pr(>F)"]][1]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(is.na(pval), "P = NA", paste0("P = ", signif(pval, 3))))


annot_df <- outs %>%
  group_by(sp) %>%
  summarise(ypos = max(anet.21p, na.rm = TRUE) + 1) %>%
  left_join(pval_df, by = "sp")

plot1 <- ggplot(outs, aes(x = tleaf, y = anet.21p)) + 
  geom_boxplot(aes(group = factor(tleaf)), fill = "#6666ff", alpha = .4, 
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
  filter(n_distinct(tleaf) >= 2) %>%  # Ensure valid ANOVA
  summarise(
    pval = tryCatch({
      summary(aov(pr.real ~ factor(tleaf)))[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA)
  ) %>%
  mutate(label = paste0("P = ", signif(pval, 3)))

annot_df <- outs %>%
  group_by(sp) %>%
  summarise(ypos = max(pr.real, na.rm = TRUE) + 1) %>%
  left_join(pval_df, by = "sp")


plot2 <- ggplot(outs, aes(x = tleaf, y = pr.real)) + 
  geom_boxplot(aes(group = factor(tleaf)), fill = "#FF9900", alpha = 0.7, 
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
  filter(n_distinct(tleaf) >= 2) %>%  # Ensure valid ANOVA
  summarise(
    pval = tryCatch({
      summary(aov(pr.percent ~ factor(tleaf)))[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA)
  ) %>%
  mutate(label = paste0("P = ", signif(pval, 3)))

annot_df <- outs %>%
  group_by(sp) %>%
  summarise(ypos = max(pr.percent, na.rm = TRUE) + 1) %>%
  left_join(pval_df, by = "sp")


plot3 <- ggplot(outs, aes(x = tleaf, y = pr.percent)) + 
  geom_boxplot(aes(group = factor(tleaf)), fill = "#FFFF99", alpha = 0.5, 
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

svg(filename = "MAIN-Figure1c.svg", width = 10, height = 10, bg = "transparent")

(plot2 / plot1 / plot3) + plot_layout(nrow = 3)

dev.off()


###------- Mixed effects models: Photorespiration ----

mod4 <- nlme::lme(pr.real ~  sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(pr.real ~  tleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

mod4 <- nlme::lme(pr.real ~  sp * tleaf, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)

###------- Mixed effects models: Rp/Anet (Phi) ----


mod6 <- nlme::lme(pr.percent ~  sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)

mod6 <- nlme::lme(pr.percent ~  tleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


mod6 <- nlme::lme(pr.percent ~  sp * tleaf, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)

###------- Mixed effects models: Photosynthesis (Anet) ----

mod5 <- nlme::lme(anet.21p ~  sp, data = outs, 
                  random = ~1| treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)

mod5 <- nlme::lme(anet.21p ~  tleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)


mod6 <- nlme::lme(anet.21p ~  sp * tleaf, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


## exclude B. pubsecense and SCAINT 
range(outs$anet.21p)
outsx <- subset(outs, sp !="Scandosorbus intermedia")
outsx <- subset(outsx, sp !="Betula pubescens")


mod5 <- nlme::lme(pr.real ~  sp * tleaf, data = outsx, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)

###------- Mixed effects models: Photosynthesis and Photorespiration ----

mod5 <- nlme::lme(pr.real ~   anet.21p, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod5)

mod4a <- nlme::lme(pr.real ~   anet.21p * tleaf, data = outs, 
                   random = ~1|sp, 
                   method = "REML", 
                   na.action=na.omit) ; anova(mod4a)


mod4 <- nlme::lme(pr.real ~   anet.21p * tleaf * sp, data = outs, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)


anova(mod4, mod4a)

### plotting averages ----

pt <- outs %>%
  group_by(sp, tleaf) %>%
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


pt25 <- subset(outs, tleaf == 25)
pt30 <- subset(outs, tleaf == 30)
pt35 <- subset(outs, tleaf == 35)

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

model <- lmer(pr.real ~ anet.21p + tleaf + (1 | sp), data = outs)
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
  do(tidy(lm(pr.real ~ tleaf, data = .)))

results2 <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(pr.percent ~ tleaf, data = .)))


results3 <- outs %>%
  group_by(sp) %>%
  do(tidy(lm(anet.21p ~ tleaf, data = .)))


outs_summary <- outs %>%
  group_by(sp, tleaf) %>%
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
  group_by(sp, tleaf) %>%
  summarise(ETR.21p.avg = mean(ETR.21p),
            se = sd(ETR.21p)/ sqrt(length(ETR.21p)))


ggplot(ETR.table, aes(y = ETR.21p.avg, x=tleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = tleaf, y = ETR.21p, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=ETR.21p.avg-se, ymax=ETR.21p.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "limegreen", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = tleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
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
  group_by(sp, tleaf) %>%
  summarise(ETR.0p.avg = mean(ETR.0p),
            se = sd(ETR.0p)/ sqrt(length(ETR.0p)))

ggplot(ETRRp.table, aes(y = ETR.0p.avg, x=tleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = tleaf, y = gsw.0p, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=ETR.0p.avg-se, ymax=ETR.0p.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "cornflowerblue", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = tleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
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
  group_by(sp, tleaf) %>%
  summarise(ETR.percent.avg = mean(ETR.percent),
            se = sd(ETR.percent)/ sqrt(length(ETR.percent)))
ggplot(ETRpercent.table, aes(y = ETR.percent.avg, x=tleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = tleaf, y = ETR.percent, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=ETR.percent.avg-se, ymax=ETR.percent.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "red3", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = tleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
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
  group_by(sp, tleaf) %>%
  summarise(pr.real.avg = mean(pr.real),
            se = sd(pr.real)/ sqrt(length(pr.real)))

ggplot(absolute.table, aes(y = pr.real.avg, x=tleaf, group=sp))+
  geom_errorbar(aes(ymin=pr.real.avg-se, ymax=pr.real.avg+se),col ="grey70", width= 2.5)+ #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  geom_point(col="grey30", size = 3, pch = 21, fill = "cornflowerblue", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = tleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
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
  group_by(sp, tleaf) %>%
  summarise(anet.21p.avg = mean(anet.21p),
            se = sd(anet.21p) / sqrt(length(anet.21p)))

ggplot(anet.table, aes(y = anet.21p.avg, x=tleaf, group=sp))+
  geom_smooth(se = F, method="lm", col = "grey80")+
  geom_errorbar(aes(ymin=anet.21p.avg-se, ymax=anet.21p.avg+se),col ="grey70", width= 2.5) +
  geom_point(col="grey30", size = 3, pch = 21, fill = "limegreen", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = tleaf, y=anet.21p, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
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
  group_by(sp, tleaf) %>%
  summarise(pr.percent.avg = mean(pr.percent),
            se = sd(pr.percent)/ sqrt(length(pr.percent)))

ggplot(percent.table, aes(y = pr.percent.avg, x=tleaf, group=sp))+
  #pmin can be used to cap the ymax. (e.g., pmin(pr.real.avg+se, 1.0))
  geom_smooth(se = F, method="lm", col = "grey80")+
  # ggpmisc::stat_poly_eq(
  #   data = outs,
  #   formula = y ~ x,
  #   aes(x = tleaf, y = pr.percent, label = paste(after_stat(p.value.label), sep = "~~~")),
  #   label.x = 22,
  #   label.y = 1.5,
  # )+
  geom_errorbar(aes(ymin=pr.percent.avg-se, ymax=pr.percent.avg+se),col ="grey70", width= 2.5)+
  geom_point(col="grey30", size = 3, pch = 21, fill = "red3", stroke = 0.85)+
  # stat_regline_equation(data = outs, aes(x = tleaf, y=pr.real, label = paste(..eq.label.., ..adj.rr.label.., paste("p = ", ..p.value..), sep = "~~~")),
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

outs$tl <- as.factor(outs$tleaf)

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
  geom_jitter(aes(size = tl), color = as.factor(outs$tleaf), alpha = 0.5) +
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
  geom_jitter(aes(size = tl), color = as.factor(outs$tleaf), alpha = 0.5) +
  xlab(" ") +
  scale_y_continuous(limits = c(0,25), 
                     name = expression(paste(italic(A)[net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')')))+
  ggthemes::theme_base()+
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        panel.border = element_rect(color = "grey70")) +
  coord_flip()
dev.off()



ggplot(pt, aes(y = pr, x=tleaf)) +
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
  facet_wrap(~tleaf)
  geom_errorbar(aes(ymin=pr.real-prse, ymax=pr.real+prse)) +
  geom_errorbarh(aes(xmin=pr.real-psse, xmax=pr.real+psse)) + 

dev.off()

  
  ggplot(outs, aes(y = pr.real, x=tleaf)) +
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

palette_colors <- c("#000000", "#E69F00", "blue", "darkred", "skyblue", "#0072B2", "#D55E00")

  
svg(filename = "Figure-2new.svg", width = 9, height = 3, bg = "transparent")
ggplot(outs_summary, aes(x = anet_mean, y = pr_mean, color = sp)) +
  geom_smooth(se = F, method = "lm", col = "grey30") +  
  geom_errorbar(data = outs_summary, aes(ymin = pr_mean - pr_se, ymax = pr_mean + pr_se), width = 0, color = "grey60") +
  geom_errorbarh(data = outs_summary, aes(xmin = anet_mean - anet_se, xmax = anet_mean + anet_se), height = 0, color = "grey60") +
  geom_point(size = 2, stroke = 0.5, aes(color = sp, shape = sp)) +
  geom_point(data = outs_summary, aes(x = anet_mean, y = pr_mean, shape = sp), 
             size = 2, stroke = 1) +
  scale_x_continuous(limits = c(0, 20), 
                     name = expression(paste(italic(A)[Net], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  scale_y_continuous(limits = c(0, 20), 
                     name = expression(paste(italic(R)[p], ' (', mu * ~'mol'~ " CO"[2]~' m'^{-2}*' s'^{-1}*')'))) +
  ggthemes::theme_base() +
  facet_wrap(~tleaf) +
  
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        panel.border = element_rect(color = "grey70")) +
  scale_shape_manual(values = c(21,22,23,24,25,1,2))+
  scale_color_manual(values = palette_colors)+
  geom_abline(slope = 1, linetype = "dashed", color = "grey50")

dev.off()

#--------------------------------------------------------------------#

### ETR analysis starts here



### Photorespiration ratio plotted over ratio of ETR used in each environment

svg(filename = "Phi-ETRdelta.svg", width = 10, height = 3.3, bg = "transparent")

ggplot(pt, aes(y = ETR.delta, x=pr, color = sp)) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) +
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
  facet_wrap(~tleaf)

dev.off()



#--------------------- test ETR over ETR
svg(filename = "ETR.percent-over-Rp.svg", width = 7, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = ETR.percent, x=pr, color = sp)) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) +
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
  facet_wrap(~tleaf)
dev.off()


##Anet

svg(filename = "ETR.21p-tleaf.svg", width = 7, height = 4.5, bg = "transparent")

ggplot(pt, aes(y = ETR.21p, x=tleaf, color = sp)) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) +
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
ggplot(pt, aes(y = ETR.0p, x=tleaf, color = sp)) +
  geom_boxplot(aes(group = cut_width(tleaf, 1))) +
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
outs$st <- as.factor(outs$tleaf)
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
outs$st <- as.factor(outs$tleaf)
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
outs$st <- as.factor(outs$tleaf)
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



mod6 <- nlme::lme(ETR.percent ~  phi * tleaf, data = pt, 
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
outs$st <- as.factor(outs$tleaf)
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
  facet_wrap(~tleaf)
dev.off()


######  Creating the dataframes for each temperature for the models ----

outs25 <- subset(outs, tleaf == 25)
outs30 <- subset(outs, tleaf == 30)
outs35 <- subset(outs, tleaf == 35)


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

# Checking the individual significance of species temperature response.

pt <- outs %>%
  group_by(sp, treeid, tleaf) %>%
  summarise(pr = mean(pr.real),
            prse = sd(pr.real)/sqrt(length(pr.real)),
            ps = mean(anet.21p),
            psse = sd(anet.21p)/ sqrt(length(anet.21p)),
            phi = mean(pr.percent),
            phise = sd(pr.percent)/sqrt(length(pr.percent)))


###PHOTORESPIRATION RATES / TEMP

Betulapen <- subset(pt, sp == "Betula pendula")

mod6 <- nlme::lme(pr ~  tleaf, data = Betulapen, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Fagus <- subset(pt, sp == "Fagus sylvatica")

mod6 <- nlme::lme(pr ~  tleaf, data = Fagus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Acer <- subset(pt, sp == "Acer platanoides")

mod6 <- nlme::lme(pr ~  tleaf, data = Acer, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Batulapub <- subset(pt, sp == "Betula pubescens")

mod6 <- nlme::lme(pr ~  tleaf, data = Batulapub, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Tilia <- subset(pt, sp == "Tilia cordata")

mod6 <- nlme::lme(pr ~  tleaf, data = Tilia, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Corylus <- subset(pt, sp == "Corylus avellana")

mod6 <- nlme::lme(pr ~  tleaf, data = Corylus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Sorbus <- subset(pt, sp == "Scandosorbus intermedia")


mod6 <- nlme::lme(pr ~  tleaf, data = Sorbus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)




###PHOTOSYNTHESIS RATES / TEMP

Betulapen <- subset(pt, sp == "Betula pendula")

mod6 <- nlme::lme(ps ~  tleaf, data = Betulapen, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Fagus <- subset(pt, sp == "Fagus sylvatica")

mod6 <- nlme::lme(ps ~  tleaf, data = Fagus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Acer <- subset(pt, sp == "Acer platanoides")

mod6 <- nlme::lme(ps ~  tleaf, data = Acer, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Batulapub <- subset(pt, sp == "Betula pubescens")

mod6 <- nlme::lme(ps ~  tleaf, data = Batulapub, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Tilia <- subset(pt, sp == "Tilia cordata")

mod6 <- nlme::lme(ps ~  tleaf, data = Tilia, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Corylus <- subset(pt, sp == "Corylus avellana")

mod6 <- nlme::lme(ps ~  tleaf, data = Corylus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


Sorbus <- subset(pt, sp == "Scandosorbus intermedia")


mod6 <- nlme::lme(ps ~  tleaf, data = Sorbus, 
                  random = ~1|treeid, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod6)


# Creating a sheet with each temperature point containing one compiled value for each species. This table is used when gathering the highest/lowest mean differences in Anet and Rp across Tleaf. 

pt <- outs %>%
  group_by(sp, tleaf) %>%
  summarise(pr = mean(pr.real),
            prse = sd(pr.real)/sqrt(length(pr.real)),
            ps = mean(anet.21p),
            psse = sd(anet.21p)/ sqrt(length(anet.21p)))

Sorbus <- subset(pt, sp == "Scandosorbus intermedia")


mod6 <- nlme::lme(pr ~  tleaf, data = Sorbus, 
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
  facet_wrap(~tleaf)
dev.off()













#-------------------------------#



#Box plots of pr.percent across the leaf temperatures

ggplot(outs, aes(y = pr.percent, x=tl, color=tl)) +
  geom_boxplot(aes(x = tl, y = pr.percent, fill = tleaf)) +
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

# Here we are running a for loop testing     the relationship between species, tleaf and the combined effect of them (I.e the temperature response of species), against pr.real (Absolute photorespiration rates ("10")), pr.percent (Proportional photorespiration rates ("11")), anet.21p (Ambient photosynthesis rates ("4")), anet.0p (O2-free photosynthesis rates ("5")) and JO2.percent (proportional oxygenation rates calculated using ETR only ("31)).


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
  
  mod1 <- nlme::lme(independentVariable ~  tleaf , data = outs, 
                    random = ~1|sp, 
                    method = "REML", 
                    na.action=na.omit) ; anova(mod1)
  
  print(paste("This analysis is based on :", names(outs)[colindeces[index]]))
  print(summary(mod1))
  
  mod1 <- nlme::lme(independentVariable ~  sp * tleaf, 
                    data = outs, 
                    random = ~1|treeid, 
                    method = "REML", 
                    na.action=na.omit) ; anova(mod1)
  
  print(paste("This analysis is based on :", names(outs)[colindeces[index]]))
  print(summary(mod1))
  
  
  print("----------------------------------------------------------------")
}


outs25 <- subset(outs, tleaf == 25)
outs30 <- subset(outs, tleaf == 30)
outs35 <- subset(outs, tleaf == 35)



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
 

mod4 <- nlme::lme(pr.real ~  anet.21p * tleaf, data = outs, 
                  random = ~1|sp, 
                  method = "REML", 
                  na.action=na.omit) ; anova(mod4)




long_data <- outs %>%
  gather(key = "method", value = "y", pr.percent, JO2.percent)

# Fit a model with an interaction term between method and tleaf
long_model <- lm(y ~ sp * tleaf * method, data = long_data)

# Check the model summary
summary(long_model)


# Fit the model with and without the interaction term between method and tleaf
model_interaction <- lm(y ~ sp * tleaf * method, data = long_data)
model_no_interaction <- lm(y ~ sp * tleaf + method, data = long_data)

# Compare the models using a likelihood ratio test
anova(model_no_interaction, model_interaction)

# Visualize the comparison

ggplot(long_data, aes(x = tleaf, y = y, color = method)) +
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





