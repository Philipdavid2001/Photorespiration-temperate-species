install.packages("GIFT", dependencies = T)
library("GIFT")

library("GIFT")
library("dplyr")
library("knitr")
library("kableExtra")
library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("tidyr")
library("ape")
library("patchwork")



all_sp <- GIFT_species()

species_lookup <- GIFT_species_lookup(genus = "Acer", epithet = "platanoides")

trait_meta <- GIFT_traits_meta()

trait_data <- GIFT_traits(trait_IDs = c("4.1.3"), agreement = 0.66, GIFT_version = "latest")

# traits ID of potential interest
tid   <- c("1.6.2", "1.6.3", "2.2.1", "2.4.1", "4.1.1", "4.1.2", "4.1.3", "4.2.1")
tname <- c("Plant_height_max", "Plant_height_mean", "Lifespan_1", "Deciduousness_1", "SLA_min", "SLA_max", "SLA_mean", "Photosynthetic_pathway") 
tunit <- c("m", "m", "years", "categorical", "cm2/g","cm2/g","cm2/g","categorical")
         



prsp <- c("Fagus sylvatica", "Scandosorbus intermedia", "Betula pubescens", 
                         "Betula pendula", "Acer platanoides", "Tilia cordata", "Corylus avellana")

tid <- c("1.6.2", "1.6.3", "2.2.1",  "4.1.1", "4.1.2", "4.1.3", "4.2.1")
tname <- c("Plant_height_max", "Plant_height_mean", "Lifespan_1", "Deciduousness_1", "SLA_min", "SLA_max", "SLA_mean", "Photosynthetic_pathway") 


tid <- c("1.6.2", "1.6.3",  "4.1.3")
tname <- c("Plant_height_max", "Plant_height_mean","SLA_mean") 
# std <- std[ , -4]

colnames(std) <- c("sp", "Plant_height_max", "Plant_height_mean","SLA_mean") 




