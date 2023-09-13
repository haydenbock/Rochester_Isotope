# Package Management -----
# install.packages(c("dplyr", "tidyverse", "vegan", "SIBER"))

library(dplyr)
library(tidyverse)
library(vegan)
library(SIBER)
library(readxl)
library(readr)

#data -----
Env <- read_excel("Environmental_Data.xlsx") 
Isotope <- read_excel("Isotope_Data.xlsx")
Joined_DF <- left_join(Env, Isotope)

#Basic Visualization ----
Full.scatter <- Joined_DF %>% ggplot(aes(x = DELTA_13C_vs_IntStandard, y = DELTA_15N_vs_air, color = Taxa, shape = Urban_Kmeans_Cluster)) +
                  geom_point()
Full.scatter

#Oribatid - omnivore/opportunist
Oribatid.scatter <- Joined_DF %>% filter(Taxa == "Oribatid") %>% ggplot(aes(x = DELTA_13C_vs_IntStandard, y = DELTA_15N_vs_air, color = Urban_Kmeans_Cluster)) +
  geom_point()

Oribatid.scatter 

#Collembola - Fungivore
Collembola.scatter <- Joined_DF %>% filter(Taxa == "Collembola") %>% ggplot(aes(x = DELTA_13C_vs_IntStandard, y = DELTA_15N_vs_air, color = Urban_Kmeans_Cluster)) +
  geom_point()

Collembola.scatter

#Mesostigmata - predator
Mesostigmata.scatter <- Joined_DF %>% filter(Taxa == "Mesostigmata") %>% ggplot(aes(x = DELTA_13C_vs_IntStandard, y = DELTA_15N_vs_air, color = Urban_Kmeans_Cluster)) +
  geom_point()

Mesostigmata.scatter






#Analysis -----

view(Joined_DF)
#plotGroupEllipses

#laymanB - calculates foodweb structure - 
# Layman 2007 - MEAN DISTANCE TO CENTROID HIGHLY IMPORTANT
# TA value indicates tight niche space; narrow diet breadth
#interesting- compare your ellipse to the max/min bayesian generated plot

