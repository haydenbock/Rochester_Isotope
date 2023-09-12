# Package Management -----
install.packages(c("dplyr", "tidyverse", "vegan", "SIBER"))

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

#Basic Plots ----


