# Rochester Isotope 
# New script - 2.26.24

#Package Management ----
# install.packages(c("dplyr", "tidyverse", "vegan", "SIBER"))

#install.packages("rjags")
#install.packages("rcartocolor")
#install.packages("ggdist")
#install.packages("ggridges")
#install.packages("ggbeeswarm")
#install.packages("gghalves")
#install.packages("agricolae")
#install.packages("bayesanova")
#install.packages("Bolstad")
#install.packages("BayesFactor")
#install.packages("coda")
#install.packages("bayesplot")

library(brms)
library(bayesplot)
library(coda)
library(ggpubr)
library(rjags)
library(dplyr)
library(vegan)
library(SIBER)
library(readxl)
library(readr)
library(tidyverse)
library(viridis)
library(tidyverse)     
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)
library(agricolae)
library(bayesanova)
library(Bolstad)
library(BayesFactor) 
library(bridgesampling)

#Import Data ----
Env <- read_excel("Environmental_Data.xlsx") 
Isotope <- read_excel("Isotope_Data.xlsx")
Joined_DF <- left_join(Env, Isotope)
My_SIBER_Data <- read_csv("My_SIBER_Data.csv") #iso1 == "∆13C", iso2 = "∆15N", group = "taxa" (1=oribatid, 2=collembola, 3=mesostigmata), community = "urbanization" (1=High, 2=Low),


#Generate Null Model ----
# Load necessary library
library(dplyr)

# Read the original CSV file
df <- My_SIBER_Data

# Find the range for iso1 and iso2
iso1_range <- range(df$iso1, na.rm = TRUE)
iso2_range <- range(df$iso2, na.rm = TRUE)

# Create a null data frame with the same structure
null_df <- df %>%
  mutate(iso1 = runif(n(), iso1_range[1], iso1_range[2]),
         iso2 = runif(n(), iso2_range[1], iso2_range[2]))

# Optionally, you can reset the group and community if needed, but as per instruction, it seems you want to keep them as is but just randomize iso1 and iso2.

# If you want to write the null data frame to a new CSV file
write.csv(null_df, "Null_My_SIBER_Data.csv", row.names = FALSE)

# SIBER ANALYSIS ----

#RESOURCE: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html

# remove previously loaded items from the current environment and remove previous graphics.
rm(list=ls())
graphics.off()

# Here, I set the seed each time so that the results are comparable. 
# This is useful as it means that anyone that runs your code, *should*
# get the same results as you, although random number generators change 
# from time to time.
set.seed(1)


# load in the included demonstration dataset
My_SIBER_Data <- read_csv("My_SIBER_Data.csv")
My_SIBER_Data <- as_tibble(My_SIBER_Data)
My_SIBER_Data <- as.data.frame(My_SIBER_Data)



# create the siber object
siber.example <- createSiberObject(My_SIBER_Data)



# Or if working with your own data read in from a *.csv file, you would use
# This *.csv file is included with this package. To find its location
# type
# fname <- system.file("extdata", "demo.siber.data.csv", package = "SIBER")
# in your command window. You could load it directly by using the
# returned path, or perhaps better, you could navigate to this folder
# and copy this file to a folder of your own choice, and create a 
# script from this vingette to analyse it. This *.csv file provides
# a template for how your own files should be formatted.

# mydata <- read.csv(fname, header=T)
# siber.example <- createSiberObject(mydata)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)  
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")



par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)


par(mfrow=c(1,1))

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5
)



# Calculate summary statistics for each group: Total Area (TA), 
#Standard Ellipse Area (SEA), and corrected Standard Ellipse Area for uneven sampling (SEAc)

# see Layman et al. 2007 Ecology for interpretation of metrics;
# NR/Y-range: ∆15N Range; Distance between the two species with the most enriched and most depleted ∆15N values (i.e., max-min)
# CR/X-range: ∆13C Range; Distance between the two species with the most enriched and most depleted ∆13C values (i.e. max-min)
# TA: Total area; Convex hull area encompassed by all species in d13C–d15N bi-plot space. This represents a measure of the total amount of niche space occupied, and thus a proxy for the total extent of trophic diversity within a food web.
# CD: Mean distance to centroid; Average Euclidean distance of each species to the d13C–d15N centroid, where the centroid is the mean d13C and d15N value for all species in the food web. This metric provides a measure of the average degree of trophic diversity within a food web.
# NND: Mean nearest neighbor distance; Mean of the Euclidean distances to each species’ nearest neighbor in bi-plot space, and thus a measure of the overall density of species packing. Food webs with a large proportion of species characterized by similar trophic ecologies will exhibit a smaller NND (increased trophic redundancy) than a web in which species are, on average, more divergent in terms of their trophic niche.
# SDNND: Standard deviation of nearest neighbor distance;  A measure of the evenness of species packing in bi-plot space that is less influenced than NND by sample size. Low SDNND values suggest more even distribution of trophic niches.

group.ML <- groupMetricsML(siber.example)
print(group.ML)
summary(group.ML)

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 10000, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 10000, p.interval = 0.95, ci.mean = T,
                  lty = 1, lwd = 2)


# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = T, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5
)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2) 



# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siber.example) 

#print out of layman metrics
print(community.ML)


# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^4 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.


ellipses.posterior <- siberMVN(siber.example, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.

SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c( 0.975) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

SEA.B.median <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.medians = T)

print(SEA.B.median)


#Here, we first calculate the proportion, and hence probability, of the SEA.B for group 1 being smaller than the SEA.B for group 2.

#High urban vs. low urban collembola; 83.55% likelihood of high urban SEA being smaller.
hc.vs.lc <- sum( SEA.B[,1] < SEA.B[,6] ) / nrow(SEA.B)
print(hc.vs.lc)

#High urban vs. low urban oribatida; 86.54% likelihood of high urban SEA being smaller.
ho.vs.lo <- sum( SEA.B[,2] < SEA.B[,4] ) / nrow(SEA.B)
print(ho.vs.lo)

#High urban vs. low urban Mesostigmata; 19.82% likelihood of high urban SEA being smaller.
hm.vs.lm <- sum( SEA.B[,3] < SEA.B[,5] ) / nrow(SEA.B)
print(hm.vs.lm)









# Traditional Stats for isotopes----


#Both isotopes are non-normal, so we should use a Kruskal-wallis test  
shapiro.test(My_SIBER_Data$iso1)
shapiro.test(My_SIBER_Data$iso2)



#make a combined treatment column if necessary
# My_SIBER_Data <- My_SIBER_Data %>% unite(Trt, c("group", "community"))
# My_SIBER_Data$Trt <- as.factor(My_SIBER_Data$Trt)


#Urban by itself: Things are significant!
UrbanKW_C.Community <- with(My_SIBER_Data,kruskal(iso1,community,group=TRUE,console=TRUE))
C_SE.Community <- UrbanKW_C.Community$means$std/sqrt(UrbanKW_C.Community$means$r)

UrbanKW_N.Community <- with(My_SIBER_Data,kruskal(iso2,community,group=TRUE,console=TRUE))
N_SE.Community <- UrbanKW_N.Community$means$std/sqrt(UrbanKW_N.Community$means$r)


#INTERACTION: Things are significant!
UrbanKW_C <- with(My_SIBER_Data,kruskal(iso1,interaction(group,community),group=TRUE,console=TRUE))
C_SE <- UrbanKW_C$means$std/sqrt(UrbanKW_C$means$r)

UrbanKW_N <- with(My_SIBER_Data,kruskal(iso2,interaction(group,community),group=TRUE,console=TRUE))
N_SE <- UrbanKW_N$means$std/sqrt(UrbanKW_N$means$r)  




# Traditional stats plot ---- 
#Plot All ∆C in one plot
TOTALC <-My_SIBER_Data %>%
  ggplot(aes(x=group, y=iso1, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 0.8) +
  
  ggdist::stat_halfeye(
    adjust = .33,## bandwidth
    alpha = 0.5,
    width = .3, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆13C (‰)") + xlab("")

#Plot All ∆N in one plot
TOTALN <-My_SIBER_Data %>%
  ggplot(aes(x=group, y=iso2, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 0.8) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    alpha = 0.5,
    width = .3, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆15N (‰)") + xlab("")


Figure1_Alt <- ggarrange(TOTALC, TOTALN,
                         labels = "A", "B")
#nrow = 1, ncol = 2,
#common.legend = TRUE, 
#legend = "bottom")

ggsave(filename = "Figures/figureone_Alt.jpg", width = 10, height = 5, device='jpeg', dpi=400)



#
#
#
#
#
#
#
#Make individual plots for each taxa.
#
# Plot isotope breadth in collembola
Collembola_Data <-My_SIBER_Data %>% 
  filter(group == "Collembola") %>% 
  select(-group) 

Urban_Color_Palette <- c("#BABABA", "#B8E186")

Col_C <- 
  Collembola_Data %>%
  ggplot(aes(x=community, y=iso1, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + 
  
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆13C (‰)") + xlab("")

Col_N <- 
  Collembola_Data %>%
  ggplot(aes(x=community, y=iso2, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + 
  
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆15N (‰)") + xlab("") + ylim(0,8)



# Plot isotope breadth in Oribatida
Oribatida_Data <-My_SIBER_Data %>% 
  filter(group == "Oribatid") %>% 
  select(-group) 

Urban_Color_Palette <- c("#BABABA", "#B8E186")

Orib_C <- 
  Oribatida_Data %>%
  ggplot(aes(x=community, y=iso1, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + 
  
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆13C (‰)") + xlab("")

Orib_N <- 
  Oribatida_Data %>%
  ggplot(aes(x=community, y=iso2, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + 
  
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆15N (‰)") + xlab("") + ylim(0,8)


# Plot isotope breadth in Mesostigmata
Mesostigmata_Data <-My_SIBER_Data %>% 
  filter(group == "Mesostigmata") %>% 
  select(-group) 

Urban_Color_Palette <- c("#BABABA", "#B8E186")

Meso_C <- 
  Mesostigmata_Data %>%
  ggplot(aes(x=community, y=iso1, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + 
  
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆13C (‰)") + xlab("")

Meso_N <- 
  Mesostigmata_Data %>%
  ggplot(aes(x=community, y=iso2, fill = community, color = community)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA) +
  
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + 
  
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14))+
  scale_color_manual(values = Urban_Color_Palette)+ 
  scale_fill_manual(values = Urban_Color_Palette) + 
  ylab("∆15N (‰)") + xlab("") + ylim(0,12)



#
# Combine Individual plots to show entire picture.
#

Figure1 <- ggarrange(Col_C,
                     Orib_C,
                     Meso_C,
                     Col_N,
                     Orib_N,
                     Meso_N,
                     labels = c("1", "2", "3", "4", "5", "6"),
                     common.legend = TRUE, legend = "bottom")

ggsave(filename = "Figures/figureone.jpg", width = 12, height = 8, device='jpeg', dpi=400)











#community wide ∆C
My_SIBER_Data %>%
  ggplot(aes(x=group, y=iso1, fill = group, color = group)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + facet_wrap(~community)


# Community wide ∆N
My_SIBER_Data %>%
  ggplot(aes(x=group, y=iso2, fill = group, color = group)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3) + facet_wrap(~community)


###
#
# How to report bayesian stats correctly : https://www.nature.com/articles/s41562-021-01177-7/tables/1
#
##
