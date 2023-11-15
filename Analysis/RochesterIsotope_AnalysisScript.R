# Project: Rochester Isotope
#    Name: Hayden Bock
#    Date: 11/15/2023
#    
#
# 
# Packages -------
# install.packages(c("dplyr", "tidyverse", "vegan", "SIBER"))

#install.packages("rjags")
library(rjags)
library(dplyr)
library(tidyverse)
library(vegan)
library(SIBER)
library(readxl)
library(readr)

# Data Import -----
  Env <- read_excel("Environmental_Data.xlsx") 
  Isotope <- read_excel("Isotope_Data.xlsx")
  Joined_DF <- left_join(Env, Isotope)
  My_SIBER_Data <- read_csv("My_SIBER_Data.csv") #iso1 == "∆13C", iso2 = "∆15N", group = "taxa" (1=oribatid, 2=collembola, 3=mesostigmata), community = "urbanization" (1=High, 2=Low),

  
# Traditional Stats for isotopes----

#analyze by urban group only  
  #significant difference in ∆C -> p = 0.0131
  UrbanAOV_C <- aov(iso1~community, data = My_SIBER_Data) #iso1 = ∆C 
  summary(UrbanAOV_C)
  
  #no significant difference in ∆N -> p = 0.821
  UrbanAOV_N <- aov(iso2~community, data = My_SIBER_Data) #iso1 = ∆C 
  summary(UrbanAOV_N)

  
#analyze by urban group and taxa 
  #urban community is significant, but no difference in taxa or interaction.
  InteractionAOV_C <- aov(iso1~community*group, data = My_SIBER_Data) #iso1 = ∆C 
  summary(InteractionAOV_C)
  
  #no difference in urbanization, but significantly different by taxa. no interaction.
  InteractionAOV_N <- aov(iso2~community*group, data = My_SIBER_Data) #iso1 = ∆C 
  summary(InteractionAOV_N)
  
  
  
  
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
  
  
  
  # Calculate summary statistics for each group: TA, SEA and SEAc
  # see Layman et al. 2007 Ecology for interpretation of metrics;
  # NR/Y-range: ∆15N Range; Distance between the two species with the most enriched and most depleted ∆15N values (i.e., max-min)
  # CR/X-range: ∆13C Range; Distance between the two species with the most enriched and most depleted ∆13C values (i.e. max-min)
  # TA: Total area; Convex hull area encompassed by all species in d13C–d15N bi-plot space. This represents a measure of the total amount of niche space occupied, and thus a proxy for the total extent of trophic diversity within a food web.
  # CD: Mean distance to centroid; Average Euclidean distance of each species to the d13C–d15N centroid, where the centroid is the mean d13C and d15N value for all species in the food web. This metric provides a measure of the average degree of trophic diversity within a food web.
  # NND: Mean nearest neighbor distance; Mean of the Euclidean distances to each species’ nearest neighbor in bi-plot space, and thus a measure of the overall density of species packing. Food webs with a large proportion of species characterized by similar trophic ecologies will exhibit a smaller NND (increased trophic redundancy) than a web in which species are, on average, more divergent in terms of their trophic niche.
  # SDNND: Standard deviation of nearest neighbor distance;  A measure of the evenness of species packing in bi-plot space that is less influenced than NND by sample size. Low SDNND values suggest more even distribution of trophic niches.
  
  group.ML <- groupMetricsML(siber.example)
  print(group.ML)
  
  
  # You can add more ellipses by directly calling plot.group.ellipses()
  # Add an additional p.interval % prediction ellilpse
  plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                    lty = 1, lwd = 2)
  
  # or you can add the XX% confidence interval around the bivariate means
  # by specifying ci.mean = T along with whatever p.interval you want.
  plotGroupEllipses(siber.example, n = 100, p.interval = 0.95, ci.mean = T,
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
  cr.p <- c(0.95, 0.99) # vector of quantiles
  
  # call to hdrcde:hdr using lapply()
  SEA.B.credibles <- lapply(
    as.data.frame(SEA.B), 
    function(x,...){tmp<-hdrcde::hdr(x)$hdr},
    prob = cr.p)
  
  # do similar to get the modes, taking care to pick up multimodal posterior
  # distributions if present
  SEA.B.modes <- lapply(
    as.data.frame(SEA.B), 
    function(x,...){tmp<-hdrcde::hdr(x)$mode},
    prob = cr.p, all.modes=T)
  
  
  # extract the posterior means
  mu.post <- extractPosteriorMeans(siber.example, ellipses.posterior)
  
  # calculate the corresponding distribution of layman metrics
  layman.B <- bayesianLayman(mu.post)
  
  
# SIBER ANALYSIS VISUALIZATION -----

  siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                   bty="L", ylim = c(0,6))
  
  # add the ML estimates (if you want). Extract the correct means 
  # from the appropriate array held within the overall array of means.
  comm1.layman.ml <- laymanMetrics(siber.example$ML.mu[[1]][1,1,],
                                   siber.example$ML.mu[[1]][1,2,]
  )
  points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)
  
  
  

  # Visualise the second community

  siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                   bty="L", ylim = c(0,6))
  
  # add the ML estimates. (if you want) Extract the correct means 
  # from the appropriate array held within the overall array of means.
  comm2.layman.ml <- laymanMetrics(siber.example$ML.mu[[2]][1,1,],
                                   siber.example$ML.mu[[2]][1,2,]
  )
  points(1:6, comm2.layman.ml$metrics, col = "red", pch = "x", lwd = 2)
  
  

  # Alternatively, pull out TA from both and aggregate them into a 
  # single matrix using cbind() and plot them together on one graph.

  
  # go back to a 1x1 panel plot
  par(mfrow=c(1,1))
  
  siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
                   xticklabels = c("Community 1", "Community 2"), 
                   bty="L", ylim = c(0,3),
                   las = 1,
                   ylab = "TA - Convex Hull Area",
                   xlab = "")
  
  