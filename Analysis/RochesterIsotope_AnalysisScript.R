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
library(tidyverse)
library(viridis)


# Data Import -----
  Env <- read_excel("Analysis/Environmental_Data.xlsx") 
  Isotope <- read_excel("Analysis/Isotope_Data.xlsx")
  Joined_DF <- left_join(Env, Isotope)
  My_SIBER_Data <- read_csv("Analysis/My_SIBER_Data.csv") #iso1 == "∆13C", iso2 = "∆15N", group = "taxa" (1=oribatid, 2=collembola, 3=mesostigmata), community = "urbanization" (1=High, 2=Low),

# Scale/Standardize data based on Cucherousset & Villeger ----
  Scaled_Data <- My_SIBER_Data %>% 
                    mutate(N_Scaled = (iso2-min(iso2))/(max(iso2)-min(iso2))) %>% 
                    mutate(C_Scaled = (iso1-min(iso1))/(max(iso1)-min(iso1)))
  
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
  
  
# Tradittional stats plot ---- 
  
  # Plot ∆C
  My_SIBER_Data %>%
    ggplot(aes(x=community, y=iso1, fill=group)) +
    geom_boxplot() +
    geom_point(position=position_dodge(width=0.75),aes(group=group))
  
  # Plot ∆N
  My_SIBER_Data %>%
    ggplot(aes(x=community, y=iso2, fill=group)) +
    geom_boxplot() +
    geom_point(position=position_dodge(width=0.75),aes(group=group))
  
  
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
  My_SIBER_Data <- read_csv("Analysis/My_SIBER_Data.csv")
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
  
  
  
  
  #work in progress
  
  
# [WORK IN PROGRESS] Calculate Isotope Diversity Indices based on Cucherousset & Villeger script/Paper ----

  #load libraries
  require(geometry)
  require(ape)
  require(rcdd)
  
  
  #import necessary objects/vectors/etc.
  nm_si<-c("d13C","d15N","dD","d34S")
  
  #Data prep for function
  My_SIBER_Data <- My_SIBER_Data %>% mutate(indiv_ID = paste(community, group, sep = "_"))
  IDiv_Data <- My_SIBER_Data %>% rename(d13C = iso1,
                           d15N = iso2,
                           Species_Code = group) %>% select(-community)
  
  IDiv_Data <- IDiv_Data %>% group_by(indiv_ID) %>% mutate(d13C_sd = sd(d13C),
                                              d15N_sd = sd(d15N))
  
  IDiv_Data_HighUrban <- IDiv_Data %>% filter(stringr::str_detect(indiv_ID, 'High'))
  IDiv_Data_LowUrban <- IDiv_Data %>% filter(stringr::str_detect(indiv_ID, 'Low'))
  
  
  #function
  IDiversity<-function(cons, weight=rep(1,nrow(cons)), nm_plot=NA, col="#477D00", transp=50, scaled=TRUE) {
    
    # names of elements used to describe consumers
    nmel<-colnames(cons)[which(colnames(cons) %in% nm_si)]
    nbel<-length(nmel)
    
    # stable isotope signature
    si<-as.matrix(cons[,nmel])
    
    # checking weighting for all individuals
    if(length(weight) != nrow(cons)) stop(paste(" error: weight does not have the same length than number of consumers"))
    
    # relative weight
    rel_weight<-weight/sum(weight)
    
    # checking number of consumers is higher than number of elements
    if (nrow(cons)<(nbel+1)) stop(paste(" error: computing indices using",nbel,"elements requires at least",nbel+1," consumers"))
    
    
    # vector to store results
    ID<-rep(NA,nbel*4+5) ; names(ID)<-c(paste("min",nmel,sep="_"), paste("max",nmel,sep="_"), paste("range",nmel,sep="_"), paste("IPos",nmel,sep="_"), c("IRic","IDiv","IDis","IEve","IUni") )
    
    ###########################################################
    # computing indices values on each axis
    
    # range of traits values
    ID[paste("min",nmel,sep="_")]<-apply(si,2,min)
    ID[paste("max",nmel,sep="_")]<-apply(si,2,max)
    ID[paste("range",nmel,sep="_")]<-ID[paste("max",nmel,sep="_")]-ID[paste("min",nmel,sep="_")]
    
    # abundance-weighted mean values
    ID[paste("IPos",nmel,sep="_")]<-rel_weight%*%si
    
    ###############################################################################################################################
    
    # generic functions for computing multidimensional diversity indices
    
    I_RED<-function(coord,relab )  {
      
      # number of species
      S<-nrow(coord) 
      
      ###########################################################
      # Richness
      IRic<-round(convhulln(coord,"FA")$vol,6)
      
      # identity of vertices
      vert0<-convhulln(coord,"Fx TO 'vert.txt'")
      vert1<-scan("vert.txt",quiet=T)
      vertices<-(vert1+1)[-1]
      
      ###########################################################
      # Evenness
      
      # inter-species Euclidean distance
      distT<-dist(coord, method="euclidian")
      
      # topology of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
      linkmst<-mst(distT)
      mstvect<-as.dist(linkmst)
      
      # pairwise cumulative relative abundances and conversion into 'dist' class
      ab2<-matrix(0,nrow=S,ncol=S)
      for (q in 1:S)
        for (r in 1:S)
          ab2[q,r]<-relab[q]+relab[r] # end of q,r
      ab2vect<-as.dist(ab2)
      
      # EW index for the (S-1) segments
      EW<-rep(0,S-1)
      flag<-1
      for (m in 1:((S-1)*S/2))
      {if (mstvect[m]!=0) {EW[flag]<-distT[m]/(ab2vect[m]) ; flag<-flag+1}}  # end of m
      
      # PEW index and comparison with 1/S-1
      minPEW<-rep(0,S-1)  ;  OdSmO<-1/(S-1)
      for (l in 1:(S-1))
        minPEW[l]<-min( (EW[l]/sum(EW)) , OdSmO )  # end of l
      
      # IEve
      IEve<-round( ( (sum(minPEW))- OdSmO) / (1-OdSmO ) ,6)
      
      ###############################################################
      # Divergence
      
      # coordinates of vertices
      coordvertices<-coord[vertices,]
      
      # coordinates of the center of gravity of the vertices (B)
      B<-apply(coordvertices,2,mean)
      
      # Euclidean dstance to B (dB)
      dB<-apply(coord, 1, function(x) { (sum((x-B)^2) )^0.5} )
      
      # mean of dB values and deviations to mean 
      meandB<-mean(dB)
      devdB<-dB-meandB
      
      # abundance-weighted mean deviation
      abdev<-relab*devdB
      ababsdev<-relab*abs(devdB)
      
      # IDiv
      IDiv<-round( (sum(abdev)+meandB) / (sum(ababsdev)+meandB) ,6)
      
      ####################################################################
      # results
      indices<-c(IRic,IEve,IDiv) ; names(indices)<-c("IRic","IEve","IDiv")
      detailsRED<-list(vertices=vertices, mst=linkmst, B=B, meandB=meandB)
      I_RED<-list(indices=indices, details=detailsRED )
      invisible(I_RED)
    } # end of function I_RED
    ########################################################################################################################################
    
    # multivariate indices from Villeger et al 2008
    ired<-I_RED(si,rel_weight)
    ID[c("IRic","IEve","IDiv")]<-ired$indices
    
    # Isotopic dispersion: scaled abundance-weighted mean distance to abundance-weighted centroid 
    dist_centr<-apply(si, 1, function(x) { (sum((x-ID[paste("IPos",nmel,sep="_")])^2) )^0.5} ) # distance to abundance-weighted centroid 
    ID["IDis"]<-(rel_weight %*% dist_centr)/ max(dist_centr) # scaling between 0(=all biomass on the centroid) and 1(=all biomass on the most extreme point)
    
    # Isotopic originality : scaled abundance weighted mean distance to nearest neighbour
    # for each organism distance to, and identity of, nearest neighbour
    dist_T<-as.matrix(dist(si,method="euclidean")) ; dist_T[which(dist_T==0)]<-NA
    oriT<-apply(dist_T, 1, min, na.rm=T )
    NN<-dist_T ; NN<-NN-apply(NN,1,min,na.rm=T) ; NN[which(NN!=0)]<-NA   ; NN[which(NN==0)]<-1
    ID["IUni"]<-(oriT %*% rel_weight) / max(oriT)   # abundance weighting and scaling by maximal distance between 2 points
    
    ########################################################################################################################################
    ########################################################################################################################################
    # graphical output
    if( is.na(nm_plot)==FALSE) {
      
      # setting axes limits given consumers signature for all elements
      nmelsd<-paste("sd",nmel,sep="_")
      
      min_axes<- apply(cons[,nmel], 2, min, na.rm=T) ; max_axes<- apply(cons[,nmel], 2, max, na.rm=T) # limits of each axis
      
      # limits of each axis given sd
      if(sum(nmelsd %in% colnames(cons) )==nbel) 
      {min_axes<-apply(cons[,nmel]-cons[,nmelsd], 2, min, na.rm=T)
      max_axes<-apply(cons[,nmel]+cons[,nmelsd], 2, max, na.rm=T) }
      
      # same range on each axis for graphics: min=observed minimal value - 5%  maximal range ; max=observed minimal value + maximal range + 5% maximal range
      rge_axes<-max_axes-min_axes # range on each axis
      
      newlim_axes<-matrix(0,length(nmel),2, dimnames=list(nmel,c("min","max") ) )
      newlim_axes[,"min"]<-min_axes-max(rge_axes)*0.05
      newlim_axes[,"max"]<-min_axes+max(rge_axes)*1.05
      
      rge_plot<-max(rge_axes)*1.1
      
      
      # one jpeg file per pair of elements with 6 panels
      
      for (e1 in 1:(nbel-1))
        for (e2 in (e1+1):nbel) 
        {
          # names of elements
          nmel1<-nmel[e1] ; eval(parse(text=paste("tit1<-tit_",nmel1,sep="") ) )
          nmel2<-nmel[e2] ; eval(parse(text=paste("tit2<-tit_",nmel2,sep="") ) )
          nmel12<-c(nmel1,nmel2)
          
          # creating jpeg file
          nmjpeg<-paste(nm_plot,"_",nmel1,"_",nmel2,".jpeg",sep="")
          jpeg(file=nmjpeg, res=150, width=1200, height=1800)
          layout(matrix(c(1:6),3,2,T)) ; layout.show(6)
          
          # limits of axes
          lim_1<-newlim_axes[nmel1,]
          lim_2<-newlim_axes[nmel2,]
          
          # if axes are for scaled isotope values, "Scaled" in axis title and range is at least from 0 to 1
          if (scaled==TRUE) {
            tit1<-eval(parse(text=paste("tit1<-scl_tit_",nmel1,sep="") ) ); tit2<-eval(parse(text=paste("tit2<-scl_tit_",nmel2,sep="") ) )
            lim_1<-c( min(-0.05, lim_1[1]) , max(1.05,lim_1[2])   )
            lim_2<-c( min(-0.05, lim_2[1]) , max(1.05,lim_2[2])   )
          }  # end of if scaled axes
          
          # labels  on axes
          lab_1<-pretty(lim_1,n=5,min.n=4) ; lab_1<-lab_1[which(lab_1>=lim_1[1] & lab_1<=lim_1[2])]
          lab_2<-pretty(lim_2,n=5,min.n=4) ; lab_2<-lab_2[which(lab_2>=lim_2[1] & lab_2<=lim_2[2])]
          
          ###############################################################
          # Isotopic Position
          
          # Isotopic space given axes limits set using consumers signature
          isotopic_space(nmX=tit1,nmY=tit2, limX=lim_1, limY=lim_2, labX=lab_1,labY=lab_2 )
          
          # mean value
          segments(ID[paste("IPos_",nmel1,sep="")],ID[paste("IPos_",nmel2,sep="")], ID[paste("IPos_",nmel1,sep="")], min(lim_2), lwd=1.5, col=col, lty=2)
          segments(ID[paste("IPos_",nmel1,sep="")],ID[paste("IPos_",nmel2,sep="")], min(lim_1) ,ID[paste("IPos_",nmel2,sep="")], lwd=1.5, col=col, lty=2)
          points( ID[paste("IPos_",nmel1,sep="")],ID[paste("IPos_",nmel2,sep="")], pch=22, bg="white", col=col,cex=2.5)
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab<-sqrt(rel_weight)*0.075*rge_plot  
          symbols(si[,nmel1],si[,nmel2], circles=sizeab, inches=FALSE, bg=col, fg=col, add=TRUE)
          
          # legend for abundance 
          rect(max(lim_1)-0.25*rge_plot, min(lim_2), max(lim_1), min(lim_2)+0.12*rge_plot)
          symbols(max(lim_1)-0.19*rge_plot, min(lim_2)+0.06*rge_plot, circles=sqrt(0.1)*0.075*rge_plot, inches=FALSE, bg="black", fg="black", add=TRUE, lwd=1.5)
          text(max(lim_1)-0.15*rge_plot, min(lim_2)+0.06*rge_plot,"10%", adj=c(0,0.5) ) 
          
          # error bars if any
          if(sum(nmelsd %in% colnames(cons) )==nbel) { meansexy(meanxy=si[,nmel12], sexy=cons[,paste("sd",nmel12,sep="_")],colb=col,lg=0.01*rge_plot ) }# sd
          
          # index
          mtext(side=3, tit1, at=min(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=1)
          mtext(side=3, paste(": IPos=",round(ID[paste("IPos_",nmel1,sep="")],3),sep=""), at=min(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=0)     
          
          mtext(side=3, tit2, at=mean(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=1)
          mtext(side=3, paste(": IPos=",round(ID[paste("IPos_",nmel2,sep="")],3),sep=""), at=mean(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=0)     
          
          mtext(side=3, "Isotopic Position", at=mean(lim_1), line=1.1, cex=0.8,adj=0.5, font=2)                
          
          ###############################################################
          # Isotopic Richness
          
          # Isotopic space given axes limits set using consumers signature
          isotopic_space(nmX=tit1,nmY=tit2, limX=lim_1, limY=lim_2, labX=lab_1,labY=lab_2 )
          
          # range on each axis
          dec1<-rge_plot*0.02
          segments( ID[paste("min_",nmel1,sep="")], min(lim_2)-dec1, ID[paste("min_",nmel1,sep="")], min(lim_2)+dec1, col=col , lwd=3) # min x
          segments( ID[paste("max_",nmel1,sep="")], min(lim_2)-dec1, ID[paste("max_",nmel1,sep="")], min(lim_2)+dec1, col=col , lwd=3) # max x
          
          segments( min(lim_1)-dec1, ID[paste("min_",nmel2,sep="")], min(lim_1)+dec1, ID[paste("min_",nmel2,sep="")],  col=col , lwd=3) # min y
          segments( min(lim_1)-dec1, ID[paste("max_",nmel2,sep="")], min(lim_1)+dec1, ID[paste("max_",nmel2,sep="")],  col=col , lwd=3) # max y
          
          # projected convex hull in 2D
          vert0<-convhulln(si[,nmel12],"Fx TO 'vert.txt'")
          vert1<-scan("vert.txt",quiet=T) ; vertices2D<-(vert1+1)[-1]
          polygon(si[vertices2D,nmel12],border=NA,col=paste(col,transp,sep=""))
          
          # all points (empty) then filling points being vertices in nD
          points(si[,nmel12], pch=21,bg=NA, col=col[1],cex=2)
          points(si[ired$details$vertices,nmel12], pch=21,bg=col, col=col,cex=2)
          
          # error bars if any
          if(sum(nmelsd %in% colnames(cons) )==nbel) { meansexy(meanxy=si[,nmel12], sexy=cons[,paste("sd",nmel12,sep="_")],colb=col,lg=0.01*rge_plot ) }# sd
          
          # index    
          mtext(side=3, tit1, at=min(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=1)
          mtext(side=3, paste(": ",round(ID[paste("range_",nmel1,sep="")],1), " [",round(ID[paste("min_",nmel1,sep="")],1),";",round(ID[paste("max_",nmel1,sep="")],1),"]",sep=""), at=min(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=0) 
          
          mtext(side=3, tit2, at=mean(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=1)
          mtext(side=3, paste(": ",round(ID[paste("range_",nmel2,sep="")],1), " [",round(ID[paste("min_",nmel2,sep="")],1),";",round(ID[paste("max_",nmel2,sep="")],1),"]",sep=""), at=mean(lim_1)+rge_plot*0.1, line=-0.4, cex=0.7,adj=0) 
          
          mtext(side=3, paste("Isotopic Richness=",round(ID['IRic'],3),sep=""), at=mean(lim_1), line=1.1, cex=0.8,adj=0.5, font=2)                
          
          ###############################################################
          # Isotopic Divergence
          
          # Isotopic space given axes limits set using consumers signature
          isotopic_space(nmX=tit1,nmY=tit2, limX=lim_1, limY=lim_2, labX=lab_1,labY=lab_2 )
          
          # projected convex hull in 2D
          vert0<-convhulln(si[,nmel12],"Fx TO 'vert.txt'")
          vert1<-scan("vert.txt",quiet=T) ; vertices2D<-(vert1+1)[-1]
          polygon(si[vertices2D,nmel12],border=col,col=NA, lwd=1 )
          
          segments(ired$details$B[nmel1], ired$details$B[nmel2], si[,nmel1],si[,nmel2],col=col, lty=2, lwd=2)
          points(ired$details$B[nmel1], ired$details$B[nmel2], pch=23,col=col,bg="white",cex=2.5)
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab<-sqrt(rel_weight)*0.075*rge_plot  
          symbols(si[,nmel1],si[,nmel2], circles=sizeab, inches=FALSE, bg=col, fg=col, add=TRUE)
          
          # legend for abundance 
          rect(max(lim_1)-0.25*rge_plot, min(lim_2), max(lim_1), min(lim_2)+0.12*rge_plot)
          symbols(max(lim_1)-0.19*rge_plot, min(lim_2)+0.06*rge_plot, circles=sqrt(0.1)*0.075*rge_plot, inches=FALSE, bg="black", fg="black", add=TRUE, lwd=1.5)
          text(max(lim_1)-0.15*rge_plot, min(lim_2)+0.06*rge_plot,"10%", adj=c(0,0.5) ) 
          
          # error bars if any
          if(sum(nmelsd %in% colnames(cons) )==nbel) { meansexy(meanxy=si[,nmel12], sexy=cons[,paste("sd",nmel12,sep="_")],colb=col,lg=0.01*rge_plot ) }# sd
          
          # index    
          mtext(side=3, paste("Isotopic Divergence=",round(ID['IDiv'],3),sep=""), at=mean(lim_1), line=0.5, cex=0.8,adj=0.5, font=2)               
          
          ###############################################################
          # Isotopic Dispersion
          
          # Isotopic space given axes limits set using consumers signature
          isotopic_space(nmX=tit1,nmY=tit2, limX=lim_1, limY=lim_2, labX=lab_1,labY=lab_2 )
          
          # distance to abundance weighted centroid
          segments(ID[paste("IPos_",nmel1,sep="")],ID[paste("IPos_",nmel2,sep="")], si[,nmel1],si[,nmel2],col=col, lty=3, lwd=2)
          points( ID[paste("IPos_",nmel1,sep="")],ID[paste("IPos_",nmel2,sep="")], pch=22, bg="white", col=col,cex=2.5)
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab<-sqrt(rel_weight)*0.075*rge_plot  
          symbols(si[,nmel1],si[,nmel2], circles=sizeab, inches=FALSE, bg=col, fg=col, add=TRUE)
          
          # legend for abundance 
          rect(max(lim_1)-0.25*rge_plot, min(lim_2), max(lim_1), min(lim_2)+0.12*rge_plot)
          symbols(max(lim_1)-0.19*rge_plot, min(lim_2)+0.06*rge_plot, circles=sqrt(0.1)*0.075*rge_plot, inches=FALSE, bg="black", fg="black", add=TRUE, lwd=1.5)
          text(max(lim_1)-0.15*rge_plot, min(lim_2)+0.06*rge_plot,"10%", adj=c(0,0.5) ) 
          
          # error bars if any
          if(sum(nmelsd %in% colnames(cons) )==nbel) { meansexy(meanxy=si[,nmel12], sexy=cons[,paste("sd",nmel12,sep="_")],colb=col,lg=0.01*rge_plot ) }# sd
          
          # index    
          mtext(side=3, paste("Isotopic Dispersion=",round(ID['IDis'],3),sep=""), at=mean(lim_1), line=0.5, cex=0.8,adj=0.5, font=2)               
          
          ###############################################################
          # Isotopic Evenness
          
          # Isotopic space given axes limits set using consumers signature
          isotopic_space(nmX=tit1,nmY=tit2, limX=lim_1, limY=lim_2, labX=lab_1,labY=lab_2 )
          
          # MST
          for (j in 1:nrow(ired$details$mst))
            for (i in 1:nrow(ired$details$mst))
              if (ired$details$mst[j,i]==1 & j>i) segments(si[,nmel1][j], si[,nmel2][j], si[,nmel1][i], si[,nmel2][i], col=col, lwd=1.5)
          
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab<-sqrt(rel_weight)*0.075*rge_plot  
          symbols(si[,nmel1],si[,nmel2], circles=sizeab, inches=FALSE, bg=col, fg=col, add=TRUE)
          
          # legend for abundance 
          rect(max(lim_1)-0.25*rge_plot, min(lim_2), max(lim_1), min(lim_2)+0.12*rge_plot)
          symbols(max(lim_1)-0.19*rge_plot, min(lim_2)+0.06*rge_plot, circles=sqrt(0.1)*0.075*rge_plot, inches=FALSE, bg="black", fg="black", add=TRUE, lwd=1.5)
          text(max(lim_1)-0.15*rge_plot, min(lim_2)+0.06*rge_plot,"10%", adj=c(0,0.5) ) 
          
          # error bars if any
          if(sum(nmelsd %in% colnames(cons) )==nbel) { meansexy(meanxy=si[,nmel12], sexy=cons[,paste("sd",nmel12,sep="_")],colb=col,lg=0.01*rge_plot ) }# sd
          
          # index    
          mtext(side=3, paste("Isotopic Evenness=",round(ID['IEve'],3),sep=""), at=mean(lim_1), line=0.5, cex=0.8,adj=0.5, font=2)     
          
          ###############################################################
          # Isotopic Uniqueness
          
          # isotopic space given axes limits set using consumers signature
          isotopic_space(nmX=tit1,nmY=tit2, limX=lim_1, limY=lim_2, labX=lab_1,labY=lab_2 )
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab<-sqrt(rel_weight)*0.075*rge_plot  
          symbols(si[,nmel1],si[,nmel2], circles=sizeab, inches=FALSE, bg=col, fg=col, add=TRUE)
          
          # distance to nearest neighbour
          for (k in 1:nrow(NN))
          {
            arrows( si[k,nmel1],si[k,nmel2], si[which(NN[k,]==1)[1],nmel1], si[which(NN[k,]==1)[1],nmel2], col="black", lwd=1.8, length=0.1, angle=20)
          } # end of k
          
          # legend for abundance 
          rect(max(lim_1)-0.25*rge_plot, min(lim_2), max(lim_1), min(lim_2)+0.12*rge_plot)
          symbols(max(lim_1)-0.19*rge_plot, min(lim_2)+0.06*rge_plot, circles=sqrt(0.1)*0.075*rge_plot, inches=FALSE, bg="black", fg="black", add=TRUE, lwd=1.5)
          text(max(lim_1)-0.15*rge_plot, min(lim_2)+0.06*rge_plot,"10%", adj=c(0,0.5) ) 
          
          # error bars if any
          if(sum(nmelsd %in% colnames(cons) )==nbel) { meansexy(meanxy=si[,nmel12], sexy=cons[,paste("sd",nmel12,sep="_")],colb=col,lg=0.01*rge_plot ) }# sd
          
          # index    
          mtext(side=3, paste("Isotopic Uniqueness=",round(ID['IUni'],3),sep=""), at=mean(lim_1), line=0.5, cex=0.8,adj=0.5, font=2)               
          
          ###############################################
          graphics.off()	
        } # end of e1, e2  
      
      
    } # end of plot
    
    ##############################################################################################
    # returning results	
    return(ID)	
    
  } # end of function IDiversity
  
  
  
  #Output, separated by urban class.
  High_IDiv_Metric <- IDiversity(cons=IDiv_Data_HighUrban)
  Low_IDiv_Metric <- IDiversity(cons=IDiv_Data_LowUrban)
  