#----------------------------------------------------------------
# Kéry et al. 2022: Integrated distance sampling models for simple point counts 
# Case study
# Analyse Distance Sampling (DS) data from Oregon2020 together with 
# EBird Point Count (PC) data and detection/nondetection (DND) data 
# Marc Kéry, Nicolas Strebel, Tyler Hallman, Kenneth F. Kellner
# March-November 2022 
#---------------------------------------------------------------

# Load data ----
# Distance Sampling (DS) data from Oregon2020
head(DS <- read.table("data/amro_obsdetection.csv",header=T,sep=";",as.is = T))
DS$Distance <- DS$Distance/1000 # Convert to km

# Table including covariates and Oregon2020 counts (particularly the zeros) 
head(sitecovs <- read.table("data/amro_sitecovs_simplified.csv",header=T,sep=";",as.is = T))

# eBird
head(PC <- read.table("data/amro_ebird_simplified.csv",header=T,sep=";",as.is = T))

# No DND data available --> Convert half of the ebird PC data to DND and keep the other half as PC  
set.seed(123)
PC$index_PC_DND <- sample(x = c(1,2),prob = c(0.5,0.5),size = nrow(PC),replace = T)
plot(x = PC$longitude.x,y = PC$latitude.x,asp=1,col=PC$index_PC_DND)
nrow(DND <- PC[PC$index_PC_DND==2,])
nrow(PC <- PC[PC$index_PC_DND==1,])

# Prepare data ----
# Prepare eBird data 
counts <- PC$count

# Link DS and zero observations (from sitecovs)
names(sitecovs)[names(sitecovs) == "Latitude.x"] <- "Latitude" # Homogenize names
names(sitecovs)[names(sitecovs) == "Longitude.x"] <- "Longitude"
nondetections <- sitecovs[sitecovs$Count==0,] # Select zero counts
nondetections <- nondetections[,names(nondetections) %in% names(DS)] # Columns from nondetection data that are not in DS can be removed
nondetections[,names(DS)[!names(DS) %in% names(nondetections)]] <- NA # Add columns to nondetecion data that are in DS only 
DS <- rbind(DS,nondetections) # combine

# Site as factor level
DS$site <- as.numeric(as.factor(DS$Location_ID))

# Define new truncation distance (B) for the DS data set 
newB <- 0.3 # in km

# Replace observations with greater distances than truncation distance with NA
DS$Distance[DS$Distance>newB] <- NA
sum(is.na(DS$Distance))

# Number of individuals detected (or 'captured') per site
ncap <- tapply(DS$Distance,DS$site,function(x) sum(!is.na(x)))
ncap <- as.vector(ncap) # Number of individuals detected per site

# Site ID of each non-missing observation
siteDS <- DS$site[!is.na(DS$Distance)]

# Either the truncation distance for eBird, 
# or large enough distance to be unlimited for species of interest (like maximum detection distance) 
fullDistance <- 0.5 # in km

# Distance bin width for rectangle approximation of integral
# approximation to the detection probability at each point
delta <- 0.05 # in km

# Make mid-points and chop up distance data
# Both for newB (in data set 1) and the fullDistance (PC and DND data)
midpt <- seq(delta/2, newB, delta)
midptFull <- seq(delta/2, fullDistance, delta)

# Convert distance data to distance category, or distance class data
dclass <- DS$Distance %/% delta + 1

# Number of distance intervals
nD <- length(midpt)            # In DS data
nDfull <- length(midptFull)    # In PC and DND data

# Observed categorical observations
# (Make restriction to non-missing distances for newB = 0.3)
dclass <- dclass[!is.na(dclass)]   # toss out missing dclasses

# Total number of individuals detected (after applying new B = 0.3)
# This is in data set 1, i.e., in the proper DS survey
nind <- length(dclass)

# Compute the areas associated with the count data for each site
(A_DS <- pi * newB^2 )
(A_PC_DND <- pi * fullDistance^2 )

# Number of sites per data set
nsites_DS <- length(unique(DS$site))      # Data set 1
nsites_PC <- length(unique(apply(PC[,c("latitude.x","longitude.x")],1,paste,collapse="_")))      # PC data
nsites_DND <- length(unique(apply(DND[,c("latitude.x","longitude.x")],1,paste,collapse="_")))      # DND data

# Attribute site number (according to DS) to each row in sitecovs
sitecovs$sitenr <- DS$site[match(sitecovs$Location_ID,DS$Location_ID)]

# Prepare objects for analysis in BUGS ----
# values to scale covariate data
mean.elev <- mean(c(sitecovs$Elevation_mean315, PC$Elevation_mean315, DND$Elevation_mean315))
sd.elev <- sd(c(sitecovs$Elevation_mean315, PC$Elevation_mean315, DND$Elevation_mean315))
mean.minutes.since.dawn <- mean(c(sitecovs$MinutesSinceDawn, PC$minutes_since_dawn, DND$minutes_since_dawn))
sd.minutes.since.dawn <- sd(c(sitecovs$MinutesSinceDawn, PC$minutes_since_dawn, DND$minutes_since_dawn))
mean.jdate <- mean(c(sitecovs$JulianDate, PC$julian_date, DND$julian_date))
sd.jdate <- sd(c(sitecovs$JulianDate, PC$julian_date, DND$julian_date))
  
# Data set to be fed into JAGS
str(bdata <- list(
  
  ### DS data set:
  nsites_DS = nsites_DS, # Number of sites  
  nind = nind, # Number of individuals detected 
  newB = newB, # Truncation distance
  nD = nD, # Number of distance bins
  midpt = midpt, # Midpoints of distance bins
  ncap = ncap, # Number of detected individuals per site
  dclass = dclass, # Distance category for each observation
  siteDS = siteDS, # Site number of each observation  
  A_DS = A_DS, # Area associated with the DS data
  delta = delta, # Bin width
  DSduration = rep(5, nsites_DS), # Survey duration 
  # Covariates on abundance
  year_DS = sitecovs$Year-10, # year number
  habitat_DS = ((sitecovs$CanopyCover_315[order(sitecovs$sitenr)])-50)/100, # Habitat covariate; scale so that 0 = -0.5, 50 = 0, 100 = 0.5 ### Changed, NS 23.12.20
  elev_DS = (sitecovs$Elevation_mean315[order(sitecovs$sitenr)]-mean.elev)/sd.elev,
  # Covariates on detectability
  cancovdetect_DS = ((sitecovs$CanopyCover_165[order(sitecovs$sitenr)])-50)/100, # Canopy Cover 
  urbandetect_DS = ((sitecovs$DevelopedMediumIntensity_165[order(sitecovs$sitenr)] + sitecovs$DevelopedHighIntensity_165[order(sitecovs$sitenr)])-50)/100,
  # Covariates on availability
  day_DS = (sitecovs$JulianDate-mean.jdate)/sd.jdate, # 
  time_DS = (sitecovs$MinutesSinceDawn-mean.minutes.since.dawn)/sd.minutes.since.dawn,

  ### PC data set:
  nsites_PC = nsites_PC, # Number of sites
  counts = counts, # Counts
  ebirdDuration = PC$duration,
  # Covariates on abundance
  year_PC = PC$year-10, # year number
  habitat_PC = ((PC$CanopyCover_315)-50)/100, # Habitat covariate
  elev_PC = (PC$Elevation_mean315-mean.elev)/sd.elev,
  # Covariates on detectability
  cancovdetect_PC = ((PC$CanopyCover_165)-50)/100, # Detectability  covariate
  urbandetect_PC = ((PC$DevelopedMediumIntensity_165 + PC$DevelopedHighIntensity_165)-50)/100,
  # Covariates on availability
  day_PC = (PC$julian_date-mean.jdate)/sd.jdate,
  time_PC = (PC$minutes_since_dawn-mean.minutes.since.dawn)/sd.minutes.since.dawn,

  ### DND data set:
  nsites_DND = nsites_DND, # Number of sites
  y = 1*(DND$count>0), # DND
  DNDduration = DND$duration,  
  # Covariates on abundance
  year_DND = DND$year-10, # year number
  habitat_DND = ((DND$CanopyCover_315)-50)/100, # Habitat covariate
  elev_DND = (DND$Elevation_mean315-mean.elev)/sd.elev,
  # Covariates on detectability
  cancovdetect_DND = ((DND$CanopyCover_165)-50)/100, # Detectability  covariate
  urbandetect_DND = ((DND$DevelopedMediumIntensity_165 + DND$DevelopedHighIntensity_165)-50)/100,
  # Covariates on availability
  day_DND = (DND$julian_date-mean.jdate)/sd.jdate,
  time_DND = (DND$minutes_since_dawn-mean.minutes.since.dawn)/sd.minutes.since.dawn,

  ### Common data for PC and DND
  fullDistance = fullDistance, # Assumed radius for PC
  nDfull = nDfull, # length of (latent) distance bins in PC data
  midptFull = midptFull, # Midpoints of (latent) distance bins
  A_PC_DND = A_PC_DND, # Area associated with the PC and DND data
  nyear = length(unique(c(DND$year,PC$year,DS$Year))) # number of years in data

))
 
 

# In this model we add random noise in the detection function of the PC/DND data
# Model
cat(file="IDS3_CaseStudy_5.txt","
    model{
      
      # Priors for all parameters in IDS model
      # -------------------------------------------------------------------
      # Partly different parameters in the detectabiliperceptibility component 
	  # of the detection model for the DS and the PC/DND portions of the data
      for(d in 1:2){           # d indexes the two data types: DS vs. PC/DND
        alpha0[d] <- log(mean.sigma[d])  # sigma intercept on log scale and ...
        mean.sigma[d] ~ dunif(0, 0.1) # ... on natural scale
        tau.eps[d] <- pow(sd.eps[d], -2)
		sd.eps[d] ~ dunif(0, 1)
      }
      alpha1 ~ dnorm(0, 1)     # Canopy cover
      alpha2 ~ dnorm(0, 1)     # Urban area

      # Shared parameters in the availability component of the detection model
      gamma0 <- log(mean.phi)  # phi intercept on log scale and ...
      mean.phi ~ dunif(0, 1)   # ... on natural scale
      gamma1 ~ dnorm(0, 0.5)   # Effect of day of year on singing rate
      gamma2 ~ dnorm(0, 0.5)   # Effect of day of year on singing rate (quadratic)
      gamma3 ~ dnorm(0, 1)     # Effect of time of day on singing rate
      gamma4 ~ dnorm(0, 1)     # Effect of time of day on singing rate (quadratic)

      # Shared parameters in the abundance model
      # Random intercept:
      for (i in 1:nyear) {
        ann.beta0[i] ~ dnorm(beta0, tau.beta0)
      }
      beta0 <- log(mean.lambda) # lambda intercept on log scale and ...
      mean.lambda ~ dunif(0, 80) # ... on natural scale
      tau.beta0 <- pow(sd.beta0,-2) 
      sd.beta0 ~ dunif(0, 2)
      # Covariates:
      beta1 ~ dnorm(0, 1)  # Effect of habitat (canopy cover) on abundance
      beta2 ~ dnorm(0, 0.1)  # Effect of habitat (canopy cover) on abundance (quadratic)
      beta3 ~ dnorm(0, 1)  # Effect of elevation on abundance
      beta4 ~ dnorm(0, 1)  # Effect of elevation on abundance (quadratic)
      
      
      # Submodel for the DS data
      # -------------------------------------------------------------------
      # Hierarchical construction of the likelihood
      # Model for binned distance observations of every detected individual
      for(i in 1:nind){       # Loop over all detected individuals
        dclass[i] ~ dcat(fc[siteDS[i],])               # Part 1 of HM
      }
      
      # Construction of the cell probabilities for the nD distance bands
      # This is for the truncation distance for the DS data (here, newB = 0.3 km)

      for(s in 1:nsites_DS){    # Loop over all sites in data set 1
        for(g in 1:nD){       # midpt = mid-point of each distance band
          log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s]^2)
          pi[s,g] <- ((2 * midpt[g] ) / newB^2) * delta # prob. per interval
          f[s,g] <- p[s,g] * pi[s,g]
          fc[s,g] <- f[s,g] / pcap[s]
        }
        # Rectangular approx. of integral that yields the Pr(capture)
        pcap[s] <- sum(f[s,])
        
        ### Log-linear models on abundance, detectability, and availability
        # Abundance (lambda)
        log(lambda1[s]) <- ann.beta0[year_DS[s]] + beta1 * habitat_DS[s] + beta2 * pow(habitat_DS[s],2) + beta3 * elev_DS[s] + beta4 * pow(elev_DS[s],2) # Log-linear model for abundance 
        
        # Detectability/perceptability (sigma)
        log(sigma[s]) <- alpha0[1] + alpha1 * cancovdetect_DS[s] + alpha2 * urbandetect_DS[s] + eps1[s]  # Log-Linear model for detection probability
		    eps1[s] ~ dnorm(0, tau.eps[1])

        # Availability (phi)
        log(phi[s]) <- gamma0 + gamma1 * day_DS[s] + gamma2 * pow(day_DS[s],2) + gamma3 * time_DS[s] + gamma4 * pow(time_DS[s],2)  # Log-linear model for availability
        theta[s] <- 1-exp(-DSduration[s]*phi[s])  # Effect of duration on availability

        # Multiply availability with detection probability
        pDS[s] <- pcap[s] * theta[s]

        ### Binomial mixture part (in conditional multinomial specification)
        ncap[s] ~ dbin(pDS[s], N1[s])  # Part 2 of HM: number captured
        N1[s] ~ dpois(A_DS * lambda1[s])  # Note use of A_DS as an offset

		# Assess model fit: compute Bayesian p-value for Freeman-Tukey discrepancy
        # Compute fit statistic for observed data (DS data portion)
        evalDS[s] <- pDS[s] * N1[s]
        EDS[s] <- pow((sqrt(ncap[s]) - sqrt(evalDS[s])), 2)

        # Generate replicate DS count data and compute same fit stats for them
        ncap.new[s] ~ dbin(pDS[s], N1[s])
        EDS.new[s] <- pow((sqrt(ncap.new[s]) - sqrt(evalDS[s])), 2)
	   }
      # Add up fit stats across sites for DS data
      fitDS <- sum(EDS[])
      fitDS.new <- sum(EDS.new[])
      
      # Submodel for the PC data
      # -----------------------------------------
      # Parameters on abundance, availability and detectability 
	  # are al shared among the submodels for the different data sets
      # For spatial reconciliation of the data, we use A as an offset
      
      # Likelihood for the PC data
      # Parameters on abundance are shared among all three submodels (DS, PC, DND),
      #	parameters on sigma are shared among PC and DND
      
      for(s in 1:nsites_PC){
        
        # Now we compute the average p over all distance bands for the unlimited-distance surveys
        
        for(g in 1:nDfull){       # midpt = mid-point of each distance band
          log(p2[s,g]) <- -midptFull[g] * midptFull[g] / (2 * sigmaPC[s]^2)
          pi2[s,g] <- ((2 * midptFull[g] ) / fullDistance^2) * delta # prob. per interval
          f2[s,g] <- p2[s,g] * pi2[s,g]
        }
        # Rectangular approx. of integral that yields the Pr(capture)
        pcap2[s] <- sum(f2[s,])
        
        ### Log-linear models on abundance, detectability and availability
        # Abundance
        log(lambda2[s]) <- ann.beta0[year_PC[s]] + beta1 * habitat_PC[s] + beta2 * pow(habitat_PC[s],2) + beta3 * elev_PC[s] + beta4 * pow(elev_PC[s],2) # Log-linear model on abundance 
    
        # Detectability
        log(sigmaPC[s]) <- alpha0[2] + alpha1 * cancovdetect_PC[s] + alpha2 * urbandetect_PC[s] +eps2[s]# Log-Linear model for detection probability
        eps2[s] ~ dnorm(0, tau.eps[2])

        # Availability
        log(phi2[s]) <- gamma0 + gamma1 * day_PC[s] + gamma2 * pow(day_PC[s],2) + gamma3 * time_PC[s] + gamma4 * pow(time_PC[s],2) 
        theta2[s] <- 1-exp(-ebirdDuration[s]*phi2[s]) # Effect of duration on availability

        # Multiply availability with detection probability
        pPC[s] <- pcap2[s] * theta2[s] 

        ### Binomial mixture part 
        N2[s] ~ dpois(A_PC_DND * lambda2[s])
        counts[s] ~ dbinom(pPC[s], N2[s]) 
		
		# Assess model fit: compute Bayesian p-value for Freeman-Tukey discrepancy
        # Compute fit statistic for observed data (PC portion)
        evalPC[s] <- pPC[s] * N2[s]
        EPC[s] <- pow((sqrt(counts[s]) - sqrt(evalPC[s])), 2)

        # Generate replicate point count data and compute same fit stats for them
        counts.new[s] ~ dbin(pPC[s], N2[s])
        EPC.new[s] <- pow((sqrt(counts.new[s]) - sqrt(evalPC[s])), 2)
	   }
      # Add up fit stats across sites for PC data
      fitPC <- sum(EPC[])
      fitPC.new <- sum(EPC.new[])

      # Submodel for the DND data
      # -----------------------------------------
      for(s in 1:nsites_DND){
        # Again, we compute the average p over all distance bands for the unlimited-distance surveys

        for(g in 1:nDfull){       # midpt = mid-point of each distance band
          log(p3[s,g]) <- -midptFull[g] * midptFull[g] / (2 * sigmaDND[s]^2)
          pi3[s,g] <- ((2 * midptFull[g] ) / fullDistance^2) * delta # prob. per interval
          f3[s,g] <- p3[s,g] * pi3[s,g]
        }
        # Rectangular approx. of integral that yields the Pr(capture)
        pcap3[s] <- sum(f3[s,])
        
        ### Log-linear models on abundance, detectability and availability
        # Abundance
        log(lambda3[s]) <- ann.beta0[year_DND[s]] + beta1 * habitat_DND[s] + beta2 * pow(habitat_DND[s],2) + beta3 * elev_DND[s] + beta4 * pow(elev_DND[s],2) # Log-linear model on abundance 
        
        # Detectability
        log(sigmaDND[s]) <- alpha0[2] + alpha1 * cancovdetect_DND[s] + alpha2 * urbandetect_DND[s] + eps3[s]# Log-Linear model for detection probability
        eps3[s] ~ dnorm(0, tau.eps[2])
        
        # Availability
        log(phi3[s]) <- gamma0 + gamma1 * day_DND[s] + gamma2 * pow(day_DND[s],2) + gamma3 * time_DND[s] + gamma4 * pow(time_DND[s],2) 
        theta3[s] <- 1-exp(-DNDduration[s]*phi3[s]) # Effect of duration on availability

        # Temporal scaling based on duration
        pDND[s] <- pcap3[s] * theta3[s] 

        ### Royle-Nichols part
        N3[s] ~ dpois(A_PC_DND * lambda3[s])
        y[s] ~ dbern(1 - (1 - pDND[s])^N3[s])

		# Assess model fit: DND are binary data, so have to aggregate
		# Compute the observed and expected number of sites with detections
		# and treat that as a fit statistics for the DND portion of the data
        y.new[s] ~ dbern(1 - (1 - pDND[s])^N3[s])    # Create replicate data sets
      }
      # Add up fit stats across sites for DND data
      fitDND <- sum(y[])
      fitDND.new <- sum(y.new[])

      # Compute Bayesian p-value for all three portions of the data
      bpvDS <- step(fitDS.new - fitDS)
      bpvPC <- step(fitPC.new - fitPC)
      bpvDND <- step(fitDND.new - fitDND)
    }
    ")

# Inits
Nst1 <- ncap + 1
Nst2 <- counts + 1
Nst3 <- 1*(DND$count>0) + 1
inits <- function(){list(mean.sigma = runif(2, 0, 0.1), alpha1 = rnorm(1, 0, 0.1), alpha2 = rnorm(1, 0, 0.1), 
                         mean.phi = runif(1, 0, 1), gamma1 = rnorm(1, 0, 0.1), gamma2 = rnorm(1, 0, 0.1), 
                         gamma3 = rnorm(1, 0, 0.1), gamma4 = rnorm(1, 0, 0.1),
                         mean.lambda = runif(1, 1, 60), beta1 = rnorm(1, 0, 0.1), beta2 = rnorm(1, 0, 0.1), 
                         beta3 = rnorm(1, 0, 0.1), beta4 = rnorm(1, 0, 0.1), 
						 sd.eps = runif(2, 1, 2), N1 = Nst1, N2 = Nst2, N3 = Nst3)}  

# Params to save
params <- c("mean.sigma", "alpha0", "alpha1", "alpha2", "sd.eps",
            "mean.phi", "gamma0", "gamma1", "gamma2", "gamma3", "gamma4", 
            "mean.lambda", "beta0", "beta1", "beta2", "beta3", "beta4", "sd.beta0", "ann.beta0",
			"fitDS", "fitDS.new", "fitPC", "fitPC.new", "fitDND", "fitDND.new", 
			"bpvDS", "bpvPC", "bpvDND")

# MCMC settings
na <- 10  ;  nc <- 3  ;  ni <- 12  ;  nb <- 2  ;  nt <- 2 # test, 30 sec

na <- 3000;  nc <- 3;  ni <- 15000;  nb <- 10000;  nt <- 5 # 393 min
na <- 5000;  nc <- 3;  ni <- 60000;  nb <- 20000;  nt <- 40 # 660 min
na <- 10000;  nc <- 4;  ni <- 120000;  nb <- 20000;  nt <- 100 # 1500 min

na <- 10000;  nc <- 4;  ni <- 120000;  nb <- 60000;  nt <- 60


# Launch JAGS, check convergence and summarize posteriors ----
# Launch
library(jagsUI)
start <- Sys.time()
set.seed(123)
out5XX <- jags(bdata, inits, params, "IDS3_CaseStudy_5.txt", n.adapt = na,
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
difftime(Sys.time(),start)
traceplot(out5XX)
print(out5XX, 2) 

# Visualize and summarize results  ----
library(sf); library(terra)

# Load environmental data
head(centroids <- read.table("data/env_centroids.csv",header=T,sep=";",as.is = T))

#load in file with desired crs
polygon.grid <- st_read("data/geodata/grid_1km.shp")

# mask for study area
mask <- st_read("data/geodata/StudyArea.shp")

# adjust variables as they were in the model data
centroids$cancov.scaled <- (centroids$cancov - 50)/100
centroids$elev.scaled <- (centroids$elev - mean.elev)/sd.elev

# Create a function to make predictions based on model output
predict.abundances <- function(environmental = environmental, parameters = parameters){
  prediction <- exp(parameters[1] + parameters[2] * environmental$cancov.scaled + parameters[3] * I(environmental$cancov.scaled^2) + 
                      parameters[4] * environmental$elev.scaled + parameters[5] * I(environmental$elev.scaled^2))
  return(prediction)
}

# Extract model output samples as a matrix
samples.matrix <- as.matrix(out5XX$samples)
test.raster <- rast(nlyrs=1, crs = crs(vect(polygon.grid)), extent = ext(polygon.grid), resolution = 1000)

### Plot map
centroids.predictions <- apply(samples.matrix[,16:20], 1, FUN = predict.abundances, environmental = centroids) # Make predictions for each km-square and iteration
centroids.predictions <- cbind(centroids, med=apply(centroids.predictions,1,median))
centroids.predictions.output_sf = st_as_sf(centroids.predictions, coords = c("X", "Y"), crs = st_crs(polygon.grid))
centroids.predictions.output_rast <- rasterize(vect(centroids.predictions.output_sf), test.raster, field = "med")
centroids.predictions.output_rast <- mask(centroids.predictions.output_rast, vect(mask)) 

# Plot density map
plot(centroids.predictions.output_rast,main="Estimates based on jags analysis")

### Estimate abundance throughout study area for each sample
# Select km-squares within study area 
within.study.area <- !(is.na(values(centroids.predictions.output_rast)))

# Make predictions for each km-square and iteration
centroids.predictions <- apply(samples.matrix[,16:20], 1, FUN = predict.abundances, environmental = centroids[within.study.area,])

# Summarize population estimates
pop.estimates <- apply(centroids.predictions,2,sum)

(mean.pop.estimate <- mean(pop.estimates)) # Mean for the entire study area
(uci.pop.estimate <- quantile(pop.estimates, probs = 0.975)) # upper
(lci.pop.estimate <- quantile(pop.estimates, probs = 0.025)) # lower

round(range(apply(centroids.predictions,1,median)),1) # range of median estimates per km-square

# Analysis in unmarked  ----
# Difference of the model fitted in unmarked compared to IDS3_CaseStudy_5.txt:
#  - no annual random effects on abundance, 
#  - no different levels of heterogeneity in the detection functions estimated for different portions of the data 
# replace "IDS3_CaseStudy_5.txt" by "IDS3_CaseStudy_simple.txt" on line 403 to fit a (jags) model equivalent to the one in unmarked
# Install dev version of unmarked that has IDS() function
# this only needs to be run once
# this requires Rtools to be installed if on Windows, and Xcode if on Mac
# First install dependencies. If you have problems with the code below,
# make sure all these R packages are the latest CRAN versions, *especially* TMB
depends <- tools::package_dependencies('unmarked')[[1]]
depends <- depends[! depends %in% c("graphics", "methods", "parallel", "stats", "utils")]
install.packages(depends)

# Install from Github
remotes::install_github("kenkellner/unmarked", ref="IDS")

# Load package
library(unmarked)

# Format JAGS input data for use in unmarked

#### Distance data ####

# Covariates
site_covs_DS <- data.frame(elev=bdata$elev_DS, habitat=bdata$habitat_DS,
                           cancov=bdata$cancovdetect_DS, urban=bdata$urbandetect_DS,
                           day=bdata$day_DS, time=bdata$time_DS)

# y matrix
y_DS <- matrix(0, nrow=bdata$nsites_DS, ncol=max(bdata$dclass))
dclass_fac <- factor(bdata$dclass, levels=1:6)
site_idx <- sort(unique(bdata$siteDS))
for (i in site_idx){
  y_DS[i,] <- table(dclass_fac[bdata$siteDS==i])
}

# Distance breaks
db <- c(bdata$midpt - 0.025, max(bdata$midpt) + 0.025)

# Survey durations
dur_DS <- bdata$DSduration

# unmarked frame
umf_DS <- unmarkedFrameDS(y=y_DS, siteCovs=site_covs_DS, dist.breaks=db,
                          survey="point", unitsIn='km')

#### Point count data ####

# Covariates
site_covs_PC <- data.frame(elev=bdata$elev_PC, habitat=bdata$habitat_PC,
                           cancov=bdata$cancovdetect_PC, urban=bdata$urbandetect_PC,
                           day=bdata$day_PC, time=bdata$time_PC)

# y matrix
y_PC <- matrix(bdata$counts, ncol=1)

# Survey durations
dur_PC <- bdata$ebirdDuration

# unmarked frame
umf_PC <- unmarkedFramePCount(y=y_PC, siteCovs=site_covs_PC)

#### Detection-nondetection data ####

# Covariates
site_covs_DND <- data.frame(elev=bdata$elev_DND, habitat=bdata$habitat_DND,
                            cancov=bdata$cancovdetect_DND, urban=bdata$urbandetect_DND,
                            day=bdata$day_DND, time=bdata$time_DND)

# y matrix
y_DND <- matrix(bdata$y, ncol=1)

# Survey durations
dur_DND <- bdata$DNDduration

# unmarked frame
umf_DND <- unmarkedFrameOccu(y=y_DND, siteCovs=site_covs_DND)

#### Fit model ####

(mod <- IDS(lambdaformula = ~habitat + I(habitat^2) + elev + I(elev^2),
           detformulaDS = ~cancov + urban,
           detformulaPC = NULL,
           detformulaOC = NULL,
           dataDS = umf_DS, dataPC = umf_PC, dataOC = umf_DND,
           availformula = ~day + I(day^2) + time + I(time^2),
           durationDS = dur_DS, durationPC = dur_PC, durationOC = dur_DND,
           maxDistPC = bdata$fullDistance, maxDistOC = bdata$fullDistance,
           K=300, unitsOut='kmsq'))

# Generate abundance map
# requires 'plot predictions per km2' section from JASG analysis earlier to have been run

# Newdata
nd <- data.frame(habitat=centroids$cancov.scaled, elev=centroids$elev.scaled)

# Generate predictions usingn newdata
pr <- predict(mod, type='lam', newdata=nd)

# Format for terra 'xyz'
xyz <- cbind(centroids[,c("X","Y")], Z=pr$Predicted)

# Create raster and mask
r <- rast(xyz, type="xyz")
rmask <- mask(r, vect(mask))

plot(rmask,main="Estimates based on unmarked analysis")