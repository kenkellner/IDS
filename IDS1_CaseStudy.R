#----------------------------------------------------------------
# Kéry et al. 2024: Integrated distance sampling models for simple point counts 
# Case study
# Analyse Distance Sampling (DS) data from Oregon2020 together with 
# EBird Point Count (PC) data 
# Marc Kéry, Nicolas Strebel, Tyler Hallman, Kenneth F. Kellner
# March 2022 - June 2023 
#---------------------------------------------------------------

# Load data ----
# Distance Sampling (DS) data from Oregon2020
head(DS <- read.table("data/amro_obsdetection.csv",header=T,sep=";",as.is = T))
DS$Distance <- DS$Distance/1000 # Convert to km

# Table including covariates and Oregon2020 counts (particularly the zeros) 
head(sitecovs <- read.table("data/amro_sitecovs_simplified.csv",header=T,sep=";",as.is = T))

# eBird
head(PC <- read.table("data/amro_ebird_simplified.csv",header=T,sep=";",as.is = T))

# Prepare data ----
# Prepare eBird data 
counts <- PC$count

# Link DS and zero observations (from sitecovs)
names(sitecovs)[names(sitecovs) == "Latitude.x"] <- "Latitude" # Homogenize names
names(sitecovs)[names(sitecovs) == "Longitude.x"] <- "Longitude"
nondetections <- sitecovs[sitecovs$Count==0,] # Select zero counts
nondetections <- nondetections[,names(nondetections) %in% names(DS)] # Columns from nondetection data that are not in DS can be removed
nondetections[,names(DS)[!names(DS) %in% names(nondetections)]] <- NA # Add columns to nondetection data that are in DS only 
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
(A_PC <- pi * fullDistance^2 )

# Number of sites per data set
nsites_DS <- length(unique(DS$site))      # Data set 1
nsites_PC <- length(unique(apply(PC[,c("latitude.x","longitude.x")],1,paste,collapse="_")))      # PC data

# Attribute site number (according to DS) to each row in sitecovs
sitecovs$sitenr <- DS$site[match(sitecovs$Location_ID,DS$Location_ID)]

# Prepare objects for analysis in BUGS ----
# values to scale covariate data
mean.elev <- mean(c(sitecovs$Elevation_mean315, PC$Elevation_mean315))
sd.elev <- sd(c(sitecovs$Elevation_mean315, PC$Elevation_mean315))
mean.minutes.since.dawn <- mean(c(sitecovs$MinutesSinceDawn, PC$minutes_since_dawn))
sd.minutes.since.dawn <- sd(c(sitecovs$MinutesSinceDawn, PC$minutes_since_dawn))
mean.jdate <- mean(c(sitecovs$JulianDate, PC$julian_date))
sd.jdate <- sd(c(sitecovs$JulianDate, PC$julian_date))
  
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
  fullDistance = fullDistance, # Assumed radius for PC
  nDfull = nDfull,             # length of (latent) distance bins in PC data
  midptFull = midptFull,       # Midpoints of (latent) distance bins
  A_PC = A_PC,                 # Area associated with the PC data
  nyear = length(unique(c(PC$year,DS$Year))) # number of years in data

))
 
 

# Model definition in the BUGS language
# Add random noise in the detection function intercepts and 
# allow the variance to be different by design (DS vs. PC)
# Note all distances are in units of 1km (and area in 1km2 units)
cat(file="IDS1_CaseStudy_model_5.txt","
    model{
      
      # Priors for all parameters in IDS model
      # -------------------------------------------------------------------
      # Partly different parameters in the detectability-perceptibility component 
	    # of the detection model for the DS and the PC portions of the data
      for(d in 1:2){           # d indexes the two data types: DS vs. PC
        alpha0[d] <- log(mean.sigma[d])  # sigma intercept on log scale and ...
        mean.sigma[d] ~ dunif(0, 1)    # ... on the natural scale (0 - 1 km)
        tau.eps[d] <- pow(sd.eps[d], -2)
		    sd.eps[d] ~ dt(0, 1, 1)I(0,)     # Magnitude of that noise
      }
      # Covariates:
      alpha1 ~ dnorm(0, 1)     # Canopy cover
      alpha2 ~ dnorm(0, 1)     # Urban area

      # Shared parameters in the availability component of the detection model
      gamma0 <- log(mean.phi)  # phi intercept on log scale and ...
      mean.phi ~ dunif(0, 1)   # ... on natural scale
	    # mean.phi ~ dbeta(2, 2)   # ... on natural scale
      gamma1 ~ dnorm(0, 1)     # Effect of day of year on singing rate
      gamma2 ~ dnorm(0, 1)     # Effect of day of year on singing rate (quadratic)
      gamma3 ~ dnorm(0, 1)     # Effect of time of day on singing rate
      gamma4 ~ dnorm(0, 1)     # Effect of time of day on singing rate (quadratic)

      # Shared parameters in the abundance model
      # Random intercept for each year
      for (i in 1:nyear) {     # Loop over 7 years
        ann.beta0[i] ~ dnorm(beta0, tau.beta0)
      }
      beta0 <- log(mean.lambda)  # lambda intercept on log scale and ...
      mean.lambda ~ dnorm(0, 0.001)I(0,) # ... on natural scale
      tau.beta0 <- pow(sd.beta0,-2) 
      sd.beta0 ~ dt(0, 1, 2)I(0,) # Magnitude of noise in lambda intercept
      # Covariates:
      beta1 ~ dnorm(0, 1)      # Effect of habitat (canopy cover) on abundance
      beta2 ~ dnorm(0, 0.1)    # Effect of habitat (canopy cover) on abundance (quadratic)
      beta3 ~ dnorm(0, 1)      # Effect of elevation on abundance
      beta4 ~ dnorm(0, 1)      # Effect of elevation on abundance (quadratic)
      
      # Submodel for the DS data
      # -------------------------------------------------------------------
      # Hierarchical construction of the likelihood
      # Model for binned distance observations of every detected individual
      for(i in 1:nind){       # Loop over all detected individuals
        dclass[i] ~ dcat(fc[siteDS[i],])               # Part 1 of HM
      }
      
      # Construction of the cell probabilities for the nD distance bands
      # This is for the truncation distance for the DS data (here, newB = 0.3 km)

      for(s in 1:nsites_DS){  # Loop over all sites in data set 1
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
		      eps1[s] ~ dnorm(0, tau.eps[1])        # Note here eps1 has one precision and below another

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
      # ------------------------
      # Parameters on abundance shared among both submodels (DS, PC),
	    # parameters on availability are mostly shared (same coefficients, different ran-ef variance),
      #	parameters on sigma are different among DS and PC
      # For spatial reconciliation of the data, we use A as an offset
      
      # Likelihood for the PC data
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
        log(sigmaPC[s]) <- alpha0[2] + alpha1 * cancovdetect_PC[s] + alpha2 * urbandetect_PC[s] + eps2[s] # Log-Linear model for detection probability
        eps2[s] ~ dnorm(0, tau.eps[2])

        # Availability
        log(phi2[s]) <- gamma0 + gamma1 * day_PC[s] + gamma2 * pow(day_PC[s],2) + gamma3 * time_PC[s] + gamma4 * pow(time_PC[s],2) 
        theta2[s] <- 1-exp(-ebirdDuration[s]*phi2[s]) # Effect of duration on availability

        # Multiply availability with detection probability
        pPC[s] <- pcap2[s] * theta2[s] 

        ### Binomial mixture part 
        N2[s] ~ dpois(A_PC * lambda2[s])
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

      # Compute Bayesian p-value for both portions of the data
      bpvDS <- step(fitDS.new - fitDS)
      bpvPC <- step(fitPC.new - fitPC)
    }
    ")

# Inits
Nst1 <- ncap + 1
Nst2 <- counts + 1
inits <- function(){list(mean.sigma = runif(2, 0, 1), alpha1 = rnorm(1, 0, 0.1), alpha2 = rnorm(1, 0, 0.1), 
                         mean.phi = runif(1, 0, 1), gamma1 = rnorm(1, 0, 0.1), gamma2 = rnorm(1, 0, 0.1), 
                         gamma3 = rnorm(1, 0, 0.1), gamma4 = rnorm(1, 0, 0.1),
                         mean.lambda = runif(1, 1, 60), beta1 = rnorm(1, 0, 0.1), beta2 = rnorm(1, 0, 0.1), 
                         beta3 = rnorm(1, 0, 0.1), beta4 = rnorm(1, 0, 0.1), 
						 sd.eps = runif(2, 0.1, 1), N1 = Nst1, N2 = Nst2)}  

# Params to save
params <- c("mean.sigma", "alpha0", "alpha1", "alpha2", "sd.eps",
            "mean.phi", "gamma0", "gamma1", "gamma2", "gamma3", "gamma4", 
            "mean.lambda", "beta0", "beta1", "beta2", "beta3", "beta4", "sd.beta0", "ann.beta0",
			"fitDS", "fitDS.new", "fitPC", "fitPC.new", "bpvDS", "bpvPC")

# MCMC settings
na <- 10  ;  nc <- 4  ;  ni <- 12  ;  nb <- 2  ;  nt <- 2 # test, 30 sec
na <- 10000;  nc <- 4;  ni <- 120000;  nb <- 60000;  nt <- 60   # As for the paper
na <- 10000;  nc <- 10;  ni <- 400000;  nb <- 200000;  nt <- 200 # takes a while...


# Launch JAGS, check convergence and summarize posteriors ----
library(jagsUI)
start <- Sys.time()
set.seed(123)
out5A <- jags(bdata, inits, params, "IDS1_CaseStudy_model_5.txt", n.adapt = na,
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
difftime(Sys.time(),start)
traceplot(out5A)
print(out5A, 2) 



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
samples.matrix <- as.matrix(out5A$samples)
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

# Write raster
writeRaster(centroids.predictions.output_rast, "IDS_AMRO_PredictedAbundance.tif", filetype = "GTiff", overwrite = TRUE)



# Analysis in unmarked  ----
# Because unmarked IDS() function doesn't support random effects we fit a different model
# Difference of the model fitted in unmarked compared to IDS1_CaseStudy_5.txt:
#  - no annual random effects on abundance, 
#  - no different levels of heterogeneity in the detection functions estimated for different portions of the data 

# Load package
# Needs a version of unmarked with the IDS() function; see the README file in this folder.
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


#### Fit model ####

mod <- IDS(lambdaformula = ~habitat + I(habitat^2) + elev + I(elev^2),
           detformulaDS = ~cancov + urban,
           detformulaPC = NULL,
           dataDS = umf_DS, dataPC = umf_PC,
           availformula = ~day + I(day^2) + time + I(time^2),
           durationDS = dur_DS, durationPC = dur_PC, 
           maxDistPC = bdata$fullDistance, 
           K=300, unitsOut='kmsq')

summary(mod)

# Generate abundance map
# requires 'plot predictions per km2' section from JAGS analysis earlier to have been run

# Newdata
nd <- data.frame(habitat=centroids$cancov.scaled, elev=centroids$elev.scaled)

# Generate predictions usingn newdata
pr <- predict(mod, type='lam', newdata=nd)

# Format for terra 'xyz'
xyz <- cbind(centroids[,c("X","Y")], Z=pr$Predicted)

# Create raster and mask
r <- rast(xyz, type="xyz", crs=crs(mask))
rmask <- mask(r, vect(mask))

plot(rmask,main="Estimates based on unmarked analysis")
