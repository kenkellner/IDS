#----------------------------------------------------------------
# Kéry et al. 2022: Integrated distance sampling models for simple point counts 
# Case study
# Analyse Distance Sampling (DS) data from Oregon2020 together with EBird Point Count (PC) data and detection/nondetection (DND) data 
# Marc Kéry, Nicolas Strebel, Tyler Hallman, Kenneth F. Kellner
# March 2022 
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
  A_PC_DND = A_PC_DND # Area associated with the PC and DND data

))
 

# Model
cat(file="IDS3_CaseStudy.txt","
    model{
      
      # Priors for parameters
      # -------------------------------------------------------------------
      # Separate parameters in the perceptibility component of the detection model
      alpha0 ~ dnorm(0, 0.01)    # For proper DS data set
      alpha1 ~ dnorm(0, 0.01)    # Canopy cover
      alpha2 ~ dnorm(0, 0.01)    # Urban area

      # Shared parameters in the availability component of the detection model
      gamma0 ~ dnorm(0, 0.1)  # Singing rate intercept
      gamma1 ~ dnorm(0, 0.1)  # Effect of day of year on singing rate
      gamma2 ~ dnorm(0, 0.1)  # Effect of day of year on singing rate (quadratic)
      gamma3 ~ dnorm(0, 0.1)  # Effect of time of day on singing rate
      gamma4 ~ dnorm(0, 0.1)  # Effect of time of day on singing rate (quadratic)

      # Shared parameters in the abundance model
      beta0 ~ dnorm(0, 0.01)  # Abundance intercept
      beta1 ~ dnorm(0, 0.01)  # Effect of habitat (canopy cover) on abundance
      beta2 ~ dnorm(0, 0.01)  # Effect of habitat (canopy cover) on abundance (quadratic)
      beta3 ~ dnorm(0, 0.01)  # Effect of elevation on abundance
      beta4 ~ dnorm(0, 0.01)  # Effect of elevation on abundance (quadratic)
      
      
      # Submodel for the DS data
      # -------------------------------------------------------------------
      # Hierarchical construction of the likelihood
      # Model for binned distance observations of every detected individual
      for(i in 1:nind){       # Loop over all detected individuals
        dclass[i] ~ dcat(fc[siteDS[i],])               # Part 1 of HM
      }
      
      # Construction of the cell probabilities for the nD distance bands
      # This is for the truncation distance for the DS data (here, 0.3)

      for(s in 1:nsites_DS){    # Loop over all sites in data set 1
        for(g in 1:nD){       # midpt = mid-point of each distance band
          log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s]^2)
          pi[s,g] <- ((2 * midpt[g] ) / newB^2) * delta # prob. per interval
          f[s,g] <- p[s,g] * pi[s,g]
          fc[s,g] <- f[s,g] / pcap[s]
        }
        # Rectangular approx. of integral that yields the Pr(capture)
        pcap[s] <- sum(f[s,])
        
        ### Log-linear models on abundance, detectability and availability
        # Abundance
        log(lambda1[s]) <- beta0 + beta1 * habitat_DS[s] + beta2 * pow(habitat_DS[s],2) + beta3 * elev_DS[s] + beta4 * pow(elev_DS[s],2) # Log-linear model for abundance 
        
        # Detectability
        log(sigma[s]) <- alpha0 + alpha1 * cancovdetect_DS[s] + alpha2 * urbandetect_DS[s]   # Log-Linear model for detection probability

        # Availability
        log(phi[s]) <- gamma0 + gamma1 * day_DS[s] + gamma2 * pow(day_DS[s],2) + gamma3 * time_DS[s] + gamma4 * pow(time_DS[s],2)  # Log-linear model for availability
        theta[s] <- 1-exp(-DSduration[s]*phi[s])  # Effect of duration on availability

        # Multiply availability with detection probability
        pDS[s] <- pcap[s] * theta[s]

        ### Binomial mixture part (in conditional multinomial specification)
        ncap[s] ~ dbin(pDS[s], N1[s])  # Part 2 of HM: number captured
        N1[s] ~ dpois(A_DS * lambda1[s]) 

      }
      
      # Submodel for the PC data
      # -----------------------------------------
      # Parameters on abundance, availability and detectability are al shared among the submodels for the different data sets
      # Note, though, that since we're modeling ABUNDANCE as basic parameter in
      # the state model, and NOT density, we must accommodate
      # the different areas to which the parameters of the abundance model
      # in the two data sets pertain. 
      # We do this with the 'N.scaling' proportionality (this is an offset).
      
      # Likelihood for the PC data
      # Note the exp(offset) equal to 'N.scaling', which is the ratio 
      # A_PC(eBird) / A_DS
      # Parameters on abundance are shared among all three submodels (DS, PC, DND), parameters on sigma are shared among PC and DND
      
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
        log(lambda2[s]) <- beta0 + beta1 * habitat_PC[s] + beta2 * pow(habitat_PC[s],2) + beta3 * elev_PC[s] + beta4 * pow(elev_PC[s],2) # Log-linear model on abundance 
    
        # Detectability
        log(sigmaPC[s]) <- alpha0 + alpha1 * cancovdetect_PC[s] + alpha2 * urbandetect_PC[s] # Log-Linear model for detection probability
        
        # Availability
        log(phi2[s]) <- gamma0 + gamma1 * day_PC[s] + gamma2 * pow(day_PC[s],2) + gamma3 * time_PC[s] + gamma4 * pow(time_PC[s],2) 
        theta2[s] <- 1-exp(-ebirdDuration[s]*phi2[s]) # Effect of duration on availability

        # Multiply availability with detection probability
        pPC[s] <- pcap2[s] * theta2[s] 

        ### Binomial mixture part 
        N2[s] ~ dpois(A_PC_DND * lambda2[s])
        counts[s] ~ dbinom(pPC[s], N2[s]) 
      }

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
        log(lambda3[s]) <- beta0 + beta1 * habitat_DND[s] + beta2 * pow(habitat_DND[s],2) + beta3 * elev_DND[s] + beta4 * pow(elev_DND[s],2) # Log-linear model on abundance 
        
        # Detectability
        log(sigmaDND[s]) <- alpha0 + alpha1 * cancovdetect_DND[s] + alpha2 * urbandetect_DND[s] # Log-Linear model for detection probability
        
        # Availability
        log(phi3[s]) <- gamma0 + gamma1 * day_DND[s] + gamma2 * pow(day_DND[s],2) + gamma3 * time_DND[s] + gamma4 * pow(time_DND[s],2) 
        theta3[s] <- 1-exp(-DNDduration[s]*phi3[s]) # Effect of duration on availability

        # Temporal scaling based on duration
        pDND[s] <- pcap3[s] * theta3[s] 

        ### Royle-Nichols part
        N3[s] ~ dpois(A_PC_DND * lambda3[s])
        y[s] ~ dbern(1 - (1 - pDND[s])^N3[s])
      }
    }
    ")

# Inits
Nst1 <- ncap + 1
Nst2 <- counts + 1
Nst3 <- 1*(DND$count>0) + 1
inits <- function(){list(alpha0 = rnorm(1, 0, 0.1), alpha1 = rnorm(1, 0, 0.1), alpha2 = rnorm(1, 0, 0.1), 
                         gamma0 = rnorm(1, 0, 0.1), gamma1 = rnorm(1, 0, 0.1), gamma2 = rnorm(1, 0, 0.1), 
                         gamma3 = rnorm(1, 0, 0.1), gamma4 = rnorm(1, 0, 0.1),
                         beta0 = rnorm(1, 0, 0.1), beta1 = rnorm(1, 0, 0.1), beta2 = rnorm(1, 0, 0.1), 
                         beta3 = rnorm(1, 0, 0.1), beta4 = rnorm(1, 0, 0.1), 
                         N1 = Nst1, N2 = Nst2, N3 = Nst3)}  

# Params to save
params <- c("alpha0", "alpha1", "alpha2", 
            "gamma0", "gamma1", "gamma2", "gamma3", "gamma4", 
            "beta0", "beta1", "beta2", "beta3", "beta4")

# MCMC settings
na <- 10  ;  nc <- 3  ;  ni <- 12  ;  nb <- 2  ;  nt <- 2 # xx test
#na <- 1000;  nc <- 3;  ni <- 2000;  nb <- 1000;  nt <- 3 # 54 min
#na <- 1000;  nc <- 3;  ni <- 20000;  nb <- 10000;  nt <- 5 # 4.8 hrs

# Launch JAGS, check convergence and summarize posteriors ----
# Launch
library(jagsUI)
start <- Sys.time()
out <- jags(bdata, inits, params, "IDS3_CaseStudy.txt", n.adapt = na,
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,parallel = TRUE)
difftime(Sys.time(),start)
traceplot(out)
print(out, 2) 

# Estimate abundance throughout study area for each sample ----
# Load environmental data
head(centroids <- read.table("data/env_centroids.csv",header=T,sep=";",as.is = T))

# adjust variables as they were in the model data
centroids$cancov.scaled <- (centroids$cancov - 50)/100
centroids$elev.scaled <- (centroids$elev - mean.elev)/sd.elev

# Create a function to do this
predict.abundances <- function(environmental = environmental, parameters = parameters){

  prediction <- exp(parameters[1] + parameters[2] * environmental$cancov.scaled + parameters[3] * I(environmental$cancov.scaled^2) + 
                                    parameters[4] * environmental$elev.scaled + parameters[5] * I(environmental$elev.scaled^2))
  return(prediction)
  
}

# Extract model output samples as a matrix
samples.matrix <- as.matrix(out$samples)

# Predict
centroids.predictions <- apply(samples.matrix[,9:13], 1, FUN = predict.abundances, environmental = centroids)

# Select median
centroids.predictions.median <- apply(centroids.predictions, 1, median)
centroids.predictions.output <- cbind(centroids, centroids.predictions.median)

# Calculate population estimate with credible intervals
pop.estimates <- apply(centroids.predictions,2, sum)

(mean.pop.estimate <- mean(pop.estimates))
(uci.pop.estimate <- quantile(pop.estimates, probs = 0.975))
(lci.pop.estimate <- quantile(pop.estimates, probs = 0.025))

# Plot predictions per km2 ----
library(sf); library(terra)

#load in file with desired crs
polygon.grid <- st_read("data/geodata/grid_1km.shp")

mask <- st_read("data/geodata/StudyArea.shp")

centroids.predictions.output_sf = st_as_sf(centroids.predictions.output, coords = c("X", "Y"), crs = st_crs(polygon.grid))

test.raster <- rast(nlyrs=1, crs = crs(vect(polygon.grid)), extent = ext(polygon.grid), resolution = 1000)

centroids.predictions.output_rast <- rasterize(vect(centroids.predictions.output_sf), test.raster, field = "centroids.predictions.median")

centroids.predictions.output_rast <- mask(centroids.predictions.output_rast, vect(mask)) # mask study area

# Plot density map
plot(centroids.predictions.output_rast,main="Estimates based on jags analysis")

# Analysis in unmarked ----

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