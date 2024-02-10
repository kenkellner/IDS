
# Supporting Information for:
# Kéry, M., Royle, J.A., Hallman, T., Robinson, W.D., Strebel, N., Kellner, K.F., 2024:
# Integrated distance sampling models for simple point counts, Ecology.


# Written by the authors, Nov 2022


# Simulation 1: Identifiability of separate observation process parameters in IDS1 and IDS2
# ----------------------------------------------------------------------------------------
#
# In this first simulation we test for identifiability of the new IDS models in their
# most basic set-up: without any covariates. We just have intercepts for density (shared)
# between two data types, and two constants for the detection function scale sigma in 
# the distance sampling (DS) data and in the point count (PC) or detection/nondetection (DND)
# data. In the two latter, we vary sigma randomly between 0.1 and 1.3, and the identifiability
# test then amounts (for the most part) to comparing the true and the estimated values of \
# sigma(PC) or sigma(DND). 
# Distance sampling code in JAGS is for binned distances and taken from Chapter 8
# in Kéry & Royle, Academic Press, 2016.
# We test identifiability separately for IDS1 and IDS2.
# This is the code for Simulation 1 in the paper.



# (1a) Combination of Distance-sampling and point count data (= IDS1 model)
# -------------------------------------------------------------------------

library(AHMbook)
library(jagsUI)

# Number of simulation reps
simrep <- 3          # for a test run
simrep <- 1000       # for real

# Create R object to save results (posterior summary table from JAGS)
estimates1 <- array(NA, dim = c(1014, 11, simrep))

# Vector to save the true values of sigmaEB
true.dataDS1 <- list()
true.dataPC <- list()
true.sigmaPC <- numeric(simrep)


# Launch the simulation
for(k in 1:simrep){        # Loop over simreps

  cat(paste('\n*** IDS1: Sim Number', k, '***\n'))

  # Draw a new value of sigPC from a U(0.1, 1.3) and save it
  sigPC <- runif(1, 0.1, 1.3)
  true.sigmaPC[k] <- sigPC
  cat(paste('Current value of sigmaPC is', round(sigPC, 3), '\n '))

  # Simulate a distance sampling data set of size 250 sites and save it
  dat1 <- simHDS(type="point", nsites = 250,                         # WAS 250 sites
    mean.lambda = 100, beta.lam = 0, mean.sigma = 1, beta.sig = 0, 
    B = 5, discard0 = FALSE, show.plot = TRUE)
  true.dataDS1[[k]] <- dat1

  # Simulate a simple point count data set of size 1,000 sites and save it
  # First, we simulate a DS data set and then 'degrade' it by 'losing' the distance data
  dat2 <- simHDS(type="point", nsites = 1000,                          # WAS 1000 sites                   
    mean.lambda = 100, beta.lam = 0, mean.sigma = sigPC, beta.sig = 0, 
    B = 5, discard0 = FALSE, show.plot = TRUE)
  true.dataPC[[k]] <- dat2

  # Create the PC data subset (call this data set 2)
  data2 <- dat2$data
  ncap2 <- tapply(data2[,5], data2[,1], function(x) sum(!is.na(x)))
  C2 <- as.vector(ncap2)  # Number of individuals detected per site

  # Create the DS data subset (call this data set 1)
  data1 <- dat1$data

  # Subset the DS data to a truncation distance (B) of 2
  data1[,'d'][data1[,'d'] > 2] <- NA
  sum(!is.na(data1[,'d'])) # print out number of detections inside B

  # Prepare the DS data and summarize some too
  # Number of individuals detected (or 'captured') per site
  ncap <- tapply(data1[,5], data1[,1], function(x) sum(!is.na(x)))
  ncap <- as.vector(ncap)  # Number of individuals detected per site

  # Prepare other data
  # Site ID of each non-missing observation
  site1 <- data1[!is.na(data1[,'d']),1]

  # Define new B for the DS data set (data 1)
  newB <- 2

  # The old B for getting the unlimited distance pbar 
  fullDistance <- 5

  # Distance bin width for rectangle approximation of integral
  # over detection function, to get average detection probability
  delta <- 0.1

  # Make mid-points and chop up distance data
  # Both for newB and the fullDistance
  midpt <- seq(delta/2, newB, delta)
  midptFull <- seq(delta/2, fullDistance, delta)

  # Bin the data: convert distance data to distance category
  dclass <- data1[,'d'] %/% delta + 1

  # Number of distance intervals
  nD <- length(midpt)
  nDfull <- length(midptFull)

  # Observed categorical observations
  # (Make restriction to non-missing distances for newB = 2)
  dclass <- dclass[!is.na(dclass)]   # toss out missing dclasses

  # Total number of individuals detected (after applying new B = 2)
  nind <- length(dclass)

  # Compute the areas associated with the data for each site
  A1 <- pi * newB^2         # [1] 12.56637, area associated with DS surveys
  A2 <- pi * fullDistance^2 # [1] 78.53982, area associated with PC surveys

  # Bundle and summarize data set to be fed into JAGS
  nsites1 <- dat1$nsites   # Data set 1
  nsites2 <- length(C2)    # Data set 2
  str(bdata <- list(nsites1 = nsites1, nind = nind, newB = newB, nD = nD,
    midpt = midpt, fullDistance = fullDistance, nDfull = nDfull, 
    midptFull = midptFull, delta = delta, ncap = ncap, dclass = dclass,
    site1 = site1, A1 = A1, nsites2 = nsites2, C2 = C2, A2 = A2))



# Specify model in BUGS language
cat(file="modelIDS1.txt","
model{

# Priors for parameters
# Separate parameters in the detection model
alpha0DS ~ dnorm(0, 0.01)    # For proper DS data set
alpha0PC ~ dnorm(0, 0.01)    # For simple point count data set
# Shared parameters in the abundance model
beta0 ~ dnorm(0, 0.01)       # Abundance intercept

# Submodel for the Distance sampling data (= data set 1)
# -------------------------------------------------------------------
# Hierarchical construction of the likelihood; see Chap. 8 in AHM1 book
# Model for binned distance observations of every detected individual
for(i in 1:nind){       # Loop over all detected individuals
  dclass[i] ~ dcat(fc[site1[i],])               # Part 1 of HM
}

# Construction of the cell probabilities for the nD distance bands
# This is for the truncation distance for the DS data (here, 2)
for(s in 1:nsites1){    # Loop over all sites in data set 1
  for(g in 1:nD){       # midpt = mid-point of each distance band
    log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s]^2)
    pi[s,g] <- ((2 * midpt[g] ) / newB^2) * delta # prob. per interval
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g] / pcap[s]
  }
  # Rectangular approx. of integral that yields the Pr(capture)
  pcap[s] <- sum(f[s,])

  # Remaining part of model in conditional multinomial specification
  ncap[s] ~ dbin(pcap[s], N1[s])  # Part 2 of HM: number captured
  # Note A1 is offset (survey area) and lambda1 becomes per-unit abundance, i.e., density
  N1[s] ~ dpois(A1 * lambda1[s])  # Part 3 of HM: true number there

  # Log-linear models for density and detection in DS data set
  log(lambda1[s]) <- beta0    # Just a constant now
  log(sigma[s]) <- alpha0DS   # Model for sigma just a constant now
}

# Now we compute the average p over all distance bands, pbar, 
# for the unlimited-distance surveys (i.e., for B = 5).
# It is this pcap2 (or here, it's average pbar2) which we will use 
# below in the submodel for the simple point count data (= data set 2)
# Note that for now pcap2 is identical for all PC sites
# since we don't have any detection covariates in the model
for(s in 1:nsites2){        # Loop over all sites in data set 2
  # Linear model for detection in PC data set
  log(sigmaPC[s]) <- alpha0PC      # just a constant for now

  for(g in 1:nDfull){       # midpt = mid-point of each distance band
    log(p2[s,g]) <- -midptFull[g] * midptFull[g] / (2 * sigmaPC[s]^2)
    pi2[s,g] <- ((2 * midptFull[g] ) / fullDistance^2) * delta # prob. per interval
    f2[s,g] <- p2[s,g] * pi2[s,g]
  }
  # Rectangular approx. of integral that yields the Pr(capture)
  pcap2[s] <- sum(f2[s,])
}
pbar2 <- mean(pcap2)  # Detection prob. out to fullDistance (here, B = 5)


# Submodel for the simple point count data (= data set 2)
# -------------------------------------------------------
# We have defined the priors already above:
#  - sigma in the detection function is *NOT* shared between data sets 1 and 2
#  - in contrast, the abundance parameter (beta0) is shared
# Note, though, that since we're modeling ABUNDANCE as basic parameter in
# the state model, and NOT density, we must accommodate the different 
# areas with which the parameter in the abundance model is associated.
# We do this by using A2 as an offset below.

# Likelihood for the simple PC data
# This is for the Poisson/Binomial variant of this submodel
# Remember that beta0 is a shared param, but not sigma
for(s in 1:nsites2){
  # Note A2 is again an offset (survey area) and 
  # lambda2 then becomes per-unit abundance, i.e., density
  N2[s] ~ dpois(A2 * lambda2[s])
  log(lambda2[s]) <- beta0
  C2[s] ~ dbinom(pbar2, N2[s])
}

# Derived parameters
# For data set 1 (with proper DS surveys)
Ntotal1 <- sum(N1[])            # Total pop. size over all DS sites (within newB = 2)
tot.area1 <- nsites1 * A1       # Total surveyed area within newB for all DS sites
D1 <- Ntotal1 / tot.area1       # Average density over all DS sites (within newB)
mean.lambda <- exp(beta0)       # Average density in both data sets
pbar1 <- mean(pcap)             # Average detection prob in DS surveys out to newB
mean.sigmaDS <- mean(sigma)     # Average sigma in DS surveys

# For data set 2 (with simple point counts)
Ntotal2 <- sum(N2[])            # Total pop. size over all PC survey areas (within B=5)
tot.area2 <- nsites2 * A2       # Total surveyed area for all PC sites
D2 <- Ntotal2 / tot.area2       # Average density over all sites
Ntotal <- Ntotal1 + Ntotal2     # Ntotal within the B(DS)=2 and B(PC)=5 circles
D <- Ntotal / (tot.area1 + tot.area2)# Best (combined) estimate of density
mean.sigmaPC <- mean(sigmaPC)   # Average sigma in eBird surveys
}
")

  # Inits
  Nst1 <- ncap + 1
  Nst2 <- C2 + 1
  inits <- function(){list(alpha0DS = rnorm(1, 0, 0.1), 
    alpha0PC = rnorm(1, -1, 0.1), beta0 = rnorm(1, 0, 0.1), 
     N1 = Nst1, N2 = Nst2)}
  
  
  # Params to save
  params <- c("alpha0DS", "alpha0PC", "beta0", "mean.sigmaDS",
    "mean.sigmaPC", "Ntotal1", "D1", "mean.lambda", "pbar1", "pbar2",
	"Ntotal2", "D2", "D", "N2")
  
  # MCMC settings
  na <- 30;  nc <- 3;  ni <- 150;  nb <- 50;  nt <- 2 # for test only
  na <- 3000;  nc <- 3;  ni <- 15000;  nb <- 5000;  nt <- 10 # for real
  
  # Launch JAGS (ART xx min), check convergence and summarize posteriors
  out1 <- jags(bdata, inits, params, "modelIDS1.txt", n.adapt = na,
    n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
    parallel = TRUE)
  print(out1, 2)
  
  # Save posterior summary
  estimates1[,,k] <- out1$summary

}            # End of loop over simreps



# Compute estimation errors for the site-level abundance N
# in the PC part of the data (used later)

# Choose number of simreps to summarize
Kmax <- 3                 # For test run
Kmax <- 1000              # For real

# Compute errors of the N2 estimates (abundance in PC data)
# Get point estimates
tmp_esti <- estimates1[14:1013,1, 1:Kmax]

# Get true values
tmp_truth <- array(NA, dim = dim(tmp_esti))
for(k in 1:Kmax){
  tmp_truth[,k] <- true.dataPC[[k]]$N.true
}

# Get absolute difference
absDiff1 <- tmp_esti - tmp_truth

# Tests and replacement of Inf with NA
sum(is.na(absDiff1)) ; sum(absDiff1 == 'Inf') ; sum(absDiff1 == '-Inf')

# Numerical summary
mean(absDiff1) ; sd(absDiff1)  ;  summary(c(absDiff1))





###########################################################################



# (1b) Combination of Distance-sampling and detection/nondetection data (= IDS2 model)
# ------------------------------------------------------------------------------------

# Number of simulation reps
simrep <- 3      # for test only
simrep <- 1000   # for real

# Create R object to save results (posterior summary table)
 estimates2 <- array(NA, dim = c(1014, 11, simrep))

# Lists to save the data and vector to save the true values of sigmaEB
 true.dataDS2 <- list()
 true.dataDND <- list()
 true.sigmaDND <- numeric(simrep)


# Launch the simulation
for(k in 1:simrep){        # Loop over simreps

  cat(paste('\n*** IDS2: Sim Number', k, '***\n'))

  # Draw a new value of simDND from a U(0.1, 1.3) and save it now
  sigDND <- runif(1, 0.1, 1.3)
  true.sigmaDND[k] <- sigDND
  cat(paste('\nCurrent value of sigmaDND is', round(sigDND, 3), '\n\n'))

  # Simulate another DS data set of size 250 sites and save it
  dat1 <- simHDS(type="point", nsites = 250, 
    mean.lambda = 100, beta.lam = 0, mean.sigma = 1, beta.sig = 0, 
    B = 5, discard0 = FALSE, show.plot = TRUE)
  true.dataDS2[[k]] <- dat1

  # Simulate a DND data set of size 1,000 sites and save
  # First, we simulate a DS data set and then 'degrade' it by 'losing' the distance data
  # and turning the counts into binary DND data
  dat2 <- simHDS(type="point", nsites = 1000, 
    mean.lambda = 100, beta.lam = 0, mean.sigma = sigDND, beta.sig = 0, 
    B = 5, discard0 = FALSE, show.plot = TRUE)
  true.dataDND[[k]] <- dat2

  # New computation of C2(analogous to what we do below for data set 1)
  data2 <- dat2$data
  ncap2 <- tapply(data2[,5], data2[,1], function(x) sum(!is.na(x)))
  C2 <- as.vector(ncap2)  # Number of individuals detected per site

  # Turn the DND data into detection/nondetection (which we call y)
  y <- as.numeric(C2 > 0)

  # Create the DS data subset (call this data set 1)
  data1 <- dat1$data

  # Subset the DS data to a truncation distance (B) of 2
  data1[,'d'][data1[,'d'] > 2] <- NA
  sum(!is.na(data1[,'d'])) # print out number of detections inside B

  # Prepare the DS data and summarize some too
  # Number of individuals detected (or 'captured') per site
  ncap <- tapply(data1[,5], data1[,1], function(x) sum(!is.na(x)))
  ncap <- as.vector(ncap)  # Number of individuals detected per site

  # Prepare other data
  # Site ID of each non-missing observation
  site1 <- data1[!is.na(data1[,'d']),1]

  # Define new B for the DS data set (data 1)
  newB <- 2

  # The old B for getting the unlimited distance pbar 
  fullDistance <- 5

  # Distance bin width for rectangle approximation of integral
  delta <- 0.1

  # Make mid-points and chop up distance data
  # Both for newB and the fullDistance
  midpt <- seq(delta/2, newB, delta)
  midptFull <- seq(delta/2, fullDistance, delta)

  # Convert distance to distance category
  dclass <- data1[,'d'] %/% delta + 1

  # Number of distance intervals
  nD <- length(midpt)
  nDfull <- length(midptFull)

  # Observed categorical observations
  # (Make restriction to non-missing distances for newB = 2)
  dclass <- dclass[!is.na(dclass)]   # toss out missing dclasses

  # Total number of individuals detected (after applying new B = 2)
  nind <- length(dclass)

  # Compute the areas associated with the data for each site
  A1 <- pi * newB^2         # [1] 12.56637, area associated with DS surveys
  A2 <- pi * fullDistance^2 # [1] 78.53982, area associated with DND surveys

  # Bundle and summarize data set to be fed into JAGS
  nsites1 <- dat1$nsites   # Data set 1 (with true DS surveys)
  nsites2 <- length(y)     # Data set 2 (with simple DND surveys)
  str(bdata <- list(nsites1 = nsites1, nind = nind, newB = newB, nD = nD,
    midpt = midpt, fullDistance = fullDistance, nDfull = nDfull, 
    midptFull = midptFull, delta = delta, ncap = ncap, 
    dclass = dclass, site1 = site1, A1 = A1, 
    nsites2 = nsites2, y = y, A2 = A2))


# Specify model in BUGS language
cat(file="modelIDS2.txt","
model{

# Priors for parameters
# Separate parameters in the detection model
alpha0DS ~ dnorm(0, 0.01)    # For proper DS data set
alpha0DND ~ dnorm(0, 0.01)   # For DND data set
# Shared parameters in the abundance/density model
beta0 ~ dnorm(0, 0.01)       # Abundance/density intercept

# Submodel for the proper distance sampling data (= data set 1)
# ------------------------------------------------------------
# Hierarchical construction of the likelihood
# Model for binned distance observations of every detected individual
for(i in 1:nind){       # Loop over all detected individuals
  dclass[i] ~ dcat(fc[site1[i],])               # Part 1 of HM
}

# Construction of the cell probabilities for the nD distance bands
# This is for the truncation distance for the DS data (here, 2)
for(s in 1:nsites1){    # Loop over all sites in data set 1
  for(g in 1:nD){       # midpt = mid-point of each distance band
    log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s]^2)
    pi[s,g] <- ((2 * midpt[g] ) / newB^2) * delta # prob. per interval
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g] / pcap[s]
  }
  # Rectangular approx. of integral that yields the Pr(capture)
  pcap[s] <- sum(f[s,])

  # Remaining part of model in conditional multinomial specification
  ncap[s] ~ dbin(pcap[s], N1[s])  # Part 2 of HM: number captured
  # Note A1 is offset (survey area 1): lambda1 becomes per-unit abundance, i.e., density
  N1[s] ~ dpois(A1*lambda1[s])       # Part 3 of HM: true number there

  # Log-linear models for abundance/density and detection in DS data set
  log(lambda1[s]) <- beta0  # Model for density also just a constant
  log(sigma[s]) <- alpha0DS # Model for sigma just a constant here
}

# Now we compute the average p over all distance bands, pbar, 
# for the unlimited-distance DND surveys (i.e., for B = 5).
# It is this pcap2 (or here, it's average pbar2) which we will use 
# below in the submodel for the DND data (= data set 2)
# Note that for now pcap2 is identical for all DND sites
# since we don't have any detection covariates in the model
for(s in 1:nsites2){        # Loop over all sites in data set 2
  # Linear model for detection in DND data set
  log(sigmaDND[s]) <- alpha0DND      # just a constant for now

  for(g in 1:nDfull){       # midpt = mid-point of each distance band
    log(p2[s,g]) <- -midptFull[g] * midptFull[g] / (2 * sigmaDND[s]^2)
    pi2[s,g] <- ((2 * midptFull[g] ) / fullDistance^2) * delta # prob. per interval
    f2[s,g] <- p2[s,g] * pi2[s,g]
  }
  # Rectangular approx. of integral that yields the Pr(capture)
  pcap2[s] <- sum(f2[s,])
}
pbar2 <- mean(pcap2)  # Detection prob. out to fullDistance (here, B = 5)


# Submodel for the DND data (= data set 2)
# -----------------------------------------
# We have defined the priors already above:
#  - sigma in the detection function is not shared between the data sets
#  - but the abundance/density parameter (beta0) is shared
# Note, though, that since we're modeling ABUNDANCE as basic parameter in
# the state model, and NOT density, we must accommodate the different 
# areas with which the parameter in the abundance model is associated.
# We do this again by using A2 as an offset.

# Likelihood for the DND data
# (has the form of a Royle-Nichols (RN) model)
# Remember that beta0 is a shared param, but not sigma
for(s in 1:nsites2){
  # Note A2 is again an offset (survey area) and 
  # lambda2 then becomes per-unit abundance, i.e., density
  N2[s] ~ dpois(A2 * lambda2[s])
  log(lambda2[s]) <- beta0       # Just a constant here

  # Crucially, here now is the RN observation model
  y[s] ~ dbern(1 - (1 - pbar2)^N2[s])
}

# Derived parameters
# For data set 1 (with proper DS surveys)
Ntotal1 <- sum(N1[])            # Total pop. size over all DS sites (within newB = 2)
tot.area1 <- nsites1 * A1       # Total surveyed area within newB for all DS sites
D1 <- Ntotal1 / tot.area1       # Average density over all DS sites (within newB)
mean.lambda <- exp(beta0)       # Average density in both data sets
pbar1 <- mean(pcap)             # Average detection prob in DS surveys out to newB
mean.sigmaDS <- mean(sigma)     # Average sigma in DS surveys

# For data set 2 (with DND surveys)
Ntotal2 <- sum(N2[])            # Total pop. size over all DND survey areas (within B=5)
tot.area2 <- nsites2 * A2       # Total surveyed area for all DND sites
D2 <- Ntotal2 / tot.area2       # Average density over all sites with DND surveys
Ntotal <- Ntotal1 + Ntotal2     # Ntotal within the B(DS)=2 and B(DND)=5 circles
D <- Ntotal / (tot.area1 + tot.area2)# Best (combined) estimate of density
mean.sigmaDND <- mean(sigmaDND) # Average sigma in DND surveys

}
")


  # Inits
  Nst1 <- ncap + 1
  Nst2 <- C2 + 1
  inits <- function(){list(alpha0DS = rnorm(1, 0, 0.1), 
    alpha0DND = rnorm(1, -1, 0.1), beta0 = rnorm(1, 0, 0.1), 
    N1 = Nst1, N2 = Nst2)}
  
  
  # Params to save
  params <- c("alpha0DS", "alpha0DND", "beta0", "mean.sigmaDS",
    "mean.sigmaDND", "Ntotal1", "D1", "mean.lambda",
    "pbar1", "pbar2", "Ntotal2", "D2", "D", "N2")
  
  # MCMC settings
  na <- 30;  nc <- 3;  ni <- 150;  nb <- 50;  nt <- 2 # for test only
  na <- 3000;  nc <- 3;  ni <- 15000;  nb <- 5000;  nt <- 10 # for real

  # Launch JAGS (ART xx min), check convergence and summarize posteriors
  out2 <- jags(bdata, inits, params, "modelIDS2.txt", n.adapt = na,
    n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
    parallel = TRUE)
  print(out2, 2)
  
  # Save posterior summary
  estimates2[,,k] <- out2$summary

}            # End of loop over simreps



# NOTE: Before summarization of results, you may have to toss out 'invalid'
# solutions. That is, model runs that have not converged.



# Compute errors of the N2 estimates
# (site-level abundance in sites with DND surveys only, used below)

# Choose number of simreps to summarize
Kmax <- 3                 # For test run
Kmax <- 1000              # For real

# Get point estimates
tmp_esti <- estimates2[14:1013,1, 1:Kmax]

# Get true values
tmp_truth <- array(NA, dim = dim(tmp_esti))
for(k in 1:Kmax){
  tmp_truth[,k] <- true.dataDND[[k]]$N.true
}

# Get absolute and relative differences
absDiff2 <- tmp_esti - tmp_truth

# Tests and replacement of Inf with NA
sum(is.na(absDiff2)) ; sum(absDiff2 == 'Inf') ; sum(absDiff2 == '-Inf')

# Numerical summary
mean(absDiff2) ; sd(absDiff2)  ;  summary(c(absDiff2))



# ------------------------------------ #
# Next produces Fig. 1 in main paper   #
# ------------------------------------ #

# Double-triple plot of sigma1, sigma2 and N2
par(mfrow = c(3, 2), mar = c(5, 6, 5, 2), cex.lab = 1.8, cex.axis = 1.8, cex.main = 1.8)

# IDS1: sigma1 (sigma in true DS data)
hist(100*estimates1[4, 1, ]-100, xlab = "Absolute error [m]", col = "grey", main = "sigma in DS data (DS + PC)", cex = 2, freq = TRUE, breaks = 20, xlim = c(-12, 12))
abline(v = 0, col = 'red', lwd = 3)
abline(v = mean(100*estimates1[4, 1, ]-100, na.rm = TRUE), col = 'blue', lty = 2, lwd = 3)

# IDS2: sigma1 (sigma in true DS data)
hist(100*estimates2[4, 1, ]-100, xlab = "Absolute error [m]", col = "grey", main = "sigma in DS data (DS + DND)", cex = 2, freq = TRUE, breaks = 20, xlim = c(-12, 12))
abline(v = 0, col = 'red', lwd = 3)
abline(v = mean(100*estimates2[4, 1, ]-100, na.rm = TRUE), col = 'blue', lty = 2, lwd = 3)

# IDS1: Estimated vs. true value of sigma2 (for the simple PC data)
plot(100*true.sigmaPC[1:Kmax], 100*estimates1[5,1,1:Kmax], xlab = "True value", ylab = "Estimated value", frame = FALSE, col = rgb(0,0,0,0.3), pch = 16, main = "sigma in PC data", cex = 1.5, xlim = c(0, 140), ylim = c(0, 140))
segments(100*true.sigmaPC[1:Kmax], 100*estimates1[5,3, 1:Kmax], 100*true.sigmaPC[1:Kmax], 100*estimates1[5,7, 1:Kmax])
abline(0, 1, col = 'red', lwd = 2)
abline(lm(I(100*estimates1[5,1,1:Kmax]) ~ I(100*true.sigmaPC[1:Kmax])), col = 'blue', lty = 2, lwd = 2)

# IDS2: Estimated vs. true value of sigma2 (for the simple DND data)
plot(100*true.sigmaDND[1:Kmax], 100*estimates2[5,1,1:Kmax], xlab = "True value", ylab = "Estimated value", frame = FALSE, col = rgb(0,0,0,0.3), pch = 16, main = "sigma in DND data", cex = 1.3, xlim = c(0, 140), ylim = c(0, 140))
segments(100*true.sigmaDND[1:Kmax], 100*estimates2[5,3,1:Kmax], 100*true.sigmaDND[1:Kmax], 100*estimates2[5,7, 1:Kmax])
abline(0, 1, col = 'red', lwd = 2)
abline(lm(I(100*estimates2[5,1,1:Kmax]) ~ I(100*true.sigmaDND[1:Kmax])), col = 'blue', lty = 2, lwd = 2)

# IDS1: Estimation error of site-level abundance N2 (sites with simple PC data)
hist(absDiff1, xlab = "Absolute error", col = "grey", main = "Site abundance in PC data", cex = 2, freq = TRUE, breaks = 100, xlim = c(-50, 50))
abline(v = 0, col = 'red', lwd = 3)
abline(v = mean(absDiff1), col = 'blue', lty = 2, lwd = 3)

# IDS2: Estimation error of site-level abundance N2 (sites with simple DND data)
hist(absDiff2, xlab = "Absolute error", col = "grey", main = "Site abundance in DND data", cex = 2, freq = TRUE, breaks = 100, xlim = c(-50, 50))
abline(v = 0, col = 'red', lwd = 3)
abline(v = mean(absDiff2), col = 'blue', lty = 2, lwd = 3)
