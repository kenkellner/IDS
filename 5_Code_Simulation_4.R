
# Supporting Information for:
# Kéry, M., Royle, J.A., Hallman, T., Robinson, W.D., Strebel, N., Kellner, K.F., 2024:
# Integrated distance sampling models for simple point counts, Ecology.


# Written by the authors, Nov 2022



# Simulation 4: How well can availability be estimated in an IDS model?
# ---------------------------------------------------------------------

# In this simulation we use the outcomes from running function simHDS() to also simulate
# an availability process in the data and then see how well we can recover estimates of it
# in an IDS model that accommodates availability in addition to perceptibility.
# Note that the information about the former will come from variable survey duration in
# the PC data, while most of the information about the latter will come from the actual 
# distance sampling data.
# We will run 3 variants of this simulation: with PC data from 1000, 3000, and from 6000 sites.
#

# Load package
# Needs a version of unmarked with the IDS() function; see README file in this folder
library(unmarked)

# Setup
library(AHMbook)

# We will simulate sample sizes that resemble somewhat those in the Robin case study 
# in the paper, i.e., about 3000 sites with HDS surveys and about 1000 sites with PC.
# Then, below we run two variants then where instead we have 3000 or 6000 sites with PC data.
# We choose a similar distribution of survey durations as in the Oregon data.

# Try and draw survey durations in the simulated PC data from such a much more skewed distribution and see how good our estimates then become.

# Generate a skewed distribution for survey duration in the PC data
xxx <- round(exp(rnorm(10000, mean = 1, sd = 1)), 0)
xxx[xxx > 30] <- 30
xxx <- xxx + 3
plot(table(xxx), type = 'h', lend = 'butt', col = 'grey', lwd = 10, xlab = 'Original duration data', frame = FALSE)

# Then, we squeeze it into the range of 2 – 5.

# Squeeze between 2 and 5
old_range <- max(xxx) - min(xxx)
new_range <- 5 - 2
new_xxx <- (((xxx - min(xxx)) * new_range) / old_range) + 2
range(new_xxx)
plot(table(new_xxx), type = 'h', lend = 'butt', col = 'grey', lwd = 10, xlab = 'Squeezed duration data', frame = FALSE)

# Below, we will use this code to generate such a skewed distribution for the survey durations.



# Simulation 4a: 1000 sites with PC data with variable duration
# -------------------------------------------------------------

# For the HDS data we keep their duration at 1/6 of the maximum duration, i.e., at 5 / 6.

# Number of simulation reps
simrep <- 100               # For testing
simrep <- 1000              # For real

# Create new R object to save true values of phi
true.phiA <- numeric(simrep)

# Create new R object to save the estimates
estiA <- array(NA, dim = c(5, 4, simrep))
dimnames(estiA) <- list(c("log mean_lam", "beta_lam", "log sigmaDS", 
  "log sigmaPC", "log phi"), c("MLE", "SE", "LCL", "UCL"), NULL)


# Choose param values in state and in detection all at once
# Shared stuff
mean.lambda <- 100 # For B = 500, should yield density of 1 per ha
beta.lam <- 1
B <- 500

# PC data
nsites_pc <- 1000  # Start with only 1k sites as in the Oregon Case study
sigPC <- 70        # Note this is in meters now

# HDS data
nsites_hds <- 3000 # But now 3k HDS sites
sigHDS <- 100      # meters !
J <- 4             # Number of distance bins
maxDist_hds <- 200 # Maximum distance of observation in the HDS protocol,
                   # Also this is in meters
db <- seq(0, maxDist_hds, length.out=(J+1)) # Distance bins




# Launch the sim
system.time(
for(k in 1:simrep){        # Loop over all simreps
  cat("\n\n*** Simulation Number", k, "***\n")

  ## Draw some value of phi from -2.3 to 0.7
  # Corresponds to rates of about 0.1 to 2
  phi <- runif(1, -2.3, 0.7)
  true.phiA[k] <- phi

  ## Simulate Point Count data set (raw, without availability so far)
  dat2.raw <- simHDS(type="point", nsites = nsites_pc, 
    mean.lambda = mean.lambda, beta.lam = beta.lam, mean.sigma = sigPC, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)

  # Re-format to put in a unmarkedFramePCount object
  y_pc <- matrix(0, nrow=nsites_pc, ncol=1)
  for (i in 1:nsites_pc){
    dsub <- dat2.raw$data[dat2.raw$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    y_pc[i,1] <- length(dsub)
  }
  # Site covariates
  sc_pc <- data.frame(habitat=dat2.raw$habitat)
  # Create unmarked frame for PC data set
  umf_pc <- unmarkedFramePCount(y=y_pc, siteCovs=sc_pc)

  ## Simulate a regular distance sampling data set
  dat1.raw <- simHDS(type="point", nsites = nsites_hds, 
    mean.lambda = mean.lambda, beta.lam = beta.lam, mean.sigma = sigHDS, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)
  # Re-format to put into an `unmarkedFrameDS` object:
  # Convert long-form to y matrix
  y_hds <- matrix(0, nrow=nsites_hds, ncol=J)
  for (i in 1:nsites_hds){
    dsub <- dat1.raw$data[dat1.raw$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    for (j in 1:J){
      y_hds[i,j] <- sum(dsub >= db[j] & dsub < db[j+1])
    }
  }
  # Site covariates
  sc_hds <- data.frame(habitat=dat1.raw$habitat)
  # Create unmarked frame
  umf_hds <- unmarkedFrameDS(y=y_hds, siteCovs=sc_hds, survey="point",
              dist.breaks=db, unitsIn="m")

  # Simulate availability process on existing dataset
  # HDS: constant duration, and therefore, constant availability
  dds <- rep(5/6, nsites_hds)
  pdds <- 1-exp(-1*dds*exp(phi))
  umf_hds2 <- umf_hds
  for (i in 1:nrow(umf_hds2@y)){
    for (j in 1:ncol(umf_hds2@y)){
      umf_hds2@y[i,j] <- rbinom(1, umf_hds2@y[i,j], pdds[i])
    }
  }


  # PC data
  # here now make duration skewed as shown above
  dpc_temp <- round(exp(rnorm(nsites_pc, mean = 1, sd = 1)), 0)
  dpc_temp [dpc_temp > 30] <- 30
  dpc_temp <- dpc_temp + 3
  old_range <- max(dpc_temp) - min(dpc_temp)
  new_range <- 5 - 2
  dpc <- (((dpc_temp - min(dpc_temp)) * new_range) / old_range) + 2
  pdpc <- 1-exp(-1*dpc*exp(phi))
  umf_pc2 <- umf_pc
  for (i in 1:nrow(umf_pc2@y)){
    for (j in 1:ncol(umf_pc2@y)){
      umf_pc2@y[i,j] <- rbinom(1, umf_pc2@y[i,j], pdpc[i])
    }
  }

  # Fit the IDS with HDS + PC data

  # Availability formula argument is availformula
  # Duration data provided to durationDS, durationPC, durationOC
  # respectively
  # must each be vectors of length equal to the number of sites 
  # for each data type
  # We add an "error-catch-saver" using try()
  out <- try(IDS(lambdaformula = ~habitat, 
    detformulaDS = ~1, detformulaPC = ~1, 
    availformula=~1, durationDS=dds, durationPC=dpc,
    dataDS = umf_hds2, dataPC = umf_pc2, maxDistPC = 500,
    K = 1000, control = list(trace=1, REPORT = 10)) )
  if(class(out) == "try-error") {damn <- TRUE}
  else{
  tmp <- summary(out)
  ses <- c(tmp$lam$SE, tmp$ds$SE, tmp$pc$SE, tmp$phi$SE) 
  cat("* True phi (log-scale)", phi, "*\n")

  # Extract and save param estimates
  estiA[, 1, k] <- coef(out)  # Save MLEs
  estiA[, 2, k] <- ses        # Save SEs
  estiA[1:2, 3:4,  k] <- confint(out, type = 'lam') # CIs
  estiA[3, 3:4, k] <- confint(out, type = 'ds')
  estiA[4, 3:4, k] <- confint(out, type = 'pc')
  estiA[5, 3:4, k] <- confint(out, type = 'phi')
  }
}
) # End of system.time()



# Weed out numerical failures 
# Check for 'NaN's for standard errors

# Count number of estimates with NaN SEs per model fit
tmp <- apply(estiA[,2,], 2, function(x) sum(x == 'NaN'))
table(tmp)    # Count of cases

# Count number model fits with any NaN standard errors
tmp <- apply(estiA[,2,], 2, function(x) sum(x == 'NaN')>0)
table(tmp)    # Count of cases

# Large SEs > 5
tmp <- apply(estiA[,2,], 2, function(x) sum(x > 5))
table(tmp)

# Freakish high point estimates (10 times larger than truth)
# Flag freak point estimates and then tally them up
tempo1 <- exp(estiA[1,1,]) > 10
tempo3 <- exp(estiA[3,1,]) > 1000
tempo4 <- exp(estiA[4,1,]) > 1000
tempo5 <- exp(estiA[5,1,]) > 10
freaks <- (tempo1 + tempo3 + tempo4 + tempo5) > 0
sum(freaks)

# Flag cases with bad SEs: either NaN or large (>5) SEs and tally up
tempo1 <- apply(estiA[,2,], 2, function(x) as.numeric(sum(x == 'NaN')>0))
tempo2 <- apply(estiA[,2,], 2, function(x) sum(x > 5, na.rm = TRUE))
bad.SE <- (tempo1 + tempo2) > 0
sum(bad.SE)

# Combine the two exclusion criteria
really.bad <- (freaks + bad.SE) > 0
sum(really.bad)

# Make a copy and set to NA all 'really bad' cases
estiAA <- estiA
estiAA[,,really.bad] <- NA

# Make a plot (which, however, isn't yet the one in the paper)
par(mfrow = c(2,2), mar = c(6,5,6,2), cex.lab = 1.5, cex.axis = 1.5,
cex.main = 1.6)
# lambda intercept on natural scale
hist(exp(estiAA[1,1,]), col = 'grey', main = 'lambda intercept', xlab = 'MLEs', breaks = 100, freq = F, ylab = '', las = 1, xlim = c(0, 10))
abline(v = 1, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiAA[1,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Intercept of availability (phi)
plot(exp(true.phiA), exp(estiAA[5,1,]), main = 'availability (phi)', xlab = 'True phi', ylab = 'MLE of phi', pch = 16, col = rgb(0,0,0,0.4), frame = F, xlim = c(0.1, 2), ylim = c(0, 3.6))
abline(0, 1, col = 'red', lwd = 3, lty = 1)
abline(lm(exp(estiAA[5,1,]) ~ exp(true.phiA)), col = 'blue', lwd = 3, lty = 2)

# Sigma of HN detection function for HDS data
hist(exp(estiAA[3,1,]), col = 'grey', main = 'sigma (DS)', xlab = 'MLEs', breaks = 40, freq = F, ylab = '', las = 1)
abline(v = 100, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiAA[3,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Sigma of HN detection function for PC data
hist(exp(estiAA[4,1,]), col = 'grey', main = 'sigma (PC)', xlab = 'MLEs', breaks = 30, freq = F, ylab = '', las = 1)
abline(v = 70, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiAA[4,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)







# Simulation 4b: 3000 sites with PC data with variable duration
# -------------------------------------------------------------

# Number of simulation reps
simrep <- 100               # For testing
simrep <- 1000              # For real

# Create new R object to save true values of phi
true.phiB <- numeric(simrep)

# Create new R object to save the estimates
estiB <- array(NA, dim = c(5, 4, simrep))
dimnames(estiB) <- list(c("log mean_lam", "beta_lam", "log sigmaDS", 
  "log sigmaPC", "log phi"), c("MLE", "SE", "LCL", "UCL"), NULL)

# Choose param values in state and in detection all at once
# ALL SAME AS ABOVE, EXCEPT FOR NUMBER OF PC SITES
nsites_pc <- 3000  # Now 3k PC sites

# Launch the sim
system.time(
for(k in 1:simrep){        # Loop over all simreps
  cat("\n\n*** Simulation Number", k, "***\n")

  ## Draw some value of phi from -2.3 to 0.7
  # Corresponds to rates of about 0.1 to 2
  phi <- runif(1, -2.3, 0.7)
  true.phiB[k] <- phi

  ## Simulate Point Count data set (raw, without availability so far)
  dat2.raw <- simHDS(type="point", nsites = nsites_pc, 
    mean.lambda = mean.lambda, beta.lam = beta.lam, mean.sigma = sigPC, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)

  # Re-format to put in a unmarkedFramePCount object
  y_pc <- matrix(0, nrow=nsites_pc, ncol=1)
  for (i in 1:nsites_pc){
    dsub <- dat2.raw$data[dat2.raw$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    y_pc[i,1] <- length(dsub)
  }
  # Site covariates
  sc_pc <- data.frame(habitat=dat2.raw$habitat)
  # Create unmarked frame for PC data set
  umf_pc <- unmarkedFramePCount(y=y_pc, siteCovs=sc_pc)

  ## Simulate a regular distance sampling data set
  dat1.raw <- simHDS(type="point", nsites = nsites_hds, 
    mean.lambda = mean.lambda, beta.lam = beta.lam, mean.sigma = sigHDS, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)
  # Re-format to put into an `unmarkedFrameDS` object:
  # Convert long-form to y matrix
  y_hds <- matrix(0, nrow=nsites_hds, ncol=J)
  for (i in 1:nsites_hds){
    dsub <- dat1.raw$data[dat1.raw$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    for (j in 1:J){
      y_hds[i,j] <- sum(dsub >= db[j] & dsub < db[j+1])
    }
  }
  # Site covariates
  sc_hds <- data.frame(habitat=dat1.raw$habitat)
  # Create unmarked frame
  umf_hds <- unmarkedFrameDS(y=y_hds, siteCovs=sc_hds, survey="point",
              dist.breaks=db, unitsIn="m")

  # Simulate availability process on existing dataset
  # HDS: constant duration, and therefore, constant availability
  dds <- rep(5/6, nsites_hds)
  pdds <- 1-exp(-1*dds*exp(phi))
  umf_hds2 <- umf_hds
  for (i in 1:nrow(umf_hds2@y)){
    for (j in 1:ncol(umf_hds2@y)){
      umf_hds2@y[i,j] <- rbinom(1, umf_hds2@y[i,j], pdds[i])
    }
  }

  # PC data
  # here now make duration skewed as shown above
  dpc_temp <- round(exp(rnorm(nsites_pc, mean = 1, sd = 1)), 0)
  dpc_temp [dpc_temp > 30] <- 30
  dpc_temp <- dpc_temp + 3
  old_range <- max(dpc_temp) - min(dpc_temp)
  new_range <- 5 - 2
  dpc <- (((dpc_temp - min(dpc_temp)) * new_range) / old_range) + 2
  pdpc <- 1-exp(-1*dpc*exp(phi))
  umf_pc2 <- umf_pc
  for (i in 1:nrow(umf_pc2@y)){
    for (j in 1:ncol(umf_pc2@y)){
      umf_pc2@y[i,j] <- rbinom(1, umf_pc2@y[i,j], pdpc[i])
    }
  }

  # Fit the IDS with HDS + PC data

  # Availability formula argument is availformula
  # Duration data provided to durationDS, durationPC, durationOC
  # respectively
  # must each be vectors of length equal to the number of sites 
  # for each data type
  # We add an "error-catch-saver" using try()
  out <- try(IDS(lambdaformula = ~habitat, 
    detformulaDS = ~1, detformulaPC = ~1, 
    availformula=~1, durationDS=dds, durationPC=dpc,
    dataDS = umf_hds2, dataPC = umf_pc2, maxDistPC = 500,
    K = 1000, control = list(trace=1, REPORT = 10)) )
  if(class(out) == "try-error") {damn <- TRUE}
  else{
  tmp <- summary(out)
  ses <- c(tmp$lam$SE, tmp$ds$SE, tmp$pc$SE, tmp$phi$SE) 
  cat("* True phi (log-scale)", phi, "*\n")

  # Extract and save param estimates
  estiB[, 1, k] <- coef(out)  # Save MLEs
  estiB[, 2, k] <- ses        # Save SEs
  estiB[1:2, 3:4,  k] <- confint(out, type = 'lam') # CIs
  estiB[3, 3:4, k] <- confint(out, type = 'ds')
  estiB[4, 3:4, k] <- confint(out, type = 'pc')
  estiB[5, 3:4, k] <- confint(out, type = 'phi')
  }
}
) # End of system.time()


# Weed out numerical failures 
# 'NaN's for standard errors ?

# Count number of estimates with NaN SEs per model fit
tmp <- apply(estiB[,2,], 2, function(x) sum(x == 'NaN'))
table(tmp)    # Count of cases

# Count number model fits with any NaN standard errors
tmp <- apply(estiB[,2,], 2, function(x) sum(x == 'NaN')>0)
table(tmp)    # Count of cases

# Large SEs > 5
tmp <- apply(estiB[,2,], 2, function(x) sum(x > 5))
table(tmp)

# Freakish high point estimates (10 times larger than truth)
# Flag freak point estimates and then tally them up
tempo1 <- exp(estiB[1,1,]) > 10
tempo3 <- exp(estiB[3,1,]) > 1000
tempo4 <- exp(estiB[4,1,]) > 1000
tempo5 <- exp(estiB[5,1,]) > 10
freaks <- (tempo1 + tempo3 + tempo4 + tempo5) > 0
sum(freaks)

# Flag cases with bad SEs: either NaN or large (>5) SEs and tally up
tempo1 <- apply(estiB[,2,], 2, function(x) as.numeric(sum(x == 'NaN')>0))
tempo2 <- apply(estiB[,2,], 2, function(x) sum(x > 5, na.rm = TRUE))
bad.SE <- (tempo1 + tempo2) > 0
sum(bad.SE)

# Combine the two exclusion criteria
really.bad <- (freaks + bad.SE) > 0
sum(really.bad)

# Make a copy and set to NA all 'really bad' cases
estiBB <- estiB
estiBB[,,really.bad] <- NA

# Make a plot (which, however, isn't yet the one in the paper)
par(mfrow = c(2,2), mar = c(6,5,6,2), cex.lab = 1.5, cex.axis = 1.5,
cex.main = 1.6)
# lambda intercept on natural scale
hist(exp(estiBB[1,1,]), col = 'grey', main = 'lambda intercept', xlab = 'MLEs', breaks = 100, freq = F, ylab = '', las = 1, xlim = c(0, 10))
abline(v = 1, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiBB[1,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Intercept of availability (phi)
plot(exp(true.phiB), exp(estiBB[5,1,]), main = 'availability (phi)', xlab = 'True phi', ylab = 'MLE of phi', pch = 16, col = rgb(0,0,0,0.4), frame = F, xlim = c(0.1, 2), ylim = c(0, 3.6))
abline(0, 1, col = 'red', lwd = 3, lty = 1)
abline(lm(exp(estiBB[5,1,]) ~ exp(true.phiB)), col = 'blue', lwd = 3, lty = 2)

# Sigma of HN detection function for HDS data
hist(exp(estiBB[3,1,]), col = 'grey', main = 'sigma (DS)', xlab = 'MLEs', breaks = 40, freq = F, ylab = '', las = 1)
abline(v = 100, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiBB[3,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Sigma of HN detection function for PC data
hist(exp(estiBB[4,1,]), col = 'grey', main = 'sigma (PC)', xlab = 'MLEs', breaks = 30, freq = F, ylab = '', las = 1)
abline(v = 70, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiBB[4,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)






# Simulation 4c: 6000 sites with PC data with variable duration
# -------------------------------------------------------------

# Number of simulation reps
simrep <- 100               # For testing
simrep <- 1000              # For real

# Create new R object to save true values of phi
true.phiC <- numeric(simrep)

# Create new R object to save the estimates
estiC <- array(NA, dim = c(5, 4, simrep))
dimnames(estiC) <- list(c("log mean_lam", "beta_lam", "log sigmaDS", 
  "log sigmaPC", "log phi"), c("MLE", "SE", "LCL", "UCL"), NULL)

# Choose param values in state and in detection all at once
# ALL SAME AS ABOVE, EXCEPT FOR NUMBER OF PC SITES
nsites_pc <- 6000  # Now 6k PC sites

# Launch the sim
system.time(
for(k in 1:simrep){        # Loop over all simreps
  cat("\n\n*** Simulation Number", k, "***\n")

  ## Draw some value of phi from -2.3 to 0.7
  # Corresponds to rates of about 0.1 to 2
  phi <- runif(1, -2.3, 0.7)
  true.phiC[k] <- phi

  ## Simulate Point Count data set (raw, without availability so far)
  dat2.raw <- simHDS(type="point", nsites = nsites_pc, 
    mean.lambda = mean.lambda, beta.lam = beta.lam, mean.sigma = sigPC, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)

  # Re-format to put in a unmarkedFramePCount object
  y_pc <- matrix(0, nrow=nsites_pc, ncol=1)
  for (i in 1:nsites_pc){
    dsub <- dat2.raw$data[dat2.raw$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    y_pc[i,1] <- length(dsub)
  }
  # Site covariates
  sc_pc <- data.frame(habitat=dat2.raw$habitat)
  # Create unmarked frame for PC data set
  umf_pc <- unmarkedFramePCount(y=y_pc, siteCovs=sc_pc)

  ## Simulate a regular distance sampling data set
  dat1.raw <- simHDS(type="point", nsites = nsites_hds, 
    mean.lambda = mean.lambda, beta.lam = beta.lam, mean.sigma = sigHDS, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)
  # Re-format to put into an `unmarkedFrameDS` object:
  # Convert long-form to y matrix
  y_hds <- matrix(0, nrow=nsites_hds, ncol=J)
  for (i in 1:nsites_hds){
    dsub <- dat1.raw$data[dat1.raw$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    for (j in 1:J){
      y_hds[i,j] <- sum(dsub >= db[j] & dsub < db[j+1])
    }
  }
  # Site covariates
  sc_hds <- data.frame(habitat=dat1.raw$habitat)
  # Create unmarked frame
  umf_hds <- unmarkedFrameDS(y=y_hds, siteCovs=sc_hds, survey="point",
              dist.breaks=db, unitsIn="m")

  # Simulate availability process on existing dataset
  # HDS: constant duration, and therefore, constant availability
  dds <- rep(5/6, nsites_hds)
  pdds <- 1-exp(-1*dds*exp(phi))
  umf_hds2 <- umf_hds
  for (i in 1:nrow(umf_hds2@y)){
    for (j in 1:ncol(umf_hds2@y)){
      umf_hds2@y[i,j] <- rbinom(1, umf_hds2@y[i,j], pdds[i])
    }
  }

  # PC data
  # here now make duration skewed as shown above
  dpc_temp <- round(exp(rnorm(nsites_pc, mean = 1, sd = 1)), 0)
  dpc_temp [dpc_temp > 30] <- 30
  dpc_temp <- dpc_temp + 3
  old_range <- max(dpc_temp) - min(dpc_temp)
  new_range <- 5 - 2
  dpc <- (((dpc_temp - min(dpc_temp)) * new_range) / old_range) + 2
  pdpc <- 1-exp(-1*dpc*exp(phi))
  umf_pc2 <- umf_pc
  for (i in 1:nrow(umf_pc2@y)){
    for (j in 1:ncol(umf_pc2@y)){
      umf_pc2@y[i,j] <- rbinom(1, umf_pc2@y[i,j], pdpc[i])
    }
  }

  # Fit the IDS with HDS + PC data

  # Availability formula argument is availformula
  # Duration data provided to durationDS, durationPC, durationOC
  # respectively
  # must each be vectors of length equal to the number of sites 
  # for each data type
  # We add an "error-catch-saver" using try()
  out <- try(IDS(lambdaformula = ~habitat, 
    detformulaDS = ~1, detformulaPC = ~1, 
    availformula=~1, durationDS=dds, durationPC=dpc,
    dataDS = umf_hds2, dataPC = umf_pc2, maxDistPC = 500,
    K = 1000, control = list(trace=1, REPORT = 10)) )
  if(class(out) == "try-error") {damn <- TRUE}
  else{
  tmp <- summary(out)
  ses <- c(tmp$lam$SE, tmp$ds$SE, tmp$pc$SE, tmp$phi$SE) 
  cat("* True phi (log-scale)", phi, "*\n")

  # Extract and save param estimates
  estiC[, 1, k] <- coef(out)  # Save MLEs
  estiC[, 2, k] <- ses        # Save SEs
  estiC[1:2, 3:4,  k] <- confint(out, type = 'lam') # CIs
  estiC[3, 3:4, k] <- confint(out, type = 'ds')
  estiC[4, 3:4, k] <- confint(out, type = 'pc')
  estiC[5, 3:4, k] <- confint(out, type = 'phi')
  }
}
) # End of system.time()



# Weed out numerical failures 
# 'NaN's for standard errors ?

# Count number of estimates with NaN SEs per model fit
tmp <- apply(estiC[,2,], 2, function(x) sum(x == 'NaN'))
table(tmp)    # Count of cases

# Count number model fits with any NaN standard errors
tmp <- apply(estiC[,2,], 2, function(x) sum(x == 'NaN')>0)
table(tmp)    # Count of cases

# Large SEs > 5
tmp <- apply(estiC[,2,], 2, function(x) sum(x > 5))
table(tmp)

# Freakish high point estimates (10 times larger than truth)
# Flag freak point estimates and then tally them up
tempo1 <- exp(estiC[1,1,]) > 10
tempo3 <- exp(estiC[3,1,]) > 1000
tempo4 <- exp(estiC[4,1,]) > 1000
tempo5 <- exp(estiC[5,1,]) > 10
freaks <- (tempo1 + tempo3 + tempo4 + tempo5) > 0
sum(freaks)

# Flag cases with bad SEs: either NaN or large (>5) SEs and tally up
tempo1 <- apply(estiC[,2,], 2, function(x) as.numeric(sum(x == 'NaN')>0))
tempo2 <- apply(estiC[,2,], 2, function(x) sum(x > 5, na.rm = TRUE))
bad.SE <- (tempo1 + tempo2) > 0
sum(bad.SE)

# Combine the two exclusion criteria
really.bad <- (freaks + bad.SE) > 0
sum(really.bad)

# Make a copy and set to NA all 'really bad' cases
estiCC <- estiC
estiCC[,,really.bad] <- NA

# Make a plot (which, however, isn't yet the one in the paper)
par(mfrow = c(2,2), mar = c(6,5,6,2), cex.lab = 1.5, cex.axis = 1.5,
cex.main = 1.6)
# lambda intercept on natural scale
hist(exp(estiCC[1,1,]), col = 'grey', main = 'lambda intercept', xlab = 'MLEs', breaks = 100, freq = F, ylab = '', las = 1, xlim = c(0, 10))
abline(v = 1, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiCC[1,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Intercept of availability (phi)
plot(exp(true.phiC), exp(estiCC[5,1,]), main = 'availability (phi)', xlab = 'True phi', ylab = 'MLE of phi', pch = 16, col = rgb(0,0,0,0.4), frame = F, xlim = c(0.1, 2), ylim = c(0, 3.6))
abline(0, 1, col = 'red', lwd = 3, lty = 1)
abline(lm(exp(estiCC[5,1,]) ~ exp(true.phiC)), col = 'blue', lwd = 3, lty = 2)

# Sigma of HN detection function for HDS data
hist(exp(estiCC[3,1,]), col = 'grey', main = 'sigma (DS)', xlab = 'MLEs', breaks = 40, freq = F, ylab = '', las = 1)
abline(v = 100, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiCC[3,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Sigma of HN detection function for PC data
hist(exp(estiCC[4,1,]), col = 'grey', main = 'sigma (PC)', xlab = 'MLEs', breaks = 30, freq = F, ylab = '', las = 1)
abline(v = 70, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiCC[4,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)





# --------------------------------- #
# Next produces Figure 4 in paper   #
# --------------------------------- #

# Note use estiAA, estiBB, and estiCC to toss out numerical failures

# Plots
par(mfrow = c(2,3), mar = c(6,6,4,3), cex.lab = 1.8, cex.axis = 1.8,
cex.main = 2.2)
xlim <- c(0, 2.5)
# lambda intercept on natural scale: 1000 PC sites
hist(exp(estiAA[1,1,]), col = 'grey', main = '1000 PC sites', xlab = 'lambda', breaks = 100, freq = F, ylab = '', las = 1, xlim = xlim)
abline(v = 1, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiAA[1,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# lambda intercept on natural scale: 3000 PC sites
hist(exp(estiBB[1,1,]), col = 'grey', main = '3000 PC sites', xlab = 'lambda', breaks = 100, freq = F, ylab = '', las = 1, xlim = xlim)
abline(v = 1, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiBB[1,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# lambda intercept on natural scale: 6000 PC sites
hist(exp(estiCC[1,1,]), col = 'grey', main = '6000 PC sites', xlab = 'lambda', breaks = 100, freq = F, ylab = '', las = 1, xlim = xlim)
abline(v = 1, col = 'red', lwd = 3, lty = 1)
abline(v = mean(exp(estiCC[1,1,]), na.rm = T), col = 'blue', lwd = 3, lty = 2)

# Intercept of availability (phi): 1000 PC sites
plot(exp(true.phiA), exp(estiAA[5,1,]), main = '', xlab = 'True phi', ylab = 'MLE of phi', pch = 16, col = rgb(0,0,0,0.4), frame = F, xlim = c(0.1, 2), ylim = c(0, 4))
abline(0, 1, col = 'red', lwd = 3, lty = 1)
abline(lm(exp(estiAA[5,1,]) ~ exp(true.phiA)), col = 'blue', lwd = 3, lty = 2)

# Intercept of availability (phi): 3000 PC sites
plot(exp(true.phiB), exp(estiBB[5,1,]), main = '', xlab = 'True phi', ylab = 'MLE of phi', pch = 16, col = rgb(0,0,0,0.4), frame = F, xlim = c(0.1, 2), ylim = c(0, 4))
abline(0, 1, col = 'red', lwd = 3, lty = 1)
abline(lm(exp(estiBB[5,1,]) ~ exp(true.phiB)), col = 'blue', lwd = 3, lty = 2)

# Intercept of availability (phi): 6000 PC sites
plot(exp(true.phiC), exp(estiCC[5,1,]), main = '', xlab = 'True phi', ylab = 'MLE of phi', pch = 16, col = rgb(0,0,0,0.4), frame = F, xlim = c(0.1, 2), ylim = c(0, 4))
abline(0, 1, col = 'red', lwd = 3, lty = 1)
abline(lm(exp(estiCC[5,1,]) ~ exp(true.phiC)), col = 'blue', lwd = 3, lty = 2)
