
# Supporting Information for:
# Kéry, M., Royle, J.A., Hallman, T., Robinson, W.D., Strebel, N., Kellner, K.F., 2024:
# Integrated distance sampling models for simple point counts, Ecology.


# Written by the authors, Nov 2022


# Simulation 3: How many DS sites are required to obtain adequate estimates of density?
# -------------------------------------------------------------------------------------
# We conduct an "admixture" simulation study, where we add varying numbers of sites 
# with the more information-rich HDS data to a larger and fixed number of sites with
# simple point counts. 

# Load package
# Needs a version of unmarked with the IDS() function; see README file in this folder
library(unmarked)

# Setup
library(AHMbook)


# Number of added distance-sampling point counts
nlevels <- 6     # Number of levels of admixture factor
Nadmix <- c(1, 20, 40, 60, 80 , 100)

# Number of simulation reps at each level of the admix factor
simrep <- 100           # For test run
simrep <- 1000         # For real

# Create R object to save the true parameter values
true.vals <- array(NA, dim = c(4, simrep))
rownames(true.vals) <- c("mean_lam", "beta_lam", "log sigmaDS", "log sigmaPC")

# Create R object to save the estimates
esti1 <- array(NA, dim = c(4, 4, nlevels, simrep))
dimnames(esti1) <- list(c("log mean_lam", "beta_lam", "log sigmaDS", 
  "log sigmaPC"), c("MLE", "SE", "LCL", "UCL"), Nadmix, NULL)


# Choose param values in state and in detection all at once
# Shared stuff
mean.lambda <- 100 # For B = 500, should yield density of 1 per ha
B <- 500

# PC data
nsites_pc <- 200
sigPC <- 70        # Note this is in meters now

# HDS data
#nsites_hds <- 20  # This is variable now ! see below
sigHDS <- 100      # meters !
J <- 4             # Number of distance bins
maxDist_hds <- 200 # Maximum distance of observation in the HDS protocol,
                   # Also this is in meters
db <- seq(0, maxDist_hds, length.out=(J+1)) # Distance bins



# Launch the sim
system.time(
for(v in 1:nlevels){  # Loop over all 6 levels of the admixture factor
  for(k in 1:simrep){        # Loop over all simreps
  cat("\n\n*** Simulation Number", k, "***\n")
  cat("\n* Admixture", Nadmix[v], "*\n")

  # Determine the number of HDS sites admixed
  nsites_hds <- Nadmix[v]


  ## Simulate Point Count data set
  dat2 <- simHDS(type="point", nsites = nsites_pc, 
    mean.lambda = mean.lambda, beta.lam = 1, mean.sigma = sigPC, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)
  # Re-format to put in a unmarkedFramePCount object
  y_pc <- matrix(0, nrow=nsites_pc, ncol=1)
  for (i in 1:nsites_pc){
    dsub <- dat2$data[dat2$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    y_pc[i,1] <- length(dsub)
  }
  # Site covariates
  sc_pc <- data.frame(habitat=dat2$habitat)
  # Create unmarked frame for PC data set
  umf_pc <- unmarkedFramePCount(y=y_pc, siteCovs=sc_pc)

  ## Simulate a regular distance sampling data set
  dat1 <- simHDS(type="point", nsites = nsites_hds, 
    mean.lambda = mean.lambda, beta.lam = 1, mean.sigma = sigHDS, 
    beta.sig = 0, B = B, discard0 = FALSE, show.plot = FALSE)
  # Re-format to put into an `unmarkedFrameDS` object:
  # Convert long-form to y matrix
  y_hds <- matrix(0, nrow=nsites_hds, ncol=J)
  for (i in 1:nsites_hds){
    dsub <- dat1$data[dat1$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    for (j in 1:J){
      y_hds[i,j] <- sum(dsub >= db[j] & dsub < db[j+1])
    }
  }
  # Site covariates
  sc_hds <- data.frame(habitat=dat1$habitat)
  # Create unmarked frame
  umf_hds <- unmarkedFrameDS(y=y_hds, siteCovs=sc_hds, survey="point",
              dist.breaks=db, unitsIn="m")

  # Fit the IDS with HDS + PC data
  # We add an "error-catch-saver" using try()
  out <- try(IDS(lambdaformula = ~habitat, 
    detformulaDS = ~1, detformulaPC = ~1, 
    dataDS = umf_hds, dataPC = umf_pc, maxDistPC = 500,
    K = 1000, control = list(trace=1, REPORT = 10)) )
  if(class(out) == "try-error") {fuck <- TRUE}
  else{
  tmp <- summary(out)
  ses <- c(tmp$lam$SE, tmp$ds$SE, tmp$pc$SE) 


  # Extract and save param estimates
  esti1[, 1, v, k] <- coef(out)  # Save MLEs
  esti1[, 2, v, k] <- ses        # Save SEs
  esti1[1:2, 3:4, v, k] <- confint(out, type = 'lam') # CIs
  esti1[3, 3:4, v, k] <- confint(out, type = 'ds')
  esti1[4, 3:4, v, k] <- confint(out, type = 'pc')
  }
}
}
) # End of system.time()



# Weed out numerical failures
# Flag cases with bad SEs: either NaN or large (>5) SEs and tally up
tempo1 <- apply(esti1[,2,,], c(2,3), function(x) as.numeric(sum(x == 'NaN')>0))
tempo2 <- apply(esti1[,2,,], c(2,3), function(x) sum(x > 5, na.rm = TRUE))
bad.SE <- (tempo1 + tempo2) > 0
apply(bad.SE, 1, mean, na.rm = TRUE)

# Freakish high point estimates 
# Flag freak point estimates and then tally them up
tempo1 <- exp(esti1[1,1,,]) > 10
tempo3 <- exp(esti1[3,1,,]) > 1000
tempo4 <- exp(esti1[4,1,,]) > 1000
freaks <- (tempo1 + tempo3 + tempo4) > 0
apply(freaks, 1, mean, na.rm = TRUE)

# Proportion of any badness as a function of admixture
really.bad <- (freaks + bad.SE) > 0
apply(really.bad, 1, mean, na.rm = TRUE)

# Set to NA all 'really bad' cases
for(k in 1:6){
#  for(l in 1:1000){
  for(l in 1:100){
    esti1[,,k,l][really.bad[k,l]] <- NA
  }
}



# --------------------------------- #
# Next produces Figure 3 in paper   #
# --------------------------------- #


# lambda intercept on 'real' scale
par(mfrow = c(2,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5,
cex.main = 1.5)
boxplot(exp(esti1[1,1,,]) ~ Nadmix[1:6], outline = TRUE, frame = FALSE, col = 'grey', main = 'lambda intercept', xlab = 'Number of DS sites added', ylab = '', las = 1, ylim = c(0, 8), na.action = na.pass)
abline(h = 1, col = 'red', lwd = 3, lty = 2)

# lambda slope on link scale
boxplot(esti1[2,1,,] ~ Nadmix, outline = TRUE, frame = FALSE, col = 'grey', xlab = 'Number of DS sites added', main = 'lambda slope', ylab = "", las = 1, ylim = c(0.9, 1.1), na.action = na.pass)
abline(h = 1, col = 'red', lwd = 3, lty = 2)

# Sigma of HN detection function for HDS data
boxplot(exp(esti1[3,1,,]) ~ Nadmix, outline = TRUE, frame = FALSE, col = 'grey', xlab = 'Number of DS sites added', main= 'sigma (DS)', ylab = "", las = 1, ylim = c(0, 240), na.action = na.pass)
abline(h = 100, col = 'red', lwd = 3, lty = 2)

# Sigma of HN detection function for PC data
boxplot(exp(esti1[4,1,,]) ~ Nadmix, outline = TRUE, frame = FALSE, col = 'grey', xlab = 'Number of DS sites added', main = 'sigma (PC)', ylab = "", las = 1, ylim = c(0, 150), na.action = na.pass)
abline(h = 70, col = 'red', lwd = 3, lty = 2)
