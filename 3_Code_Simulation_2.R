
# Supporting Information for:
# Kéry, M., Royle, J.A., Hallman, T., Robinson, W.D., Strebel, N., Kellner, K.F., 2024:
# Integrated distance sampling models for simple point counts, Ecology.


# Written by the authors, Nov 2022


# Simulation 2: Identifiability of IDS1 with distinct covariate effects in the observation model
# -----------------------------------------------------------------------------------------------
# Here we conduct two sets of simulations to answer the following related questions:
# (i) Does the IDS model allow DS and PC detection to have different covariate relationships?
# (ii) Are relationships still identifiable if the same covariates are related to both detection
# and density?
# We have two very similar simulations. In both we simulate data with a model that 
# has one covariate in abundance and one in detection, and we analyse the resulting data 
# with the data-generating model. However, in the first case, we have two independent
# covariates, which we call habitat (for abundance) and wind (for detection), 
# while in the second case, we have only a single covariate (habitat), for which we 
# simulate an effect in abundance and two different effects in detection.
#
# In contrast to simulations 1 and 1B, we now use the functionality of the new IDS() function 
# in R package unmarked for model fitting, which is much faster than MCMC with JAGS.
# Needs a version of unmarked with the IDS() function; see the README file in this folder.

# Load package
library(unmarked)

# Setup
library(AHMbook)



# Simulation 2a: Different covariates in abundance and detection
# --------------------------------------------------------------
# Note that in the data simulation function, setting mean.lambda <- 100 will 
# yield a density of 1 per ha, if we set B = 500

# Number of simulation reps
simrep <- 100                   # For test run only ... this code runs very swift !
simrep <- 1000                  # For real

# Create R object to save the true parameter values
true.vals1 <- array(NA, dim = c(6, simrep))
rownames(true.vals1) <- c("mean_lam", "beta_lam", "log sigmaDS", "beta sigmaDS", "log sigmaPC", "beta sigmaPC")

# Create R object to save the estimates
esti1 <- array(NA, dim = c(6, 4, simrep))
dimnames(esti1) <- list(c("mean_lam", "beta_lam", "log sigmaDS", "beta sigmaDS", "log sigmaPC", "beta sigmaPC"), c("MLE", "SE", "LCL", "UCL"), NULL)


# Choose param values in state and in detection all at once
# Shared stuff
mean.lambda <- 100 # For B = 500, should yield density of 1 per ha
beta.lambda <- 1
B <- 500

# HDS data
nsites_hds <- 200
sigHDS <- 100      # meters !
J <- 4             # Number of distance bins
maxDist_hds <- 200 # Maximum distance of observation in the HDS protocol,
                   # Also this is in meters
db <- seq(0, maxDist_hds, length.out=(J+1)) # Distance bins

# PC data
nsites_pc <- 1000
sigPC <- 150       # Note this is in meters now


# Launch the sim
system.time(
  for(k in 1:simrep){        # Loop over all simreps

  cat("\n\n*** Simulation Number", k, "***\n")

  ## Pick a value for the two slopes of the detection covariate
  # randomly from a uniform(-0.5, 0.5)
  slope.sigHDS <- runif(1, -0.5, 0.5)
  slope.sigPC <- runif(1, -0.5, 0.5)

  # Save all the true values (some of this pretty redundant....)
  true.vals1[1, k] <- mean.lambda
  true.vals1[2, k] <- beta.lambda
  true.vals1[3, k] <- log(sigHDS)
  true.vals1[4, k] <- slope.sigHDS
  true.vals1[5, k] <- log(sigPC)
  true.vals1[6, k] <- slope.sigPC

  ## Simulate Point Count data set
  dat2 <- simHDS(type="point", nsites = nsites_pc, 
    mean.lambda = mean.lambda, beta.lam = beta.lambda, 
    mean.sigma = sigPC, beta.sig = slope.sigPC, B = B, 
    discard0 = FALSE, show.plot = FALSE)
  # Re-format to put in a unmarkedFramePCount object
  y_pc <- matrix(0, nrow=nsites_pc, ncol=1)
  for (i in 1:nsites_pc){
    dsub <- dat2$data[dat2$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    y_pc[i,1] <- length(dsub)
  }
  # Two site covariates: habitat (on lambda) and wind (on sigma)
  sc_pc1 <- data.frame(habitat=dat2$habitat)
  sc_pc2 <- data.frame(wind=dat2$wind)

  # Create unmarked frame for PC data set
  umf_pc <- unmarkedFramePCount(y=y_pc,
    siteCovs=data.frame(habitat=sc_pc1, wind=sc_pc2))
#  summary(umf_pc)

  ## Simulate a regular distance sampling data set
  dat1 <- simHDS(type="point", nsites = nsites_hds, 
    mean.lambda = mean.lambda, beta.lam = beta.lambda, 
    mean.sigma = sigHDS, beta.sig = slope.sigHDS, B = B, 
    discard0 = FALSE, show.plot = FALSE)
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
  sc_hds1 <- data.frame(habitat=dat1$habitat)
  sc_hds2 <- data.frame(wind=dat1$wind)
  # Create unmarked frame
  umf_hds <- unmarkedFrameDS(y=y_hds, 
    siteCovs=data.frame(habitat=sc_hds1, wind=sc_hds2), 
    survey="point", dist.breaks=db, unitsIn="m")

  # Fit the IDS with HDS + PC data
  # We add an "error-catch-saver" using try()
  out <- try(IDS(lambdaformula = ~habitat, 
    detformulaDS = ~wind, detformulaPC = ~wind, 
    dataDS = umf_hds, dataPC = umf_pc, maxDistPC = 500,
    K = 1000, control = list(trace=1, REPORT = 10)) )
  if(class(out) == "try-error") {fuck <- TRUE}
  else{
  tmp <- summary(out)
  ses <- c(tmp$lam$SE, tmp$ds$SE, tmp$pc$SE) 


  # Extract and save param estimates
  esti1[, 1, k] <- coef(out)  # Save MLEs
  esti1[, 2, k] <- ses        # Save SEs
  esti1[1:2, 3:4, k] <- confint(out, type = 'lam') # CIs
  esti1[3:4, 3:4, k] <- confint(out, type = 'ds')
  esti1[5:6, 3:4, k] <- confint(out, type = 'pc')
  }
}
) # End of system.time()



# Weed out the numerical failures, e.g., when NAs
sum(is.na(esti1))
bad1 <- rep(0, simrep)
for(k in 1:simrep){
  bad1[k] <- ifelse(sum(is.na(esti1[,,k])) > 0, 1, 0)
}
sum(bad1)
toss1 <- which(bad1 == 1)
# We also take a sigma intercept estimate >1000 on the real scale as indicating a problem.
bad.sigma.DS1 <- rep(0, simrep)
for(k in 1:simrep){
  bad.sigma.DS1[k] <- ifelse(exp(esti1[3,1,k]) > 1000, 1, 0)
}
sum(bad.sigma.DS1)
toss1 <- which(rowSums(cbind(bad1, bad.sigma.DS1)) > 0)
length(toss1)


# Compare true values and estimated values
# (This figure is not in the paper)
# Note that when nobody must be tossed out, then have to 
# eliminate the '-toss1' in next block of code, since code will break if toss1 is empty

par(mfrow = c(2,3), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
nb <- 50
# Abundance model
hist(exp(esti1[1,1, -toss1]), breaks = nb, col = 'grey', freq = FALSE, xlab = 'Abundance intercept', main = 'mean.lambda')
abline(v = 1, col = 'red', lwd = 3, lty = 3)
hist(esti1[2,1, -toss1], breaks = nb, col = 'grey', freq = FALSE, xlab = 'Abundance slope', main = 'Slope of abundance \non habitat covariate (beta.lambda)')
abline(v = 1, col = 'red', lwd = 3, lty = 3)
# Detection model in DS data
hist(exp(esti1[3,1, -toss1]), breaks = nb, col = 'grey', freq = FALSE, xlab = 'Sigma intercept (DS)', main = 'mean.sigma.DS')
abline(v = 100, col = 'red', lwd = 3, lty = 3)
plot(true.vals1[4,-toss1], esti1[4,1, -toss1], xlab = 'True sigma slope (DS)', ylab = 'Estimated sigma slope (DS)', main = 'Slope of sigma(DS)\non wind covariate', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 3, lty = 3)
# Detection model in PC data
hist(exp(esti1[5,1, -toss1]), breaks = nb, col = 'grey', freq = FALSE, xlab = 'Sigma intercept (PC)', main = 'mean.sigma.PC')
abline(v = 150, col = 'red', lwd = 3, lty = 3)
plot(true.vals1[6, -toss1], esti1[6,1, -toss1], xlab = 'True sigma slope (PC)', ylab = 'Estimated sigma slope (PC)', main = 'Slope of sigma(PC)\non wind covariate', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 3, lty = 3)
# This looks just wonderful :)




# Simulation 2b: Use the same abundance and detection covariate for both HDS and PC data
# --------------------------------------------------------------------------------------
# This must be the harder case. For it, we first need to write a variant of the simHDS() function
# that enables one to generate such data where a single covariate affects 
# abundance and detection. Here is a new version of the function, adapting the version 
# that is in the AHMbook package. We have to execute the next block of code to 
# define that function and then later use it.

# --------------- Start of function definition ----------------------
simHDSX <- function (type = c("line", "point"), nsites = 100, mean.lambda = 2, 
    beta.lam = 1, mean.sigma = 1, beta.sig.hab = -0.5, beta.sig.wind = 0, B = 3, discard0 = TRUE, 
    show.plot = TRUE) 
{
    nsites <- round(nsites[1])
    type <- match.arg(type)
    habitat <- rnorm(nsites)  # Covariate 1
    wind <- runif(nsites, -2, 2) # covariate
    lambda <- exp(log(mean.lambda) + beta.lam * habitat)
    N <- rpois(nsites, lambda)
    N.true <- N
    sigma <- exp(log(mean.sigma) + beta.sig.hab * habitat + beta.sig.wind * wind)
    data <- NULL
    for (i in 1:nsites) {
        if (N[i] == 0) {
            data <- rbind(data, c(i, NA, NA, NA, NA))
            next
        }
        if (type == "line") {
            d <- runif(N[i], 0, B)
            p <- exp(-d * d/(2 * (sigma[i]^2)))
            y <- rbinom(N[i], 1, p)
            u <- v <- rep(NA, N[i])
            d <- d[y == 1]
            u <- u[y == 1]
            v <- v[y == 1]
            y <- y[y == 1]
        }
        if (type == "point") {
            u <- runif(N[i], 0, 2 * B)
            v <- runif(N[i], 0, 2 * B)
            d <- sqrt((u - B)^2 + (v - B)^2)
            N.true[i] <- sum(d <= B)
            p <- exp(-d * d/(2 * (sigma[i]^2)))
            pp <- ifelse(d <= B, 1, 0) * p
            y <- rbinom(N[i], 1, pp)
            u <- u[y == 1]
            v <- v[y == 1]
            d <- d[y == 1]
            y <- y[y == 1]
        }
        if (sum(y) > 0) 
            data <- rbind(data, cbind(rep(i, sum(y)), y, u, v, 
                d))
        else data <- rbind(data, c(i, NA, NA, NA, NA))
    }
    colnames(data) <- c("site", "y", "u", "v", "d")
    if (discard0) 
        data <- data[!is.na(data[, 2]), ]
    if (show.plot) {
        if (type == "line") {
            op <- par(mfrow = c(1, 3))
            on.exit(par(op))
            tryPlot <- try({
                hist(data[, "d"], col = "lightblue", breaks = 20, 
                  main = "Frequency of distances", xlab = "Distance")
                ttt <- table(data[, 1])
                n <- rep(0, nsites)
                n[as.numeric(rownames(ttt))] <- ttt
                plot(habitat, n, main = "Observed counts (n) vs. habitat")
                plot(wind, n, main = "Observed counts (n) vs. wind speed")
            }, silent = TRUE)
            if (inherits(tryPlot, "try-error")) 
                tryPlotError(tryPlot)
        }
        if (type == "point") {
            op <- par(mfrow = c(2, 2))
            on.exit(par(op))
            tryPlot <- try({
                plot(data[, "u"], data[, "v"], pch = 16, main = "Located individuals in point transects", 
                  xlim = c(0, 2 * B), ylim = c(0, 2 * B), col = data[, 
                    1], asp = 1)
                points(B, B, pch = "+", cex = 3, col = "black")
                plotrix::draw.circle(B, B, B)
                hist(data[, "d"], col = "lightblue", breaks = 20, 
                  main = "Frequency of distances", xlab = "Distance")
                ttt <- table(data[, 1])
                n <- rep(0, nsites)
                n[as.numeric(rownames(ttt))] <- ttt
                plot(habitat, n, main = "Observed counts (n) vs. habitat")
                plot(wind, n, main = "Observed counts (n) vs. wind speed")
            }, silent = TRUE)
            if (inherits(tryPlot, "try-error")) 
                tryPlotError(tryPlot)
        }
    }
    list(type = type, nsites = nsites, mean.lambda = mean.lambda, 
        beta.lam = beta.lam, mean.sigma = mean.sigma, 
        beta.sig.hab = beta.sig.hab, beta.sig.wind = beta.sig.wind, 
        B = B, data = data, habitat = habitat, wind = wind, 
        N = N, N.true = N.true)
}
# --------------- End of function definition ----------------------

# We specify one effect of habitat on abundance and two different ones on the 
# two detection functions (for the HDS data and the PC data).

library(unmarked)

# Number of simulation reps
simrep <- 100        # For test case
simrep <- 1000

# Create R object to save the true parameter values
true.vals2 <- array(NA, dim = c(6, simrep))
rownames(true.vals2) <- c("mean_lam", "beta_lam", "log sigmaDS", "beta sigmaDS", "log sigmaPC", "beta sigmaPC")

# Create R object to save the estimates
esti2 <- array(NA, dim = c(6, 4, simrep))
dimnames(esti2) <- list(c("mean_lam", "beta_lam", "log sigmaDS", "beta sigmaDS", "log sigmaPC", "beta sigmaPC"), c("MLE", "SE", "LCL", "UCL"), NULL)


# Choose param values in state and in detection all at once
# Shared stuff
mean.lambda <- 100 # For B = 500, should yield density of 1 per ha
beta.lambda <- 1
B <- 500

# HDS data
nsites_hds <- 200
sigHDS <- 100      # meters !
J <- 4             # Number of distance bins
maxDist_hds <- 200 # Maximum distance of observation in the HDS protocol,
                   # Also this is in meters
db <- seq(0, maxDist_hds, length.out=(J+1)) # Distance bins

# PC data
nsites_pc <- 1000
sigPC <- 150       # Note this is in meters now


# Launch the sim
system.time(
  for(k in 1:simrep){        # Loop over all simreps

  ## Pick a value for the two slopes of the detection covariate
  # randomly from a uniform(-0.5, 0.5)
  slope.sigHDS <- runif(1, -0.5, 0.5)
  slope.sigPC <- runif(1, -0.5, 0.5)

  # Save all the true values (some of this pretty redundant....)
  true.vals2[1, k] <- mean.lambda
  true.vals2[2, k] <- beta.lambda
  true.vals2[3, k] <- log(sigHDS)
  true.vals2[4, k] <- slope.sigHDS
  true.vals2[5, k] <- log(sigPC)
  true.vals2[6, k] <- slope.sigPC

  ## Simulate Point Count data set
  dat2 <- simHDSX(type="point", nsites = nsites_pc, 
    mean.lambda = mean.lambda, beta.lam = beta.lambda, 
    mean.sigma = sigPC, beta.sig.hab = slope.sigPC, 
    beta.sig.wind = 0, B = B, discard0 = FALSE, show.plot = FALSE)
  # Re-format to put in a unmarkedFramePCount object
  y_pc <- matrix(0, nrow=nsites_pc, ncol=1)
  for (i in 1:nsites_pc){
    dsub <- dat2$data[dat2$data[,"site"]==i,,drop=FALSE][,"d"]
    if(any(is.na(dsub))) next # Skip if no detections
    y_pc[i,1] <- length(dsub)
  }
  # Two site covariates: habitat (on lambda) and wind (on sigma)
  sc_pc1 <- data.frame(habitat=dat2$habitat)

  # Create unmarked frame for PC data set
  umf_pc <- unmarkedFramePCount(y=y_pc,
    siteCovs=data.frame(habitat=sc_pc1))

  ## Simulate a regular distance sampling data set
  dat1 <- simHDSX(type="point", nsites = nsites_hds, 
    mean.lambda = mean.lambda, beta.lam = beta.lambda, 
    mean.sigma = sigHDS, beta.sig.hab = slope.sigHDS, 
    beta.sig.wind = 0,B = B, discard0 = FALSE, show.plot = FALSE)
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
  sc_hds1 <- data.frame(habitat=dat1$habitat)
  # Create unmarked frame
  umf_hds <- unmarkedFrameDS(y=y_hds, 
    siteCovs=data.frame(habitat=sc_hds1, wind=sc_hds2), 
    survey="point", dist.breaks=db, unitsIn="m")

  # Fit the IDS with HDS + PC data
  # We add an "error-catch-saver" using try()
  out <- try(IDS(lambdaformula = ~habitat, 
    detformulaDS = ~habitat, detformulaPC = ~habitat, 
    dataDS = umf_hds, dataPC = umf_pc, maxDistPC = 500,
    K = 1000, control = list(trace=1, REPORT = 10)) )
  if(class(out) == "try-error") {fuck <- TRUE}
  else{
  tmp <- summary(out)
  ses <- c(tmp$lam$SE, tmp$ds$SE, tmp$pc$SE) 


  # Extract and save param estimates
  esti2[, 1, k] <- coef(out)  # Save MLEs
  esti2[, 2, k] <- ses        # Save SEs
  esti2[1:2, 3:4, k] <- confint(out, type = 'lam') # CIs
  esti2[3:4, 3:4, k] <- confint(out, type = 'ds')
  esti2[5:6, 3:4, k] <- confint(out, type = 'pc')
  }
}
) # End of system.time()


# Weed out numerical failures:
sum(is.na(esti2))
bad2 <- rep(0, simrep)
for(k in 1:simrep){
  bad2[k] <- ifelse(sum(is.na(esti2[,,k])) > 0, 1, 0)
}
sum(bad2)
toss2 <- which(bad2 == 1)
bad.sigma.DS2 <- rep(0, simrep)

for(k in 1:simrep){
  bad.sigma.DS2[k] <- ifelse(exp(esti2[3,1,k]) > 1000, 1, 0)
}
sum(bad.sigma.DS2)
toss2 <- which(rowSums(cbind(bad2, bad.sigma.DS2)) > 0)
length(toss2)

# Compare true values and estimated values
# (This figure is not in the paper)
# NOTE: Next block of code will again break when toss2 is empty set

par(mfrow = c(2,3), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
nb <- 50
# Abundance model
hist(exp(esti2[1,1, -toss2]), breaks = nb, col = 'grey', freq = FALSE, xlab = 'Abundance intercept', main = 'mean.lambda')
abline(v = 1, col = 'red', lwd = 3, lty = 3)
hist(esti2[2,1, -toss2], breaks = nb, col = 'grey', freq = FALSE, xlab = 'Abundance slope', main = 'Slope of abundance \non habitat covariate (beta.lambda)')
abline(v = 1, col = 'red', lwd = 3, lty = 3)
# Detection model in DS data
hist(exp(esti2[3,1, -toss2]), breaks = nb, col = 'grey', freq = FALSE, xlab = 'Sigma intercept (DS)', main = 'mean.sigma.DS')
abline(v = 100, col = 'red', lwd = 3, lty = 3)
plot(true.vals2[4,-toss2], esti2[4,1, -toss2], xlab = 'True sigma slope (DS)', ylab = 'Estimated sigma slope (DS)', main = 'Slope of sigma(DS)\non habitat covariate', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 3, lty = 3)
# Detection model in PC data
hist(exp(esti2[5,1, -toss2]), breaks = nb, col = 'grey', freq = FALSE, xlab = 'Sigma intercept (PC)', main = 'mean.sigma.PC')
abline(v = 150, col = 'red', lwd = 3, lty = 3)
plot(true.vals2[6, -toss2], esti2[6,1, -toss2], xlab = 'True sigma slope (PC)', ylab = 'Estimated sigma slope (PC)', main = 'Slope of sigma(PC)\non habitat covariate', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 3, lty = 3)





# --------------------------------- #
# Next produces Figure 2 in paper   #
# --------------------------------- #

# Compare true values and estimated values
par(mfrow = c(4, 2), mar = c(4, 5, 5, 2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
nb <- 30

# Detection model in DS data
# Intercept (Sim Set 1)
hist(exp(esti1[3,1, -toss1]), breaks = nb, col = 'grey', freq = TRUE, xlab = '', main = 'Simulation 2a\nIntercept (DS)')
abline(v = 100, col = 'red', lwd = 2, lty = 1)
abline(v = mean(exp(esti1[3,1, -toss1])), col = 'blue', lwd = 2, lty = 2)
# Intercept (Sim Set 2)
hist(exp(esti2[3,1, -toss2]), breaks = nb, col = 'grey', freq = TRUE, xlab = '', main = 'Simulation 2b\nIntercept (DS)')
abline(v = 100, col = 'red', lwd = 2, lty = 1)
abline(v = mean(exp(esti2[3,1, -toss2])), col = 'blue', lwd = 2, lty = 2)

# Slope (Sim Set 1)
plot(true.vals1[4,-toss1], esti1[4,1, -toss1], xlab = 'Truth', ylab = 'Estimate', main = 'Slope (DS)', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 2, lty = 1)
abline(lm(esti1[4,1, -toss1]~ true.vals1[4,-toss1]), col = 'blue', lwd = 2, lty = 2)

# Slope (Sim Set 2)
plot(true.vals2[4,-toss2], esti2[4,1, -toss2], xlab = 'Truth', ylab = 'Estimate', main = 'Slope (DS)', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 2, lty = 1)
abline(lm(esti2[4,1, -toss2]~ true.vals2[4,-toss2]), col = 'blue', lwd = 2, lty = 2)


# Detection model in PC data
# Intercept (Sim Set 1)
hist(exp(esti1[5,1, -toss1]), breaks = nb, col = 'grey', freq = TRUE, xlab = '', main = 'Intercept (PC)')
abline(v = 150, col = 'red', lwd = 2, lty = 1)
abline(v = mean(exp(esti1[5,1, -toss1])), col = 'blue', lwd = 2, lty = 2)

# Intercept (Sim Set 2)
hist(exp(esti2[5,1, -toss2]), breaks = nb, col = 'grey', freq = TRUE, xlab = '', main = 'Intercept (PC)')
abline(v = 150, col = 'red', lwd = 2, lty = 1)
abline(v = mean(exp(esti2[5,1, -toss2])), col = 'blue', lwd = 2, lty = 2)

# Slope (Sim Set 1)
plot(true.vals1[6, -toss1], esti1[6,1, -toss1], xlab = 'Truth', ylab = 'Estimate', main = 'Slope (PC)', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 2, lty = 1)
abline(lm(esti1[6,1, -toss1]~ true.vals1[6,-toss1]), col = 'blue', lwd = 2, lty = 2)

# Slope (Sim Set 2)
plot(true.vals2[6, -toss2], esti2[6,1, -toss2], xlab = 'Truth', ylab = 'Estimate', main = 'Slope (PC)', pch = 16, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, col = 'red', lwd = 2, lty = 1)
abline(lm(esti2[6,1, -toss2]~ true.vals2[6,-toss2]), col = 'blue', lwd = 2, lty = 2)
