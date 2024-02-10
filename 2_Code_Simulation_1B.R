
# Supporting Information for:
# Kéry, M., Royle, J.A., Hallman, T., Robinson, W.D., Strebel, N., Kellner, K.F., 2024:
# Integrated distance sampling models for simple point counts, Ecology.


# Written by the authors, Nov 2022


# Simulation 1B: Can we estimate separate detection functions in IDS1 under widely varying
# values of the input parameters for the data ?
#
#
# This builds on Simulation 1, where we had constants for all parameters and simply varied 
# the value of the detection function sigma in the PC or the DND portions of the data
# among simulation reps.
# Here now, for each simrep we change the values of lambda (shared for both data sets), 
# sigma in both data sets, and of the value of the truncation distance in the DS part of the data.
# In addition we have an effect of a single habitat covariate in density, for greater realism.
#
# We work with model IDS1 only, i.e., combine proper DS data with simple PC data.
# We will vary the following inputs in the data simulation within the following ranges 
# (by picking a random number from uniform distributions):
#    lambda: runif(1, 10, 500) (Corresponding to densities of 0.1 to 5)
#    sigma(DS): runif(1, 0.2, 1.2)
#    sigma(PC): runif(1, 0.2, 1.2)
#    B (in DS data set): runif(1, 1, 3), rounded to 1 digit

# We vary these by a "response surface design" and use the same sample sizes 
# as in Simulation 1: 250 sites with proper DS surveys and 1000 sites with simple PC surveys.

library(AHMbook)
library(jagsUI)

# Number of simulation reps
simrep <- 3         # For a test run
simrep <- 1000      # For real

# Create R object to save results (posterior summary table)
estimates <- array(NA, dim = c(1016, 11, simrep))

# Vector to save the true values of everything and the two data sets
true.dataDS <- list()
true.dataPC <- list()
true.vals <- array(NA, dim = c(4, simrep))
rownames(true.vals) <- c('lambda', 'sigma_DS', 'sigma_PC', 'B_DS')

# Launch simulation
for(k in 1:simrep){      # Loop over simreps

  cat(paste('\n*** Model IDS1: Sim Number', k, '***\n'))

  # Draw new values of lambda, 
  # sigma(DS) and sigma(PC) [both from a U(0.2, 1.2)]
  # and of B(DS) [from U(1, 3)] and save them
  lam <- runif(1, 10, 500)
  sigDS <- runif(1, 0.2, 1.2)
  sigPC <- runif(1, 0.2, 1.2)
  newB <- round(runif(1, 1, 3), 1)   # Note rounded
  true.vals[1, k] <- lam
  true.vals[2, k] <- sigDS
  true.vals[3, k] <- sigPC
  true.vals[4, k] <- newB
  cat(paste('True lambda:', round(lam, 3), '\nTrue sigmaDS:', round(sigDS, 3), '\nTrue sigmaPC:', round(sigPC, 3), '\nB in DS data:', round(newB,3), '\n'))

  # Simulate another Distance-sampling data set of size 250 sites and save
  dat1 <- simHDS(type="point", nsites = 250, 
    mean.lambda = lam, beta.lam = 1, mean.sigma = sigDS, beta.sig = 0, 
    B = 5, discard0 = FALSE, show.plot = TRUE)
  true.dataDS[[k]] <- dat1

  # Simulate a simple PC data set of size 1,000 sites and save
  dat2 <- simHDS(type="point", nsites = 1000, 
    mean.lambda = lam, beta.lam = 1, mean.sigma = sigPC, beta.sig = 0, 
    B = 5, discard0 = FALSE, show.plot = TRUE)
  true.dataPC[[k]] <- dat2

  # Data wrangling similar to Simulation 1
  # Get C2(counts in data set 2)
  data2 <- dat2$data
  habitat2 <- dat2$habitat
  ncap2 <- tapply(data2[,5], data2[,1], function(x) sum(!is.na(x)))
  C2 <- as.vector(ncap2)  # Number of individuals detected per site

  # Create the Distance sampling data subset (call this data set 1)
  data1 <- dat1$data
  habitat1 <- dat1$habitat

  # Subset the DS data to a truncation distance (B) of newB
  data1[,'d'][data1[,'d'] > newB] <- NA
  sum(!is.na(data1[,'d']))

  # Prepare the DS data and summarize some too
  # Number of individuals detected (or 'captured') per site
  ncap <- tapply(data1[,5], data1[,1], function(x) sum(!is.na(x)))
  ncap <- as.vector(ncap)  # Number of individuals detected per site

  # Prepare other data
  # Site ID of each non-missing observation
  site1 <- data1[!is.na(data1[,'d']),1]

  # Define new B for the DS data set (data 1)
  #  newB <- 2

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
  dclass <- dclass[!is.na(dclass)]   # toss out missing dclasses

  # Total number of individuals detected (after applying new B = 2)
  nind <- length(dclass)

  # Compute the areas associated with the data for each site
  A1 <- pi * newB^2         # Area associated with DS surveys: variable now
  A2 <- pi * fullDistance^2 # [1] 78.53982, area associated with PC surveys

  # Bundle and summarize data set to be fed into JAGS
  nsites1 <- dat1$nsites   # Data set 1
  nsites2 <- length(C2)    # Data set 2
  bdata <- list(nsites1 = nsites1, nind = nind, newB = newB, nD = nD,
    midpt = midpt, fullDistance = fullDistance, nDfull = nDfull, 
    midptFull = midptFull, delta = delta, ncap = ncap, 
    habitat1 = habitat1, dclass = dclass, site1 = site1, A1 = A1, 
    nsites2 = nsites2, C2 = C2, habitat2 = habitat2, A2 = A2)


# Specify model in BUGS language
cat(file="modelIDS1.txt","
model{

# Priors for parameters
# Separate parameters in the detection model
alpha0DS ~ dnorm(0, 0.01)    # For proper DS data set
alpha0PC ~ dnorm(0, 0.01)    # For simple PC data set
# Shared parameters in the abundance/density model
beta0 ~ dnorm(0, 0.01)       # Abundance intercept
beta1 ~ dnorm(0, 0.01)       # Abundance slope on habitat

# Submodel for the Distance sampling data (= data set 1)
# -------------------------------------------------------------------
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
  # Note A1 is offset (survey area) and lambda1 becomes per-unit abundance, i.e., density
  N1[s] ~ dpois(A1 * lambda1[s])  # Part 3 of HM: true number there

  # Log-linear models for abundance and detection in DS data set
  log(lambda1[s]) <- beta0 + beta1 * habitat1[s] # Linear model abundance
  log(sigma[s]) <- alpha0DS   # Model for sigma just a constant for now
}

# Now we compute the average p over all distance bands, pbar, 
# for the unlimited-distance surveys (i.e., for B = 5).
# It is this pcap2 (or here, it's average pbar2) which we will use 
# below in the submodel for the simple PC data (= data set 2)
# Note that for now pcap2 is identical for all simple PC sites
# since we don't have any detection covariates in the model
for(s in 1:nsites2){        # Loop over all sites in data set 2
  # Linear model for detection in simple PC data set
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


# Submodel for the simple PC data (= data set 2)
# ----------------------------------------------
# We have defined the priors already above:
#  - sigma in the detection function is not shared between the data sets
#  - but the two abundance parameters (beta0, beta1) are shared
#
# Likelihood for the simple PC data
# Note the exp(offset) equal to 'N.scaling', which is the ratio 
# A2(simple PC) / A1(DS with B = 2)
# Remember that beta0 and beta1 are shared params, but not sigma
for(s in 1:nsites2){
  # Note A2 is again an offset (survey area) and 
  # lambda2 then becomes per-unit abundance, i.e., density
  N2[s] ~ dpois(A2 * lambda2[s])
  log(lambda2[s]) <- beta0 + beta1 * habitat2[s]
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
    beta1 = rnorm(1, 0, 0.1), N1 = Nst1, N2 = Nst2)}

# Params to save
  params <- c("alpha0DS", "alpha0PC", "beta0", "beta1", "mean.sigmaDS",
    "mean.sigmaPC", "Ntotal1", "D1", "mean.lambda",
    "pbar1", "pbar2", "Ntotal2", "D2", "Ntotal", "D", "N2")

# MCMC settings (Choose one)
#   na <- 100;  nc <- 3;  ni <- 200;  nb <- 100;  nt <- 2  # for test run
  na <- 3000;  nc <- 3;  ni <- 25000;  nb <- 5000;  nt <- 2 # for real

# Launch JAGS (ART 1 min), check convergence and summarize posteriors
out <- jags(bdata, inits, params, "modelIDS1.txt", n.adapt = na,
  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
  parallel = TRUE)
# traceplot(out)     # not shown
print(out, 2)

# Save posterior summary
estimates[,,k] <- out$summary

}            # End of loop over simreps


# NOTE: Before summarization of results, you may have to toss out 'invalid'
# solutions. That is, model runs that have not converged.


# Summarize in plots

# Choose how many simreps to summarize
maxK <- 3         # For test run
maxK <- 1000      # for real

# Note we will again express all sigma in units of metres, 
# which means we here multiply relevant things by 100.

# Compute estimation error for all estimates of N2, the population sizes in the simulated 
# simple PC data, expressed by the difference between the estimates and the truth. 

# Compute errors of the N2 estimates
# Get point estimates
tmp_esti <- estimates[16:1015, 1, ]

# Get true values
tmp_truth <- array(NA, dim = dim(tmp_esti))
for(k in 1:maxK){
  tmp_truth[, k] <- true.dataPC[[k]]$N.true
}

# Get absolute and relative differences
absDiff <- tmp_esti - tmp_truth

# Tests and replacement of Inf with NA
sum(is.na(absDiff)) ; sum(absDiff == 'Inf') ; sum(absDiff == '-Inf')

# Also compute the estimation error between the true and the estimated density 
# in all surveyed sites (i.e., 250 DS and 1000 PC sites)

# Get true values of sumN inside of the 5 unit circle
trueNtotalDS <- trueNtotalPC <- numeric(maxK)
for(k in 1:maxK){
  trueNtotalDS[k] <- sum(true.dataDS[[k]]$N.true)
  trueNtotalPC[k] <- sum(true.dataPC[[k]]$N.true)
}
trueNtotal <- trueNtotalDS + trueNtotalPC
trueD <- trueNtotal / (1250 * pi * 5^2)

# Get point estimates and CRI of mean density
Desti <- estimates[15, 1, ]

# Compute absolute and relative differences of density
absDiffDens <- Desti - trueD

# Tests and replacement of Inf with NA
sum(is.na(absDiffDens)) ; sum(absDiffDens == 'Inf') ; sum(absDiffDens == '-Inf')




# ---------------------------------------- #
# Next produces Figure S1 in Appendix S2   #
# ---------------------------------------- #



# Estimated versus true values of sigma in the two data sets
xlim <- ylim <- c(15, 150)
par(mfrow = c(2, 2), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# Estimated versus true values of sigma in the two data sets
plot(100*true.vals[2, ], 100*estimates[5, 1, ], xlab = "True sigma (DS)", ylab = "Estimate", frame = FALSE, col = rgb(0,0,0,0.3), pch = 16, main = "sigma at DS locations", cex = 2, xlim = xlim, ylim = ylim)
segments(100*true.vals[2, ], 100*estimates[5, 3, ], 100*true.vals[2, ], 100*estimates[5, 7, ])
abline(0, 1, col = 'red', lwd = 3)
abline(lm(I(100*estimates[5, 1, ]) ~ I(100*true.vals[2, ])), col = 'blue', lty = 2, lwd = 3)

plot(100*true.vals[3, ], 100*estimates[6, 1, ], xlab = "True sigma (PC)", ylab = "Estimate", frame = FALSE, col = rgb(0,0,0,0.3), pch = 16, main = "sigma at PC locations", cex = 2, xlim = xlim, ylim = ylim)
segments(100*true.vals[3, ], 100*estimates[6, 3, ], 100*true.vals[3, ], 100*estimates[6, 7, ])
abline(0, 1, col = 'red', lwd = 3)
abline(lm(I(100*estimates[6, 1, ]) ~ I(100*true.vals[3, ])), col = 'blue', lty = 2, lwd = 3)

# Absolute error of site-level abundance
hist(absDiff, xlab = "Absolute error in N", col = "grey", main = "Site-level N at PC locations", cex = 2, freq = TRUE, breaks = 700, xlim = c(-100, 100))
abline(v = 0, col = 'red', lwd = 2)
abline(v = mean(absDiff), col = 'blue', lty = 2, lwd = 2)

# Absolute error of average density (all plots confounded)
hist(absDiffDens, xlab = "Absolute error in D", col = "grey", main = "Mean density at all locations", cex = 2, freq = TRUE, breaks = 50, xlim = c(-0.8, 0.8))
abline(v = 0, col = 'red', lwd = 2)
abline(v = mean(absDiffDens), col = 'blue', lty = 2, lwd = 2)
