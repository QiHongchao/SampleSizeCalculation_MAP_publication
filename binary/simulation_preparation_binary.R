##Package
library(rjags)
library(RBesT)
library(ggplot2)
library(ggpubr)

##Summary data for historical controls, data from liver transplant clinical trials (Friede et al., 2017, RSM)
event_historical <- c(4, 6, 4)
nsub_historical <- c(20, 54, 34)
J <- length(nsub_historical)
##Precision for the between-study SD
prec_between_sd <- 1 
##New trial, OR
OR_candidate <- seq(2, 6, 2)/10
##Setting of simulation
num_iter <- 20000
perc_burnin <- 0.2
num_chains <- 4
num_simulation <- 5000 ##Number of simulation for each sample size
##Random seeds for Step 1
set.seed(1234)
seeds_jags <- sample(.Machine$integer.max, num_chains)

##Derive the MAP prior of the control proportion based on the historical data
jags_data_map <- list(J = J, prec_between_sd = prec_between_sd,
                      nsub_historical = nsub_historical, event_historical = event_historical)

##MAP analysis
map_jags_model <- jags.model(file = "./jags/map_binary.txt", data = jags_data_map, 
                             n.chains = num_chains, n.adapt = 1000, quiet = T,
                             inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[1]),
                                          list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[2]),
                                          list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[3]),
                                          list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[4])))

##Burnin stage
update(map_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")

##Sampling after burnin
params <- c("logit_p0", "logit_p0_new", "between_sd")
map_jags <- coda.samples(map_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                         progress.bar = "none")
##The samples, derive the MAP prior from historical controls
map_sample <- as.data.frame(do.call(rbind, map_jags))
##The prior for p0_new
map_sample$p0_new <- 1/(1 + exp(-map_sample$logit_p0_new))

##Propose a list of candidate sample sizes
samplesize_num <- 6
samplesize_candidate <- c(seq(20, 120, length.out = samplesize_num), 
                          seq(30, 280, length.out = samplesize_num),
                          seq(350, 750, length.out = samplesize_num))
interval <- c(20, 50, 80)
hypothesis_testing <- matrix(NA, num_simulation, samplesize_num)
##Random seeds
set.seed(2021)
seeds_simulation <- matrix(sample(.Machine$integer.max, num_simulation*samplesize_num), 
                           nrow = num_simulation)
seeds_jags_simulation <- matrix(sample(.Machine$integer.max, num_simulation*num_chains*samplesize_num), 
                                nrow = num_simulation*num_chains)

##Data set to save power values
power_res <- data.frame(num_simulation = num_simulation, samplesize = samplesize_candidate, 
                        OR = rep(OR_candidate, each = samplesize_num), 
                        power = NA)

##ESS by RBesT
beta_mix <- mixfit(map_sample$p0_new, type = "beta", Nc = 2)
RBesT::ess(beta_mix, method = "moment")
RBesT::ess(beta_mix, method = "morita")
RBesT::ess(beta_mix, method = "elir")
