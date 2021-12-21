##Package
library(rjags)

##Historical data, data from liver transplant clinical trial (Friede et al., 2017, RSM)
event_historical <- c(4, 6, 4)
nsub_historical <- c(20, 54, 34)
J <- length(nsub_historical)
##This precision can be justified, sd = 2 for logitp0 is quite uninformative
##Halfnormal(1) has its 95th percentile at approx. 2
prec_between_sd <- 1 
##New trial, OR
OR_candidate <- seq(2, 8, 2)/10
##Setting of simulation
num_iter <- 5000
perc_burnin <- 0.2
num_chains <- 4
set.seed(1234)
seeds_jags <- sample(.Machine$integer.max, num_chains)
num_simulation <- 5000 ##Number of simulation for each sample size

##Step 1: Get the MAP prior of the control proportion based on the historical data
jags_data_map <- list(J = J, prec_between_sd = prec_between_sd,
                      nsub_historical = nsub_historical, event_historical = event_historical)

##MAP
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
##The samples, derive the MAP prior from the data set
map_sample <- as.data.frame(do.call(rbind, map_jags))
##The prior for p0_new
map_sample$p0_new <- 1/(1 + exp(-map_sample$logit_p0_new))


##step 7: a list of candidate sample sizes
samplesize_num <- 11
samplesize_candidate <- c(seq(20, 120, length.out = samplesize_num), 
                          seq(30, 280, length.out = samplesize_num),
                          seq(50, 750, length.out = samplesize_num),
                          seq(500, 3500, length.out = samplesize_num))
interval <- c(10, 25, 70, 300)
hypothesis_testing <- matrix(NA, num_simulation, samplesize_num)

##Random seeds for JAGS
set.seed(2021)
seeds_simulation <- matrix(sample(.Machine$integer.max, num_simulation*samplesize_num), 
                           nrow = num_simulation)
seeds_jags_simulation <- matrix(sample(.Machine$integer.max, num_simulation*num_chains*samplesize_num), 
                                nrow = num_simulation*num_chains)

##Data set to save power values
power_res <- data.frame(num_simulation = num_simulation, samplesize = samplesize_candidate, 
                        OR = rep(OR_candidate, each = samplesize_num), 
                        power = NA)

