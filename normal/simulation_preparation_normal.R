##Package
library(rjags)
library(RBesT)
library(ggplot2)
library(ggpubr)

##Data preparation
J <- 5
estimate_historical <- c(6.1, 7.3, 5.7, 7.3, 7.5)
se_historical <- c(0.7, 0.66, 0.81, 0.63, 0.64)
var_historical <- se_historical^2
nsub_historical <- c(301, 275, 155, 304, 291)
sd_historical <- se_historical*sqrt(nsub_historical)
sd <- sd_current <- ceiling(mean(sd_historical))
##Simulation settings
num_chains <- 4
num_iter <- 20000
perc_burnin <- 0.2
num_simulation <- 5000
set.seed(1234)
seeds_jags <- sample(.Machine$integer.max, num_chains)
trteff_candidate <- 2:4

##MAP data
jags_data_map <- list(J = J, estimate = estimate_historical, se = se_historical, sd = sd)

##MAP
bayesianma <- function(file, jags_data_map, params) {
  ##MAP model
  map_jags_model <- jags.model(file = file, data = jags_data_map, 
                               n.chains = num_chains, n.adapt = 1000, quiet = T,
                               inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[1]),
                                            list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[2]),
                                            list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[3]),
                                            list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags[4])))
  
  ##Burnin stage
  update(map_jags_model, n.iter = num_iter*perc_burnin, progress.bar = "none")
  
  ##Sampling after burnin
  params <- params
  
  map_jags <- coda.samples(map_jags_model, params, n.iter = num_iter*(1-perc_burnin), 
                           progress.bar = "none")
  
  map_sample <- as.data.frame(do.call(rbind, map_jags))
  
  return(list(map_jags = map_jags, map_sample = map_sample))
}
##Random effects meta analysis
bma_random <- bayesianma(file = "./jags/map_normal.txt", 
                         jags_data_map = jags_data_map,
                         params = c("mu", "sigma", "theta_new"))
summary(bma_random$map_jags)
plot(bma_random$map_jags)
hist(bma_random$map_sample$mu, breaks = 50)
hist(bma_random$map_sample$theta_new, breaks = 50)

##ESS calculation
normal_mix <- mixfit(bma_random$map_sample$theta_new, type = "norm", Nc = 1)
RBesT::ess(normal_mix, method = "moment", sigma = sd_current)
RBesT::ess(normal_mix, method = "morita", sigma = sd_current)
RBesT::ess(normal_mix, method = "elir", sigma = sd_current)

##Propose a grid of candidate sample sizes
samplesize_num <- 6
samplesize_candidate <- c(seq(200, 600, length.out = samplesize_num), seq(80, 280, length.out = samplesize_num),
                          seq(10, 150, length.out = samplesize_num))
interval <- c(80, 40, 28)
hypothesis_testing <- matrix(NA, num_simulation, samplesize_num)
##Random seeds
set.seed(2021)
seeds_simulation <- matrix(sample(.Machine$integer.max, num_simulation*samplesize_num), 
                           nrow = num_simulation)
seeds_jags_simulation <- matrix(sample(.Machine$integer.max, num_simulation*num_chains*samplesize_num), 
                                nrow = num_simulation*num_chains)

##Data set to save power values
power_res <- data.frame(num_simulation = num_simulation, samplesize = samplesize_candidate, 
                        trteff = rep(trteff_candidate, each = samplesize_num), 
                        power = NA)
