rm(list=ls())

##Simulation preparation
source("simulation_preparation_normal.R")

start <- Sys.time()
for (k in 1:length(trteff_candidate)) {
  trteff <- trteff_candidate[k]
  for (j in 1:samplesize_num) {
    samplesize <- samplesize_candidate[(k-1)*samplesize_num+j]
    
    ##Simulate a new trial data set; analyze it using the MAC approach; hypothesis testing; do the above many times
    for (i in 1:num_simulation) {
      set.seed(seeds_simulation[i, j])
      control_mean <- sample(bma_random$map_sample$theta_new, 1)
      control_data <- rnorm(samplesize, control_mean, sd_current)
      treatment_data <- rnorm(samplesize, control_mean+trteff, sd_current)
      
      ##MAC data
      jags_data_mac <- list(J = J, estimate = estimate_historical, se = se_historical, sd = sd,
                            mean_control = mean(control_data), tau_new = samplesize/var(control_data),
                            mean_treatment = mean(treatment_data), tau1_new = samplesize/var(treatment_data))
      
      ##MAC
      mac_jags_model <- jags.model(file = "./jags/mac_normal.txt", data = jags_data_mac, 
                                   n.chains = num_chains, n.adapt = 500, quiet = T,
                                   inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+1, j]),
                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+2, j]),
                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+3, j]),
                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[i*num_chains, j])))
      
      ##Burnin stage
      update(mac_jags_model, n.iter = num_iter*perc_burnin, progress.bar = "none")
      
      ##Sampling after burnin
      params <- c("theta_new", "theta1_new", "trteff", "sigma")
      mac_jags <- coda.samples(mac_jags_model, params, n.iter = num_iter*(1 - perc_burnin), 
                               progress.bar = "none")
      mac_sample <- as.data.frame(do.call(rbind, mac_jags))
      
      ##Hypothesis testing
      hypothesis_testing[i, j] <- (quantile(mac_sample$trteff, 0.025) > 0 | quantile(mac_sample$trteff, 0.975) < 0)
      
      ##Progress monitoring
      if (i%%(num_simulation/10) == 0) {
        print(paste0("treatment effect: ", trteff, ", sample size: ", samplesize, ", simulation: ", i))
      }
    }
    ##Calculate the power based on the chosen sample size
    power_res$power[(k-1)*samplesize_num+j] <- mean(hypothesis_testing[,j])
  }
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_mac_", k, ".csv"), row.names = F)
}
Sys.time() - start
