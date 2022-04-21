rm(list=ls())

##Simulation preparation
source("simulation_preparation_binary.R")

start <- Sys.time()
for (k in 1:length(OR_candidate)) {
  OR <- OR_candidate[k]
for (j in 1:samplesize_num) {
  ##Specify sample size per arm
  samplesize <- samplesize_candidate[(k-1)*samplesize_num+j]
  
  ##Simulate a new trial data set; analyze it using the MAP approach; hypothesis testing; do the above many times
  for (i in 1:num_simulation) {
    set.seed(seeds_simulation[i, j])
    p0_new_sim <- sample(map_sample$p0_new, 1)
    p1_new_sim <- OR*p0_new_sim/(1+(OR-1)*p0_new_sim)
    event_control <- rbinom(1, samplesize, p0_new_sim)
    event_treatment <- rbinom(1, samplesize, p1_new_sim)
    
    ##MAC analysis
    jags_data_mac <- list(J = J, prec_between_sd = prec_between_sd,
                          nsub_historical = nsub_historical, event_historical = event_historical,
                          nsub_control_new = samplesize, nsub_treatment_new = samplesize,
                          event_control_new = event_control, event_treatment_new = event_treatment)
    
    mac_jags_model <- jags.model(file = "./jags/mac_binary.txt", data = jags_data_mac, 
                                 n.chains = num_chains, n.adapt = 1000, quiet = T,
                                 inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+1, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+2, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+3, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[i*num_chains, j])))
    
    ##Burnin stage
    update(mac_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("p0_new", "p1", "diff")
    mac_jags <- coda.samples(mac_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                             progress.bar = "none")
    mac_sample <- as.data.frame(do.call(rbind, mac_jags))
    
    ##Hypothesis testing
    hypothesis_testing[i, j] <- (quantile(mac_sample$diff, 0.025) > 0 | quantile(mac_sample$diff, 0.975) < 0)
    
    ##Progress monitoring
    if (i%%(num_simulation/10) == 0) {
      print(paste0("OR: ", OR, " sample size: ", samplesize, " simulation: ", i))
    }
  }
  ##step 6: calculate the power based on the specified sample size
  power_res$power[(k-1)*samplesize_num+j] <- mean(hypothesis_testing[,j])
}
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_binary_mac_", k, ".csv"), row.names = F)
}
Sys.time() - start
