rm(list=ls())

##Set working directory
if (Sys.getenv("USERNAME") == "043712") {
  setwd("V:\\Users\\043712(H. Qi)\\Documents\\PhD projects\\Repo_publications\\SampleSizeCalculation_MAP_publication\\binary")
  .libPaths("V:\\Users\\043712(H. Qi)\\Documents\\Rlib")
} else {
  setwd("C:\\EMC\\Research\\PhD projects\\Repo_publications\\SampleSizeCalculation_MAP_publication\\binary")
}

##Simulation preparation, step 1 conducted in the file
source("simulation_preparation_binary.R")

start <- Sys.time()
for (k in 1:length(OR_candidate)) {
  OR <- OR_candidate[k]
for (j in 1:samplesize_num) {
  ##choose sample size per arm
  samplesize <- samplesize_candidate[(k-1)*samplesize_num+j]
  
  ##step 2 - step 5: sample a new trial data set; analyze it using the pooling approach; hypothesis testing; do the above N times
  for (i in 1:num_simulation) {
    set.seed(seeds_simulation[i, j])
    ##sample p0_new from the MAP prior
    p0_new_sim <- sample(map_sample$p0_new, 1)
    p1_new_sim <- OR*p0_new_sim/(1+(OR-1)*p0_new_sim)
    event_control <- rbinom(1, samplesize, p0_new_sim)
    event_treatment <- rbinom(1, samplesize, p1_new_sim)
    
    ##Pooling analysis, pool the historical controls and the current control
    jags_data_pool <- list(nsub_control_all = sum(nsub_historical)+samplesize, 
                           event_control_all = sum(event_historical)+event_control, 
                           nsub_treatment = samplesize,
                           event_treatment = event_treatment)
    
    pool_jags_model <- jags.model(file = "./jags/pool_binary.txt", data = jags_data_pool, 
                                 n.chains = num_chains, n.adapt = 1000, quiet = T,
                                 inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+1, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+2, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+3, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[i*num_chains, j])))
    
    ##Burnin stage
    update(pool_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("p0_new", "p1", "diff")
    
    pool_jags <- coda.samples(pool_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                             progress.bar = "none")
    
    pool_sample <- as.data.frame(do.call(rbind, pool_jags))
    
    hypothesis_testing[i, j] <- (quantile(pool_sample$diff, 0.025) > 0 | quantile(pool_sample$diff, 0.975) < 0)
    
    ##Progress
    if (i%%(num_simulation/2) == 0) {
      print(paste0("OR: ", OR, " sample size: ", samplesize, " simulation: ", i))
    }
  }
  ##step 6: calculate the power based on the chosen sample size
  power_res$power[(k-1)*samplesize_num+j] <- mean(hypothesis_testing[,j])
}
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_binary_pool_", k, ".csv"), row.names = F)
}
Sys.time() - start
