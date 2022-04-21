rm(list=ls())

##Simulation preparation
source("simulation_preparation_binary.R")

start <- Sys.time()
for (k in 1:length(OR_candidate)) {
  OR <- OR_candidate[k]
for (j in 1:samplesize_num) {
  ##Specify sample size per arm
  samplesize <- samplesize_candidate[(k-1)*samplesize_num+j]
  
  ##Simulate a new trial data set; analyze it using the NB approach; hypothesis testing; do the above many times
  for (i in 1:num_simulation) {
    set.seed(seeds_simulation[i, j])
    p0_new_sim <- sample(map_sample$p0_new, 1)
    p1_new_sim <- OR*p0_new_sim/(1+(OR-1)*p0_new_sim)
    event_control <- rbinom(1, samplesize, p0_new_sim)
    event_treatment <- rbinom(1, samplesize, p1_new_sim)
    
    ##No borrowing analysis
    jags_data_nb <- list(nsub_control_new = samplesize, nsub_treatment_new = samplesize,
                          event_control_new = event_control, event_treatment_new = event_treatment)
    
    nb_jags_model <- jags.model(file = "./jags/nb_binary.txt", data = jags_data_nb,
                                 n.chains = num_chains, n.adapt = 1000, quiet = T,
                                 inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+1, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+2, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+3, j]),
                                              list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[i*num_chains, j])))

    ##Burnin stage
    update(nb_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")

    ##Sampling after burnin
    params <- c("p0", "p1", "diff")
    nb_jags <- coda.samples(nb_jags_model, params, n.iter = num_iter * (1 - perc_burnin),
                             progress.bar = "none")
    nb_sample <- as.data.frame(do.call(rbind, nb_jags))
    
    ##Hypothesis testing
    hypothesis_testing[i, j] <- (quantile(nb_sample$diff, 0.025) > 0 | quantile(nb_sample$diff, 0.975) < 0)
    
    ##Progress monitoring
    if (i%%(num_simulation/10) == 0) {
      print(paste0("OR: ", OR, " sample size: ", samplesize, " simulation: ", i))
    }
  }
  ##Calculate the power based on the specified sample size
  power_res$power[(k-1)*samplesize_num+j] <- mean(hypothesis_testing[,j])
}
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_binary_nb_", k, ".csv"), row.names = F)
}
Sys.time() - start

##Sample size for the treatment arm, calculate the exact sample size that achieves 80% power using linear interpolation
for (i in 1:length(OR_candidate)) {
  res <- read.csv(paste0("./results/hypothesis_testing_binary_nb_", i, ".csv"))
  power_res$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(res, 2, mean)
  nb <- power_res[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  samplesize_up_nb <- nb$samplesize[min(which(nb$power>=0.8))]
  power_low_nb <- nb$power[min(which(nb$power>=0.8))-1]
  power_up_nb <- nb$power[min(which(nb$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  print(paste0("OR: ", OR_candidate[i], ", sample size for the treatment arm: ", ceiling(samplesize_nb)))
}
write.csv(power_res, "./results/power_res_binary_nb.csv", row.names = F)
