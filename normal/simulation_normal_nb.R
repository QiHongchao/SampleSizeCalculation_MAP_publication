rm(list=ls())

##Simulation preparation
source("simulation_preparation_normal.R")

start <- Sys.time()
for (k in 1:length(trteff_candidate)) {
  trteff <- trteff_candidate[k]
  for (j in 1:samplesize_num) {
    samplesize <- samplesize_candidate[(k-1)*samplesize_num+j]
    
    ##Simulate a new trial data set; analyze it using the NB approach; hypothesis testing; do the above many times
    for (i in 1:num_simulation) {
      set.seed(seeds_simulation[i, j])
      control_mean <- sample(bma_random$map_sample$theta_new, 1)
      control_data <- rnorm(samplesize, control_mean, sd_current)
      treatment_data <- rnorm(samplesize, control_mean+trteff, sd_current)
      
      ##NB data
      jags_data_nb <- list(mean_control = mean(control_data), tau_new = samplesize/var(control_data),
                            mean_treatment = mean(treatment_data), tau1_new = samplesize/var(treatment_data))
      
      ##NB
      nb_jags_model <- jags.model(file = "./jags/nb_normal.txt", data = jags_data_nb, 
                                   n.chains = num_chains, n.adapt = 500, quiet = T,
                                   inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+1, j]),
                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+2, j]),
                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[(i-1)*num_chains+3, j]),
                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seeds_jags_simulation[i*num_chains, j])))
      
      ##Burnin stage
      update(nb_jags_model, n.iter = num_iter*perc_burnin, progress.bar = "none")
      
      ##Sampling after burnin
      params <- c("theta_new", "theta1_new", "trteff")
      nb_jags <- coda.samples(nb_jags_model, params, n.iter = num_iter*(1 - perc_burnin), 
                               progress.bar = "none")
      nb_sample <- as.data.frame(do.call(rbind, nb_jags))
      
      ##Hypothesis testing
      hypothesis_testing[i, j] <- (quantile(nb_sample$trteff, 0.025) > 0 | quantile(nb_sample$trteff, 0.975) < 0)
      
      ##Progress
      if (i%%(num_simulation/10) == 0) {
        print(paste0("trteff: ", trteff, " sample size: ", samplesize, " simulation: ", i))
      }
    }
  }
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_nb_", k, ".csv"), row.names = F)
}
Sys.time() - start

##Required sample size
for (i in 1:length(trteff_candidate)) {
  power_res$power[power_res$trteff==trteff_candidate[i]] <- apply(read.csv(paste0("./results/hypothesis_testing_nb_", i, ".csv")), 2, mean)
  nb <- power_res[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  ##NB
  samplesize_up_nb <- nb$samplesize[min(which(nb$power>=0.8))]
  power_low_nb <- nb$power[min(which(nb$power>=0.8))-1]
  power_up_nb <- nb$power[min(which(nb$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  print(paste0("trteff: ", trteff_candidate[i], ", sample size for the treatment arm: ", ceiling(samplesize_nb)))
}
write.csv(power_res, "./results/power_res_nb.csv", row.names = F)
