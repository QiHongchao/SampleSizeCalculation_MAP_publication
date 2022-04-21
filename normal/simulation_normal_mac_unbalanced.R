##The script must be executed after the execution of "simulation_normal_nb.R"
rm(list=ls())

##Simulation preparation
source("simulation_preparation_normal.R")


##The purpose of the unbalanced trial is to make sure that the sample size in the treatment arm is not reduced
power_res_nb <- read.csv("./results/power_res_nb.csv")
samplesize_treatment_candidate <- NULL
##Sample size for the treatment arm
for (i in 1:length(trteff_candidate)) {
  nb <- power_res_nb[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  ##NB
  samplesize_up_nb <- nb$samplesize[min(which(nb$power>=0.8))]
  power_low_nb <- nb$power[min(which(nb$power>=0.8))-1]
  power_up_nb <- nb$power[min(which(nb$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  print(paste0("trteff: ", trteff_candidate[i], ", sample size for the treatment arm: ", ceiling(samplesize_nb)))
  samplesize_treatment_candidate[i] <- ceiling(samplesize_nb)
}

start <- Sys.time()
for (k in 1:length(trteff_candidate)) {
  trteff <- trteff_candidate[k]
  samplesize_treatment <- samplesize_treatment_candidate[k]
  for (j in 1:samplesize_num) {
    samplesize_control <- samplesize_candidate[(k-1)*samplesize_num+j]
    
    ##Simulate a new trial data set; analyze it using the MAC approach; hypothesis testing; do the above many times
    for (i in 1:num_simulation) {
      set.seed(seeds_simulation[i, j])
      control_mean <- sample(bma_random$map_sample$theta_new, 1)
      control_data <- rnorm(samplesize_control, control_mean, sd_current)
      treatment_data <- rnorm(samplesize_treatment, control_mean+trteff, sd_current)
      
      ##MAC data
      jags_data_mac <- list(J = J, estimate = estimate_historical, se = se_historical, sd = sd, 
                            mean_control = mean(control_data), tau_new = samplesize_control/var(control_data),
                            mean_treatment = mean(treatment_data), tau1_new = samplesize_treatment/var(treatment_data))
      
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
      params <- c("theta_new", "trteff")
      mac_jags <- coda.samples(mac_jags_model, params, n.iter = num_iter*(1 - perc_burnin), 
                               progress.bar = "none")
      mac_sample <- as.data.frame(do.call(rbind, mac_jags))
      
      ##Hypothesis testing
      hypothesis_testing[i, j] <- (quantile(mac_sample$trteff, 0.025) > 0 | quantile(mac_sample$trteff, 0.975) < 0)
      
      ##Progress monitoring
      if (i%%(num_simulation/10) == 0) {
        print(paste0("treatment effect: ", trteff, ", sample size: ", samplesize_control, ", simulation: ", i))
      }
    }
    ##Calculate the power based on the chosen sample size
    power_res$power[(k-1)*samplesize_num+j] <- mean(hypothesis_testing[,j])
  }
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_mac_unbalanced_", k, ".csv"), row.names = F)
}
Sys.time() - start
