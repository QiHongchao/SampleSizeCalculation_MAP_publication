##The script should be executed after "samplesize_map_binary_nb.R"
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

##The purpose of the unbalanced trial is to make sure that the sample size in the treatment arm is not reduced
power_res_nb <- read.csv("./results/power_res_binary_nb.csv")
samplesize_treatment_candidate <- NULL
##Sample size for the treatment arm
for (i in 1:length(OR_candidate)) {
  nb <- power_res_nb[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  ##NB
  samplesize_up_nb <- nb$samplesize[min(which(nb$power>=0.8))]
  power_low_nb <- nb$power[min(which(nb$power>=0.8))-1]
  power_up_nb <- nb$power[min(which(nb$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  print(paste0("trteff: ", OR_candidate[i], " sampleszie for the treatment arm: ", ceiling(samplesize_nb)))
  samplesize_treatment_candidate[i] <- ceiling(samplesize_nb)
}

start <- Sys.time()
for (k in 1:length(OR_candidate)) {
  OR <- OR_candidate[k]
  ##prespecify the sample size of the treatment arm
  samplesize_treatment <- samplesize_treatment_candidate[k]
for (j in 1:samplesize_num) {
  ##choose the sample size of the control arm
  samplesize_control <- samplesize_candidate[(k-1)*samplesize_num+j]
  
  ##step 2 - step 5: sample a new trial data set; analyze it using the MAC approach; hypothesis testing; do the above N times
  for (i in 1:num_simulation) {
    set.seed(seeds_simulation[i, j])
    p0_new_sim <- sample(map_sample$p0_new, 1)
    p1_new_sim <- OR*p0_new_sim/(1+(OR-1)*p0_new_sim)
    event_control <- rbinom(1, samplesize_control, p0_new_sim)
    event_treatment <- rbinom(1, samplesize_treatment, p1_new_sim)
    
    ##MAC analysis
    jags_data_mac <- list(J = J, prec_between_sd = prec_between_sd,
                          nsub_historical = nsub_historical, event_historical = event_historical,
                          nsub_control_new = samplesize_control, nsub_treatment_new = samplesize_treatment,
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
    
    hypothesis_testing[i, j] <- (quantile(mac_sample$diff, 0.025) > 0 | quantile(mac_sample$diff, 0.975) < 0)
    
    ##Progress
    if (i%%(num_simulation/2) == 0) {
      print(paste0("OR: ", OR, " sample size: ", samplesize_control, " simulation: ", i))
    }
  }
  ##step 6: calculate the power based on the chosen sample size
  power_res$power[(k-1)*samplesize_num+j] <- mean(hypothesis_testing[,j])
}
  ##Hypothesis testing for different candidate sample sizes
  write.csv(hypothesis_testing, paste0("./results/hypothesis_testing_binary_mac_unbalanced_", k, ".csv"), row.names = F)
}
Sys.time() - start
