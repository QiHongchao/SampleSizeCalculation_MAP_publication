rm(list = ls())

##Data preparation
source("simulation_preparation_normal.R")

##Combine the results
nb <- read.csv("./results/power_res_nb.csv")
nb$Approach <- "No borrowing"
mac <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(trteff_candidate)),
                  samplesize = samplesize_candidate, trteff = rep(trteff_candidate, each = samplesize_num),
                  power = NA, Approach = "MAP+balanced")
mac_ub <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(trteff_candidate)),
                     samplesize = samplesize_candidate, trteff = rep(trteff_candidate, each = samplesize_num),
                     power = NA, Approach = "MAP+unbalanced")
for (i in 1:length(trteff_candidate)) {
  mac$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./results/hypothesis_testing_mac_", i, ".csv")), 2, mean)
  mac_ub$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./results/hypothesis_testing_mac_unbalanced_", i, ".csv")), 2, mean)
}

res_normal <- rbind(nb, mac, mac_ub)

##Power curves
powercurves <- list()
for (i in 1:length(trteff_candidate)) {
  powercurves[[i]] <- ggplot(data = res_normal[res_normal$trteff==trteff_candidate[i],], 
                             aes(x = samplesize, y = power*100, group = Approach, 
                                 col = Approach, shape = Approach, linetype = Approach)) + 
    geom_point() + geom_line() + geom_hline(yintercept=80, linetype="dashed", color="red") +
    labs(x = "Sample size", y = "Power (%)") + labs(tag = paste0("(", LETTERS[i], ")")) +
    scale_x_continuous(breaks=res_normal$samplesize[res_normal$trteff==trteff_candidate[i]]) +
    scale_y_continuous(breaks=seq(0, 100, 10)) + 
    expand_limits(y=c(0, 100)) + 
    theme_classic() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    ggtitle(paste0("Treatment effect = ", trteff_candidate[i]))
}
ggarrange(powercurves[[1]], powercurves[[2]], powercurves[[3]],
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

##Required sample sizes and saved sample sizes
savedsamplesize <- data.frame(trteff = trteff_candidate, nb = NA, map = NA, mapub = NA, savedmap = NA, savedmapub = NA)
for (i in 1:length(trteff_candidate)) {
  nb_trteff <- nb[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  map_trteff <- mac[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  map_ub_trteff <- mac_ub[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  ##NB
  samplesize_up_nb <- nb_trteff$samplesize[min(which(nb_trteff$power>=0.8))]
  power_low_nb <- nb_trteff$power[min(which(nb_trteff$power>=0.8))-1]
  power_up_nb <- nb_trteff$power[min(which(nb_trteff$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  ##MAP
  samplesize_up_map <- map_trteff$samplesize[min(which(map_trteff$power>=0.8))]
  power_low_map <- map_trteff$power[min(which(map_trteff$power>=0.8))-1]
  power_up_map <- map_trteff$power[min(which(map_trteff$power>=0.8))]
  samplesize_map <- (0.8 - (power_up_map - (power_up_map - power_low_map)/interval[i]*samplesize_up_map))*interval[i]/(power_up_map - power_low_map)
  ##MAP UB
  samplesize_up_map_ub <- map_ub_trteff$samplesize[min(which(map_ub_trteff$power>=0.8))]
  power_low_map_ub <- map_ub_trteff$power[min(which(map_ub_trteff$power>=0.8))-1]
  power_up_map_ub <- map_ub_trteff$power[min(which(map_ub_trteff$power>=0.8))]
  samplesize_map_ub <- (0.8 - (power_up_map_ub - (power_up_map_ub - power_low_map_ub)/interval[i]*samplesize_up_map_ub))*interval[i]/(power_up_map_ub - power_low_map_ub)
  ##Data set
  savedsamplesize$nb[i] <- 2*ceiling(samplesize_nb)
  savedsamplesize$map[i] <- 2*ceiling(samplesize_map)
  savedsamplesize$mapub[i] <- ceiling(samplesize_map_ub) + ceiling(samplesize_nb)
  savedsamplesize$savedmap[i] <- savedsamplesize$nb[i] - savedsamplesize$map[i]
  savedsamplesize$savedmapub[i] <- savedsamplesize$nb[i] - savedsamplesize$mapub[i]
}
savedsamplesize
