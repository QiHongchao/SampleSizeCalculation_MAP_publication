rm(list = ls())

##Set working directory
if (Sys.getenv("USERNAME") == "043712") {
  setwd("V:\\Users\\043712(H. Qi)\\Documents\\PhD projects\\Repo_publications\\SampleSizeCalculation_MAP_publication\\binary")
  .libPaths("V:\\Users\\043712(H. Qi)\\Documents\\Rlib")
} else {
  setwd("C:\\EMC\\Research\\PhD projects\\Repo_publications\\SampleSizeCalculation_MAP_publication\\binary")
}

##Data preparation
source("simulation_preparation_binary.R")

##Package
library(ggplot2)
library(ggpubr)

##Combine the results
nb <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(OR_candidate)),
                 samplesize = samplesize_candidate, OR = rep(OR_candidate, each = samplesize_num),
                 power = NA, Approach = "No borrowing")
mac <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(OR_candidate)),
                  samplesize = samplesize_candidate, OR = rep(OR_candidate, each = samplesize_num),
                  power = NA, Approach = "MAP+balanced")
mac_ub <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(OR_candidate)),
                     samplesize = samplesize_candidate, OR = rep(OR_candidate, each = samplesize_num),
                     power = NA, Approach = "MAP+unbalanced")
for (i in 1:length(OR_candidate)) {
      nb$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./results/hypothesis_testing_binary_nb_", i, ".csv")), 2, mean)
      mac$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./results/hypothesis_testing_binary_mac_", i, ".csv")), 2, mean)
      mac_ub$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./results/hypothesis_testing_binary_mac_unbalanced_", i, ".csv")), 2, mean) 
}

res_binary <- rbind(nb, mac, mac_ub)

##Power curves
powercurves <- list()
for (i in 1:length(OR_candidate)) {
  powercurves[[i]] <- ggplot(data = res_binary[res_binary$OR==OR_candidate[i],], aes(x = samplesize, y = power*100, group = Approach, col = Approach)) + 
    geom_point() + geom_line() + geom_hline(yintercept=80, linetype="dashed", color="red") +
    labs(x = "Sample size", y = "Power (%)") + labs(tag = paste0("(", LETTERS[i], ")")) +
    scale_x_continuous(breaks=res_binary$samplesize[res_binary$OR==OR_candidate[i]]) +
    scale_y_continuous(breaks=seq(floor(min(res_binary$power)*100), ceiling(max(res_binary$power)*100), 4)) + 
    theme_classic() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("OR = ", OR_candidate[i]))
}
ggarrange(powercurves[[1]], powercurves[[2]], powercurves[[3]], powercurves[[4]],
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

##Saved sample size
savedsamplesize <- data.frame(OR = OR_candidate, nb = NA, map = NA, mapub = NA, savedmap = NA, savedmapub = NA)
for (i in 1:length(OR_candidate)) {
  nb_OR <- nb[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  map_OR <- mac[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  map_ub_OR <- mac_ub[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  ##NB
  samplesize_up_nb <- nb_OR$samplesize[min(which(nb_OR$power>=0.8))]
  power_low_nb <- nb_OR$power[min(which(nb_OR$power>=0.8))-1]
  power_up_nb <- nb_OR$power[min(which(nb_OR$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  ##MAP
  samplesize_up_map <- map_OR$samplesize[min(which(map_OR$power>=0.8))]
  power_low_map <- map_OR$power[min(which(map_OR$power>=0.8))-1]
  power_up_map <- map_OR$power[min(which(map_OR$power>=0.8))]
  samplesize_map <- (0.8 - (power_up_map - (power_up_map - power_low_map)/interval[i]*samplesize_up_map))*interval[i]/(power_up_map - power_low_map)
  ##MAP UB
  samplesize_up_map_ub <- map_ub_OR$samplesize[min(which(map_ub_OR$power>=0.8))]
  power_low_map_ub <- map_ub_OR$power[min(which(map_ub_OR$power>=0.8))-1]
  power_up_map_ub <- map_ub_OR$power[min(which(map_ub_OR$power>=0.8))]
  samplesize_map_ub <- (0.8 - (power_up_map_ub - (power_up_map_ub - power_low_map_ub)/interval[i]*samplesize_up_map_ub))*interval[i]/(power_up_map_ub - power_low_map_ub)
  
  ##Saved sample size
  print(paste("OR:", OR_candidate[i], "saved sample size, balanced:", round(2*(samplesize_map - samplesize_nb), 2),
              "saved sample size, unbalanced:", round(samplesize_map_ub - samplesize_nb, 2)))
  savedsamplesize$nb[i] <- round(samplesize_nb, 2)
  savedsamplesize$map[i] <- round(samplesize_map, 2)
  savedsamplesize$mapub[i] <- round(samplesize_map_ub, 2)
  savedsamplesize$savedmap[i] <- round(2*(samplesize_map - samplesize_nb), 2)
  savedsamplesize$savedmapub[i] <- round(samplesize_map_ub - samplesize_nb, 2)
}
savedsamplesize
