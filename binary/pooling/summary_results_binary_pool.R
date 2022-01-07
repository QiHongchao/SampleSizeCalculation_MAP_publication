rm(list = ls())

##Set working directory
if (Sys.getenv("USERNAME") == "043712") {
  setwd("V:\\Users\\043712(H. Qi)\\Documents\\PhD projects\\Repo_publications\\SampleSizeCalculation_MAP_publication\\binary")
  .libPaths("V:\\Users\\043712(H. Qi)\\Documents\\Rlib")
} else {
  setwd("C:\\EMC\\Research\\PhD projects\\Repo_publications\\SampleSizeCalculation_MAP_publication\\binary")
}

##Data preparation
source("pooling/simulation_preparation_binary_pool.R")

##Package
library(ggplot2)
library(ggpubr)

##Combine the results
nb <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(OR_candidate)),
                 samplesize = samplesize_candidate, OR = rep(OR_candidate, each = samplesize_num),
                 power = NA, Approach = "No borrowing")
pool <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(OR_candidate)),
                  samplesize = samplesize_candidate, OR = rep(OR_candidate, each = samplesize_num),
                  power = NA, Approach = "Pool+balanced")
pool_ub <- data.frame(num_simulation = rep(num_simulation, samplesize_num*length(OR_candidate)),
                     samplesize = samplesize_candidate, OR = rep(OR_candidate, each = samplesize_num),
                     power = NA, Approach = "Pool+unbalanced")
for (i in 1:length(OR_candidate)) {
  nb$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./pooling/hypothesis_testing_binary_nb_", i, ".csv")), 2, mean)
  pool$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./pooling/hypothesis_testing_binary_pool_", i, ".csv")), 2, mean)
  pool_ub$power[((i-1)*samplesize_num + 1):(i*samplesize_num)] <- apply(read.csv(paste0("./pooling/hypothesis_testing_binary_pool_unbalanced_", i, ".csv")), 2, mean) 
}

res_binary <- rbind(nb, pool, pool_ub)

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
ggarrange(powercurves[[1]], powercurves[[2]], powercurves[[3]],
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

##Saved sample size
savedsamplesize <- data.frame(OR = OR_candidate, nb = NA, pool = NA, poolub = NA, savedpool = NA, savedpoolub = NA)
for (i in 1:length(OR_candidate)) {
  nb_OR <- nb[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  pool_OR <- pool[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  pool_ub_OR <- pool_ub[((i-1)*samplesize_num + 1):(i*samplesize_num),]
  ##NB
  samplesize_up_nb <- nb_OR$samplesize[min(which(nb_OR$power>=0.8))]
  power_low_nb <- nb_OR$power[min(which(nb_OR$power>=0.8))-1]
  power_up_nb <- nb_OR$power[min(which(nb_OR$power>=0.8))]
  samplesize_nb <- (0.8 - (power_up_nb - (power_up_nb - power_low_nb)/interval[i]*samplesize_up_nb))*interval[i]/(power_up_nb - power_low_nb)
  ##pool
  samplesize_up_pool <- pool_OR$samplesize[min(which(pool_OR$power>=0.8))]
  power_low_pool <- pool_OR$power[min(which(pool_OR$power>=0.8))-1]
  power_up_pool <- pool_OR$power[min(which(pool_OR$power>=0.8))]
  samplesize_pool <- (0.8 - (power_up_pool - (power_up_pool - power_low_pool)/interval[i]*samplesize_up_pool))*interval[i]/(power_up_pool - power_low_pool)
  ##pool UB
  samplesize_up_pool_ub <- pool_ub_OR$samplesize[min(which(pool_ub_OR$power>=0.8))]
  power_low_pool_ub <- pool_ub_OR$power[min(which(pool_ub_OR$power>=0.8))-1]
  power_up_pool_ub <- pool_ub_OR$power[min(which(pool_ub_OR$power>=0.8))]
  samplesize_pool_ub <- (0.8 - (power_up_pool_ub - (power_up_pool_ub - power_low_pool_ub)/interval[i]*samplesize_up_pool_ub))*interval[i]/(power_up_pool_ub - power_low_pool_ub)
  
  ##Saved sample size
  print(paste("OR:", OR_candidate[i], "saved sample size, balanced:", round(2*(samplesize_pool - samplesize_nb), 2),
              "saved sample size, unbalanced:", round(samplesize_pool_ub - samplesize_nb, 2)))
  savedsamplesize$nb[i] <- round(samplesize_nb, 2)
  savedsamplesize$pool[i] <- round(samplesize_pool, 2)
  savedsamplesize$poolub[i] <- round(samplesize_pool_ub, 2)
  savedsamplesize$savedpool[i] <- round(2*(samplesize_pool - samplesize_nb), 2)
  savedsamplesize$savedpoolub[i] <- round(samplesize_pool_ub - samplesize_nb, 2)
}
savedsamplesize
