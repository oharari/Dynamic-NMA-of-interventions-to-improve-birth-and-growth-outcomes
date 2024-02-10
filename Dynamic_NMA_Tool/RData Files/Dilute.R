wd = "C:/Users/oharari/Desktop/RData Files"
setwd(wd)

#outcome = "0-6_Stunting"
#outcome = "LAZ"
outcome = "HAZ"

fname = paste0(outcome, "_NMA_Object.RData")
load(fname)

monkey = 1

if(monkey){
  M = my_NMA2$meta.sim$sims.matrix
  my_NMA2$meta.sim$sims.matrix = M[seq(1, nrow(M), 10), ]
  
  my_NMA2$meta.sim$sims.list = NULL
  my_NMA2$meta.sim$sims.array = NULL
  my_NMA2$meta.sim$mean = my_NMA2$meta.sim$median = my_NMA2$meta.sim$sd = NULL
  my_NMA2$meta.sim$last.values = NULL
  
  save(my_NMA2, file = fname)
}


# my_NMA = my_NMA2
# baseline = my_NMA$trt_names[1]
# events_are_bad = T
# base_vs_all = T
# M = post_superiority_samples(my_NMA, baseline, events_are_bad, base_vs_all)
