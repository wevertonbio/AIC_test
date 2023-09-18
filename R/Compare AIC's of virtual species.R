#### Test AIC's of virtual species ####
#Load the functions...

#Load packages
library(dplyr) #Sorry...
library(future.apply)

#### Linear ####
#Load data
load("C:/Users/wever/Desktop/AIC_test/Data/data_linear.RData")

#Get grid with formulas
variables <- colnames(data_l1)[2:5]
var_comb <- kuenm::all_var_comb(variables, 2)
formula_x <- lapply(seq_along(var_comb), function(i){
  f_l <- get_formulas_maxnet(independent = var_comb[[1]], type = "l")
  f_q <- get_formulas_maxnet(independent = var_comb[[i]], type = "q")
  f_lq <- get_formulas_maxnet(independent = var_comb[[i]], type = "lq")
  f_lqp <- get_formulas_maxnet(independent = var_comb[[i]], type = "lqp")
  f_final <- c(f_l, f_q, f_lq, f_lqp)
  names(f_final) <- c("l", "q", "lq", "lqp")
  return(f_final)
}) %>% unlist()
#Store formulas in dataframe
formula_d <- data.frame(formula = paste(formula_x, -1),
                        features = names(formula_x))
#Expand grid
regm <- c(0.1, 0.5, 1, 2, 3, 5)
f_grid <- expand.grid(formula_d$formula, regm)
f_grid$features <- formula_d$features
colnames(f_grid) <- c("Formulas", "regm", "Features")
#Factor to character
f_grid$Formulas <- as.character(f_grid$Formulas)

#Data l1
m_l1 <- compare_aic(data = data_l1, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
dir.create("../ResultsAIC")
write.csv(m_l1, "../ResultsAIC/m_l1.csv", row.names = F)

#Data l2
m_l2 <- compare_aic(data = data_l2, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
write.csv(m_l2, "../ResultsAIC/m_l2.csv", row.names = F)

#Data l3
m_l3 <- compare_aic(data = data_l3, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
write.csv(m_l3, "../ResultsAIC/m_l3.csv", row.names = F)

#Data l4
m_l4 <- compare_aic(data = data_l4, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
write.csv(m_l4, "../ResultsAIC/m_l4.csv", row.names = F)

####Get summary of statistics####
#Import results
res_l <- list.files("../ResultsAIC", pattern = "m_l", full.names = T)
res_l <- lapply(seq_along(res_l), function(i){
  res_i <- read.csv(res_l[i])
  sum_i <- eval_stats(res_i)
  return(sum_i)
})

#Explore Plot ####
plot(res[[1]]$npar, res[[1]]$AIC)
plot(res[[1]]$npar, res[[1]]$AIC_warren)
plot(res[[2]]$npar, res[[2]]$AIC)
plot(res[[2]]$npar, res[[2]]$AIC_warren)
plot(res[[3]]$npar, res[[3]]$AIC)
plot(res[[3]]$npar, res[[3]]$AIC_warren)
plot(res[[4]]$npar, res[[4]]$AIC)
plot(res[[4]]$npar, res[[4]]$AIC_warren)

plot(res[[1]]$Omission_rate_at_5.mean, res[[1]]$AIC)
plot(res[[1]]$Omission_rate_at_5.mean, res[[1]]$AIC_warren)

plot(res[[2]]$Omission_rate_at_5.mean, res[[2]]$AIC)
plot(res[[2]]$Omission_rate_at_5.mean, res[[2]]$AIC_warren)

plot(res[[3]]$Omission_rate_at_5.mean, res[[3]]$AIC)
plot(res[[3]]$Omission_rate_at_5.mean, res[[3]]$AIC_warren)

plot(res[[4]]$Omission_rate_at_5.mean, res[[4]]$AIC)
plot(res[[4]]$Omission_rate_at_5.mean, res[[4]]$AIC_warren)


res[[1]] %>% filter(AIC == min(AIC))
res[[1]] %>% filter(AIC_warren == min(AIC_warren))

res[[2]] %>% filter(AIC == min(AIC))
res[[2]] %>% filter(AIC_warren == min(AIC_warren))

res[[3]] %>% filter(AIC == min(AIC))
res[[3]] %>% filter(AIC_warren == min(AIC_warren))

res[[4]] %>% filter(AIC == min(AIC))
res[[4]] %>% filter(AIC_warren == min(AIC_warren))

hist(data_l1[, 2])

b <- data_l1 %>% filter(pres_back == 1) %>% dplyr::select(2)
hist(b$bio_1)


#### Quadratic ####
#Load data
load("C:/Users/wever/Desktop/AIC_test/Data/data_quadratic.RData")

#Data q2
m_q2 <- compare_aic(data = data_q2, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
write.csv(m_q2, "../ResultsAIC/m_q2.csv", row.names = F)

#Data q3
m_q3 <- compare_aic(data = data_q3, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
write.csv(m_q3, "../ResultsAIC/m_q3.csv", row.names = F)

#Data q4
m_q4 <- compare_aic(data = data_q4, pr_bg = "pres_back",
                    formula_grid = f_grid, n_folds = 4, parallel = T,
                    ncores = 8, progress_bar = T)
#Save results
write.csv(m_q4, "../ResultsAIC/m_q4.csv", row.names = F)

####Get summary of statistics####
#Import results
res_q <- list.files("../ResultsAIC", pattern = "m_q", full.names = T)
res_q <- lapply(seq_along(res_q), function(i){
  res_i <- read.csv(res_q[i])
  sum_i <- eval_stats(res_i)
  return(sum_i)
})

res_q[[1]] %>% filter(AIC == min(AIC, na.rm = T))
res_q[[1]] %>% filter(AIC_warren == min(AIC_warren, na.rm = T))

res_q[[2]] %>% filter(AIC == min(AIC, na.rm = T))
res_q[[2]] %>% filter(AIC_warren == min(AIC_warren, na.rm = T))

res_q[[3]] %>% filter(AIC == min(AIC, na.rm = T))
res_q[[3]] %>% filter(AIC_warren == min(AIC_warren, na.rm = T))


par(mfrow = c(1, 2))
plot(res_q[[3]]$npar, res_q[[3]]$AIC)
plot(res_q[[3]]$npar, res_q[[3]]$AIC_warren)

cov(data_q4[data_q4$pres_back == 1, -1])
