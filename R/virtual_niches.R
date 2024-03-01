library(geodata)
library(evniche)

# function for niches with linear responses
linear_respone <- function(variables, slopes) {
  l <- sapply(1:ncol(variables), function(x) {
    variables[, x] * slopes[x]
  })
  if (ncol(variables) != 1) {
    l <- apply(l, 1, sum)
  } else {
    l <- l[, 1]
  }
  return((l - min(l)) / (max(l) - min(l)))
}

# directory for Data
dir.create("Data")


# downloading environmental variables
var <- worldclim_global(var = "bio", res = 10, path = "Data")

## renaming variables
names(var) <- gsub("wc2.1_10m_", "", names(var))
names(var)

## selecting a subset of variables
selected <- c(1, 7, 12, 15)

var <- var[[selected]]

# uniform distribution of environments according to ranges of environments
## ranges
ranges <- minmax(var)

## sets of environments
set.seed(1)
data_env <- sapply(1:ncol(ranges), function (x) {
  runif(100000, min = ranges[1, x], max = ranges[2, x])
})

colnames(data_env) <- colnames(ranges)

head(data_env)


# presences denoting distinct types of response to variables
## descriptive stats of example data (background)
summary(data_env)

## sample of background
set.seed(1)
background <- data_env[sample(100000, 10000), ]

#background <- data_l1[data_l1$pres_back == 0, -1]
save(background, file = "Data/background.RData")

## combinations
comb <- list(v1_2 = paste0("bio_", c(1, 12)),
             v1_3 = paste0("bio_", c(1, 7, 12)),
             v1_4 = paste0("bio_", c(1, 7, 12, 15)))

### ranges for niches
qranges <- sapply(1:ncol(data_env), function(x) {
  quantile(data_env[, x])[c(2, 4)]
})

colnames(qranges) <- colnames(ranges)

## uniform responses of variables
set.seed(20)
bio7 <- runif(1000, ranges[1, 2], ranges[2, 2])
bio12 <- runif(1000, ranges[1, 3], ranges[2, 3])
bio15 <- runif(1000, ranges[1, 4], ranges[2, 4])

save(bio7, bio12, bio15, data_env, file = "Data/example_data.RData")

## niches with linear responses
### defining probabilities for niches with linear responses
lin_1var <- linear_respone(variables = matrix(data_env[, 1], ncol = 1),
                           slopes = 3)
lin_2var <- linear_respone(variables = data_env[, comb[[1]]],
                           slopes = c(500, 5))
lin_3var <- linear_respone(variables = data_env[, comb[[2]]],
                           slopes = c(500, 1000, 5))
lin_4var <- linear_respone(variables = data_env[, comb[[3]]],
                           slopes = c(500, 1000, 5, 300))

plot(data_env[, 1], lin_1var)

### rows to select data with increasing probabilities towards the high end
set.seed(12)
bl1 <- sort(sample(nrow(data_env), 1000, prob = lin_1var))
bl2 <- sort(sample(nrow(data_env), 1000, prob = lin_2var))
bl3 <- sort(sample(nrow(data_env), 1000, prob = lin_3var))
bl4 <- sort(sample(nrow(data_env), 1000, prob = lin_4var))

### data with linear responses
data_l1 <- cbind(bio_1 = data_env[bl1, 1],
                 bio_7 = bio7, bio_12 = bio12, bio_15 = bio15)
data_l2 <- cbind(bio_1 = data_env[bl2, 1], bio_7 = bio7,
                 bio_12 = data_env[bl2, 3], bio_15 = bio15)
data_l3 <- cbind(bio_1 = data_env[bl3, 1], bio_7 = data_env[bl3, 2],
                 bio_12 = data_env[bl3, 3], bio_15 = bio15)
data_l4 <- cbind(bio_1 = data_env[bl4, 1], bio_7 = data_env[bl4, 2],
                 bio_12 = data_env[bl4, 3], bio_15 = data_env[bl4, 4])

### preparing data for models (linear responses)
data_l1 <- rbind(data.frame(pres_back = 1, data_l1),
                 data.frame(pres_back = 0, background))

data_l2 <- rbind(data.frame(pres_back = 1, data_l2),
                 data.frame(pres_back = 0, background))

data_l3 <- rbind(data.frame(pres_back = 1, data_l3),
                 data.frame(pres_back = 0, background))

data_l4 <- rbind(data.frame(pres_back = 1, data_l4),
                 data.frame(pres_back = 0, background))

### save the data
save(data_l1, data_l2, data_l3, data_l4, file = "Data/data_linear.RData")

## niches with quadratic responses (linear as well)
### virtual niches
v_niches <- lapply(1:length(comb), function(x) {
  ## variances
  vars <- evniche:::var_from_range(range = qranges[, comb[[x]]])

  ## variance-covariance matrix
  cov <- 0
  varcov <- evniche:::var_cov_matrix(variances = vars, covariances = cov)

  ## centroid
  cent <- evniche:::centroid(range = qranges[, comb[[x]]])

  ## ellipsoid characteristics (virtual niche)
  ell_features(centroid = cent, covariance_matrix = varcov, level = 0.99)
})

names(v_niches) <- names(comb)

### generating data from niches and background
vd_pre_data <- lapply(1:length(v_niches), function(x) {
  set.seed(1)
  virtual_data(features = v_niches[[x]], from = "ellipsoid", n = 1000)
})

names(vd_pre_data) <- names(v_niches)

### complete missing columns
vd_pre_data[[1]] <- cbind(bio_1 = vd_pre_data[[1]][, 1], bio_7 = bio7,
                          bio_12 = vd_pre_data[[1]][, 2], bio_15 = bio15)
vd_pre_data[[2]] <- cbind(vd_pre_data[[2]], bio_15 = bio15)

### quadratic for one variable
quad_1var <- cbind(bio_1 = vd_pre_data[[1]][, 1], bio_7 = bio7, bio_12 = bio12,
                   bio_15 = bio15)

### preparing data for models (linear and quadratic responses)
data_q1 <- rbind(data.frame(pres_back = 1, quad_1var),
                 data.frame(pres_back = 0, background))

data_q2 <- rbind(data.frame(pres_back = 1, vd_pre_data[[1]]),
                 data.frame(pres_back = 0, background))

data_q3 <- rbind(data.frame(pres_back = 1, vd_pre_data[[2]]),
                 data.frame(pres_back = 0, background))

data_q4 <- rbind(data.frame(pres_back = 1, vd_pre_data[[3]]),
                 data.frame(pres_back = 0, background))

### save the data
save(data_q1, data_q2, data_q3, data_q4, file = "Data/data_quadratic.RData")


## niches with product responses (linear and quadratic as well)
### virtual niches
v_nichesp <- lapply(1:length(comb), function(x) {
  ## variances
  vars <- evniche:::var_from_range(range = qranges[, comb[[x]]])

  ## covariance limits
  cov_lim <- covariance_limits1(range = qranges[, comb[[x]]], prop_reg = 0.001)
  ## variance-covariance matrix
  covm <- cov_lim$max_covariance
  varcov <- evniche:::var_cov_matrix(variances = vars, covariances = covm)

  ## centroid
  cent <- evniche:::centroid(range = qranges[, comb[[x]]])

  ## ellipsoid characteristics (virtual niche)
  ell_features(centroid = cent, covariance_matrix = varcov, level = 0.99)
})

names(v_nichesp) <- names(comb)

### generating data from niches and background
vd_pre_datap <- lapply(1:length(v_nichesp), function(x) {
  set.seed(1)
  virtual_data(features = v_nichesp[[x]], from = "ellipsoid", n = 1000)
})

names(vd_pre_datap) <- names(v_nichesp)

### complete missing columns
vd_pre_datap[[1]] <- cbind(bio_1 = vd_pre_datap[[1]][, 1], bio_7 = bio7,
                          bio_12 = vd_pre_datap[[1]][, 2], bio_15 = bio15)
vd_pre_datap[[2]] <- cbind(vd_pre_datap[[2]], bio_15 = bio15)

### preparing data for models (linear and quadratic responses)
data_p2 <- rbind(data.frame(pres_back = 1, vd_pre_datap[[1]]),
                 data.frame(pres_back = 0, background))

data_p3 <- rbind(data.frame(pres_back = 1, vd_pre_datap[[2]]),
                 data.frame(pres_back = 0, background))

data_p4 <- rbind(data.frame(pres_back = 1, vd_pre_datap[[3]]),
                 data.frame(pres_back = 0, background))

### save the data
save(data_p2, data_p3, data_p4, file = "Data/data_product.RData")
