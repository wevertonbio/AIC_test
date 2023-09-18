library(geodata)
library(evniche)

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

# uniform distribution of environments accornig to ranges of environments
## ranges
ranges <- minmax(var)

## sets of environments
set.seed(1)
data_env <- sapply(1:ncol(ranges), function (x) {
  runif(100000, min = ranges[1, x], max = ranges[2, x])
})

colnames(data_env) <- colnames(ranges)

head(data_env)
save(data_env, file = "Data/example_data.RData")



# presences denoting sitinct types of response to variables
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

## niches with linear responses
### uniform values along the range of variables
nls <- 44

set.seed(40)
bu1 <- sort(runif(nls, ranges[1, 1], ranges[2, 1]))
bu7 <- sort(runif(nls, ranges[1, 2], ranges[2, 2]))
bu12 <- sort(runif(nls, ranges[1, 3], ranges[2, 3]))
bu15 <- sort(runif(nls, ranges[1, 4], ranges[2, 4]))

### values with increasing probabilities towards high
reps <- 1:nls
reps[(nls-3):nls] <- reps[(nls-3):nls] + 1:4

bl1 <- unlist(lapply(1:nls, function(x) {rep(bu1[x], reps[x])}))
bl7 <- unlist(lapply(1:nls, function(x) {rep(bu7[x], reps[x])}))
bl12 <- unlist(lapply(1:nls, function(x) {rep(bu12[x], reps[x])}))
bl15 <- unlist(lapply(1:nls, function(x) {rep(bu15[x], reps[x])}))

### data with linear responses
data_l1 <- cbind(bio_1 = bl1, bio_7 = bio7, bio_12 = bio12, bio_15 = bio15)
data_l2 <- cbind(bio_1 = bl1, bio_7 = bio7, bio_12 = bl12, bio_15 = bio15)
data_l3 <- cbind(bio_1 = bl1, bio_7 = bl7, bio_12 = bl12, bio_15 = bio15)
data_l4 <- cbind(bio_1 = bl1, bio_7 = bl7, bio_12 = bl12, bio_15 = bl15)

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

  ## covariance limits
  #cov_lim <- evniche:::covariance_limits(range = qranges[, comb[[x]]])
  ## variance-covariance matrix
  cov <- 0  # cov_lim$max_covariance * ifelse(x == 6, 0.5, 0.2) # covariance selected
  varcov <- evniche:::var_cov_matrix(variances = vars, covariances = cov)

  ## centroid
  cent <- evniche:::centroid(range = qranges[, comb[[x]]])

  ## ellipsoid characteristics (virtual niche)
  ell_features(centroid = cent, covariance_matrix = varcov, level = 0.99)
})

names(v_niches) <- names(comb)

### generating data from niches and background
vd_pre_data <- lapply(1:length(v_niches), function(x) {
  ## predict suitability
  predp <- evniche:::ell_predict(data = data_env[, comb[[x]]],
                                 features = v_niches[[x]])

  ## generate new data
  ## based on the ellipsoid and available conditions
  set.seed(1)
  virtual_data(features = v_niches[[x]], from = "ellipsoid", n = 1000)
})

names(vd_pre_data) <- names(v_niches)

### complete missing rows
vd_pre_data[[1]] <- cbind(bio_1 = vd_pre_data[[1]][, 1], bio_7 = bio7,
                          bio_12 = vd_pre_data[[1]][, 1], bio_15 = bio15)
vd_pre_data[[2]] <- cbind(vd_pre_data[[2]], bio_15 = bio15)

### quadratic fro one variable
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

