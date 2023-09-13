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
  runif(1000000, min = ranges[1, x], max = ranges[2, x])
})

colnames(data_env) <- colnames(ranges)

head(data_env)
save(data_env, file = "Data/example_data.RData")



# virtual niches
## descriptive stats of example data (background)
summary(data_env)

## combinations
comb <- list(v1_2 = paste0("bio_", c(1, 12)),
             v1_3 = paste0("bio_", c(1, 7, 12)),
             v1_4 = paste0("bio_", c(1, 7, 12, 15)))

## niches with linear responses
### linear probabilities
lprobs <- 1:1000000 / 1000000
sprobs <- 1:100000 / 100000
s1probs <- 1:10000 / 10000
s2probs <- 1:1000 / 1000

### creating positive records: linear response to one variable
data_envs <- data_env[order(data_env[, "bio_1"]), ]
data_pres <- data_envs[sample(1000000, 100000, prob = lprobs), ]
data_l1 <- data_pres[sample(100000, 100), ]

### creating positive records: linear response to two variables
data_pres <- data_pres[order(data_pres[, "bio_12"]), ]
data_pres <- data_pres[sample(100000, 10000, prob = sprobs), ]
data_l2 <- data_pres[sample(10000, 100), ]

### creating positive records: linear response to three variables
data_pres <- data_pres[order(data_pres[, "bio_7"]), ]
data_pres <- data_pres[sample(10000, 1000, prob = s1probs), ]
data_l3 <- data_pres[sample(1000, 100), ]

### creating positive records: linear response to four variables
data_pres <- data_pres[order(data_pres[, "bio_15"]), ]
data_l4 <- data_pres[sample(1000, 100, prob = s2probs), ]

save(data_l1, data_l2, data_l3, data_l4, file = "Data/data_linear.RData")

## niches with quadratic responses (linear as well)
### ranges for niches
qranges <- lapply(1:ncol(data_env), function(x) {
  quantile(data_env[, x])[c(2, 4)]
})

names(qranges) <- colnames(ranges)


## niches with product responses (linear and quadratic as well)

