#### Functions to compare AICs ####

compare_aic <- function(data, #Data in SWD format
                        pr_bg, #Column name with presence (1) or background (0)
                        formula_grid, #Grid with formulas
                        n_folds = 4, #Columns name with k_folds or vector indicating k_folds
                        parallel = TRUE,
                        ncores = 6,
                        progress_bar = TRUE) { #Show progress bar?
  ####Get k-folds####
  pre <- which(data[, pr_bg] == 1)
  aus <- which(data[, pr_bg] == 0)
  foldp <- sample(cut(seq(1, length(pre)), breaks = n_folds, labels = FALSE))
  folda <- sample(cut(seq(1, length(aus)), breaks = n_folds, labels = FALSE))
  #Create columns to save pr_bg
  data$folds <- NA
  data$folds[which(data[,pr_bg] == 1)] <- foldp
  data$folds[which(data[,pr_bg] == 0)] <- folda


  #Make and register cluster
  cl <- parallel::makeCluster(ncores)
  #Register objectes in cluster
  plan(cluster, workers=cl)
  #Get bar progress
  handlers(handler_progress(format="[:bar] :percent :eta :message"))

  #Start looping
  with_progress({
    pg <- progressor(along=1:nrow(formula_grid))
    eval_m <- future_lapply(1:nrow(formula_grid), function(x){
      tryCatch(
        {#Get grid x
          grid_x <- formula_grid[x,] #Get i candidate model
          formula_x <- as.formula(grid_x$Formulas) #Get formula from grid x
          reg_x <- grid_x$regm #Get regularization multiplier from grid x

          #Complete model with AIC do japones
          m_aic <- glmnet_mx(p = occ_bg[, pr_bg], data = occ_bg,
                             f = formula_x, regmult = reg_x, calculate_AIC = T)

          #AIC from Warren (Kuenm)
          prediction <- predict.glmnet_mx(m_aic, newdata = occ_bg, type = "exponential")
          vals <- predict.glmnet_mx(m_aic, occ_bg[occ_bg$pr_ab == 1, ],
                                    type = "exponential")
          LL <- sum(log(vals + .Machine$double.eps))
          npar <- length(m_aic$betas)
          AICc <- ((2 * npar) - (2 * LL)) + (2 * npar * (npar +
                                                           1)/(length(vals) - npar - 1))
          #Model with kfolds
          k_fold <- occ_bg$fold
          mods <- lapply(1:length(unique(k_fold)), function(i) {
            occ_bg_i <- subset(occ_bg, occ_bg$fold != i) #Select i k-fold
            #Run model
            mod_i <- glmnet_mx(p = occ_bg_i[,pr_bg], data = occ_bg_i,
                               f = formula_x, regmult = reg_x,
                               calculate_AIC = FALSE)

            #Predict model only to background
            pred_i <- as.numeric(predict(object = mod_i, newdata = occ_bg, clamp = FALSE,
                                         type = "cloglog"))
            pred_i <- cbind(occ_bg[,c("pr_ab", "fold")], "suit" = pred_i)

            # #Predict model to spatraster, only to check
            # vars_r <- terra::rast("Models/Araucaria_angustifolia/PCA_variables.tiff")
            # pred_r <- terra::predict(vars_r, mod_i, type = "cloglog", na.rm = TRUE)
            # plot(pred_r)

            #Extract suitability in train and test points
            train_xy <- subset(occ_bg, occ_bg$fold != i & pr_ab == 1)[, c("x", "y")]
            test_xy <- subset(occ_bg, occ_bg$fold == i & pr_ab == 1)[, c("x", "y")]
            suit_val_cal <- subset(pred_i, occ_bg$fold != i & pr_ab == 1)[, "suit"]
            suit_val_eval <- subset(pred_i, occ_bg$fold == i & pr_ab == 1)[, "suit"]
            ####Calculate omission rate following kuenm####
            om_rate <- omrat_maxnet(threshold = c(5, 10),
                                    pred_train = suit_val_cal,
                                    pred_test = suit_val_eval)

            # ####Calculate pROC following kuenm####
            # proc_i <- proc_maxnet(pred_test = suit_val_eval, model = pred_i$suit, threshold = 5,
            #                       rand_percent = 50, iterations = 500, parallel = F)

            #Calculate pROC following enmpa
            proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                      prediction = pred_i$suit)

            ####Save metrics in a dataframe
            df_eval <- data.frame(Replicate = i,
                                  t(data.frame(om_rate)),
                                  proc_auc_ratio = proc_i$pROC_summary[1],
                                  proc_pval = proc_i$pROC_summary[2],
                                  AIC = m_aic$AIC,
                                  AIC_warren = AICc, row.names = NULL,
                                  npar = npar)
            df_eval2 <- cbind(grid_x, df_eval)
            return(df_eval2)
          })
          #Return evaluaiton final
          eval_final <- do.call("rbind", mods)
          return(eval_final)
          #Progress
          pg(sprintf("i=%g", i))}, #Print progress
        error = function(e) NULL)
    }, future.chunk.size=1)
  })
  df_replicates <- do.call("rbind", eval_m)
  return(df_replicates)
} #End of function


#Function to summarize data
eval_stats <- function(calib_results){
  toagg <- c("Omission_rate_at_5", "Omission_rate_at_10",
             "proc_auc_ratio", "proc_pval")
  agg_formula <- "~ Formulas + regm + Features"
  xy <- lapply(toagg, function(x) {
    do.call(data.frame, aggregate(as.formula(paste(x, agg_formula)),
                                  data = calib_results, FUN = function(y) c(mean = round(mean(y),
                                                                                         4), sd = round(sd(y), 4))))
  })
  #Summarize stats
  stats <- Reduce(function(x, y) merge(x, y,
                                       by = c("Formulas", "regm", "Features")),
                  xy)

  stats_AICS <- calib_results[!duplicated(calib_results[,c("Formulas", "regm", "AIC", "AIC_warren",
                                                           "npar")]),] [,c("Formulas", "regm", "AIC", "AIC_warren",
                                                                           "npar")]
  stats_final <- merge(stats, stats_AICS, by = c("Formulas", "regm"))
  return(stats_final)
}





####Omission Rate####
omrat_maxnet <- function(threshold = 5, pred_train, pred_test) {
  om_rate <- vector("numeric", length = length(threshold))
  for (i in 1:length(threshold)) {
    val <- ceiling(length(pred_train) * threshold[i]/100) + 1
    omi_val_suit <- sort(pred_train)[val]
    om_rate[i] <- as.numeric(length(pred_test[pred_train <
                                                    omi_val_suit])/length(pred_test))
  }
  names(om_rate) <- paste("Omission_rate_at_", threshold, sep = "")
  return(om_rate)
}
# t2 <- omrat_maxnet(threshold = 5,
#                    pred_train =  suit_val_cal,
#                    pred_test = suit_val_eval)


#### New functions for maxnet ####
glmnet_mx <- function(p, data, f, regmult = 1.0,
                      regfun = maxnet.default.regularization,
                      addsamplestobackground = TRUE, weights_1_0 = c(1, 100),
                      calculate_AIC = FALSE, ...) {
  if (anyNA(data)) {
    stop("NA values in data table. Please remove them and rerun.")
  }
  if (!is.vector(p)) {
    stop("p must be a vector.")
  }

  if (addsamplestobackground) {
    pdata <- data[p == 1, ]
    ndata <- data[p == 0, ]

    # add to background any presence data that isn't there already
    wadd <- !do.call(paste, pdata) %in% do.call(paste, ndata)
    p <- c(p, rep(0, sum(wadd)))
    toadd <- pdata[wadd, ]
    data <- rbind(data, toadd)
  }

  mm <- model.matrix(f, data)
  reg <- regfun(p, mm) * regmult
  weights <- ifelse(p == 1, weights_1_0[1], weights_1_0[2])
  lambdas <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) *
    sum(p) / sum(weights)

  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",
                          standardize = FALSE, penalty.factor = reg,
                          lambda = lambdas, weights = weights, ...)

  bb <- coef(model)[, 200]

  class(model) <- c("glmnet_mx", class(model))
  if (length(model$lambda) < 200) {
    msg <- "glmnet failed to complete regularization path. Model may be infeasible."
    if (!addsamplestobackground) {
      msg <- paste(msg, "\n\tTry re-running with addsamplestobackground = TRUE.")
    }
    stop(msg)
  }

  # AIC calculation
  filter <- bb[-1] != 0
  bb <- c(bb[1], bb[-1][filter])

  if (calculate_AIC) {
    model$AIC <- aic_glmnetmx(x = mm[, filter], y = p, beta = bb)
  } else {
    model$AIC <- NA
  }

  # returning other elements
  model$betas <- bb[-1]
  model$alpha <- 0
  rr <- predict.glmnet_mx(model, data[p == 0, , drop = FALSE],
                          type = "exponent", clamp = FALSE)
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw))
  model$alpha <- -log(sum(rr))

  model$penalty.factor <- reg
  model$featuremins <- apply(mm, 2, min)
  model$featuremaxs <- apply(mm, 2, max)

  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2, max)

  means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data)[!vv], function(n) {
    which.max(table(data[p == 1, n, drop = FALSE]))
  })
  names(majorities) <- names(data)[!vv]
  model$samplemeans <- unlist(c(means, majorities))

  model$levels <- lapply(data, levels)

  return(model)
}



#### maxnet default regularization ####
maxnet.default.regularization <- function (p, m) {
  isproduct <- function(x) {grepl(":", x) & !grepl("\\(", x)}
  isquadratic <- function(x) {grepl("^I\\(.*\\^2\\)", x)}
  ishinge <- function(x) {grepl("^hinge\\(", x)}
  isthreshold <- function(x) {grepl("^thresholds\\(", x)}
  iscategorical <- function(x) {grepl("^categorical\\(", x)}

  regtable <- function(name, default) {
    if (ishinge(name)) {
      return(list(c(0, 1), c(0.5, 0.5)))
    }
    if (iscategorical(name)) {
      return(list(c(0, 10, 17), c(0.65, 0.5, 0.25)))
    }
    if (isthreshold(name)) {
      return(list(c(0, 100), c(2, 1)))
    }
    default
  }

  lregtable <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
  qregtable <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5, 0.25, 0.05))
  pregtable <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9, 0.55, 0.05))

  mm <- m[p == 1, ]
  np <- nrow(mm)
  lqpreg <- lregtable

  if (sum(isquadratic(colnames(mm)))) {
    lqpreg <- qregtable
  }
  if (sum(isproduct(colnames(mm)))) {
    lqpreg <- pregtable
  }

  classregularization <- sapply(colnames(mm), function(n) {
    t <- regtable(n, lqpreg)
    approx(t[[1]], t[[2]], np, rule = 2)$y
  })
  classregularization <- classregularization / sqrt(np)


  ishinge <- grepl("^hinge\\(", colnames(mm))

  hmindev <- sapply(1:ncol(mm), function(i) {
    if (!ishinge[i]) {
      return(0)
    }
    avg <- mean(mm[, i])
    std <- max(sd(mm[, i]), 1 / sqrt(np))
    std * 0.5 / sqrt(np)
  })

  tmindev <- sapply(1:ncol(mm), function(i) {
    ifelse(isthreshold(colnames(mm)[i]) && (sum(mm[, i]) ==
                                              0 || sum(mm[, i]) == nrow(mm)), 1, 0)
  })

  pmax(0.001 * (apply(m, 2, max) - apply(m, 2, min)), hmindev,
       tmindev, apply(as.matrix(mm), 2, sd) * classregularization)
}


#### aic for glmnet like maxent ####
aic_glmnetmx <- function(x, y, beta) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.")
  }
  if (mode(x) != "numeric") {
    stop("x must be numeric.")
  }

  x_glm <- glm(y ~ x, family = "binomial")
  predict_glm <- predict(x_glm)

  pr0 <- 1 - 1 / (1 + exp(predict_glm))
  pr <- 1 - 1 / (1 + exp(cbind(1, x) %*% beta))
  cc <- as.vector(beta)
  jc <- abs(cc) > 1e-05
  xj <- cbind(1, x)[, jc]
  Pi <- diag(as.vector(pr * (1 - pr)))
  Pi0 <- diag(as.vector(pr0 * (1 - pr0)))
  j22 <- t(xj) %*% Pi %*% xj
  i22 <- t(xj) %*% Pi0 %*% xj

  aic <- -2 * sum(y * log(pr) + (1 - y) * log(1 - pr)) +
    2 * sum(diag(solve(j22) %*% i22))

  return(aic)
}


#### predict gmlm like maxent to new data ####
predict.glmnet_mx <- function (object, newdata, clamp = FALSE,
                               type = c("link", "exponential",
                                        "cloglog", "logistic")) {
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v]),
                           object$varmax[v])
    }
  }

  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
               names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
               terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
               terms)
  f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))

  mm <- model.matrix(f, data.frame(newdata))

  if (clamp) {
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                 object$featuremaxs[names(object$betas)]))
  }

  link <- (mm %*% object$betas) + object$alpha
  type <- match.arg(type)

  # return prediction
  if (type == "link") {
    return(link)
  }
  if (type == "exponential") {
    return(exp(link))
  }
  if (type == "cloglog") {
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic") {
    return(1/(1 + exp(-object$entropy - link)))
  }
}


####Get formulas####
get_formulas_maxnet <- function (independent, type = "lqpht",
                                 minvar = 1, maxvar = NULL)
{
  if (is.character(type)) {
    if (!all(unlist(strsplit(type, "")) %in% c("l", "p",
                                               "q", "h", "t"))) {
      stop("'type' must be: 'l', 'p', 'q', 'h', 't', or a combination of those three.")
    }
  }  else {
    stop("'type' must be of class character.")
  }
  predictors <- independent
  npred <- length(predictors)
  aux <- " "
  #Linear
  if (grepl("l", type)) {
    aux <- paste(aux, paste(predictors, collapse = " + "),
                 sep = " + ")
  }
  #Quadratic
  if (grepl("q", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("I(", predictors[v], "^2)"),
                   sep = " + ")
    }
  }
  #Product
  if (grepl("p", type)) {
    if (npred > 1) {
      inter_tab <- utils::combn(predictors, 2)
      aux_inter <- paste0(" + ", paste(apply(inter_tab,
                                             2, paste, collapse = ":"), collapse = " + "))
      if (length(aux_inter) > 0) {
        aux <- paste0(aux, aux_inter)
      }
    }
    else {
      if (grepl("l", type) | grepl("q", type)) {
        message("'p' is is only possible with 2 or more independent variables.",
                "\nReturning other combinations.")
      }
      else {
        stop("'p' is is only possible with 2 or more independent variables.",
             "\nTry other combinations of type.")
      }
    }
  }
  #Hinge
  if (grepl("h", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("hinge(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  #Threshold
  if (grepl("t", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("thresholds(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  out <- paste0("~",
                gsub(pattern = "  \\+ ", x = aux, replacement = ""))
  return(out)
}

#my_formula <- get_formulas(independent = c("PC1", "PC2"), type = "lq")


# ####PROC####
# #Not necessary, use enmpa::proc_enm
# proc_maxnet <- function(pred_test, model, threshold = 5,
#                         rand_percent = 50,
#                         iterations = 500, parallel = FALSE) {
#   min_pred <- min(model, na.rm = T)
#   max_pred <- max(model, na.rm = T)
#   vals <- na.omit(model)
#   nvals <- length(vals)
#   vals <- c(vals, pred_test)
#   vals <- as.numeric(cut(vals, 500))
#   pred_test <- vals[(nvals + 1):length(vals)]
#   vals <- vals[1:nvals]
#   classpixels <- as.data.frame(table(vals), stringsAsFactors = FALSE)
#   colnames(classpixels) <- c("value", "count")
#   if (min_pred == max_pred) {
#     warning("\nmodel has no variability, pROC will return NA.\n")
#     p_roc <- rep(NA, 2)
#     names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", threshold,
#                              "%"), "pval_pROC")
#     auc_ratios <- rep(NA, 3)
#     names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC",
#                            "AUC_ratio")
#     p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
#   } else {
#     classpixels <- classpixels %>% dplyr::mutate(value = rev(value),
#                                                  count = rev(count), totpixperclass = cumsum(count),
#                                                  percentpixels = totpixperclass/sum(count)) %>% dplyr::arrange(value)
#     error_sens <- 1 - (threshold/100)
#     prediction_errors <- classpixels[, "value"]
#     fractional_area <- classpixels[, "percentpixels"]
#     n_data <- length(pred_test)
#     n_samp <- ceiling((rand_percent/100) * n_data)
#     big_classpixels <- matrix(rep(prediction_errors, each = n_samp),
#                               ncol = length(prediction_errors))
#     partial_AUC <- 1:iterations %>% purrr::map_df(~kuenm:::calc_aucDF(big_classpixels,
#                                                                       fractional_area, pred_test, n_data, n_samp, error_sens))
#     naID <- !is.na(partial_AUC$auc_ratio)
#     nona_valproc <- partial_AUC$auc_ratio[naID]
#     mauc <- mean(nona_valproc)
#     proc <- sum(nona_valproc <= 1)/length(nona_valproc)
#     p_roc <- c(mauc, proc)
#     names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", threshold,
#                              "%"), "pval_pROC")
#     auc_ratios <- partial_AUC
#     names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC",
#                            "AUC_ratio")
#     p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
#   }
# }
# #t <- proc_maxnet(pred_test = suit_val_eval, model = pred_i)

# ####AIC de Warren####
# #It's not necessary, calculating inside looping
# aic_maxnet <- function(exp_pred_occ, ncoefs, pred_aic = NULL) {
#   #Get sum of whole area
#   p.sum <- sum(values(pred_aic, na.rm=TRUE))
#   exp_pred_occ <- exp_pred_occ/p.sum
#
#   #Calculate AIC following ENMeval::aic.maxent()
#   n.occs <- length(exp_pred_occ)
#   AIC.valid <- ncoefs < n.occs
#   for (i in 1:length(AIC.valid)) {
#     if (AIC.valid[i] == FALSE) {
#       message(paste("Warning: model", names(AIC.valid)[i],
#                     "has more non-zero coefficients (ncoef) than occurrence records for training, so AIC cannot be calculated."))
#     }
#   }
#   LL <- sum(log(exp_pred_occ), na.rm = TRUE)
#   AICc <- (2 * ncoefs - 2 * LL) + (2 * (ncoefs) * (ncoefs +
#                                                      1)/(n.occs - ncoefs - 1))
#   AICc <- sapply(1:length(AICc), function(x)
#     ifelse(AIC.valid[x] ==FALSE | is.infinite(AICc[x]), NA, AICc[x]))
#   return(AICc)
#   }
