####Compare AIC - Example####
####Import data####
occ_bg <- read.csv("Examples/SWD_file.csv")
f_grid <- read.csv("Examples/Grid.csv")

#Calibrate models and get AICs
df_aic <- compare_aic(data = occ_bg, formula_grid = f_grid,
                      pr_bg = "pr_ab", n_folds = 4, parallel = TRUE,
                      ncores = 8, progress_bar = TRUE)

#Summarize data
s_aic <- eval_stats(df_aic)
