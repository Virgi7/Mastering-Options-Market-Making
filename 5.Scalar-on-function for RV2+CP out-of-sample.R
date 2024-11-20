# Scalar-on-functon Model: RV(ti+10min)^2 = α + int(β(s) · X(s)^2) + ϵ(ti+10min)
# with the introduction of non-conformal prediction method 
#####
# Libraries
#####
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)
library(tidyr)
library(progress)
library(pbapply)
pboptions(type='none')
library(dbscan)
library(gridExtra)
Sys.setenv("ARROW_WITH_SNAPPY" = "ON")
#install.packages("arrow",force = TRUE)
library(arrow)
library("readxl")
#install.packages("caret")
library(caret)
library(ISLR2)
library(car)
library(np)
library(splines)
library(fda)
library(magrittr)
library(KernSmooth)
library(openxlsx)
library(zoo)
library(broom)
library(data.table)
library(fda.usc)
library(viridis)
library(scales)
library(patchwork)
source("fRegress std error mia versione.R")
#####
# Dataset upload
#####
# Target upload
train <- read_excel('train.xlsx')
# Stocks list
stock_id_unique <- unique(train$stock_id)
print(stock_id_unique)
# The analyzed stock was chosen randomly
stock_id_chosen <- 61 #sample(stock_id_unique, 1)

# Order Book data frame
file_pattern <- paste0('book_train.parquet/stock_id=', stock_id_chosen, '/*.parquet')
file_paths <- list.files(path = dirname(file_pattern), pattern = basename(file_pattern), full.names = TRUE)
df_book_data <- lapply(file_paths, read_parquet)
df_book_data <- do.call(rbind, df_book_data)
#####
# Utilities
#####
# MSE and RMSE functions
mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}
rmse <- function(y_true, y_pred) {
  sqrt(mse(y_true, y_pred))
}

# Log-returns calculation (with an initial NA)
log_return <- function(list_stock_prices) {
  return(c(NA, diff(log(list_stock_prices))))  
}
sqlogret_per_time_id <- function(file_path) {
  df_book_data <- read_parquet(file_path)
  
  df_book_data <- df_book_data %>%
    mutate(
      wap = (bid_price1 * ask_size1 + ask_price1 * bid_size1) / (bid_size1 + ask_size1),
      log_return = ave(wap, time_id, FUN = function(x) log_return(x))^2
    ) %>%
    filter(!is.na(log_return))
  
  return(df_book_data)
}

# RV(ti) calculation
realized_volatility <- function(series_log_return) {
  return(sqrt(sum(series_log_return^2, na.rm = TRUE))) 
}

# RV(ti) calculation for each bucket
realized_volatility_per_time_id <- function(file_path, prediction_column_name) {
  df_book_data <- as.data.table(read_parquet(file_path))
  df_book_data[, wap := (bid_price1 * ask_size1 + ask_price1 * bid_size1) / (bid_size1 + ask_size1)]
  df_book_data[, log_return := log_return(wap), by = .(time_id)]
  df_book_data <- df_book_data[!is.na(log_return)]
  
  df_realized_vol_per_stock <- df_book_data[, .(volatility = realized_volatility(log_return)), by = .(time_id)]
  
  setnames(df_realized_vol_per_stock, "volatility", prediction_column_name)
  
  stock_id <- strsplit(file_path, "=")[[1]][2]
  stock_id <- strsplit(stock_id, "/")[[1]][1]
  
  df_realized_vol_per_stock[, row_id := paste0(stock_id, "-", time_id)]
  
  return(df_realized_vol_per_stock[, .(row_id, pred=get(prediction_column_name))])
  
}

#  RV(ti) calculation for each bucket in different files
past_realized_volatility_per_stock <- function(list_file, prediction_column_name) {
  df_past_realized <- data.table()  
  
  # For each file in the list
  for (file in list_file) {
    df_past_realized <- rbind(df_past_realized,
                              realized_volatility_per_time_id(file, prediction_column_name))
  }
  
  return(df_past_realized)
}
#####
# Dataset preparation
#####
#  RV(ti) calculation for each bucket
df_past_realized_train <- past_realized_volatility_per_stock(list_file = file_paths,
                                                             prediction_column_name = 'pred')

# Selection of realized volatility for the current stock
train_stock <- train %>% filter(stock_id == stock_id_chosen)

# Matching of 'stock_id' and 'time_id' --> 'row_id'
train_stock <- train_stock %>% mutate(row_id = paste(stock_id, time_id, sep = "-"))
train_stock <- train_stock %>% dplyr::select(row_id, target)

# Merge between df_past_realized_train and train on 'row_id', in order to construct the data set needed 
# for the linear model
df_joined <- train_stock %>%
  left_join(df_past_realized_train %>% dplyr::select(row_id, pred), by = "row_id")
df_joined <- df_joined %>%
  mutate(row_id = as.integer(sapply(strsplit(row_id, "-"), `[`, 2)))

df_book_data <- sqlogret_per_time_id(file_paths)

time_buckets <- as.vector(df_joined$row_id)
df_past_realized_train$row_id <- time_buckets

# Training-Validation-Test Splitting
split_time_buckets_no_replacement <- function(time_buckets, train_size = 0.56, val_size = 0.24, test_size = 0.2, seed = 123) {
  
  # Check that consider all the dataset
  if (round(train_size + val_size + test_size, 2) != 1) {
    stop("Train, Validation and Test have to sum up at 1")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Total number of observations
  n <- length(time_buckets)
  
  # Number of observations for each set
  n_train <- round(n * train_size)
  n_val <- round(n * val_size)
  n_test <- n - n_train - n_val 
  
  # Sampling without substitution for training
  train_time_buckets <- sample(time_buckets, size = n_train, replace = FALSE)
  
  # Remove the training buckets
  remaining_buckets <- setdiff(time_buckets, train_time_buckets)
  
  val_time_buckets <- sample(remaining_buckets, size = n_val, replace = FALSE)
  
  # Remaining buckets will be assigned to test set
  test_time_buckets <- setdiff(remaining_buckets, val_time_buckets)
  
  return(list(train = train_time_buckets, validation = val_time_buckets, test = test_time_buckets))
}

split_results <- split_time_buckets_no_replacement(time_buckets)
#####
# Training
#####
# Training data: True RV(ti+10min)^2
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$train, ]
vol <- as.numeric(vol$target)^2
# Order book data
df_book_data_train <- df_book_data[df_book_data$time_id %in% split_results$train, ]

# Time series
# Column selection
column_name <- "log_return"  
# Dataframe creation
df_dynamic <- df_book_data_train[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
# Conversion of df_dynamic in data.table for a faster process
df_dynamic <- as.data.table(df_dynamic)
# Matrix creation (600 rows and n_obs columns)
time_istants <- unique(df_book_data_train$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
# Multiple NAs at the beginning and at the end of the matrix columns reach to some
# problems in the interpolation phase
# Counting NA at the beginning
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
# Counting NA at the end
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)

df_matrix <- as.matrix(df)
# Substitution of NAs with 0s which mean no-changes in prices
for (i in 1:length(na_start)){
  if(na_start[i]!=0){
    if (na_start[i]>1)
      df_matrix[1:na_start[i], i]<-0
    else
      df_matrix[1, i]<-0
  }
}
for (i in 1:length(na_end)) {
  if (na_end[i] != 0) {
    if (na_end[i] > 1) {
      df_matrix[(600 - na_end[i] + 1):600, i] <- 0
    } else {
      df_matrix[600, i] <- 0
    }
  }
}
Xobs0 <- df_matrix
# Spline interpolation
Xobs0 <- na.spline(Xobs0)

n_obs = length(vol)
abscissa = 0:599

volatility_squared <- vol

# X(t)
logret_list = vector("list", 2)
logret_list[[1]] = rep(1, n_obs)
nbasis_beta = 10
grade = 4
m = grade
degree = m-1
nbasis_covariates = 110 
logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, logret_basis)
logretfd = logretSmooth$fd
logret_list[[2]] = logretfd

# B(t)
conbasisinter = create.constant.basis(rangeval=c(0,600))
beta_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_beta, norder=m)
betalist = vector("list", 2)
betalist[[1]] = conbasisinter
betalist[[2]] = beta_basis

# Scalar-on-function Model
fRegressList = fRegress(volatility_squared, logret_list, betalist)
betaestlist_nonpen = fRegressList$betaestlist

# Target selection for the chosen stock
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
# Not squared
volatility <- as.numeric(vol[vol$time_id %in% split_results$train, ]$target)

volatility_squaredhat_nonpen = fRegressList$yhatfdobj
# Residuals evaluation
volatility_squaredres_nonpen = volatility_squared - volatility_squaredhat_nonpen

# Diagnosis
SSE1.1 = sum(volatility_squaredres_nonpen^2)
SSE0 = sum((volatility_squared - mean(volatility_squared))^2)
RSQ1 = (SSE0 - SSE1.1)/SSE0
Fratio1 = ((SSE0 - SSE1.1)/5)/(SSE1.1/29)
RMSE.1 = rmse(volatility, sqrt(volatility_squaredhat_nonpen))

# Confidence interval at 95% for B(t) function
sigmaE. = sum(volatility_squaredres_nonpen^2)/(n_obs - fRegressList$df)
sigmaE = sigmaE. * diag(rep(1, n_obs))
y2cMap = logretSmooth$y2cMap
stderrList = fRegress.stderr(fRegressList, y2cMap, sigmaE)
betafdPar = betaestlist_nonpen[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
plot(betafd, xlab = "Seconds", ylab = "Log-return2 Reg. Coeff.", lwd = 2)
lines(betafd + betastderrfd, lty=2, lwd=1)
lines(betafd - betastderrfd, lty=2, lwd=1)

# Penalized model implementation
lambda0 = 2^2
betafdPar = fdPar(beta_basis, 2, lambda0)
betalist[[2]] = betafdPar

# Cross validation
loglam = seq(-15.5,0,0.5)
nlam = length(loglam)
SSE1.CV = matrix(0, nlam, 1)
logret_betafd_CV = matrix(0, nlam, nbasis_beta)
RMSE.CV = matrix(0, nlam, 1)
GCVscore =  matrix(0, nlam, 1)
for(ilam in 1:nlam){
  lambda = 2^loglam[ilam]
  betalist_i = betalist
  betafdPar2 = betalist_i[[2]]
  betafdPar2$lambda = lambda
  betalist_i[[2]] = betafdPar2
  fReg_i = fRegress(volatility_squared, logret_list, betalist_i)
  betaestlist_i = fReg_i$betaestlist
  logret_betafd_CV[ilam] = betaestlist_i[[2]]$fd
  volatility_squaredhat_CV = fReg_i$yhatfdobj
  volatility_squaredres_CV = volatility_squared - volatility_squaredhat_CV
  sigmaE. = sum(volatility_squaredres_CV^2)/(n_obs - fReg_i$df)
  sigmaE = sigmaE. * diag(rep(1, n_obs))
  y2cMap = logretSmooth$y2cMap
  stderrList = fRegress.stderr(fReg_i, y2cMap, sigmaE)
  betafdPar = betaestlist_i[[2]]
  betafd = betafdPar$fd
  betastderrList = stderrList$betastderrlist
  betastderrfd = betastderrList[[2]]
  plot(betafd, xlab = "Seconds", ylab = "Log-return2 Reg. Coeff.", lwd = 2)
  lines(betafd + 2*betastderrfd, lty=2, lwd=1)
  lines(betafd - 2*betastderrfd, lty=2, lwd=1)
  SSE1.CV[ilam] = sum(volatility_squaredres_CV^2)
  RMSE.CV[ilam] = rmse(volatility, sqrt(volatility_squaredhat_CV))
  GCVscore[ilam] = fReg_i$GCV
}

# Minimum GCV score is for lambda = 2^-15
minGCVscore <- min(GCVscore)
min_index <- which.min(GCVscore) 
exp_lambda_opt <- loglam[min_index]
lambda_opt_GCV <- 2^exp_lambda_opt

lambda_opt = lambda_opt_GCV

betalist_pen = betalist
betafdPar2 = betalist_pen[[2]]
betafdPar2$lambda = lambda_opt
betalist_pen[[2]] = betafdPar2
fRegressList_pen = fRegress(volatility_squared, logret_list, betalist_pen)
betaestlist_pen = fRegressList_pen$betaestlist
logret_betafd_pen = betaestlist_pen[[2]]$fd
volatility_squaredhat_pen = fRegressList_pen$yhatfdobj
volatility_squaredres_pen = volatility_squared - volatility_squaredhat_pen
sigmaE. = sum(volatility_squaredres_pen^2)/(n_obs - fRegressList_pen$df)
sigmaE = sigmaE. * diag(rep(1, n_obs))
y2cMap = logretSmooth$y2cMap
stderrList = fRegress.stderr(fRegressList_pen, y2cMap, sigmaE)
betafdPar = betaestlist_pen[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
RMSE.2 = rmse(volatility, sqrt(volatility_squaredhat_pen))
SSE1.2 = sum(volatility_squaredres_pen^2)
RSQ2 = (SSE0 - SSE1.2)/SSE0
# Final plot
plot(betafd, xlab = "Seconds", ylab = "Log-return2 Reg. Coeff.", lwd = 2)
lines(betafd + betastderrfd, lty=2, lwd=1)
lines(betafd - betastderrfd, lty=2, lwd=1)

# Residuals evaluation for penalized model
volatility_res_pen = as.numeric(volatility - sqrt(volatility_squaredhat_pen))
# Scalar-on-function Model for residuals
fRegressList_pen_RESID = fRegress(log(abs(volatility_res_pen)), logret_list, betalist_pen)

#####
# Validation
#####
# Validation data: True RV(ti+10min)^2
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$validation, ]
vol <- as.numeric(vol$target)^2
# Order book data
df_book_data_validation <- df_book_data[df_book_data$time_id %in% split_results$validation, ]

# Column selection
column_name <- "log_return" 
# Dataframe creation
df_dynamic <- df_book_data_validation[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
# Conversion of df_dynamic in data.table for a faster process
df_dynamic <- as.data.table(df_dynamic)
# Matrix creation (600 rows and n_obs columns)
time_istants <- unique(df_book_data_validation$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
# Multiple NAs at the beginning and at the end of the matrix columns reach to some
# problems in the interpolation phase
# Counting NA at the beginning
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
# Counting NA at the end
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)

df_matrix <- as.matrix(df)
# Substitution of NAs with 0s which mean no-changes in prices
for (i in 1:length(na_start)){
  if(na_start[i]!=0){
    if (na_start[i]>1)
      df_matrix[1:na_start[i], i]<-0
    else
      df_matrix[1, i]<-0
  }
}
for (i in 1:length(na_end)) {
  if (na_end[i] != 0) {
    if (na_end[i] > 1) {
      df_matrix[(600 - na_end[i] + 1):600, i] <- 0
    } else {
      df_matrix[600, i] <- 0
    }
  }
}
Xobs0 <- df_matrix
# Spline interpolation
Xobs0 <- na.spline(Xobs0)

# New basis for X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

# B(t) extraction from training phase and prediction
beta_est <- fRegressList_pen$betaestlist[[2]]$fd
vol_pred_validation <- sqrt(inprod(logretfd, beta_est))

# B(t) extraction from training phase and residuals prediction
beta_est_res <- fRegressList_pen_RESID$betaestlist[[2]]$fd
res_pred_validation <- inprod(logretfd, beta_est_res)

# Target selection for the chosen stock
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$validation, ]
# Not squared
volatility <- as.numeric(vol$target)

# Non-conformity measure computation
non_conformity_score_validation <- abs(vol_pred_validation - volatility)/exp(res_pred_validation)

quan_value_of_ncscore_90 <- as.numeric(quantile(non_conformity_score_validation, 0.90))
quan_value_of_ncscore_95 <- as.numeric(quantile(non_conformity_score_validation, 0.95))
quan_value_of_ncscore_99 <- as.numeric(quantile(non_conformity_score_validation, 0.99))

#####
# Test
#####
# Test data: True RV(ti+10min)^2 
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$test, ]
vol <- as.numeric(vol$target)^2
# Order book data
df_book_data_test <- df_book_data[df_book_data$time_id %in% split_results$test, ]

# Column selection
column_name <- "log_return" 
# Dataframe creation
df_dynamic <- df_book_data_test[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
# Conversion of df_dynamic in data.table for a faster process
df_dynamic <- as.data.table(df_dynamic)
# Matrix creation (600 rows and n_obs columns)
time_istants <- unique(df_book_data_test$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
# Multiple NAs at the beginning and at the end of the matrix columns reach to some
# problems in the interpolation phase
# Counting NA at the beginning
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
# Counting NA at the end
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)

df_matrix <- as.matrix(df)
# Substitution of NAs with 0s which mean no-changes in prices
for (i in 1:length(na_start)){
  if(na_start[i]!=0){
    if (na_start[i]>1)
      df_matrix[1:na_start[i], i]<-0
    else
      df_matrix[1, i]<-0
  }
}
for (i in 1:length(na_end)) {
  if (na_end[i] != 0) {
    if (na_end[i] > 1) {
      df_matrix[(600 - na_end[i] + 1):600, i] <- 0
    } else {
      df_matrix[600, i] <- 0
    }
  }
}
Xobs0 <- df_matrix
# Spline interpolation
Xobs0 <- na.spline(Xobs0)

# New basis for X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

# Prediction based on the fitted model in the training phase 
vol_pred_test <- sqrt(inprod(logretfd, beta_est))
# Residuals prediction 
res_pred_test <- exp(inprod(logretfd, beta_est_res))

# Target selection for the chosen stock
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$test, ]
# Not squared
true_volatility_test <- as.numeric(vol$target)
# RMSE evaluation based on test data
RMSE.test <- rmse(true_volatility_test, vol_pred_test)

# Conformal intervals construction
# Conformal intervals at 90%
cp_lower_90 <-  vol_pred_test - quan_value_of_ncscore_90 * res_pred_test
cp_upper_90 <-  vol_pred_test + quan_value_of_ncscore_90 * res_pred_test
cp_lower_90[cp_lower_90<0]<-0

# Conformal intervals at 95%
cp_lower_95 <-  vol_pred_test - quan_value_of_ncscore_95 * res_pred_test
cp_upper_95 <-  vol_pred_test + quan_value_of_ncscore_95 * res_pred_test
cp_lower_95[cp_lower_95<0]<-0

# Conformal intervals at 99%
cp_lower_99 <-  vol_pred_test - quan_value_of_ncscore_99 * res_pred_test
cp_upper_99 <-  vol_pred_test + quan_value_of_ncscore_99 * res_pred_test
cp_lower_99[cp_lower_99<0]<-0

#####
# CP diagnosis and representation
#####
test_buckets <- time_buckets[time_buckets %in% split_results$test]

# Widths evaluation
width_90 <- cp_upper_90 - cp_lower_90
width_95 <- cp_upper_95 - cp_lower_95
width_99 <- cp_upper_99 - cp_lower_99

# Dataframe for ggplot, saved and exported for the representation
df_width <- data.frame(
  time = test_buckets,
  width_90 = width_90,
  width_95 = width_95,
  width_99 = width_99
)
saveRDS(df_width, file = "df_width_sf.rds")

cat(sprintf("Medium width of the 90%% conformal region: %.8f\n", mean(width_90)))
cat(sprintf("Medium width of the 95%% conformal region: %.8f\n", mean(width_95)))
cat(sprintf("Medium width of the 99%% conformal region: %.8f\n", mean(width_99)))

# Coverages evaluation
num_points_outside_0 <- sum(true_volatility_test < cp_lower_90 | true_volatility_test > cp_upper_90)
coverage_score_0 <- 1.0 - num_points_outside_0 / length(true_volatility_test)

num_points_outside_1 <- sum(true_volatility_test < cp_lower_95 | true_volatility_test > cp_upper_95)
coverage_score_1 <- 1.0 - num_points_outside_1 / length(true_volatility_test)

num_points_outside_2 <- sum(true_volatility_test < cp_lower_99 | true_volatility_test > cp_upper_99)
coverage_score_2 <- 1.0 - num_points_outside_2 / length(true_volatility_test)

cat(sprintf("Coverage Score of the 90%% conformal region: %.3f%%\n", coverage_score_0 * 100))
cat(sprintf("Coverage Score of the 95%% conformal region: %.3f%%\n", coverage_score_1 * 100))
cat(sprintf("Coverage Score of the 99%% conformal region: %.3f%%\n", coverage_score_2 * 100))

# Inside/Outside flags for the true volatility values
inside_interval_90 <- true_volatility_test >= cp_lower_90 & true_volatility_test <= cp_upper_90
inside_interval_95 <- true_volatility_test >= cp_lower_95 & true_volatility_test <= cp_upper_95
inside_interval_99 <- true_volatility_test >= cp_lower_99 & true_volatility_test <= cp_upper_99

# Dataframe for the CP plot
pred_vol_test_from_linear = df_past_realized_train[df_past_realized_train$row_id %in% split_results$test, ]$pred
df_plots <- data.frame(true_volatility_test, vol_pred_test, pred_vol_test_from_linear, cp_lower_95, cp_upper_95, inside_interval_95, cp_lower_99, cp_upper_99, inside_interval_99,  cp_lower_90, cp_upper_90, inside_interval_90)

# Sample selection of the regions
df_plots <- df_plots %>% sample_n(100)

# Global range for the representation
global_min <- min(c(df_plots$cp_lower_90, df_plots$cp_upper_90, 
                    df_plots$cp_lower_95, df_plots$cp_upper_95, 
                    df_plots$cp_lower_99, df_plots$cp_upper_99))
global_max <- max(c(df_plots$cp_lower_90, df_plots$cp_upper_90, 
                    df_plots$cp_lower_95, df_plots$cp_upper_95, 
                    df_plots$cp_lower_99, df_plots$cp_upper_99))

# Common step
common_y_step <- 0.005  

p90 <- ggplot(df_plots, aes(x = pred_vol_test_from_linear, y = true_volatility_test)) +
  geom_point(aes(color = inside_interval_90), size = 3) +
  geom_errorbar(aes(ymin = cp_lower_90, ymax = cp_upper_90), width = 0.0005, color = "#004777") +
  scale_color_manual(values = c("red", "#004777"), labels = c("Outside Interval", "Inside Interval")) +
  labs(title = "Conformal Regions 
       with a Significance Level of 90%",
       x = "RV(ti)", y = "Predicted Values of RV(ti+10min)",
       color = "") +
  scale_y_continuous(
    breaks = seq(global_min, global_max, by = common_y_step),  
    limits = c(global_min, global_max)  
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "top",  # Muovere la legenda in alto
    legend.text = element_text(size = 14)
  )

p95 <- ggplot(df_plots, aes(x = pred_vol_test_from_linear, y = true_volatility_test)) +
  geom_point(aes(color = inside_interval_95), size = 3) +
  geom_errorbar(aes(ymin = cp_lower_95, ymax = cp_upper_95), width = 0.0005, color = "#004777") +
  scale_color_manual(values = c("red", "#004777"), labels = c("Outside Interval", "Inside Interval")) +
  labs(title = "Conformal Regions with 
       a Significance Level of 95%",
       x = "RV(ti)", y = "Predicted Values of RV(ti+10min)",
       color = "") +
  scale_y_continuous(
    breaks = seq(global_min, global_max, by = common_y_step),
    limits = c(global_min, global_max)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "top",  
    legend.text = element_text(size = 14)
  )

p99 <- ggplot(df_plots, aes(x = pred_vol_test_from_linear, y = true_volatility_test)) +
  geom_point(aes(color = inside_interval_99), size = 3) +
  geom_errorbar(aes(ymin = cp_lower_99, ymax = cp_upper_99), width = 0.0005, color = "#004777") +
  scale_color_manual(values = c("red", "#004777"), labels = c("Outside Interval", "Inside Interval")) +
  labs(title = "Conformal Regions with 
       a Significance Level of 99%",
       x = "RV(ti)", y = "Predicted Values of RV(ti+10min)",
       color = "") +
  scale_y_continuous(
    breaks = seq(global_min, global_max, by = common_y_step),
    limits = c(global_min, global_max)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "top",  
    legend.text = element_text(size = 14)
  )

# Visualization of the three different regions
p90+p95+p99
