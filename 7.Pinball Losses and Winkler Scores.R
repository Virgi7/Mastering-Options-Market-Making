# Pinball Loss and Winkler Score w.r. to different coverages
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
library(qgam)
library(patchwork)
source("fRegress std error mia versione.R")
########################################
###### Scalar-on-function for RV2 ######
########################################
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
#####
# Dataset preparation
#####
df_book_data <- sqlogret_per_time_id(file_paths)

time_buckets <- train[train$stock_id == stock_id_chosen, "time_id", drop=FALSE]
time_buckets <- as.vector(time_buckets$time_id)

# Training-Validation-Test Splitting
split_time_buckets_no_replacement <- function(time_buckets, train_size = 0.5, val_size = 0.3, test_size = 0.2, seed = 123) {
  
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
# Training data
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$train, ]
vol <- as.numeric(vol$target)^2

df_book_data_train <- df_book_data[df_book_data$time_id %in% split_results$train, ]
column_name <- "log_return" 
df_dynamic <- df_book_data_train[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_train$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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

# Target volatility of the training, not squared
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
volatility <- as.numeric(vol[vol$time_id %in% split_results$train, ]$target)

volatility_squaredhat_nonpen = fRegressList$yhatfdobj
volatility_squaredres_nonpen = volatility_squared - volatility_squaredhat_nonpen

SSE1.1 = sum(volatility_squaredres_nonpen^2)
SSE0 = sum((volatility_squared - mean(volatility_squared))^2)
RSQ1 = (SSE0 - SSE1.1)/SSE0
Fratio1 = ((SSE0 - SSE1.1)/5)/(SSE1.1/29)
RMSE.1 = rmse(volatility, sqrt(volatility_squaredhat_nonpen))

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

# Penalized form of the model
lambda0 = 2^2
betafdPar = fdPar(beta_basis, 2, lambda0)
betalist[[2]] = betafdPar

loglam = seq(-15,0,0.5)
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

# Residuals fitting
volatility_res_pen = as.numeric(volatility - sqrt(volatility_squaredhat_pen))
fRegressList_pen_RESID = fRegress(log(abs(volatility_res_pen)), logret_list, betalist_pen)

#####
# Validation
#####
# Validation data
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$validation, ]
vol <- as.numeric(vol$target)^2

df_book_data_validation <- df_book_data[df_book_data$time_id %in% split_results$validation, ]
column_name <- "log_return"  
df_dynamic <- df_book_data_validation[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_validation$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

# X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

# B(t)
beta_est <- fRegressList_pen$betaestlist[[2]]$fd
vol_pred_validation <- sqrt(inprod(logretfd, beta_est))
beta_est_res <- fRegressList_pen_RESID$betaestlist[[2]]$fd
res_pred_validation <- inprod(logretfd, beta_est_res)

vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$validation, ]
volatility <- as.numeric(vol$target)

# Non-conformity measure
non_conformity_score_validation <- abs(vol_pred_validation - volatility)/exp(res_pred_validation)

coverage <- seq(0,1, 0.02)
quan_value_of_ncscore_alfa <- as.numeric(quantile(non_conformity_score_validation, coverage))

#####
# Test
#####
# Test data
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$test, ]
vol <- as.numeric(vol$target)^2

df_book_data_test <- df_book_data[df_book_data$time_id %in% split_results$test, ]
column_name <- "log_return" 
df_dynamic <- df_book_data_test[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_test$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

# X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

vol_pred_test <- sqrt(inprod(logretfd, beta_est))
res_pred_test <- exp(inprod(logretfd, beta_est_res))

vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$test, ]
true_volatility_test <- as.numeric(vol$target)

RMSE.test <- rmse(true_volatility_test, vol_pred_test)

# Pinball Loss e Winkler score
cp_lower_alfa <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_volatility_test)) 
cp_upper_alfa <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_volatility_test)) 
Pinball_Loss_alfa <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_volatility_test)) 
Winkler_Score_alfa <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_volatility_test)) 

winkler_score <- function(true_value, lower_bound, upper_bound, alpha) {
  if (true_value >= lower_bound && true_value <= upper_bound) {
    # Case 1: Real value is in the interval
    return(upper_bound - lower_bound)
  } else if (true_value < lower_bound) {
    # Case 2: Real value is lower than the interval
    return((upper_bound - lower_bound) + (2 / alpha) * (lower_bound - true_value))
  } else {
    # Case 3: Real value is higher than the interval
    return((upper_bound - lower_bound) + (2 / alpha) * (true_value - upper_bound))
  }
}

for (i in 1:length(quan_value_of_ncscore_alfa)) {
  cp_lower_alfa[i, ] <- vol_pred_test - quan_value_of_ncscore_alfa[i] * res_pred_test
  cp_upper_alfa[i, ] <- vol_pred_test + quan_value_of_ncscore_alfa[i] * res_pred_test
  
  if (coverage[i] >= 0.50) {
    Pinball_Loss_alfa[i, ] <- pinLoss(true_volatility_test, cp_upper_alfa[i, ], coverage[i])
  } else {
    Pinball_Loss_alfa[i, ] <- pinLoss(true_volatility_test, cp_lower_alfa[i, ], 1 - coverage[i])
  }
  
  alfa_wink = 1 - coverage[i]
  Winkler_Score_alfa[i, ] <- mapply(winkler_score, true_volatility_test, cp_lower_alfa[i, ], cp_upper_alfa[i, ], MoreArgs = list(alpha = 1 - alfa_wink))
}

############################
####### Linear Model #######
############################
#####
# Utilities
#####

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

#####
# Training
#####
# Training data
true_vol_train <- as.numeric(df_joined %>% filter(row_id %in% split_results$train) %>% pull(target))
pred_vol_train <- df_joined %>% filter(row_id %in% split_results$train) %>% pull(pred)

df_book_data_train <- df_book_data[df_book_data$time_id %in% split_results$train, ]
column_name <- "log_return" 
df_dynamic <- df_book_data_train[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_train$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

n_obs = length(true_vol_train)
abscissa = 0:599

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

# Linear Model
model_train <- lm(true_vol_train ~ pred_vol_train)
res_train <- resid(model_train)

# Residuals fitting
fRegressList_RESID = fRegress(log(abs(res_train)), logret_list, betalist)

#####
# Validation
#####
# Validation data
true_vol_valid <- as.numeric(df_joined %>% filter(row_id %in% split_results$validation) %>% pull(target))
pred_vol_valid <- df_joined %>% filter(row_id %in% split_results$validation) %>% pull(pred)

df_book_data_valid <- df_book_data[df_book_data$time_id %in% split_results$validation, ]
column_name <- "log_return"  
df_dynamic <- df_book_data_valid[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_valid$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

vol_predicted_valid <- predict(model_train, newdata = data.frame(pred_vol_train = pred_vol_valid))

# X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd
beta_est_res <- fRegressList_RESID$betaestlist[[2]]$fd
res_pred_validation <- inprod(logretfd, beta_est_res)
non_conformity_score_valid <- abs(vol_predicted_valid - true_vol_valid) / exp(res_pred_validation) 

coverage <- seq(0, 1, 0.02)
quan_value_of_ncscore_alfa <- as.numeric(quantile(non_conformity_score_valid, coverage))

#####
# Test
#####
# Test data
true_vol_test <- as.numeric(df_joined %>% filter(row_id %in% split_results$test) %>% pull(target))
pred_vol_test <- df_joined %>% filter(row_id %in% split_results$test) %>% pull(pred)

df_book_data_test <- df_book_data[df_book_data$time_id %in% split_results$test, ]
column_name <- "log_return"  
df_dynamic <- df_book_data_test[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_test$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

# X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

vol_predicted_test <- predict(model_train, newdata = data.frame(pred_vol_train = pred_vol_test))
res_predicted_test <- exp(inprod(logretfd, beta_est_res))
RMSE.test <- rmse(true_vol_test, vol_predicted_test)

# Pinball Loss e Winkler score
cp_lower_alfa_linear <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_vol_test)) 
cp_upper_alfa_linear <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_vol_test)) 
Pinball_Loss_alfa_linear <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_vol_test)) 
Winkler_Score_alfa_linear <- matrix(NA, length(quan_value_of_ncscore_alfa), length(true_vol_test)) 


for (i in 1:length(quan_value_of_ncscore_alfa)) {
  cp_lower_alfa_linear[i, ] <- vol_pred_test - quan_value_of_ncscore_alfa[i] * res_predicted_test
  cp_upper_alfa_linear[i, ] <- vol_pred_test + quan_value_of_ncscore_alfa[i] * res_predicted_test
  
  if (coverage[i] >= 0.50) {
    Pinball_Loss_alfa_linear[i, ] <- pinLoss(true_vol_test, cp_upper_alfa_linear[i, ], coverage[i])
  } else {
    Pinball_Loss_alfa_linear[i, ] <- pinLoss(true_vol_test, cp_lower_alfa_linear[i, ], 1 - coverage[i])
  }
  
  alfa_wink = 1 - coverage[i]
  Winkler_Score_alfa_linear[i, ] <- mapply(winkler_score, true_vol_test, cp_lower_alfa_linear[i, ], cp_upper_alfa_linear[i, ], MoreArgs = list(alpha = 1 - alfa_wink))
}
Pinball_Loss_alfa_linear[1,]<-rep(0,766)
#########################################
#### Scalar-on-function for DeltaRV #####
#########################################
#####
# Dataset preparation
#####

df_past_realized <- past_realized_volatility_per_stock(list_file = file_paths,
                                                       prediction_column_name = 'pred')

df_book_data <- sqlogret_per_time_id(file_paths)
time_buckets <- train[train$stock_id == stock_id_chosen, "time_id", drop=FALSE]
time_buckets <- as.vector(time_buckets$time_id)
df_past_realized$row_id <- time_buckets
#####
# Training
#####
# Training data
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$train, ]
vol <- as.numeric(vol$target)
df_past_realized_train <- df_past_realized[df_past_realized$row_id %in% split_results$train, ]
delta_vol_train <- vol - df_past_realized_train$pred

df_book_data_train <- df_book_data[df_book_data$time_id %in% split_results$train, ]
column_name <- "log_return" 
df_dynamic <- df_book_data_train[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_train$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

n_obs = length(vol)
abscissa = 0:599

# X(t)
logret_list = vector("list", 2)
logret_list[[1]] = rep(1, n_obs)
nbasis_beta = 10
grade = 4
m = grade
degree = m-1
nbasis_covariates = 120 
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
fRegressList = fRegress(delta_vol_train, logret_list, betalist)
betaestlist_nonpen = fRegressList$betaestlist

delta_volhat_nonpen = fRegressList$yhatfdobj
delta_volres_nonpen = delta_vol_train - delta_volhat_nonpen

SSE1.1 = sum(delta_volres_nonpen^2)
SSE0 = sum((delta_vol_train - mean(delta_vol_train))^2)
RSQ1 = (SSE0 - SSE1.1)/SSE0
Fratio1 = ((SSE0 - SSE1.1)/5)/(SSE1.1/29)
volatility_hat_train <- delta_volhat_nonpen + df_past_realized_train$pred
RMSE.1 = rmse(volatility_hat_train, vol)

sigmaE. = sum(delta_volres_nonpen^2)/(n_obs - fRegressList$df)
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

# Penalized version of the model
lambda0 = 2^2
betafdPar = fdPar(beta_basis, 2, lambda0)
betalist[[2]] = betafdPar

loglam = seq(-15,0,0.5)
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
  fReg_i = fRegress(delta_vol_train, logret_list, betalist_i)
  betaestlist_i = fReg_i$betaestlist
  logret_betafd_CV[ilam] = betaestlist_i[[2]]$fd
  delta_volhat_CV = fReg_i$yhatfdobj
  delta_volres_CV = delta_vol_train - delta_volhat_CV
  sigmaE. = sum(delta_volres_CV^2)/(n_obs - fReg_i$df)
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
  SSE1.CV[ilam] = sum(delta_volres_CV^2)
  volatility_hat_train = delta_volhat_CV + df_past_realized_train$pred
  RMSE.CV[ilam] = rmse(volatility_hat_train, vol)
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
fRegressList_pen = fRegress(delta_vol_train, logret_list, betalist_pen)
betaestlist_pen = fRegressList_pen$betaestlist
logret_betafd_pen = betaestlist_pen[[2]]$fd
delta_volhat_pen = fRegressList_pen$yhatfdobj
delta_volres_pen = delta_vol_train - delta_volhat_pen
sigmaE. = sum(delta_volres_pen^2)/(n_obs - fRegressList_pen$df)
sigmaE = sigmaE. * diag(rep(1, n_obs))
y2cMap = logretSmooth$y2cMap
stderrList = fRegress.stderr(fRegressList_pen, y2cMap, sigmaE)
betafdPar = betaestlist_pen[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
volatility_hat_train = delta_volhat_pen + df_past_realized_train$pred
RMSE.2 = rmse(volatility_hat_train, vol)
RMSE.2 = rmse(delta_vol_train, delta_volhat_pen)
SSE1.2 = sum(delta_volres_pen^2)
RSQ2 = (SSE0 - SSE1.2)/SSE0
# Final plot
plot(betafd, xlab = "Seconds", ylab = "Log-return2 Reg. Coeff.", lwd = 2)
lines(betafd + betastderrfd, lty=2, lwd=1)
lines(betafd - betastderrfd, lty=2, lwd=1)

delta_volres_pen_real = as.numeric(delta_volhat_pen + df_past_realized_train$pred - vol) 
fRegressList_pen_RESID = fRegress(log(abs(delta_volres_pen_real)), logret_list, betalist_pen)

#####
# Validation
#####
# Validation data
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$validation, ]
vol <- as.numeric(vol$target)
df_past_realized_validation <- df_past_realized[df_past_realized$row_id %in% split_results$validation, ]
delta_vol_validation <- vol - df_past_realized_validation$pred

df_book_data_validation <- df_book_data[df_book_data$time_id %in% split_results$validation, ]
column_name <- "log_return" 
df_dynamic <- df_book_data_validation[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_validation$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

# X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

# B(t)
beta_est <- fRegressList_pen$betaestlist[[2]]$fd
delta_vol_pred_validation <- inprod(logretfd, beta_est)
beta_est_res <- fRegressList_pen_RESID$betaestlist[[2]]$fd
res_pred_validation <- inprod(logretfd, beta_est_res)

# Non-conformity measure
vol_res_pred_valid <- delta_vol_pred_validation + df_past_realized_validation$pred - vol
non_conformity_score_validation <- abs(vol_res_pred_valid)/exp(res_pred_validation)

quan_value_of_ncscore_alfa <- as.numeric(quantile(non_conformity_score_validation, coverage))
#####
# Test
#####
# Test data
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
vol <- vol[vol$time_id %in% split_results$test, ]
vol <- as.numeric(vol$target)
df_past_realized_test <- df_past_realized[df_past_realized$row_id %in% split_results$test, ]
delta_vol_test <- vol - df_past_realized_test$pred

df_book_data_test <- df_book_data[df_book_data$time_id %in% split_results$test, ]
column_name <- "log_return" 
df_dynamic <- df_book_data_test[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data_test$time_id)
df <- matrix(NA, nrow = 600, ncol = length(time_istants))
colnames(df) <- time_istants
rownames(df) <- 1:600
df_dynamic[, {
  df[seconds_in_bucket + 1, as.character(time_id)] <<- value
}, by = .(time_id)]
count_na_start <- function(x) {
  sum(cumprod(is.na(x)))
}
count_na_end <- function(x) {
  sum(cumprod(rev(is.na(x))))
}
na_start <- apply(df, 2, count_na_start)
na_end <- apply(df, 2, count_na_end)
df_matrix <- as.matrix(df)
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
Xobs0 <- na.spline(Xobs0)

# X(t)
new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
logretfd = logretSmooth$fd

vol_pred_test <- inprod(logretfd, beta_est) + df_past_realized_test$pred
res_pred_test <- exp(inprod(logretfd, beta_est_res))
RMSE.test <- rmse(vol_pred_test, vol)


# Pinball Loss e Winkler score
cp_lower_alfa <- matrix(NA, length(quan_value_of_ncscore_alfa), length(vol)) 
cp_upper_alfa <- matrix(NA, length(quan_value_of_ncscore_alfa), length(vol)) 
Pinball_Loss_alfa_delta <- matrix(NA, length(quan_value_of_ncscore_alfa), length(vol)) 
Winkler_Score_alfa_delta <- matrix(NA, length(quan_value_of_ncscore_alfa), length(vol)) 

for (i in 1:length(quan_value_of_ncscore_alfa)) {
  cp_lower_alfa[i, ] <- vol_pred_test - quan_value_of_ncscore_alfa[i] * res_pred_test
  cp_upper_alfa[i, ] <- vol_pred_test + quan_value_of_ncscore_alfa[i] * res_pred_test
  
  if (coverage[i] >= 0.50) {
    Pinball_Loss_alfa_delta[i, ] <- pinLoss(true_volatility_test, cp_upper_alfa[i, ], coverage[i])
  } else {
    Pinball_Loss_alfa_delta[i, ] <- pinLoss(true_volatility_test, cp_lower_alfa[i, ], 1 - coverage[i])
  }
  
  alfa_wink = 1 - coverage[i]
  Winkler_Score_alfa_delta[i, ] <- mapply(winkler_score, true_volatility_test, cp_lower_alfa[i, ], cp_upper_alfa[i, ], MoreArgs = list(alpha = 1 - alfa_wink))
}

#####
# Representation of Pinball losses and Winkler scores across the three models
#####
# Preparation of Winkler score dataframe
df_winkler <- data.frame(
  coverage = coverage[3:length(coverage)-1],
  winkler_score_1 = rowMeans(Winkler_Score_alfa_linear[3:length(coverage)-1,]),  
  winkler_score_2 = rowMeans(Winkler_Score_alfa[3:length(coverage)-1,]),
  winkler_score_3 = rowMeans(Winkler_Score_alfa_delta[3:length(coverage)-1,])
)
df_winkler_long <- df_winkler %>%
  pivot_longer(cols = c(winkler_score_1, winkler_score_2, winkler_score_3), 
               names_to = "winkler_type", 
               values_to = "winkler_score") %>%
  mutate(winkler_type = case_when(
    winkler_type == "winkler_score_1" ~ "Winkler score Linear",
    winkler_type == "winkler_score_2" ~ "Winkler score Scalar-on-function for RV²",
    winkler_type == "winkler_score_3" ~ "Winkler score Scalar-on-function for ΔRV"
  ))

df_winkler_long <- df_winkler_long %>%
  mutate(winkler_type = factor(winkler_type, levels = c(
    "Winkler score Linear", 
    "Winkler score Scalar-on-function for RV²",
    "Winkler score Scalar-on-function for ΔRV"
  )))

# Preparation of Pinball loss dataframe
df_pinball <- data.frame(
  coverage = coverage,
  pinball_loss_1 = rowMeans(Pinball_Loss_alfa_linear),
  pinball_loss_2 = rowMeans(Pinball_Loss_alfa),
  pinball_loss_3 = rowMeans(Pinball_Loss_alfa_delta)   
)

df_pinball_long <- df_pinball %>%
  pivot_longer(cols = c(pinball_loss_1, pinball_loss_2, pinball_loss_3), 
               names_to = "pinball_type", 
               values_to = "pinball_loss") %>%
  mutate(pinball_type = case_when(
    pinball_type == "pinball_loss_1" ~ "Pinball loss Linear",
    pinball_type == "pinball_loss_2" ~ "Pinball loss Scalar-on-function for RV²",
    pinball_type == "pinball_loss_3" ~ "Pinball loss Scalar-on-function for ΔRV"
  ))

df_pinball_long <- df_pinball_long %>%
  mutate(pinball_type = factor(pinball_type, levels = c(
    "Pinball loss Linear", 
    "Pinball loss Scalar-on-function for RV²", 
    "Pinball loss Scalar-on-function for ΔRV"
  )))


color_palette <- c("#ffdc7c", "#6ba292", "#89043D")

# Combined plots for Pinball loss
ggplot(df_pinball_long, aes(x = coverage, y = pinball_loss)) +
  geom_line(aes(color = pinball_type, linetype = pinball_type), size = 1.2) +  
  geom_point(aes(fill = pinball_type, shape = pinball_type), size = 3, stroke = 1, color = "black") +  
  
  labs(
    title = "Comparison of Pinball Losses Across Coverages",
    x = "Coverage",
    y = "Pinball Loss"
  ) +
  
  theme_light(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"), 
    legend.text = element_text(size = 14),
    panel.border = element_blank()
  ) +
  
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  ylim(0, max(df_pinball_long$pinball_loss) * 1.1)


# Combined plots for Winkler score
ggplot(df_winkler_long, aes(x = coverage, y = winkler_score)) +
  geom_line(aes(color = winkler_type, linetype = winkler_type), size = 1.2) +  
  geom_point(aes(fill = winkler_type, shape = winkler_type), size = 3, stroke = 1, color = "black") +  
  
  labs(
    title = "Comparison of Winkler Scores Across Coverages",
    x = "Coverage",
    y = "Winkler Score"
  ) +
  
  theme_light(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16), 
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"), 
    legend.text = element_text(size = 14),
    panel.border = element_blank()
  ) +
  
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  ylim(0, max(df_winkler_long$winkler_score) * 1.1) +
  scale_y_log10()

