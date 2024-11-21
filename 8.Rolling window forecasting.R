# PROVO A IMPOSTARE LA MOVING WINDOW : DALL'ANALISI DELLA VOLATILITA OTTENIAMO
# CHE LE NOSTRE DATE VANNO DAL 1 GENNAIO 2020 AL 31 DICEMBRE 2020 E CHE PER
# OGNIU GIORNO ABBIAMO 12 RILEVAZIONI DI BUCKET DI 10 MIN, QUINDI QUELLO
# CHE POSSIAMO FARE È PREDIRRE UN GIORNO PER QUELLO SUCCESSIVO
# Scalar Delta(vol) on functional model with conformal prediction
## Librerie
###########
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
library(grid)
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
stock_id_chosen <- 41 #sample(stock_id_unique, 1)

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

# RV(ti) calculation for each bucket in different files
past_realized_volatility_per_stock <- function(list_file, prediction_column_name) {
  df_past_realized <- data.table()  
  
  for (file in list_file) {
    df_past_realized <- rbind(df_past_realized,
                              realized_volatility_per_time_id(file, prediction_column_name))
  }
  
  return(df_past_realized)
}

# Functions to evaluate the accuracy of the conformal regions and compute the final quantities
# threshold identifies the minimum difference between realized and implied for which it is possible
# to identify a trade
check_volatility_95 <- function(realized_vol, implied_vol, threshold = 5) {
  diff_check <- abs(realized_vol - implied_vol) > threshold
  
  interval_check <- ifelse(implied_vol >= df_test_predictions$cp_lower_95 & 
                             implied_vol <= df_test_predictions$cp_upper_95, 1, 0)
  
  return(list(difference_check = diff_check, interval_check = interval_check))
}
check_volatility_99 <- function(realized_vol, implied_vol, threshold = 5) {
  diff_check <- abs(realized_vol - implied_vol) > threshold
  
  interval_check <- ifelse(implied_vol >= df_test_predictions$cp_lower_99 & 
                             implied_vol <= df_test_predictions$cp_upper_99, 1, 0)
  
  return(list(difference_check = diff_check, interval_check = interval_check))
}
# Counting the matches of realized/non-realized trade and possible/non-possible trade
count_matches <- function(result_95) {
  match_count <- 0
  for (element in result_95) {
    list_tf <- element[[1]]
    list_01 <- element[[2]]
    
    for (i in 1:length(list_01)) {
      if (list_01[i] == 0 && list_tf[i] == TRUE) {
        match_count <- match_count + 1
      }
    }
  }
  
  return(match_count)
}
#####
# Dataset preparation
#####
#  RV(ti) calculation for each bucket
df_past_realized <- past_realized_volatility_per_stock(list_file = file_paths,
                                                       prediction_column_name = 'pred')


df_book_data <- sqlogret_per_time_id(file_paths)
time_buckets <- train[train$stock_id == stock_id_chosen, "time_id", drop=FALSE]
time_buckets <- as.vector(time_buckets$time_id)
time_buckets <- time_buckets
df_past_realized$row_id <- time_buckets

# Correct time_id order 
time_buckets_ordered <- read_excel('correct_order_time_id.xlsx')

# Implied volatilities of the 3 matched stocks
IV_stocks <- read_excel('Stocks IV.xlsx')

# Reordering of the volatilities with respect to the time order
vol_realized <- as.data.frame(train[train$stock_id == stock_id_chosen, c("time_id", "target")])
vol_realized_sorted <- vol_realized[match(time_buckets_ordered$time_id, vol_realized$time_id), ]
df_past_realized_sorted <- df_past_realized[match(time_buckets_ordered$time_id, df_past_realized$row_id), ]
date <- IV_stocks$`CSCO US Equity`

# Selection of the interesting date: from 01/01/2020 to 31/12/2020 and assignment of them
date_rep <- rep(date, each = 12)
vol_realized_sorted_imp <- vol_realized_sorted[1:3012,]
vol_realized_sorted_imp$date <- date_rep
df_past_realized_sorted_imp <- df_past_realized_sorted[1:3012,]
df_past_realized_sorted_imp$date <- date_rep
df_book_data_sorted <- df_book_data %>%
  filter(time_id %in% vol_realized_sorted_imp$time_id) %>%
  filter(time_id %in% time_buckets_ordered$time_id) %>%
  arrange(factor(time_id, levels = time_buckets_ordered$time_id))
df_dates <- data.frame(time_id=vol_realized_sorted_imp$time_id, date=vol_realized_sorted_imp$date)
df_book_data_sorted <- df_book_data_sorted %>%
  left_join(df_dates, by = "time_id")

# Rolling window selection
rolling_window_forecast <- function(time_buckets, train_size, val_size, test_size, num_windows) {
  
  n <- length(time_buckets)
  
  windows <- list()
  
  if (train_size + val_size + test_size > n) {
    stop("Sum of train_size, val_size e test_size has to be equal to the number of observations")
  }
  
  for (i in 0:(num_windows - 1)) {
    start_index <- i * test_size + 1
    end_train_index <- start_index + train_size - 1
    end_val_index <- end_train_index + val_size
    
    if (end_val_index > n) break 
    
    train_window <- time_buckets[start_index:end_train_index]
    val_window <- time_buckets[(end_train_index + 1):end_val_index]
    test_window <- time_buckets[(end_val_index + 1):(end_val_index + test_size)]
    
    if (length(test_window) != test_size) break
    
    windows[[i + 1]] <- list(train = train_window, val = val_window, test = test_window)
  }
  
  return(windows)
}
# Estraction of train, validation and test
time_buckets <- IV_stocks$`CSCO US Equity`
time_buckets_test <- time_buckets[211:length(time_buckets)]
train_size = (length(time_buckets)-length(time_buckets_test))*0.6
val_size = (length(time_buckets)-length(time_buckets_test))*0.4
train_test_splits <- rolling_window_forecast(time_buckets, train_size = train_size, val_size = val_size, test_size = 1, num_windows = length(date) - train_size - val_size)
output_directory <- "plots_csco"
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}
results_list_95 <- vector("list", length(train_test_splits))
results_list_99 <- vector("list", length(train_test_splits))


for (rolling_window in 1:length(train_test_splits)) {
  
  #####
  ## Training
  #####
  # Training data
  vol <- vol_realized_sorted_imp[vol_realized_sorted_imp$date %in% train_test_splits[[rolling_window]]$train, ]
  vol <- as.numeric(vol$target)
  df_past_realized_train <- df_past_realized_sorted_imp[df_past_realized_sorted_imp$date %in% train_test_splits[[rolling_window]]$train, ]
  delta_vol_train <- vol - df_past_realized_train$pred
  
  df_book_data_train <- df_book_data_sorted[df_book_data_sorted$date %in% train_test_splits[[rolling_window]]$train, ]
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
  fRegressList = fRegress(delta_vol_train, logret_list, betalist)
  betaestlist_nonpen = fRegressList$betaestlist
  logret_betafd_nonpen = betaestlist_nonpen[[2]]$fd

  delta_volhat_nonpen = fRegressList$yhatfdobj
  # Residuals evaluation
  delta_volres_nonpen = delta_vol_train - delta_volhat_nonpen
  
  # Diagnosis
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
  
  # Penalized form of the model
  lambda0 = 2^2
  betafdPar = fdPar(beta_basis, 2, lambda0)
  betalist[[2]] = betafdPar
  
  # Cross-validation
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
  vol <- vol_realized_sorted_imp[vol_realized_sorted_imp$date %in% train_test_splits[[rolling_window]]$val, ]
  vol <- as.numeric(vol$target)
  df_past_realized_validation <- df_past_realized_sorted_imp[df_past_realized_sorted_imp$date %in% train_test_splits[[rolling_window]]$val, ]
  delta_vol_validation <- vol - df_past_realized_validation$pred
  
  df_book_data_validation <- df_book_data_sorted[df_book_data_sorted$date %in% train_test_splits[[rolling_window]]$val, ]
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

  # New basis for X(t)
  new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
  logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
  logretfd = logretSmooth$fd
  beta_est <- fRegressList_pen$betaestlist[[2]]$fd

  delta_vol_pred_validation <- inprod(logretfd, beta_est)
  # B(t)
  beta_est_res <- fRegressList_pen_RESID$betaestlist[[2]]$fd
  res_pred_validation <- inprod(logretfd, beta_est_res)
  
  # Non-conformity measure
  vol_res_pred_valid <- delta_vol_pred_validation + df_past_realized_validation$pred - vol
  non_conformity_score_validation <- abs(vol_res_pred_valid)/exp(res_pred_validation)
  
  quan_value_of_ncscore_95 <- as.numeric(quantile(non_conformity_score_validation, 0.95))
  quan_value_of_ncscore_99 <- as.numeric(quantile(non_conformity_score_validation, 0.99))
  
  
  #####
  # Test
  #####
  # Test data
  vol <- vol_realized_sorted_imp[vol_realized_sorted_imp$date %in% train_test_splits[[rolling_window]]$test, ]
  vol <- as.numeric(vol$target)
  df_past_realized_test <- df_past_realized_sorted_imp[df_past_realized_sorted_imp$date %in% train_test_splits[[rolling_window]]$test, ]
  delta_vol_test <- vol - df_past_realized_test$pred
  
  df_book_data_test <- df_book_data_sorted[df_book_data_sorted$date %in% train_test_splits[[rolling_window]]$test, ]
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

  # New basis for X(t)
  new_logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
  logretSmooth = smooth.basis(abscissa, Xobs0, new_logret_basis)
  logretfd = logretSmooth$fd
  
  vol_pred_test <- inprod(logretfd, beta_est) + df_past_realized_test$pred
  res_pred_test <- exp(inprod(logretfd, beta_est_res))
  RMSE.test <- rmse(vol_pred_test, vol)
  
  # Conformal regions construction
  cp_lower_95 <-  vol_pred_test - quan_value_of_ncscore_95 * res_pred_test
  cp_lower_95[cp_lower_95<0]<-0
  cp_upper_95 <-  vol_pred_test + quan_value_of_ncscore_95 * res_pred_test
  
  cp_lower_99 <-  vol_pred_test - quan_value_of_ncscore_99 * res_pred_test
  cp_lower_99[cp_lower_99<0]<-0
  cp_upper_99 <-  vol_pred_test + quan_value_of_ncscore_99 * res_pred_test
  
  # Implied volatilities selection for the current stock
  CSCO <- data.frame(dates = IV_stocks$`CSCO US Equity`, implied_vol = IV_stocks$...2)
  iv <- CSCO$implied_vol
  iv_rep <- rep(iv, each = 12)
  df_past_realized_sorted_imp$impvol <- iv_rep
  implied_volatility_test <-  df_past_realized_sorted_imp[df_past_realized_sorted_imp$date %in% train_test_splits[[i]]$test, ]$impvol

  #####
  # CP representation
  #####
  pred_vol_test_from_linear = df_past_realized_test$pred
  df_test_predictions <- data.frame(implied_volatility_test, vol, vol_pred_test, pred_vol_test_from_linear, cp_lower_95, cp_upper_95, cp_lower_99, cp_upper_99)
  
  # Evaluation of annualized realized volatility
  from_bucket_to_annualized <- function(column_df) {
    column_df = ifelse(column_df < 0,
                       -sqrt((column_df^2) / (60 * 10) * 251 * 8 * 3600) * 100,  
                       sqrt((column_df^2) / (60 * 10) * 251 * 8 * 3600) * 100)   
    return(column_df)
  }
  
  # Data frame creation for the graph
  df_test_predictions$vol_pred_test <- mapply(from_bucket_to_annualized, df_test_predictions$vol_pred_test)
  df_test_predictions$cp_lower_99 <- mapply(from_bucket_to_annualized, df_test_predictions$cp_lower_99)
  df_test_predictions$cp_upper_99 <- mapply(from_bucket_to_annualized, df_test_predictions$cp_upper_99)
  df_test_predictions$cp_lower_95 <- mapply(from_bucket_to_annualized, df_test_predictions$cp_lower_95)
  df_test_predictions$cp_upper_95 <- mapply(from_bucket_to_annualized, df_test_predictions$cp_upper_95)
  df_test_predictions$pred_vol_test_from_linear <- mapply(from_bucket_to_annualized, df_test_predictions$pred_vol_test_from_linear)
  df_test_predictions$vol <- mapply(from_bucket_to_annualized, df_test_predictions$vol)
  
  df_test_predictions$inside_interval_95 <- df_test_predictions$vol >= df_test_predictions$cp_lower_95 & df_test_predictions$vol <= df_test_predictions$cp_upper_95
  df_test_predictions$inside_interval_99 <- df_test_predictions$vol >= df_test_predictions$cp_lower_99 & df_test_predictions$vol <= df_test_predictions$cp_upper_99
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  data_normale <- as.Date(train_test_splits[[rolling_window]]$test, format = "%Y-%m-%d")
  data_inglese <- format(data_normale, "%B %d, %Y")
  print(data_inglese)  
  
  Sys.setlocale("LC_TIME", "it_IT.UTF-8")
  time_sequence <- seq(from = as.POSIXct("10:30", format = "%H:%M"), 
                       to = as.POSIXct("16:00", format = "%H:%M"), 
                       by = "30 mins")
  
  df_test_predictions$time <- rep(time_sequence, length.out = nrow(df_test_predictions))
  common_y_limits <- range(
    c(df_test_predictions$vol, 
      df_test_predictions$cp_lower_99, 
      df_test_predictions$cp_upper_99)
  )
  # Here are needed two cases for a problem of representation
  if (all(df_test_predictions$inside_interval_99 == TRUE)) {
    cp_plot_99_enhanced <- ggplot(df_test_predictions, aes(x = time, y = vol)) +
      geom_point(aes(color = inside_interval_99), size = 3) +
      geom_errorbar(aes(ymin = cp_lower_99, ymax = cp_upper_99), width = 0.9, color = "#47A025") +
      scale_color_manual(values = c("#47A025", "#F71735"), labels = c("Inside Interval", "Outside Interval")) +
      ylim(common_y_limits) + 
      geom_hline(yintercept = df_test_predictions$implied_volatility_test, linetype = "dashed", color = "#5448C8", linewidth=1.2) +  
      labs(title = "Conformal Predictions with 
           a Significance Level of 99%",
           x = "Time", y = "Predicted Values",
           color = "Prediction Status") +
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
    
  } else {
    cp_plot_99_enhanced <- ggplot(df_test_predictions, aes(x = time, y = vol)) +
      geom_point(aes(color = inside_interval_99), size = 3) +
      geom_errorbar(aes(ymin = cp_lower_99, ymax = cp_upper_99), width = 0.9, color = "#47A025") +
      scale_color_manual(values = c("#F71735", "#47A025"), labels = c("Outside Interval", "Inside Interval")) +
      geom_hline(yintercept = df_test_predictions$implied_volatility_test, linetype = "dashed", color = "#5448C8", linewidth=1.2) + 
      ylim(common_y_limits) + 
      labs(title = "Conformal Predictions with 
           a Significance Level of 99%",
           x = "Time", y = "Predicted Values",
           color = "Prediction Status") +
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
  }
  
  if (all(df_test_predictions$inside_interval_95 == TRUE) ) {
    cp_plot_95_enhanced <- ggplot(df_test_predictions, aes(x = time, y = vol)) +
      geom_point(aes(color = inside_interval_95), size = 3) +
      geom_errorbar(aes(ymin = cp_lower_95, ymax = cp_upper_95), width = 0.9, color = "#47A025") +
      scale_color_manual(values = c("#47A025", "#F71735"), labels = c("Inside Interval", "Outside Interval")) +
      ylim(common_y_limits) +
      geom_hline(yintercept = df_test_predictions$implied_volatility_test, linetype = "dashed", color = "#5448C8", linewidth=1.2) +  
      labs(title = "Conformal Predictions with 
           a Significance Level of 95%",
           x = "Time", y = "Predicted Values",
           color = "Prediction Status") +
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
    
  } else {
    cp_plot_95_enhanced <- ggplot(df_test_predictions, aes(x = time, y = vol)) +
      geom_point(aes(color = inside_interval_95), size = 3) +
      geom_errorbar(aes(ymin = cp_lower_95, ymax = cp_upper_95), width = 0.9, color = "#47A025")+
      scale_color_manual(values = c("#F71735", "#47A025"), labels = c("Outside Interval", "Inside Interval")) +
      ylim(common_y_limits) + 
      geom_hline(yintercept = df_test_predictions$implied_volatility_test, linetype = "dashed", color = "#5448C8", linewidth=1.2) +  
      labs(title = "Conformal Predictions with 
           a Significance Level of 95%",
           x = "Time", y = "Predicted Values",
           color = "Prediction Status") +
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
  }
  combined_title <- textGrob(paste("Conformal Predictions of Realized Volatility on", data_inglese), 
                             gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Arrange the plots with the title
  final_plot <- grid.arrange(combined_title, cp_plot_99_enhanced, cp_plot_95_enhanced, 
                             ncol = 2, nrow = 2, heights = c(0.1, 1), 
                             layout_matrix = rbind(c(1, 1), c(2, 3)))
  
  # Combine the plots and ensure they share the same legend
  plot_file <- file.path(output_directory, paste0("combined_plot_", rolling_window, ".png"))
  ggsave(plot_file, final_plot, width = 12, height = 6)  
  results_list_95[[rolling_window]] <- check_volatility_95(df_test_predictions$vol, df_test_predictions$implied_volatility_test)
  results_list_99[[rolling_window]] <- check_volatility_99(df_test_predictions$vol, df_test_predictions$implied_volatility_test)
  
  }  

cat("Plots saved in the directory:", output_directory, "\n")

# Evaluation of rolling window forecasting metrics
count_matches(results_list_99)
count_matches(results_list_95)
