# Scalar-on-functon Model: RV(ti+10min) - RV(ti) = α + int(β(s) · X(s)^2) + ϵ(t)
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

# Log-returns calculation (with an initial NA) for each second in each time id
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
  
  return(df_realized_vol_per_stock[, .(row_id, get(prediction_column_name))])
}

#  RV(ti) calculation for each bucket in different files
past_realized_volatility_per_stock <- function(list_file, prediction_column_name) {
  df_past_realized <- data.table()  
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

# Target selection for the chosen stock
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, "target", drop=FALSE])
vol <- as.numeric(vol$target)
# Definition of the target as the difference between volatilities
delta_vol <- vol - df_past_realized_train$V2

df_book_data <- sqlogret_per_time_id(file_paths)

# Time series
# Column selection
column_name <- "log_return" 
# Dataframe creation
df_dynamic <- df_book_data[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
# Conversion of df_dynamic in data.table for a faster process
df_dynamic <- as.data.table(df_dynamic)
# Matrix creation (600 rows and n_obs columns)
time_istants <- unique(df_book_data$time_id)
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
#####
# Scalar-on-function fitting
#####
logret_list = vector("list", 2)
logret_list[[1]] = rep(1, n_obs)
# Parameters 
nbasis_beta = 10
grade = 5
m = grade
nbasis_covariates = 120
logret_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_covariates, norder=m)
logretSmooth = smooth.basis(abscissa, Xobs0, logret_basis)
logretfd = logretSmooth$fd
logret_list[[2]] = logretfd

conbasisinter = create.constant.basis(rangeval=c(0,600))
beta_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_beta, norder=m)
betalist = vector("list", 2)
betalist[[1]] = conbasisinter
betalist[[2]] = beta_basis

# Scalar-on-functional model
fRegressList = fRegress(delta_vol, logret_list, betalist)
# Coefficients extraction
betaestlist_nonpen = fRegressList$betaestlist
logret_betafd_nonpen = betaestlist_nonpen[[2]]$fd

# Prediction extraction
deltavol_hat_nonpen = fRegressList$yhatfdobj
# Residuals evaluation
deltavol_res_nonpen = delta_vol - deltavol_hat_nonpen
#####
# Diagnosis
#####
# R2
SSE1.1 = sum(deltavol_res_nonpen^2)
SSE0 = sum((delta_vol - mean(delta_vol))^2)
RSQ1 = (SSE0 - SSE1.1)/SSE0
# F-Ratio
Fratio1 = ((SSE0 - SSE1.1)/5)/(SSE1.1/29)
#RMSE
RMSE.1 = rmse(delta_vol, deltavol_hat_nonpen)

# Confidence interval at 95% for B(s) function
sigmaE. = sum(deltavol_res_nonpen^2)/(n_obs - fRegressList$df)
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

#####
# Scalar-on-function Model with penalization term
#####
lambda0 = 2^2
betafdPar = fdPar(beta_basis, 2, lambda0)
# Introduction of the penalization term in the functional object
betalist[[2]] = betafdPar

# Cross validation for the selection of the right lambda, without the command
# as the coefficients matrix is near the singularity
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
  fReg_i = fRegress(delta_vol, logret_list, betalist_i)
  betaestlist_i = fReg_i$betaestlist
  logret_betafd_CV[ilam] = betaestlist_i[[2]]$fd
  delta_volhat_CV = fReg_i$yhatfdobj
  delta_volres_CV = delta_vol - delta_volhat_CV
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
  RMSE.CV[ilam] = rmse(delta_vol, delta_volhat_CV)
  GCVscore[ilam] = fReg_i$GCV
}
# R2 vector
RSQCV = (SSE0 - SSE1.CV)/SSE0

# Minimum GCV score is for lambda = 2^-15
minGCVscore <- min(GCVscore)
min_index <- which.min(GCVscore) 
exp_lambda_opt <- loglam[min_index]
lambda_opt_GCV <- 2^exp_lambda_opt

lambda_opt = lambda_opt_GCV

# Fitting of the penalized model
betalist_pen = betalist
betafdPar2 = betalist_pen[[2]]
betafdPar2$lambda = lambda_opt
betalist_pen[[2]] = betafdPar2
fRegressList_pen = fRegress(delta_vol, logret_list, betalist_pen)
betaestlist_pen = fRegressList_pen$betaestlist
logret_betafd_pen = betaestlist_pen[[2]]$fd
delta_volhat_pen = fRegressList_pen$yhatfdobj
delta_volres_pen = delta_vol - delta_volhat_pen
sigmaE. = sum(delta_volres_pen^2)/(n_obs - fRegressList_pen$df)
sigmaE = sigmaE. * diag(rep(1, n_obs))
y2cMap = logretSmooth$y2cMap
stderrList = fRegress.stderr(fRegressList_pen, y2cMap, sigmaE)
betafdPar = betaestlist_pen[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
RMSE.2 = rmse(delta_vol,delta_volhat_pen)
SSE1.2 = sum(deltavol_res_nonpen^2)
RSQ2 = (SSE0 - SSE1.2)/SSE0
# Final plot
plot(betafd, xlab = "Seconds", ylab = "Log-return2 Reg. Coeff.", lwd = 2)
lines(betafd + betastderrfd, lty=2, lwd=1)
lines(betafd - betastderrfd, lty=2, lwd=1)
#####
# Plot of B(s)
#####
time_points <- seq(1, 600, length.out = 600)  
# Evaluate `betafd` and `betastderrfd` at these time points
beta_values <- eval.fd(time_points, betafd)
upper_values <- eval.fd(time_points, betafd + betastderrfd)
lower_values <- eval.fd(time_points, betafd - betastderrfd)
# Create a data frame for ggplot
df <- data.frame(
  Time = time_points,
  Beta = beta_values,
  Upper = upper_values,
  Lower = lower_values
)
ggplot(df, aes(x = Time, y = Beta)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#89043D", alpha = 0.3) +
  geom_line(color = "black", size = 0.9) +  
  geom_line(aes(y = Upper), color = "#89043D", linetype = "dashed", size = 0.6) +  
  geom_line(aes(y = Lower), color = "#89043D", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "solid", size = 0.5) +
  labs(
    title = expression(bold(beta~"Coefficients with Confidence Interval")),
    x = "Seconds",
    y = expression(bold(beta(s))) 
  ) +
  scale_x_continuous(breaks = seq(0, 600, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    legend.position = "top",  
    legend.text = element_text(size= 18)
  )