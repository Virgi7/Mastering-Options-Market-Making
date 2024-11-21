# Freedman and Lane permutation test
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

#####
# Dataset preparation
#####
df_book_data <- sqlogret_per_time_id(file_paths)

column_name <- "log_return"  
df_dynamic <- df_book_data[c("time_id", "seconds_in_bucket", column_name)]
colnames(df_dynamic)[ncol(df_dynamic)] <- "value"
df_dynamic <- as.data.table(df_dynamic)
time_istants <- unique(df_book_data$time_id)
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

# Fase 1: UNDER H0: B(s)'= 0, build a functional model assuming B(s) constant

volatility_squared <- vol

# X(t)
logret_list = vector("list", 2)
logret_list[[1]] = rep(1, n_obs)
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
beta_basis_H0 = create.constant.basis(rangeval=c(0,600))
betalist = vector("list", 2)
betalist[[1]] = conbasisinter
betalist[[2]] = beta_basis_H0

# Scalar-on-function Model
fRegressList_H0 = fRegress(volatility_squared, logret_list, betalist)
betaestlist_H0 = fRegressList_H0$betaestlist

volatility_squaredhat_H0 = fRegressList_H0$yhatfdobj
# Residuals evaluation
residuals_H0 = volatility_squared - volatility_squaredhat_H0
RMSE.constantbeta = rmse(sqrt(volatility_squaredhat_H0), vol_not_squared)

# Fase 2: Build the classical functional model with B(s) not constant

# X(t)
logret_list = vector("list", 2)
logret_list[[1]] = rep(1, n_obs)
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
nbasis_beta = 10
beta_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_beta, norder=m)
betalist = vector("list", 2)
betalist[[1]] = conbasisinter
betalist[[2]] = beta_basis

# Scalar-on-function Model
fRegressList_base = fRegress(volatility_squared, logret_list, betalist)
betaestlist_base = fRegressList_base$betaestlist
logret_betafd_base = betaestlist_base[[2]]$fd

# Computation of B(s) function
beta_t_base <- eval.fd(abscissa, logret_betafd_base)

# First derivative and L2 norm
dbeta_dt_base <- diff(beta_t_base) / diff(abscissa)
T0 <- sqrt(sum(dbeta_dt_base^2))

# Fase 3: Residuals permutation and construction of the model with these

B <- 1000 # Number of permutations
T_stat <- numeric(B) # Vector where it will be stored the values of T*

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = B,
                       complete = "=",   
                       incomplete = "-", 
                       current = ">",    
                       clear = FALSE,    
                       width = 100)      

for(perm in 1:B){
  pb$tick()
  
  # Permutation
  residuals_sampled = sample(residuals_H0)
  # y_fict = y_fitted + residuals_sampled
  volatility_squared_fict = fRegressList_H0$yhatfdobj + residuals_sampled
  
  # Fase 4: Build the model with a fictitious y 
  
  # X(t)
  logret_list = vector("list", 2)
  logret_list[[1]] = rep(1, n_obs)
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
  nbasis_beta = 10
  beta_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_beta, norder=m)
  betalist = vector("list", 2)
  betalist[[1]] = conbasisinter
  betalist[[2]] = beta_basis
  
  # Scalar-on-function Model
  fRegressList_perm = fRegress(as.numeric(volatility_squared_fict), logret_list, betalist)
  betaestlist_perm = fRegressList_perm$betaestlist
  logret_betafd_perm = betaestlist_perm[[2]]$fd
  plot(logret_betafd_perm, xlab = "Seconds", ylab = "Beta for log_returns")
  
  beta_t <- eval.fd(abscissa, logret_betafd_perm)
  
  # Fase 5: Compute the first derivative and L2 norm
  dbeta_dt <- diff(beta_t) / diff(abscissa)
  # Test statistic:
  T_stat[perm] <- sqrt(sum(dbeta_dt^2))
}

# p-value of the test
p_value = sum(T_stat>T0)/B

cat(sprintf("P-value of the Freedman & Lane test is: %.3f%%\n", p_value*100))

