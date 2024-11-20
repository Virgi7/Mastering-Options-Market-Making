# Simple Linear Regression Model: RV(ti+10min) = α + β · RV(ti) + ϵ(ti+10min)
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
#####
# Linear model fitting
#####
model <- lm(target ~ pred, data = df_joined)
# Summary
summary(model)
# Residuals evaluation
res <- residuals(model)
# Fitted values
fitted_values <- fitted(model)
#####
# Diagnosis
#####
par(mfrow = c(1, 2))

# 1. Residuals vs Fitted
plot(fitted_values, res, main = "Residuals vs Fitted",
     xlab = "Fitted values", ylab = "Residuals",
     cex.main = 1.5,     # Increase title font size
     cex.lab = 1.4,      # Increase axis label font size
     cex.axis = 1.2)     # Increase axis tick font size
abline(h = 0, col = "red", lty = 2)

# 2. Q-Q Plot of Residuals
qqnorm(res, main = "Q-Q Plot of Residuals", cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.4)
qqline(res, col = "red")

# Diagnosis metrics
mse_linear <- mse(as.numeric(df_joined$target), fitted_values)
rmse_linear <- rmse(as.numeric(df_joined$target), fitted_values)
cat('RMSE:', rmse_linear, '\n')

