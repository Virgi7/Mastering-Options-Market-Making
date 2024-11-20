# Best number of basis for the basis expansion, to achieve the lowest RMSE
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
source("fRegress std error mia versione.R")
#####
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
#####
# Dataset preparation
#####
# Order book calculation
df_book_data <- sqlogret_per_time_id(file_paths)
# Target selection for the chosen stock
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, "target", drop=FALSE])
# Target of this model is realized volatility squared
vol <- as.numeric(vol$target)^2

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
volatility_squared <- vol
# Parameters
grade = 4
m = grade
degree = m-1
# B-spline basis for the function B(s)
nbasis_beta = 10
conbasisinter = create.constant.basis(rangeval=c(0,600))
beta_basis = create.bspline.basis(rangeval=c(0,600), nbasis=nbasis_beta, norder=m)
betalist = vector("list", 2)
betalist[[1]] = conbasisinter
betalist[[2]] = beta_basis

# Target selection for the chosen stock
vol <- as.data.frame(train[train$stock_id == stock_id_chosen, c("target", "time_id")])
# Without the square
vol_not_squared <- as.numeric(vol$target)

# Define the range of basis functions to test
basis_values <- seq(50, 200, by = 10)  
results <- data.frame(nbasis_covariates = integer(), rmse = double())

# Loop over each value of basis functions for covariates
for (nbasis_covariates in basis_values) {
  
  # Create the basis for covariates with current nbasis_covariates
  logret_basis <- create.bspline.basis(rangeval = c(0, 600), nbasis = nbasis_covariates, norder = m)
  logretSmooth <- smooth.basis(abscissa, Xobs0, logret_basis)
  logretfd <- logretSmooth$fd
  
  # Fit the scalar-on-functional model
  logret_list <- list(rep(1, n_obs), logretfd) 
  fRegressList <- fRegress(volatility_squared, logret_list, betalist)
  
  # Calculate the RMSE value for this model
  volatility_squaredhat <- as.numeric(fRegressList$yhatfdobj)
  volatility_squaredhat[volatility_squaredhat < 0] <- 0
  current_rmse <- rmse(vol_not_squared, sqrt(volatility_squaredhat))
  
  # Store the RMSE value
  results <- rbind(results, data.frame(nbasis_covariates = nbasis_covariates, rmse = current_rmse))
  
  cat("Completed nbasis_covariates =", nbasis_covariates, ": RMSE =", current_rmse, "\n")
}

# Select the number of basis functions with the lowest RMSE
optimal_nbasis <- results$nbasis_covariates[which.min(results$rmse)]
cat("Optimal number of basis functions for covariates based on RMSE:", optimal_nbasis, "\n")

# Plot
ggplot(results, aes(x = nbasis_covariates, y = rmse)) +
  geom_line(color = "#3D348B", linewidth = 1.0) +
  geom_point(color = "#3D348B", size = 3) +  
  geom_vline(xintercept = optimal_nbasis, linetype = "dashed", color = "#DF2935", linewidth = 1.2) +
  annotate("text", x = optimal_nbasis, y = max(results$rmse), 
           label = paste("Optimal:", optimal_nbasis), color = "#DF2935", 
           vjust = -1, hjust = 0.5, fontface = "bold", size = 4) +
  labs(
    title = "RMSE for Different Numbers of Basis Functions for X(t)",
    x = "Kz",
    y = "RMSE"
  ) +
  scale_x_continuous(breaks = seq(50, 200, by = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    legend.position = "top",  
    #legend.title = element_text(size= 20, face = "bold"),
    legend.text = element_text(size= 18)
  )



