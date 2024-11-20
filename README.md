# Mastering Options Market Making: Forecasting Short-Term Realized Volatility Using Non-Parametric Models
This repository contains the code for a master thesis project aimed at developing an **Hedging Strategy based on a Realized Volatility Forecasting Model**.

**Volatility** plays a critical role in financial markets by quantifying the magnitude of price changes in a security, directly influencing trading strategies and risk management. Traders, especially market makers, rely on two primary metrics to assess volatility: **realized volatility**, which measures historical price movements and **implied volatility**, which reflects market expectations. The comparison between these two quantities is particularly important for market makers, as it helps guide the development of effective hedging strategies. In fact, by generating buy and sell signals based on this comparison, market makers understand the best position to take and can better manage risk and exploit profit opportunities. As such, it is essential to develop a reliable model for predicting realized volatility, while implied volatility is derived from market option prices.

Short-term realized volatility is of particular interest, as leading global market makers, such as Optiver, may want to update their positions continuously in response to market dynamics.

This study aims to develop an effective model for predicting realized volatility over very short intervals, leveraging high-frequency Limit Order Book data provided by Optiver. The goal is to forecast the realized volatility for the second 10-minutes sub-interval in a 20-minutes time window, using data from the first 10-minutes segment. The forecasting procedure is structured in two steps: first a point forecast is generated, which is then completed by a probabilistic forecast. The models applied start with _Simple Linear Regression_ and are then extended to _Non-Parametric Statistics_, focusing on functional models, in particular the _Scalar-on-function_ model for point forecasting and the _Conformal Prediction_ method for probabilistic forecasting.

The complete implementation of these models can be found in the project's respective _R scripts_.

### Repository Index

1. [Introduction](#introduction)
2. [Problem Setting](#problem-setting)
3. [Scalar-on-function Model](#scalar-on-function-model)
4. [Point-wise Forecasting and Conformal Prediction](#pointwise-forecasting-and-conformal-prediction)
5. [Hedging Strategy based on Conformal Prediction](#hedging-startegy-based-on-conformal-prediction)
6. [Conclusions and Future Developments](#conclusions-and-future-developments)
7. [Repository Files](#repository-files)
8. [How to Use](#how-to-use)

# Problem Setting

## Data overview
The dataset, sourced from Kaggle, consists of 112 anonymized stock, identified by a number from 1 to 126. For each stock high-frequency trade and order book data are available. It is important to note that all models and results discussed in the following chapters are based on analyses conducted on Stock 61, which was randomly selected from the dataset. However, it should also be highlighted that the observed results are consistent across other stocks in the dataset, and so, for they robustness, can be generalized.

For each stock, the trade and the order book dataset include detections related to various time intervals of 20 minutes during the trading day, defined as time buckets. These buckets, while not necessary ordered, maintain a consistent structure across all stocks and are defined as intervals $T_i = (t_i; t_{i+10_{min}})$, where $t_i$ represents the first 10 minutes and $t_{i+10_{min}}$ identifies the second 10 minutes, for all the 3830 intervals in the dataset, identified with $i$. For each bucket $T_i$, trades and order book data are provided  with a granularity at the level of the seconds in the first 10 minutes.

For each time bucket described in both the order book and trade data, the realized volatility over the last 10 minutes in the bucket is also provided; it is the target for the model (i.e. what has to be predicted), which will detailed further in the next section.

## Data pre-processing
It is necessary to highlight two important points; the first one is about missing data in the order book.
In fact, there are many missing data related to entirely time buckets, which means that there were no updates in bid and ask prices in that period.
And even in specific time slots, the seconds measured are different and in different quantities, but this does not represent a scarcity problem because when there is less variation in prices and quantities, the value of volatility changes less quickly and it is better to give more weight, and therefore more relevance, to the high variation buckets.

The second point concerns the decision to leave out trade data in this particular case, instead preferring order book data, which are constantly available, even for illiquid products. However, trade data can be incorporated into the model in a further analysis, should it be necessary to capture further aspects of market activity.

## Reverse Engineering and Embedding
It is important to note that, using techniques developed by participants in the Optiver challenge, which used the dataset in question, it has been possible to approximate the true time-order of the $T_i$'s and to establish an approximate mapping between the anonymized stocks and the actual market assets. These insights are used in the final phase of this paper, where it is simulated the application of the proposed models on the data.

## Simple Linear Regression
Starting from the current literature, the widely used model for forecasting realized volatility is the HAR model, introduced by \cite{corsi2012har}, who demonstrates that a simple, linear model, based on historical realized volatilities, yields strong forecasting performance. However, the HAR model is not suitable in this context due to the loss of temporal order within the 20-minutes intervals.

Building upon the forecasting capability of the HAR model, this study first implements a Linear Regression model for forecasting $RV_{t_{i+10_{min}}}$, using as a covariate the standard deviation of the stock log-returns over the first half of the 20-minutes interval; in fact, this measure is identified as the most effective approximation of realized volatility, by \cite{zhang2005tale}. It is identified by $RV_{t_i}$.

The baseline Linear Regression model used for comparison is expressed as:

$$
RV_{t_i+10_{min}} = \alpha + \beta \cdot RV_{t_i} + \epsilon_{t_i+10_{min}}, \quad \epsilon_{t_i+10_{min}} \sim N(0, \sigma^2).
$$

Where:
- $RV_{t_i}$: Realized volatility at time $t_i$.
- $\alpha$: Intercept term.
- $\beta$: Coefficient for $RV_{t_i}$.
- $\epsilon_{t_i+10_{min}}$: Normally distributed error term with mean $0$ and variance $\sigma^2$.

This first model seems to have a good prediction capability, as assessed by the Root Mean Squared Error (RMSE). The primary limitation of this Linear Regression model lies in the violation of two key assumptions: homoscedasticity and Gaussianity of the residuals. While these issues will not be solved in next models fitting, they will be partially addressed through the later application of the assumptions-free Conformal Prediction method for probabilistic forecasting.

# Scalar-on-function Model
To improve prediction accuracy, it is introduced a **linear functional model** that incorporates second-by-second **quadratic log-return series** from each $t_i$'s order book snapshot as a functional covariate.

The **central idea** is to expand each covariate, which has multiple values within each 10-minutes interval, by adding a new dimension representing seconds within the interval. This is achieved using carefully chosen **basis functions**.

- The coefficient of this covariate is modeled as a function of time, denoted as $\beta(s)$.
- The target variable remains scalar and corresponds to $RV_{t_{i+10_{min}}}^2$ to maintain comparability with the previous Linear Regression model.

The fitted model is called **Scalar-on-Function**, and its structure is given by:

$$
RV_{t_i+10_{min}}^2 = \alpha + \int_0^T X_{t_i}^2(s) \beta(s) \, ds + \epsilon_{t_{i+10_{min}}}
$$

Where:
- $X_{t_i}^2(s)$: Quadratic log-returns over seconds in the interval.
- $\beta(s)$: Coefficient function expanded using basis functions.
- $\epsilon_{t_{i+10_{min}}}$: Error term.

To identify the most influential seconds for predicting realized volatility, a regularization method is applied, as the initial fitting of the $\beta(s)$ 
function, due to its sinusoidal behaviour, does not provide this information directly. Specifically, a penalty term is added to impose smoothness on the estimated $\beta(s)$ by selecting, through cross-validation, an optimal value for the regularization parameter $\lambda$, which controls the roughness of the function. 

The RMSE of this Scalar-on-function model results slightly worst than the linear model’s one, so it is performed a permutation test, Freedman and Lane one, from which results to prefer the impact of the entire time series of squared log-returns to the covariate of the linear model $RV_{t_i}$, which synthesises these information.

So, the final proposed model combines the baseline Linear Regression model and the functional covariate improvement, resulting in a **Scalar-on-Function model** for $\Delta RV_t$. The response variable $\Delta RV_t$ is defined as:

$$
\Delta RV_t = RV_{t_{i+10\text{min}}} - RV_{t_i}.
$$

The model is expressed as:

$$
\Delta RV_t = \alpha + \int_0^T X_{t_i}^2(s) \beta(s) \, ds + \epsilon_t.
$$

Where:
- $\Delta RV_t$: Change in realized volatility over the 10-minute interval.
- $\alpha$: Intercept term.
- $X_{t_i}^2(s)$: Functional covariate of quadratic log-returns.
- $\beta(s)$: Time-dependent coefficient function.
- $\epsilon_t$: Error term.

This model results as the best performing one in terms of accuracy. Here the RMSE values, evaluated on in-sample data, expressed in bps.

| Model                                   | RMSE (in-sample) |
|-----------------------------------------|-----------------:|
| Simple Linear Regression                | 8.691            | 
| Scalar-on function $RV_{t_{i+10min}}^2$ | 8.901            |
| Scalar-on function $\Delta RV_t$        | 8.431            |


# Point-wise Forecasting and Conformal Prediction

In the second phase of this work, it is evaluated the out-of-sample point-wise forecasts produced by the three models introduced in the previous sections. Given that the models exhibit heteroscedasticity and non-Gaussianity in residuals, the predictions are further refined using the Conformal Prediction method. This method, as introduced by \cite{fontana2023conformal}, is a distribution-free and non-parametric approach requiring minimal assumptions, capable of producing statistically valid prediction intervals even in finite-sample scenarios.

## Conformal Prediction with Inductive Conformal Prediction (ICP)

This repository implements **Inductive Conformal Prediction (ICP)**, following the methodology presented in \cite{kato}. The goal of the project is to evaluate predictive accuracy and precision of different models for forecasting realized volatility using **Linear Regression** and **Scalar-on-function** approaches.

The workflow includes:
- Model fitting with training data
- Calculation of non-conformity scores
- Generation of conformal prediction intervals
- Performance evaluation using RMSE, coverage, and width metrics

The ICP method works as follows:
1. **Fit the selected model** on the training data.
2. **Compute residuals**: Using a chosen non-conformity measure.
3. **Calculate non-conformity scores**: The quantiles of these scores are computed using the validation set.
4. **Generate conformal prediction intervals**: The quantiles are then used on the test set to build prediction intervals.

The dataset is split into three sets: training, validation, and test, following the **80/20 rule**, where more weight is given to the validation set. This allows the validation set to better calibrate the conformal prediction intervals.

|     Training      |      Validation      |     Test       |
|------------------:|---------------------:|---------------:|
|       56%         |         24%          |         20%    |

The conformal regions obtained are of the type:

$$
\Gamma_{1 - \alpha}(x_i) = (\hat{y_i} - q_{1 - \alpha} \cdot \sigma_i, \; \hat{y_i} + q_{1 - \alpha} \cdot \sigma_i )
$$

where $q_{1-\alpha}$ is the quantile corresponding to the desired coverage level $1-\alpha$, of the non-conformity scores distribution R, defined as 

$R_i = \left|\frac{y_i - \hat{y}_i}{\sigma_i} \right|$ and $\sigma_i = e^{\mu_i}$, where $\mu_i$ is the predicted value of $\ln(|y_i - \hat{y}_i|)$.

With the application of this method, the first thing done is the evaluation of the out-of-sample forecasting, for each of the three method described above, confirming that the linear model and the Scalar-on-function $\Delta RV_t$ are the best performing ones.

| Model                                   | RMSE (out.of-sample) |
|-----------------------------------------|---------------------:|
| Simple Linear Regression                | 8.879                | 
| Scalar-on function $RV_{t_{i+10min}}^2$ | 9.604                |
| Scalar-on function $\Delta RV_t$        | 8.771                |

Additionally, also the goodness of fit of these regions can be then evaluated by computing some metrics, specifically related to prediction regions, like their medium width and total coverage and also to probabilistic forecasting, such as Pinball loss and Winkler score.

The conformal regions computed and compared are at 90%, 95% and 99%. Coverage is not very informative as the obtained values are similar across all the three models

VEDI SE INSERIRE I GRAFICI

# Hedging Strategy based on Conformal Prediction
The final stage of this study exploits the goodness of fit and the predictive power of the Scalar-on-function for $\Delta RV_{t}$, combined with the probabilistic power of the Conformal Prediction method. The model is applied directly to data for matched stocks: Stock 59 with AMGN, Stock 4 with CHTR, and Stock 41 with CSCO, using IV data from January 1 to December 31, 2020.

In this stage, RV forecasting is applied using a rolling window to predict intra-day RV values for each day in the last two months of the year. The hedging strategy leverages RV and IV differences; the strategy performed is to buy when RV is higher than IV and to sell when RV is lower than IV

The effectiveness of the forecasted buy and sell signals is summarized in the table below, using three empirical metrics:
- **Wrong Signal (WS)**: Counts instances where the model's signal opposes the actual RV.
- **Lost Opportunity (LO)**: Occurs when the conformal regions include implied volatility, suggesting no update, while the actual RV is higher or lower than the implied volatility.
- **Successful Trade (ST)**: Identifies the right buy/sell signals that lead to trades.

These metrics are evaluated at the 99% and 95% significance levels, considering a total of $M = 492$ signals, as there are 12 RV values per day for two months (only trading days are counted).

| **Stock**     | **CP at 99%**   | **CP at 95%**   |
|---------------|-----------------|-----------------|
| **AMGN-59**   | **WS:** 0.203%<br>**LO:** 53.252%<br>**ST:** 9.553% | **WS:** 0.203%<br>**LO:** 45.325%<br>**ST:** 17.479% |
| **CHTR-4**    | **WS:** 0%<br>**LO:** 62.398%<br>**ST:** 13.617% | **WS:** 0%<br>**LO:** 49.390%<br>**ST:** 26.626% |
| **CSCO-41**   | **WS:** 0%<br>**LO:** 71.341%<br>**ST:** 6.911% | **WS:** 0.203%<br>**LO:** 42.886%<br>**ST:** 35.163% |

*Table: Wrong Signal (WS), Lost Opportunity (LO), and Successful Trade (ST) percentages in conformal regions at different significance levels for stocks AMGN-59, CHTR-4, and CSCO-41.*


# Conclusions and Future Developments


# Repository Files
- _CrashBusters.ipynb_: This notebook contains the final implementation of our EWS product CrashBusters.
- _CrashBustersApp.ipynb_: This is the notebook of our CrashBusters Application.
- _DecisionTree\_RandomForest.ipynb_: This notebook details the use of a Decision Tree classifier to identify market crashes.
- _Logistic Regression.ipynb_: This notebook contains the implementation and evaluation of a Logistic Regression model for market anomaly detection.
- _NaiveBayes.ipynb_: This notebook contains the code referred to the Naive Bayes method for market crashes detection.
- _SVM.ipynb_: This notebook demonstrates the use of a Support Vector Machine (SVM) for classifying market crashes.
- _kNN.ipynb_:  This notebook contains the analysis of kNN classifier as market anomaly detection model.

# How to Use
- Clone the repository to your local machine.
- Navigate to the directory and open the Jupyter notebooks in your preferred environment.
- Run the notebooks in the sequence provided to replicate the experiments and view the results.
- The final model selection notebook (_CrashBusters.ipynb) will guide you through the process of selecting the most effective model based on our experiments.
