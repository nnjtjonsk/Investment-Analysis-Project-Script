# Investment Analysis Assignment 

# --- 0. Setup: Load Libraries ---
# install.packages(c("readxl", "dplyr", "ggplot2", "tidyr", "lubridate"))
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(broom)


excel_file_path <- "InvestmentMastersheet.xlsx"

company_sheet_name <- "company-data"
market_sheet_name <- "RF and Mkt"
book_value_sheet_name <- "Book Value"
stocks_for_q2_q3 <- c(10107, 10104, 55976, 12490, 11850, 14593) # MSFT, ORCL, WMT, IBM, XOM, AAPL

# --- Data Loading and Cleaning Function ---
load_and_prepare_data <- function(excel_path, company_sheet, market_sheet) {
  company_df <- read_excel(excel_path, sheet = company_sheet)
  market_df <- read_excel(excel_path, sheet = market_sheet)
  company_df$date <- as.Date(company_df$date)
  colnames(market_df)[1] <- "date"
  market_df$date <- as.Date(market_df$date)
  company_df_cleaned <- company_df %>%
    filter(!is.na(PRC) & !is.na(SHROUT) & !is.na(RET)) %>%
    mutate(PRC = abs(PRC), RET = as.numeric(as.character(RET))) %>%
    filter(!is.na(RET))
  colnames(market_df)[2] <- "RF"
  colnames(market_df)[3] <- "sprtrn"
  market_df_cleaned <- market_df %>%
    mutate(RF = RF / 100) %>%
    rename(Mkt_Return = sprtrn) %>%
    select(date, Mkt_Return, RF)
  full_data <- inner_join(company_df_cleaned, market_df_cleaned, by = "date")
  return(full_data)
}

# --- 1. Question 1: Market-Value Weighted Portfolio Analysis ---
cat("--- Starting Question 1: Market-Value Weighted Portfolio Analysis ---\n")
all_data <- load_and_prepare_data(excel_file_path, company_sheet_name, market_sheet_name)
portfolio_returns <- all_data %>% arrange(date) %>% group_by(date) %>% mutate(MarketValue = PRC * SHROUT, TotalMarketValue = sum(MarketValue, na.rm = TRUE)) %>% mutate(Weight = MarketValue / TotalMarketValue, WeightedReturn = Weight * RET) %>% summarise(Portfolio_Return = sum(WeightedReturn, na.rm = TRUE)) %>% ungroup()
avg_monthly_return <- mean(portfolio_returns$Portfolio_Return, na.rm = TRUE)
monthly_risk_sd <- sd(portfolio_returns$Portfolio_Return, na.rm = TRUE)
investment_eur <- 200000; usd_to_eur_rate <- 0.92; annual_return <- 12 * avg_monthly_return; annual_risk_sd <- sqrt(12) * monthly_risk_sd; z_score <- qnorm(0.01); var_eur <- investment_eur * (annual_return + z_score * annual_risk_sd)

worst_10_percent_threshold <- quantile(portfolio_returns$Portfolio_Return, probs = 0.10, na.rm = TRUE)
worst_returns <- portfolio_returns$Portfolio_Return[portfolio_returns$Portfolio_Return <= worst_10_percent_threshold]
guise_10_percent <- mean(worst_returns, na.rm = TRUE)
sink("Q1_results.txt"); cat("Question 1 Results...\n"); print(list(avg_monthly_return=avg_monthly_return, monthly_risk_sd=monthly_risk_sd, var_eur=var_eur, guise_10_percent=guise_10_percent)); sink()
cat("Question 1 results saved to Q1_results.txt\n\n")


# --- 2. Question 2: Mean-Variance Frontier (Full Sample) ---
cat("--- Starting Question 2: Mean-Variance Frontier (Full Sample) ---\n")
q2_data_filtered <- all_data %>% filter(PERMNO %in% stocks_for_q2_q3)
recent_tickers_q2 <- q2_data_filtered %>% group_by(PERMNO) %>% filter(date == max(date)) %>% summarise(stable_ticker = first(TICKER))
q2_data_consolidated <- q2_data_filtered %>% left_join(recent_tickers_q2, by = "PERMNO") %>% select(date, stable_ticker, RET)
q2_data_wide <- q2_data_consolidated %>% pivot_wider(names_from = stable_ticker, values_from = RET)
returns_matrix_q2 <- as.matrix(q2_data_wide[, -1]); mean_returns_q2 <- colMeans(returns_matrix_q2, na.rm = TRUE); cov_matrix_q2 <- cov(returns_matrix_q2, use = "pairwise.complete.obs"); cov_matrix_q2[is.na(cov_matrix_q2)] <- 0
calculate_portfolio_stats <- function(weights, mean_returns, cov_matrix) { port_return <- sum(weights * mean_returns); port_sd <- sqrt(t(weights) %*% cov_matrix %*% weights); return(c(ret = port_return, sd = port_sd)) }
target_returns_q2 <- seq(min(mean_returns_q2), max(mean_returns_q2), length.out = 100)
frontier_no_shorts_list_q2 <- lapply(target_returns_q2, function(target_ret) {
  objective_penalized <- function(weights) { penalty <- 1000; sum_constraint <- (sum(weights) - 1)^2; return_constraint <- (sum(weights * mean_returns_q2) - target_ret)^2; return(t(weights) %*% cov_matrix_q2 %*% weights + penalty * (sum_constraint + return_constraint)) }
  result <- optim(par = rep(1/length(mean_returns_q2), length(mean_returns_q2)), fn = objective_penalized, method = "L-BFGS-B", lower = 0, upper = 1)
  weights <- result$par / sum(result$par); stats <- calculate_portfolio_stats(weights, mean_returns_q2, cov_matrix_q2); return(data.frame(ret = stats['ret'], sd = stats['sd']))
})
frontier_no_shorts_q2 <- do.call(rbind, frontier_no_shorts_list_q2)
ones_q2 <- rep(1, length(mean_returns_q2)); inv_cov_q2 <- solve(cov_matrix_q2); a_q2 <- as.numeric(t(ones_q2) %*% inv_cov_q2 %*% ones_q2); b_q2 <- as.numeric(t(ones_q2) %*% inv_cov_q2 %*% mean_returns_q2); c_q2 <- as.numeric(t(mean_returns_q2) %*% inv_cov_q2 %*% mean_returns_q2)
frontier_sd_shorts_q2 <- sqrt((a_q2 * target_returns_q2^2 - 2*b_q2*target_returns_q2 + c_q2) / (a_q2*c_q2 - b_q2^2)); frontier_shorts_q2 <- data.frame(ret = target_returns_q2, sd = frontier_sd_shorts_q2)
individual_assets_q2 <- data.frame(TICKER = names(mean_returns_q2), ret = mean_returns_q2, sd = sqrt(diag(cov_matrix_q2)))
vw_portfolio_6_q2 <- all_data %>% filter(PERMNO %in% stocks_for_q2_q3) %>% group_by(date) %>% summarise(Portfolio_Return = weighted.mean(RET, PRC * SHROUT, na.rm = TRUE))
vw_port_6_stats_q2 <- data.frame(TICKER = "VW Portfolio", ret = mean(vw_portfolio_6_q2$Portfolio_Return, na.rm = TRUE), sd = sd(vw_portfolio_6_q2$Portfolio_Return, na.rm = TRUE))
sink("Q2_results.txt"); cat("Question 2 Results: Input Data (Full Sample)\n\n"); cat("Individual Asset Performance:\n"); print(individual_assets_q2); cat("\nValue-Weighted Portfolio of the 6 Assets:\n"); print(vw_port_6_stats_q2); sink()
mv_plot_q2 <- ggplot() + geom_path(data = frontier_shorts_q2, aes(x = sd, y = ret, color = "Short Sales Allowed"), size = 1) + geom_path(data = frontier_no_shorts_q2, aes(x = sd, y = ret, color = "No Short Sales"), size = 1) + geom_point(data = individual_assets_q2, aes(x = sd, y = ret), size = 3) + geom_text(data = individual_assets_q2, aes(x = sd, y = ret, label = TICKER), vjust = -1) + geom_point(data = vw_port_6_stats_q2, aes(x = sd, y = ret), size = 4, color = "darkgreen", shape = 18) + geom_text(data = vw_port_6_stats_q2, aes(x = sd, y = ret, label = TICKER), vjust = -1.5) + labs(title = "Mean-Variance Frontier (Full Sample)", x = "Risk (StDev)", y = "Expected Return", color = "Frontier Type") + theme_minimal() + theme(legend.position="top") + scale_y_continuous(labels = scales::percent) + scale_x_continuous(labels = scales::percent)
ggsave("mean_variance_frontier_full_sample.png", plot = mv_plot_q2, width = 10, height = 7)
cat("Question 2 plot and results saved.\n\n")


# --- 3. Question 3: Mean-Variance Frontier (Last 10 Years) ---
cat("--- Starting Question 3: Mean-Variance Frontier (Last 10 Years) ---\n")
most_recent_date <- max(all_data$date); ten_years_ago <- most_recent_date - years(10)
q3_data_filtered <- all_data %>% filter(date >= ten_years_ago)
q3_data_for_pivot <- q3_data_filtered %>% filter(PERMNO %in% stocks_for_q2_q3)
recent_tickers_q3 <- q3_data_for_pivot %>% group_by(PERMNO) %>% filter(date == max(date)) %>% summarise(stable_ticker = first(TICKER))
q3_data_consolidated <- q3_data_for_pivot %>% left_join(recent_tickers_q3, by = "PERMNO") %>% select(date, stable_ticker, RET)
q3_data_wide <- q3_data_consolidated %>% pivot_wider(names_from = stable_ticker, values_from = RET)
returns_matrix_q3 <- as.matrix(q3_data_wide[, -1]); mean_returns_q3 <- colMeans(returns_matrix_q3, na.rm = TRUE); cov_matrix_q3 <- cov(returns_matrix_q3, use = "pairwise.complete.obs"); cov_matrix_q3[is.na(cov_matrix_q3)] <- 0
target_returns_q3 <- seq(min(mean_returns_q3), max(mean_returns_q3), length.out = 100)
frontier_no_shorts_list_q3 <- lapply(target_returns_q3, function(target_ret) {
  objective_penalized <- function(weights) { penalty <- 1000; sum_constraint <- (sum(weights) - 1)^2; return_constraint <- (sum(weights * mean_returns_q3) - target_ret)^2; return(t(weights) %*% cov_matrix_q3 %*% weights + penalty * (sum_constraint + return_constraint)) }
  result <- optim(par = rep(1/length(mean_returns_q3), length(mean_returns_q3)), fn = objective_penalized, method = "L-BFGS-B", lower = 0, upper = 1)
  weights <- result$par / sum(result$par); stats <- calculate_portfolio_stats(weights, mean_returns_q3, cov_matrix_q3); return(data.frame(ret = stats['ret'], sd = stats['sd']))
})
frontier_no_shorts_q3 <- do.call(rbind, frontier_no_shorts_list_q3)
ones_q3 <- rep(1, length(mean_returns_q3)); inv_cov_q3 <- solve(cov_matrix_q3); a_q3 <- as.numeric(t(ones_q3) %*% inv_cov_q3 %*% ones_q3); b_q3 <- as.numeric(t(ones_q3) %*% inv_cov_q3 %*% mean_returns_q3); c_q3 <- as.numeric(t(mean_returns_q3) %*% inv_cov_q3 %*% mean_returns_q3)
frontier_sd_shorts_q3 <- sqrt((a_q3 * target_returns_q3^2 - 2*b_q3*target_returns_q3 + c_q3) / (a_q3*c_q3 - b_q3^2)); frontier_shorts_q3 <- data.frame(ret = target_returns_q3, sd = frontier_sd_shorts_q3)
individual_assets_q3 <- data.frame(TICKER = names(mean_returns_q3), ret = mean_returns_q3, sd = sqrt(diag(cov_matrix_q3)))
vw_portfolio_6_q3 <- q3_data_filtered %>% filter(PERMNO %in% stocks_for_q2_q3) %>% group_by(date) %>% summarise(Portfolio_Return = weighted.mean(RET, PRC * SHROUT, na.rm = TRUE))
vw_port_6_stats_q3 <- data.frame(TICKER = "VW Portfolio", ret = mean(vw_portfolio_6_q3$Portfolio_Return, na.rm = TRUE), sd = sd(vw_portfolio_6_q3$Portfolio_Return, na.rm = TRUE))
sink("Q3_results.txt"); cat("Question 3 Results: Input Data (Last 10 Years)\n\n"); cat("Individual Asset Performance:\n"); print(individual_assets_q3); cat("\nValue-Weighted Portfolio of the 6 Assets:\n"); print(vw_port_6_stats_q3); sink()
mv_plot_q3 <- ggplot() + geom_path(data = frontier_shorts_q3, aes(x = sd, y = ret, color = "Short Sales Allowed"), size = 1) + geom_path(data = frontier_no_shorts_q3, aes(x = sd, y = ret, color = "No Short Sales"), size = 1) + geom_point(data = individual_assets_q3, aes(x = sd, y = ret), size = 3) + geom_text(data = individual_assets_q3, aes(x = sd, y = ret, label = TICKER), vjust = -1) + geom_point(data = vw_port_6_stats_q3, aes(x = sd, y = ret), size = 4, color = "darkgreen", shape = 18) + geom_text(data = vw_port_6_stats_q3, aes(x = sd, y = ret, label = TICKER), vjust = -1.5) + labs(title = "Mean-Variance Frontier (Last 10 Years)", x = "Risk (StDev)", y = "Expected Return", color = "Frontier Type") + theme_minimal() + theme(legend.position="top") + scale_y_continuous(labels = scales::percent) + scale_x_continuous(labels = scales::percent)
ggsave("mean_variance_frontier_last_10_years.png", plot = mv_plot_q3, width = 10, height = 7)
cat("Question 3 plot and results saved.\n\n")

# --- 4. Question 4: Beta Estimation and Distribution ---
cat("--- Starting Question 4: Beta Estimation and Distribution ---\n")
cat("Using market portfolio from Q1 and estimating betas for each firm...\n")

# 4.1 Prepare data with excess returns for both stocks and the market
# Reuse the `portfolio_returns` object from Q1 as our market
market_excess_returns <- portfolio_returns %>%
  rename(Market_Return_VW = Portfolio_Return) %>%
  inner_join(all_data %>% select(date, RF) %>% distinct(), by = "date") %>%
  mutate(Mkt_Excess_Return = Market_Return_VW - RF) %>%
  select(date, Mkt_Excess_Return)

# Prepare individual stock data with excess returns
stock_excess_returns <- all_data %>%
  mutate(Stock_Excess_Return = RET - RF) %>%
  select(date, PERMNO, TICKER, Stock_Excess_Return)

# Combine into a single data frame for regression
regression_data <- inner_join(stock_excess_returns, market_excess_returns, by = "date")

# 4.2 Run CAPM regression for each firm, grouping by PERMNO only
beta_estimates_raw <- regression_data %>%
  group_by(PERMNO) %>% 
  filter(n() > 24) %>% 
  do(tidy(lm(Stock_Excess_Return ~ Mkt_Excess_Return, data = .))) %>%
  ungroup() %>%
  filter(term == "Mkt_Excess_Return") %>%
  select(PERMNO, beta = estimate)

# Create a lookup for the most recent ticker for each PERMNO for clean reporting
recent_tickers <- all_data %>%
  group_by(PERMNO) %>%
  filter(date == max(date)) %>%
  summarise(TICKER = first(TICKER))

# Join the most recent ticker to the beta estimates
beta_estimates <- beta_estimates_raw %>%
  left_join(recent_tickers, by = "PERMNO") %>%
  select(PERMNO, TICKER, beta) # Reorder columns for clarity

# 4.3 Save the beta estimates to a file
sink("Q4_beta_estimates.txt")
cat("Question 4: CAPM Beta Estimates for Each Firm (Consolidated by PERMNO)\n"); cat("--------------------------------------------------\n"); print(as.data.frame(beta_estimates)); sink()
cat("Question 4 beta estimates saved to Q4_beta_estimates.txt\n")

# 4.4 Create and save the histogram of the beta distribution
beta_histogram <- ggplot(beta_estimates, aes(x = beta)) +
  geom_histogram(binwidth = 0.2, fill = "dodgerblue", color = "white", boundary = 0) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.05, y = Inf, label = "Market Beta = 1.0", hjust = 0, vjust = 2, color = "red", fontface = "bold") +
  labs(title = "Cross-Sectional Distribution of CAPM Betas", subtitle = "Market defined as the value-weighted average of all firms in the sample", x = "Estimated Beta (β)", y = "Number of Firms") +
  scale_x_continuous(breaks = seq(floor(min(beta_estimates$beta) * 5) / 5, ceiling(max(beta_estimates$beta) * 5) / 5, by = 0.2)) +
  theme_minimal(base_size = 14)

ggsave("beta_distribution_histogram.png", plot = beta_histogram, width = 11, height = 7, dpi = 300)
cat("Question 4 plot saved as 'beta_distribution_histogram.png'.\n\n")

# --- 5. Question 5: Scatter Plot of Average Return vs. Beta ---
cat("--- Starting Question 5: Scatter Plot of Average Return vs. Beta ---\n")

# 5.1 Load the necessary data
cat("Loading data and beta estimates from Question 4...\n")


# 5.2 Calculate Average Excess Return for each PERMNO
# group by PERMNO to get one average return for each unique company.
average_excess_returns <- all_data %>%
  # First, calculate the excess return for every observation
  mutate(Excess_Return = RET - RF) %>%
  # Now, group by company and calculate the mean over its entire history
  group_by(PERMNO) %>%
  summarise(
    Avg_Excess_Return = mean(Excess_Return, na.rm = TRUE)
  ) %>%
  ungroup()

# 5.3 Merge the average returns with the beta estimates
# use an inner_join to ensure we only plot companies for which we have both a beta and an average return.
plot_data <- inner_join(beta_estimates, average_excess_returns, by = "PERMNO")

# 5.4 Create the Scatter Plot
cat("Generating the scatter plot...\n")
return_beta_scatterplot <- ggplot(plot_data, aes(x = beta, y = Avg_Excess_Return)) +
  geom_point(alpha = 0.6, color = "blue") +
   geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + # Add a zero line for reference
  labs(
    title = "Average Monthly Excess Return vs. Beta",
    subtitle = "Each point represents a single firm from the sample",
    x = "Estimated Beta (β)",
    y = "Average Monthly Excess Return"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_minimal(base_size = 14)


ggsave("return_vs_beta_scatterplot.png", plot = return_beta_scatterplot, width = 11, height = 7, dpi = 300)
cat("Question 5 plot saved as 'return_vs_beta_scatterplot.png'.\n\n")

cat("--- Question 5 finished successfully. ---\n")


  cat("--- Starting Question 6: Security Market Line (SML) Estimation ---\n")

  # --- 6. Question 6: Security Market Line (SML) Estimation ---
  cat("--- Starting Question 6: Security Market Line (SML) Estimation ---\n")
  
  # 6.1 Prepare the data for the cross-sectional regression
  cat("Loading data and preparing inputs...\n")

  average_excess_returns <- all_data %>%
    mutate(Excess_Return = RET - RF) %>%
    group_by(PERMNO) %>%
    summarise(Avg_Excess_Return = mean(Excess_Return, na.rm = TRUE)) %>%
    ungroup()
  
  sml_data <- inner_join(beta_estimates, average_excess_returns, by = "PERMNO")
  
  # 6.2 Estimate the Security Market Line
  cat("Running the cross-sectional regression to estimate the SML...\n")
  sml_model <- lm(Avg_Excess_Return ~ beta, data = sml_data)
  sml_results <- tidy(sml_model)
  
  # 6.3 Save the estimation results to a table
  cat("Saving SML estimation results to a file...\n")
  sink("Q6_sml_results.txt")
  cat("Question 6: Security Market Line (SML) Estimation Results\n")
  cat("Regression: Avg_Excess_Return ~ Intercept + beta\n")
  cat("--------------------------------------------------\n")
  print(sml_results)
  sink()
  cat("Question 6 SML results saved to Q6_sml_results.txt\n")
  
  # 6.4 Create the updated scatter plot with the estimated SML
  cat("Generating the updated scatter plot with the SML...\n")
  sml_plot <- ggplot(sml_data, aes(x = beta, y = Avg_Excess_Return)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_abline(
      intercept = sml_results$estimate[1],
      slope = sml_results$estimate[2],
      color = "red",
      size = 1
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(
      title = "Security Market Line (SML) for the Sample",
      subtitle = "Average Monthly Excess Return vs. Beta with Estimated SML",
      x = "Estimated Beta (β)",
      y = "Average Monthly Excess Return"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    theme_minimal(base_size = 14)
  
  # Save the plot to a high-quality PNG file
  ggsave("sml_scatterplot.png", plot = sml_plot, width = 11, height = 7, dpi = 300)
  cat("Question 6 plot saved as 'sml_scatterplot.png'.\n\n")
  
  cat("--- Question 6 finished successfully. ---\n")
  
  # --- 7. Question 7: Multi-Factor Model Regression ---
  cat("--- Starting Question 7: Multi-Factor Model Regression ---\n")
  book_value_df <- read_excel(excel_file_path, sheet = book_value_sheet_name)
  
  book_value_clean <- book_value_df %>%
    rename(PERMNO = LPERMNO, BV_date = datadate, BVPS = bkvlps) %>%
    mutate(fyear = year(BV_date)) %>%
    group_by(PERMNO, fyear) %>%
    slice(1) %>%
    ungroup() %>%
    select(PERMNO, fyear, BVPS)
  
  monthly_data <- all_data %>% mutate(year = year(date), Size = PRC * SHROUT)
  firm_characteristics <- monthly_data %>% mutate(fyear = year - 1) %>% left_join(book_value_clean, by = c("PERMNO", "fyear")) %>% group_by(PERMNO) %>% fill(BVPS, .direction = "down") %>% ungroup() %>% filter(BVPS > 0 & PRC > 0) %>% mutate(Book_to_Market = BVPS / PRC) %>% group_by(PERMNO) %>% summarise(Avg_Excess_Return = mean(RET - RF, na.rm = TRUE), Avg_Log_Size = mean(log(Size), na.rm = TRUE), Avg_BM = mean(Book_to_Market, na.rm = TRUE)) %>% ungroup()
  final_regression_data <- firm_characteristics %>% inner_join(beta_estimates, by = "PERMNO")
  three_factor_model <- lm(Avg_Excess_Return ~ beta + Avg_Log_Size + Avg_BM, data = final_regression_data)
  three_factor_results <- tidy(three_factor_model)
  sink("Q7_three_factor_results.txt"); print(three_factor_results); sink()
  cat("Question 7 results saved.\n\n")
  
  cat("--- All questions processed. Script finished successfully. ---\n")