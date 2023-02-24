# Last updated: Feb 21, 2023

library(tidyverse)
library(lubridate)
library(haven)
library(stargazer)
library(xtable)
library(lmtest)
library(sandwich)
library(broom)
library(whitestrap)
library(recipes)
library(ivreg)
library(car)

# ddir_json <- 
ddir_matt <- './data/'
ddir <- ddir_matt

# data preparation ----
wc_price_df <- read.csv(
  file.path(ddir, "w_retail_gas_prices.csv"))

wc_sales_df <- read.csv(
  file.path(ddir, "m_retail_gas_sales.csv"))

## aggregate price to monthly freq ----
wc_price_m <- wc_price_df
colnames(wc_price_m)[2] <- 'all_retail'
wc_price_m <- wc_price_m %>%
  select(Date,all_retail) %>%
  mutate(Date = mdy(Date),
         month = month(Date),
         year = year(Date)) %>%
  group_by(year,month) %>%
  summarise(price = mean(all_retail), .groups='keep')

## merge price and sales data ----
wc_sales_m <- wc_sales_df %>%
  select(Date,contains('West.Coast')) %>%
  rename(sales = contains('West.Coast')) %>%
  mutate(Date = my(Date),
         month = month(Date),
         year = year(Date))

merged <- left_join(wc_price_m,wc_sales_m,
                    by = c('year','month')) %>%
  drop_na(price,sales)

# create instruments + indicators ----
make_regdf <- function(lower,upper){
  analysis_df <- merged %>%
    filter((year >= lower) & (year <= upper)) %>%
    mutate(lprice = log(price),
           lsales = log(sales),
           aftexpl = ifelse((Date > my('Feb 2015')),1,0),
           t_expl = as.numeric(Date - my('Feb 2015')),
           t_expl_sq = t_expl ** 2,
           fct_month = as.factor(month),
           fct_year = as.factor(year),
           fct_quarter = as.factor(quarter(Date)),
           fct_week = as.factor(week(Date)),
           fct_wday = as.factor(wday(Date)),
           fct_qday = as.factor(qday(Date)))
  
  first_stage <- recipe(lprice ~ ., data = analysis_df)
  
  reg_df <- first_stage %>%
    step_dummy(fct_month,fct_year,fct_quarter,fct_week,
               fct_wday,fct_qday, one_hot = T) %>%
    prep(training = analysis_df) %>%
    bake(analysis_df)
  
  reg_df
}

# returns a regression model
run_reg <- function(df, w, z){
  # w - exogenous
  # z - instruments
  form_y <- 'lsales ~'
  form_endog <- '| lprice |'
  full_form <- paste(form_y,w,form_endog,z)
  
  tsls <- ivreg(formula = as.formula(full_form), 
                data = df)
  tsls
}

# preliminary regressions ----

base_models <- list()
for (x in 0:2) {
  w <- 't_expl + t_expl_sq'
  z <- 'aftexpl'
  if (x != 0){
    df_p <- make_regdf(2012+x,2018)
    df_m <- make_regdf(2012,2018-x)
    df_pm <- make_regdf(2012+x,2018-x)
    reg_p <- run_reg(df_p, w, z)
    reg_m <- run_reg(df_m, w, z)
    reg_pm <- run_reg(df_pm, w, z)
    base_models <- append(base_models,
                          list(reg_p,reg_m,reg_pm))
  } else {
    df <- make_regdf(2012,2018)
    reg <- run_reg(df,w,z)
    base_models <- append(base_models,reg)
  }
}
for (i in seq(1,40,8)) {
  if (i != 33){
    idx <- seq(i,i+7)
  } else{
    idx <- seq(i,i+6)
  }
  stargazer(base_models[idx],
            out = paste0('base',idx[1],'_',
                         idx[length(idx)],'.tex'),
            type='latex',digits=3, column.sep.width = "5pt",
            dep.var.labels.include = F,
            keep.stat = c('chi2'),
            title='Baseline Regressions (1-8)',
            notes = c('Baseline TSLS regression run with progressively 
                    smaller windows. First run for 2012-2018, then 
                    remaining with +1 lower, -1 upper, then +1 lower 
                    and -1 upper for window size bounds. '))
}



# graphs ----
ggplot(data = analysis_df, 
       mapping = aes(x = lprice, y = lsales)) + 
  geom_point()

ggplot(data = analysis_df, 
       mapping = aes(x = Date)) + 
  geom_point(mapping = aes(
    y = (analysis_df$lprice - mean(analysis_df$lprice)) / 
      sd(analysis_df$lprice)),
    color = 'red') +
  geom_point(mapping = aes(
    y = (analysis_df$lsales - mean(analysis_df$lsales)) / 
      sd(analysis_df$lsales)), 
    color = 'blue') + 
  scale_x_date(date_breaks = '1 year',
               date_labels = '%Y'
               ,
               lim = c(ymd('2014-1-1'),
                       ymd('2016-1-1'))
               )
  

ggplot(data = analysis_df, 
       mapping = aes(x = aftexpl, y = lprice)) + 
  geom_point()

