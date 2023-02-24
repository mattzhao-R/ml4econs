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
           aftexpl = as.factor(ifelse((Date > my('Feb 2015')),1,0)),
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
    step_dummy(aftexpl,fct_month,fct_year,fct_quarter,fct_week,
               fct_wday,fct_qday) %>%
    step_interact(terms = ~ contains('fct_month'):starts_with('fct_year')) %>%
    step_interact(terms = ~ contains('fct_quarter'):starts_with('fct_year')) %>%
    step_interact(terms = ~ ends_with('aftexpl_X1'):starts_with('fct_year')) %>%
    step_interact(terms = ~ ends_with('aftexpl_X1')*contains('fct_quarter')) %>%
    step_interact(terms = ~ ends_with('aftexpl_X1')*contains('fct_month')) %>%
    prep(training = analysis_df) %>%
    bake(analysis_df)
  
  reg_df
}

# returns a regression model
run_tsls <- function(df, w, z){
  # w - exogenous
  # z - instruments
  form_y <- 'lsales ~'
  form_endog <- '| lprice |'
  full_form <- paste(form_y,w,form_endog,z)
  
  tsls <- ivreg(formula = as.formula(full_form), 
                data = df)
  tsls
}

# run first stage
run_fstg <- function(df, w, z){
  # w - exogenous
  # z - instruments
  form_endog <- 'lprice ~'
  first_form <- paste(form_endog,z,'+',w)
  
  fstg <- ivreg(formula = as.formula(first_form), 
                data = df)
  fstg
}


# preliminary regressions ----

## baseline ----
w <- 't_expl + t_expl_sq'
z <- 'aftexplX1'
df_2012_2018 <- make_regdf(2012,2018)
base_reg_2012_2018 <- run_tsls(df_2012_2018, w, z)
df_2013_2017 <- make_regdf(2013,2017)
base_reg_2013_2017 <- run_tsls(df_2013_2017, w, z)
df_2014_2016 <- make_regdf(2014,2016)
base_reg_2014_2016 <- run_tsls(df_2014_2016, w, z)

stargazer(base_reg_2012_2018,base_reg_2013_2017,base_reg_2014_2016,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          column.labels = c('2012-2018','2013-2017','2014-2016','2014-2018'),
          keep.stat = c('chi2'),
          title='Baseline Regressions')

base_fstg_2012_2018 <- run_fstg(df_2012_2018, w, z)
base_fstg_2013_2017 <- run_fstg(df_2013_2017, w, z)
base_fstg_2014_2016 <- run_fstg(df_2014_2016, w, z)

stargazer(base_fstg_2012_2018,base_fstg_2013_2017,base_fstg_2014_2016,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          keep.stat = c('chi2'),
          title='Baseline First Stage')

## other specifications

df_2012_2018 <- make_regdf(2012,2018)

# redo these
vars <- colnames(df_2012_2018)
aft_all <- vars[str_detect(vars,'aftexpl')]
month_all <- vars[str_detect(vars,'_month_')]
year_all <- vars[str_detect(vars,'_year_')]
quarter_all <- vars[str_detect(vars,'_quarter_')]
year_ind <- setdiff(year_all,month_all)[1:6]
quarter.year <- setdiff(year_all,month_all)[7:24]
aftexpl.quarter <- intersect(aft_all,quarter_all)[1:3]
aftexpl.quarter.year <- setdiff(year_all,month_all)[25:42]
month.year <- setdiff(year_all,quarter_all)[7:72]
aftexpl.month <- intersect(aft_all,month_all)[1:11]
aftexpl.month.year <- setdiff(year_all,quarter_all)[73:138]
month_ind <- setdiff(month_all,year_all)[1:11]
quarter_ind <- setdiff(quarter_all,year_all)[1:3]
aftexpl.year <- intersect(aft_all,year_all)[1:11]

w <- paste('t_expl + t_expl_sq',
           paste0(month_ind,collapse = '+'),
           paste0(aftexpl.quarter.year,collapse = '+'),
           sep = ' + ')

z.yr_reg_2012_2018 <- run_tsls(df_2012_2018, w, z)

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

