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
library(fixest)
library(ivreg)
library(car)

# ddir_json <- 
ddir_matt <- './data/'
ddir <- ddir_matt

# data preparation ----

## prep functions ----
make_longdf <- function(filename){
  df <- read.csv(file.path(ddir, filename)) %>%
    mutate(Date = mdy(Date),
           month = month(Date),
           year = year(Date)) %>%
    pivot_longer(cols=!c(Date,month,year),
                 names_to = 'county',
                 values_to = 'value')
    df
}

make_covariates <- function(df){
  df <- df %>%
    mutate(lprice = log(price),
           lsales = log(sales),
           quarter = quarter(Date),
           month = month(Date),
           year = year(Date),
           aftexpl = ifelse((Date > my('Feb 2015')),1,0),
           t_expl = as.numeric(Date - my('Feb 2015')),
           t_expl_sq = t_expl ** 2)
  df
}


## Data Processing - making dataframe for regression ----

m_sales_df <- make_longdf('m_county_sales_2012_2022.csv') %>%
  rename(sales = value)
w_price_df <- make_longdf('w_county_price_idxsf_2000_2023.csv') %>%
  rename(price = value)
m_price_df <- w_price_df %>%
  group_by(year,month,county) %>%
  summarise(price = mean(price), .groups = 'keep')

temp <- left_join(m_price_df,m_sales_df,
                    by = c('year','month','county')) %>%
  drop_na(price,sales) %>%
  select(Date,county,price,sales,month,year)
m_merged_df <- make_covariates(temp)

### adding indicators ----

qr <- paste0('qr',m_merged_df$quarter)
mth <- paste0('mth',m_merged_df$month)
yr <- paste0('yr',m_merged_df$year)
aftexpl <- m_merged_df$aftexpl

analysis_df <- bind_cols(m_merged_df,
      data.frame(i(mth, ref = T)),
      data.frame(i(qr, ref = T)),
      data.frame(i(yr, ref = T)),
      data.frame(i(mth, i.yr, ref = T)),
      data.frame(i(qr, i.yr, ref = T)),
      data.frame(i(mth, aftexpl, ref = T)),
      data.frame(i(qr, aftexpl, ref = T)),
      data.frame(i(yr, aftexpl, ref = T)))

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

summary(lm(lsales ~ lprice, data = analysis_df))
summary(ivreg(lsales ~ t_expl + t_expl_sq | lprice | aftexpl,
              data = analysis_df))

summary(run_tsls(analysis_df,'t_expl + t_expl_sq','aftexpl'))

# stargazer(base_reg_2014_2016,
#           type='latex',digits=3, column.sep.width = "5pt",
#           dep.var.labels.include = F,
#           column.labels = c('2014-2018'),
#           keep.stat = c('chi2'),
#           title='Baseline Regressions')

# graphs ----
## dingel plot
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

