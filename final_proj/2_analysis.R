# Last updated: Feb 29, 2023

library(MASS)
library(data.table)
library(hdm)
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
library(plm)
source('3_estimators.R')

# ddir_json <- 
ddir_matt <- './data/'
ddir <- ddir_matt

df <- read.csv(file.path(ddir, 'final_data.csv'))[,-1] %>%
  mutate(date = ymd(date))


## variable column extraction ----

#### mthyr fe + instrument int ----
# indiv_mthyr_fe <- colnames(df)[str_detect(colnames(df),'mthyrfe_')]
# mthyr_fe <- paste(indiv_mthyr_fe,
#                 collapse = ' + ')
# indiv_mthyr.aftexpl <- colnames(df)[str_detect(colnames(df),'mthyr.aftexplmth')]
# mthyr.aftexpl <- paste(indiv_mthyr.aftexpl,
#                 collapse = ' + ')
# indiv_mthyr.aftexpl.dist_to_ref <- colnames(df)[str_detect(colnames(df),'mthyr.aftexpl.dist_to_ref')]
# mthyr.aftexpl.dist_to_ref <- 
#   paste(indiv_mthyr.aftexpl.dist_to_ref,
#         collapse = ' + ')

### time fe + instrument int ----

indiv_time_fe <- colnames(df)[str_detect(colnames(df),'timefe_')]
time_fe <- paste(indiv_time_fe,
                  collapse = ' + ')
indiv_time.aftexpl <- colnames(df)[str_detect(colnames(df),'timefe.aftexpl_')]
time.aftexpl <- paste(indiv_time.aftexpl,
                       collapse = ' + ')
indiv_time.aftexpl.dist_to_ref <- colnames(df)[str_detect(colnames(df),'timefe.aftexpl.dist_to_ref')]
time.aftexpl.dist_to_ref <- 
  paste(indiv_time.aftexpl.dist_to_ref,
        collapse = ' + ')


### county fe + instrument int ----
indiv_cnty_fe <- colnames(df)[str_detect(colnames(df), "cntyfe_")]
cnty_fe <- paste(indiv_cnty_fe,
                collapse = ' + ')
indiv_cnty.aftexpl <- colnames(df)[str_detect(colnames(df),'cntyfe.aftexpl_')]
cnty.aftexpl <- paste(indiv_cnty.aftexpl,
                     collapse = ' + ')

### aftexpl.dist dummies ----
indiv_aftexpl.dist_dummies <- colnames(df)[str_detect(colnames(df), "aftexpl.dist_to_ref\\d+")]
aftexpl.dist_dummies <- paste0(indiv_aftexpl.dist_dummies, collapse = ' + ')

## preliminary analysis ----
### regressions ----
#### see if distance to refinery has relevance on its own ----
ols_lprice_dist_to_ref_pre_expl <-
  lm(lprice ~ dist_to_refinery,
     data = df %>% filter(date > my('Feb 2015')))

stargazer(ols_lprice_dist_to_ref_pre_expl,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          title='Log Price vs Distance to Refinery')

# ols_base <- feols(lsales ~ lprice, 
#                   fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
#                   data = df)
# 
# etable(ols_base,
#        tex = T,
#        title = 'Baseline Regressions')


#### first stage ----

fstg_fm <- paste(
  'lprice ~',
  aftexpl.dist_dummies
)
lm_fstg <- feols(
  as.formula(fstg_fm),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)
# fm2 <- 'lprice ~ aftexpl + aftexpl.dist_to_ref'
# lm_new_fstg <- feols(
#   as.formula(fm2),
#   fixef = c(indiv_time_fe,indiv_cnty_fe),
#   data = df
# )
# fm3 <- paste(
#   'lprice ~ ',
#   paste(
#     time.aftexpl,
#     'aftexpl + aftexpl.dist_to_ref',
#     sep = '+')
# )
# lm_new_fstg_wtime.aftexpl <- feols(
#   as.formula(fm3),
#   fixef = c(indiv_time_fe,indiv_cnty_fe),
#   data = df
# )

fitstat(lm_fstg,'wald')
etable(lm_fstg,
       tex = T,
       title = 'First Stage',
       fitstat = c('n','wf.stat'),
       style.tex = style.tex('aer'))


#### checking for interaction significance in first stage (ovb) ----
##### regress aftexpl.county on U (fstg residuals) ----

betapvals_aftexpl.cnty_on_U <- c()
U <- lm_fstg$residuals
# beta <- c()

for (i in 1:length(indiv_cnty.aftexpl)) {
  y <- indiv_cnty.aftexpl[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  # temp_beta <- temp_lm$coefficients['U']
  # beta <- c(beta,
  #           temp_beta)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_aftexpl.cnty_on_U <- c(betapvals_aftexpl.cnty_on_U,
                                   beta_pval)
}

sum(betapvals_aftexpl.cnty_on_U < 0.01) / length(betapvals_aftexpl.cnty_on_U)


##### regress time.aftexpl, time.aftexpl.dist_to_ref on U ----

betapvals_time.aftexpl_on_U <- c()
U <- lm_fstg$residuals

for (i in 1:length(indiv_time.aftexpl)) {
  y <- indiv_time.aftexpl[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_time.aftexpl_on_U <- c(betapvals_time.aftexpl_on_U,
                                   beta_pval)
}

sum(
  betapvals_time.aftexpl_on_U[!is.na(betapvals_time.aftexpl_on_U)] < 0.01
  ) 


betapvals_time.aftexpl.dist_to_expl_on_U <- c()
U <- lm_fstg$residuals

for (i in 1:length(indiv_time.aftexpl.dist_to_ref)) {
  y <- indiv_time.aftexpl.dist_to_ref[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_time.aftexpl.dist_to_expl_on_U <- c(betapvals_time.aftexpl.dist_to_expl_on_U,
                                    beta_pval)
}

sum(
  betapvals_time.aftexpl.dist_to_expl_on_U[!is.na(betapvals_time.aftexpl.dist_to_expl_on_U)] < 0.01
  )




#### first stage (plm) 
# fstg_plm_fm <- 'lprice ~ aftexpl + aftexpl.dist_to_ref'
# plm_fstg <- plm(
#   formula = as.formula(fstg_plm_fm),
#   data = df,
#   effect = 'twoways',
#   model = 'random'
# )
# 
# stargazer(plm_fstg,
#           type='latex',digits=3, column.sep.width = "5pt",
#           title='PLM First Stage with Random Effects')


#### lprice ~ fixed effects ----
model <- feols(
  lprice ~ 1,
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)
# write.csv(data.frame(residuals=model$residuals),
#           'data/price_fe_residuals.csv')


#### lprice ~ aftexpl.dist dummies ----
model2 <- feols(
  as.formula(paste('lprice ~', aftexpl.dist_dummies)),
  data = df %>% mutate(lprice = model$residuals)
)
fitstat(model2,'f')

## Estimators (OLS, TSLS, JIVE, RJIVE, Post-Lasso) ----
### OLS ----
#### fixest 
fm_ols <- 'lsales ~ lprice'
ols_est <- feols(
  as.formula(fm_ols),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)

#### plm 
ols_plm_fm <- 'lsales ~ lprice'
plm_ols <- plm(
  formula = as.formula(ols_plm_fm),
  data = df,
  effect = 'twoways',
  model = 'random'
)

### TSLS ----
#### fixest 
fm_tsls <- paste(
  'lsales ~ ',
  '1 | lprice ~', 
  aftexpl.dist_dummies)
tsls_est <- feols(
  as.formula(fm_tsls),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)

# #### plm 
# tsls_plm_fm <- 'lsales ~ lprice | aftexpl + aftexpl.dist_to_ref'
# plm_tsls <- plm(
#   formula = as.formula(tsls_plm_fm),
#   data = df,
#   effect = 'time',
#   model = 'random'
# ) #NOTE: PLM TSLS did not run with twoway (error: twoway not supported)

### jive ----
# from estimators script
## manually demeaning variables for FE-IV
noint_df <- df %>%
  select(county,time,lprice,lsales,
         matches("aftexpl.dist_to_ref\\d+"))
cnty_means <- noint_df %>%
  group_by(time) %>%
  summarise(across(!county,mean), .groups='keep')
temp <- left_join(noint_df,cnty_means, by = 'time',
                 suffix = c('','_means'))
p_temp <- temp %>%
  select(!contains('_means'))
for (col in colnames(temp)) {
  if ((col != 'county') & (col != 'time') & (!str_detect(col,'_means'))){
    p_temp[,col] <- 
      (temp %>% select(all_of(col))) - (temp %>% select(contains(paste0(col,'_means'))))
  }
}
partialled_df <- bind_cols(
  p_temp,
  df %>% select(contains('timefe_'))
)
# write.csv(partialled_df,file.path(ddir,'partialled_noint.csv'))

y <- data.matrix(partialled_df$lsales)
n <- length(y)
X <- data.matrix(
  partialled_df %>%
    dplyr::select(lprice,contains('timefe_'))
)
Z <- data.matrix(
  partialled_df %>%
    dplyr::select(starts_with('aftexpl.dist_to_ref'),
                  contains('timefe_'))
  )

jive_est <- jive.est(y,X,Z,SE=T)
jive_pointest <- jive_est[[1]][[1]]
jive_se <- jive_est[[2]][[1]]
jive_tstat <- jive_pointest / jive_se
jive_pval <- pt(abs(jive_tstat),df=n-84-138-1,lower.tail = F)
jive_res <- y - (X %*% jive_est$est)
# write.csv(bind_cols(data.frame(residuals=jive_res),
#                     partialled_df),
#           'data/jive_residuals.csv')


### rjive ----



### Post Lasso ----
y <- data.matrix(df$lsales)
n <- length(y)
d <- data.matrix(
  df %>%
    dplyr::select(lprice)
)
x <- data.matrix(
  df %>%
    dplyr::select(contains('timefe_'),contains('cntyfe_'))
)
z <- data.matrix(
  df %>%
    dplyr::select(matches('aftexpl.dist_to_ref\\d+'))
)

# postlassoIV_est <- rlassoIV(x,d,y,z,select.Z = T,select.X = F,
#                         post = T)
PLIV_fm <- paste(
  'lsales ~',
  'lprice +',
  time_fe,'+',
  cnty_fe,'|',
  aftexpl_dummies,'+',
  time_fe,'+',
  cnty_fe
)
postlassoIV_est <- rlassoIV(as.formula(PLIV_fm),
                            data = df,
                            select.Z = T,select.X = F,
                            post = T)


## tests ----
### durbin watson AR1 ----
dwt_pvals <- c()
for (cnty in df$county) {
  t_df <- df %>%
    filter(county == cnty)
  dwt_results <- durbinWatsonTest(
    lm(lprice ~ time,data=t_df),alt='positive')
  dwt_pvals <- c(dwt_pvals,dwt_results[['p']])
}

### Tables ----
#### main results ----
etable(ols_est,tsls_est,ols_est,ols_est,
       tex = T,
       title = 'Main Results',
       fitstat = c('n','ivf1','ivwald'),
       style.tex = style.tex('aer')
       )

#### plm 
# stargazer(plm_ols,plm_tsls,
#           type='latex',digits=3, column.sep.width = "5pt",
#           title='Main Results')


### Graphs ----



