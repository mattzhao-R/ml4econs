# Last updated: Feb 27, 2023

library(MASS)
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
indiv_aftexpl_dummies <- colnames(df)[str_detect(colnames(df), "aftexpl.dist_to_ref\\d+")]
aftexpl_dummies <- paste0(indiv_aftexpl_dummies, collapse = ' + ')

## functions ----



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


#### saturated first stage (fixest) ----
fm <- paste(
  'lprice ~ ',
  paste(
    time.aftexpl,
    time.aftexpl.dist_to_ref,
    'aftexpl + aftexpl.dist_to_ref',
    sep = '+')
)
lm_old_fstg <- feols(
  as.formula(fm),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)
fm2 <- 'lprice ~ aftexpl + aftexpl.dist_to_ref'
lm_new_fstg <- feols(
  as.formula(fm2),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)
fm3 <- paste(
  'lprice ~ ',
  paste(
    time.aftexpl,
    'aftexpl + aftexpl.dist_to_ref',
    sep = '+')
)
lm_new_fstg_wtime.aftexpl <- feols(
  as.formula(fm3),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)

etable(lm_old_fstg,lm_new_fstg,lm_new_fstg_wtime.aftexpl,
       tex = T,
       title = 'First Stages',
       fitstat = 'f',
       style.tex = style.tex('aer'))


# checking for interaction significance in first stage (ovb)
## regress aftexpl.county on U (fstg residuals) 

betapvals_aftexpl.cnty_on_U <- c()
U <- lm_new_fstg$residuals
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


#### regress time.aftexpl, time.aftexpl.dist_to_ref on U ----

betapvals_time.aftexpl_on_U <- c()
U <- lm_new_fstg$residuals

for (i in 1:length(indiv_time.aftexpl)) {
  y <- indiv_time.aftexpl[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_time.aftexpl_on_U <- c(betapvals_time.aftexpl_on_U,
                                   beta_pval)
}

sum(betapvals_time.aftexpl_on_U < 0.01) / length(betapvals_time.aftexpl_on_U)


betapvals_time.aftexpl.dist_to_expl_on_U <- c()
U <- lm_new_fstg$residuals

for (i in 1:length(indiv_time.aftexpl.dist_to_ref)) {
  y <- indiv_time.aftexpl.dist_to_ref[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_time.aftexpl.dist_to_expl_on_U <- c(betapvals_time.aftexpl.dist_to_expl_on_U,
                                    beta_pval)
}

sum(betapvals_time.aftexpl.dist_to_expl_on_U < 0.01) / length(betapvals_time.aftexpl.dist_to_expl_on_U)




#### first stage (plm) ----
fstg_plm_fm <- 'lprice ~ aftexpl + aftexpl.dist_to_ref'
plm_fstg <- plm(
  formula = as.formula(fstg_plm_fm),
  data = df,
  effect = 'twoways',
  model = 'random'
)

stargazer(plm_fstg,
          type='latex',digits=3, column.sep.width = "5pt",
          title='PLM First Stage with Random Effects')


#### lprice ~ fixed effects ----
model <- feols(
  lprice ~ 1,
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)
model$residuals

## Estimators (OLS, TSLS, JIVE, RJIVE, Post-Lasso) ----

y <- data.matrix(df$lsales)
n <- length(y)
X <- cbind(rep(1,n),data.matrix(df %>%
                   select(lprice,
                          contains('timefe_'),contains('cntyfe_')))
           )
Z <- cbind(rep(1,n),data.matrix(df %>%
                   select(aftexpl,aftexpl.dist_to_ref,
                          contains('timefe_'),contains('cntyfe_')))
           )

### OLS ----
#### fixest ----
fm_ols <- 'lsales ~ lprice'
ols_est <- feols(
  as.formula(fm_ols),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)

#### plm ----
ols_plm_fm <- 'lsales ~ lprice'
plm_ols <- plm(
  formula = as.formula(ols_plm_fm),
  data = df,
  effect = 'twoways',
  model = 'random'
)

### TSLS ----
# X_hat <- lm_new_fstg$fitted.values
# fm_tsls <- paste(
#   'lsales ~ ',
#   'lprice ~',
#   paste(
#     time_fe,
#     cnty_fe,
#     'aftexpl + aftexpl.dist_to_ref',
#     sep = '+')
# )
#### fixest ----
fm_tsls <- paste(
  'lsales ~ ',
  '1 | lprice ~', 
  aftexpl_dummies)
tsls_est <- feols(
  as.formula(fm_tsls),
  fixef = c(indiv_time_fe,indiv_cnty_fe),
  data = df
)
# tsls.est(y,X,Z,SE=T)



#### plm ----
tsls_plm_fm <- 'lsales ~ lprice | aftexpl + aftexpl.dist_to_ref'
plm_tsls <- plm(
  formula = as.formula(tsls_plm_fm),
  data = df,
  effect = 'time',
  model = 'random'
) #NOTE: PLM TSLS did not run with twoway (error: twoway not supported)

### jive ----
# from estimators script
jive.est(y,X,Z,SE=T)


### rjive ----



### Post Lasso ----



### Tables ----
#### fixest ----
etable(ols_est,tsls_est,
       tex = T,
       title = 'Main Results',
       fitstat = 'f',
       style.tex = style.tex('aer')
       )

#### plm ----
stargazer(plm_ols,plm_tsls,
          type='latex',digits=3, column.sep.width = "5pt",
          title='PLM Main Results')


### graphs ----

# Dingel Plot (looks bad)
ggplot(data = df, 
       mapping = aes(
         x = aftexpl + aftexpl.dist_to_ref,
         y = lprice),
       ) +
  geom_point()
  


