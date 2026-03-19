# This script run regression models in brms and save output in folder stats_outputs

library('dplyr')
library('ggplot2')
library('patchwork')
library('brms')
library('tidybayes')
library('bayesplot')
options(mc.cores = parallel::detectCores())

source('dataGatherer.R')

# Scale pkpd params and rSUGamt
data_std = mydata |>  mutate(across(c(CLro:ks, rSUGamt), ~c(scale(.))))

# Uninformative priors for all models
priors = c(prior(student_t(3,0,2.5),class='Intercept'),
           prior(student_t(3,0,2.5),class='b'))

print(Sys.time())

#----1. Depends RNB on how patient is (pkpdParams) and what you've done(dROC,reversalStrategy)?------------

# Standardize pkpd params

fit_pkpdROCSUG = brm(bf( isRNB ~ 1 + CLro + V1ro + V2ro + CLsu + ke0 + EC50 + Hill + ks + dROC + sug), 
   data = data_std, 
   family= bernoulli(),  
   prior = priors,
   warmup=2000, chains=4, iter=7000,
   seed = 123,
   file = "./stats_outputs/fit_pkpdROCSUG.RDS",
   file_refit = "on_change")

print(tidybayes::summarise_draws(fit_pkpdROCSUG))
print('Ended fitting RNB~ pkpd + roc + sug')
print(Sys.time())


#---2.0  Depends RNB on JUST  rate dSUG/remaining roc?--------------------
# I exclude pkpd params due to colinearity with rSUGamt

fit_rSUGamt = brm(bf( isRNB ~ 1 + rSUGamt), 
   data = data_std, 
   family= bernoulli(),  
   prior = priors,
   warmup=2000, chains=4, iter=7000,
   seed = 123,
   file = "./stats_outputs/fit_rSUGamt.RDS",
   file_refit = "on_change")

print(tidybayes::summarise_draws(fit_rSUGamt))
print('Ended fitting RNB~ 1 + rSUGamt')
print(Sys.time())

# 2.1 Is rSUGamt different between isRNBevent?

fit_rSUGamtG = brm(bf(rSUGamt ~ 1 + isRNB ),
			data = mydata,
			family = lognormal(),
                        prior = priors,
			warmup = 2000, chains=4, iter=7000,
			seed = 123,
			file = "./stats_outputs/fit_rSUGamtG",
			file_refit = "on_change")

print(tidybayes::summarise_draws(fit_rSUGamt))
print('Ended fitting rSUGamt~ 1 + isRNB')
print(Sys.time())

# 3. --Depends amount roc in peripheral cmt on total ROC or administered sug?---------

fit_A2ro = brm(bf( A2ro ~ 1 + dROC + isRNB), 
                 data = data_std |> filter(sug=='adjusted'), 
                 family= lognormal(), 
                 prior = priors,
                 warmup=2000, chains=4, iter=7000,
                 seed = 123,
		 file =  "./stats_outputs/fit_A2ro.RDS",
   		 file_refit = "on_change")

print(tidybayes::summarise_draws(fit_A2ro))
print('Ended fitting A2ro ~ 1 + dROC + isRNB')
print(Sys.time())

# 4.---Depends amount roc in central cmt on total ROC or administered sug?---------

fit_A1ro = brm(bf( A1ro ~ 1 + dROC + isRNB), 
                 data = data_std |> filter(sug=='adjusted'), 
                 family= lognormal(), 
                 prior = priors,
                 warmup=2000, chains=4, iter=7000,
                 seed = 123,
		 file =  "./stats_outputs/fit_A1ro.RDS",
		 file_refit = "on_change")

print(tidybayes::summarise_draws(fit_A1ro))
print('Ended fitting A1ro ~ 1 + dROC + isRNB')
print(Sys.time())

# 5. Run same model pattern for t90 and characteristics of RNB

vars = c('t90','tIniRNB', 'durRNB', 'lowestTOF')

  for (v in vars){
    
    formula = paste0(v, ' ~ 1 + dROC + sug')
    
    if (v=='t90') {
      fdata = mydata |> filter(t90 >= 0)
    } else {fdata = mydata |> filter(isRNB == TRUE)}
    
    f = brm(formula = formula,
            data = fdata,
            family = 'lognormal',
            prior = priors,
            warmup = 2000,
            iter = 7000,
            chains = 4,
            seed = 123,
            file = paste0("./stats_outputs/fit_",v,".RDS"),
            file_refit = "on_change")
    
    asOrder = c('a_lowROCvarSUG', 
                       'b_dROC284.4',
                       'b_dROC379.2',
                       'b_sugfixed')
    asOrder = asOrder[seq(length(asOrder),1)]
    
    transf = f |> 
      tidy_draws() |> 
      select(starts_with("b_")) |> 
      mutate(across(c(2:4), ~exp(.x + b_Intercept)), 
             a_lowROCvarSUG = exp(b_Intercept)) |> 
      select(!b_Intercept)
    forPlot = transf |> tidyr::pivot_longer(cols = everything()) |>
      mutate(name = factor(name,levels= asOrder,
                           labels = substr(asOrder,3, 100)))
    
    fplot = ggplot(forPlot, aes(x=value, y=name)) + 
      tidybayes::stat_pointinterval() + 
      xlab('Estimate, minutes -Not the whole sample!-') +
      ggtitle(paste0(v, '~ 1 + dROC+ sug'), subtitle = 'exp(intercept + var)')
    
    ggsave(paste0('./stats_outputs/plt_', v,'.png'))

  }  
