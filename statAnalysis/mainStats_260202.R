# Main statistics

library('dplyr')
library('ggplot2')
library('patchwork')
library('brms')
library('tidybayes')
library('bayesplot')

source('dataGatherer.R')
# Support functions-----------------------------------------
transformMe = function(afit){
  transformed = afit |> tidy_draws() |> select(starts_with('b_')) |> 
    mutate(isRNB = exp(b_Intercept + b_isRNBTRUE),
           noRNB = exp(b_Intercept), diff = isRNB - noRNB) |>
    select(isRNB, noRNB, diff)
  
  short = summarise_draws(transformed) |> select(variable, median, q5,q95)
  return (short)
}
pretty = function(a,b,c) {                # function for median(q5-q95)
  paste0(as.character(round(a,2)),
         ' (',as.character(round(b,2)),
         '-', as.character(round(c,2)), ')') 
}

# Pooled MSP cases

# 1 How many RNB events ------------------------------------
print('n(%) events in MSP pooled cases')
nRNBcases = mydata |> group_by(isRNB) |> summarise(n=n()) |> 
  mutate(p = n/dim(mydata)[1])
print(nRNBcases)

# 2 Characteristics of RNB events --------------------------
targets = c('tIniRNB','durRNB', 'lowestTOF')
smz = summarise_draws(mydata |> 
                  filter(isRNB==TRUE) |> 
                  select(all_of(targets))) |> 
  select(median, q5, q95)
RNBeventsChars = cbind(targets,smz)
print(RNBeventsChars)

# 3. OR and OR-puc for the main model--------------------------

fit = readRDS("./stats_outputs/fit_pkpdROCSUG.RDS")

#Order in which params will be shown. And set labels
asOrder =  c('b_sugfixed','b_dROC379.2','b_dROC284.4', 
             'b_CLro', 'b_V1ro','b_V2ro', 
             'b_CLsu','b_ks','b_ke0','b_EC50', 'b_Hill'    )
pkpdROCSUG_summary = summarise_draws(fit) |> filter(variable %in% asOrder)
# OBS, in the logit space / coeffs as log(odds)
print(pkpdROCSUG_summary)

#Reverse order for easier visualization
srev =seq(length(asOrder),1)
asOrder = asOrder[srev]
#Make labels easier to read
varsLabs = substr(asOrder,3,1000)

transformed = fit |> tidy_draws() |>
  select(starts_with("b_") & !starts_with("b_I")) |> 
  mutate(across(where(is.numeric),~exp(.x)))  # Exp for log(odds) -> odds
forPlot = transformed |>
  tidyr::pivot_longer(cols=everything()) |>
  mutate(name = factor(name,levels=asOrder, labels = varsLabs),
         isFactor = ifelse(name%in%c('sugfixed','dROC379.2','dROC284.4'),
                           'isF','noF'))

g1 = ggplot(data = forPlot |> filter(isFactor=='isF'), 
            aes(x=value, y=name)) + 
  geom_vline(aes(xintercept=1), color='red') +
  tidybayes::stat_pointinterval(.width=c(.66,.90), shape = 'o', 
                                color='orange',
                                point_color='black',
                                size=10,
                                linewidth = 10 ) +
  xlab('Odds Ratio, OR') + ylab('')+ 
  theme_minimal() + theme(legend.position='none')


g2 = ggplot(data = forPlot |> filter(isFactor!='isF'), 
            aes(x=value, y=name)) + 
  geom_vline(aes(xintercept=1), color='red') +
  tidybayes::stat_pointinterval(.width=c(.66,.90), shape = 'o', 
                                color='darkolivegreen4',
                                point_color='black',
                                size=10,
                                linewidth = 10 ) +
  xlab('Odds Ratio per unit-change, OR-puc') + ylab('')+ 
  theme_minimal() + theme(legend.position='none')

gplot_pkpdROCSUG = g1/g2 + plot_layout(heights= c(1,1))

plot(gplot_pkpdROCSUG)
ggsave('./stats_outputs/pkpdROCSUG_plot.jpeg', 
       gplot_pkpdROCSUG, 
       units='in', width=7, height=6, 
       bg='white')


# get a tidy table
sumMain = summarise_draws(fit)

sumMain_short = summarise_draws(transformed) |>
  select(variable, median, q5, q95, rhat, ess_bulk, ess_tail) |>
  mutate(variable = substr(variable,3,100),
         rhat = substr(as.character(rhat),1,4),
         OR = pretty(median, q5, q95)) |>
  select(all_of(c('variable','OR','rhat','ess_bulk','ess_tail')))

print(sumMain_short)

# 4. Depends the rate SUG/remaining ROC on occurrence of RNB-?-----------------

# Descriptive
rate_describe = mydata |> 
  group_by(isRNB) |> summarize(qtMe(rSUGamt))
print(rate_describe)


f_rate = readRDS('./stats_outputs/fit_rSUGamtG.rds')
summaryRate = summarise_draws(f_rate)

summaryRate_short = transformMe(f_rate)
print(summaryRate_short)

gRate = ggplot(data = mydata , aes(x=rSUGamt, fill=isRNB)) + 
  geom_density(alpha=.4) + ggtitle('ratio SUG / remaning ROC at reversal time')
plot(gRate)

# 5. A1ro, diff between yes no RNB events-------------------------------------

fa1 = readRDS('./stats_outputs/fit_A1ro.RDS')
summaryfA1 = summarise_draws(fa1)


A1_describe = mydata |> filter (sug=='adjusted') |> 
  group_by(isRNB) |> summarize(qtMe(A1ro))

summaryA1_short = transformMe(fa1)
print(summaryA1_short)

gA1 = ggplot(data = mydata |> filter(sug=='adjusted') , aes(x=A1ro, fill=isRNB)) + 
  geom_density(alpha=.4) + xlab('Rocuronium -mg- central compartment')
#plot(gA1) 

# 6. A2ro, diff between yes no RNB events-------------------------------------

fa2 = readRDS('./stats_outputs/fit_A2ro.RDS')
summaryfA2 = summarise_draws(fa2)


A2_describe = mydata |> filter (sug=='adjusted') |> 
  group_by(isRNB) |> summarize(qtMe(A2ro))

summaryA2_short = transformMe(fa2)
print(summaryA2_short)

gA2 = ggplot(data = mydata |> filter(sug=='adjusted') , aes(x=A2ro, fill=isRNB)) + 
  geom_density(alpha=.4) + xlab('Rocuronium -mg- peripheral compartment')
#plot(gA2) 

gAroc = gA1 + gA2
plot(gAroc)

wrapup = list(nRNB = nRNBcases,
              howRNBwere = RNBeventsChars,
              mainModel_summary = sumMain_short,
              plot_mainModel = gplot_pkpdROCSUG,
              rSUGamt_describe = rate_describe,
              rSUGamt_model_short = summaryRate,
              rSUGamt_model_plot = gRate,
              fAro1_summary = summaryA1_short,
              Aro1_describe = A1_describe,
              fAro2_summary = summaryA2_short,
              Aro2_describe = A2_describe,
              Aro_plot = gAroc)

saveRDS(wrapup, './stats_outputs/mainStats.RDS')
