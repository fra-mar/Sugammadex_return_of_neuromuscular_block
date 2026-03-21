# Dataset creator forlater brms analysis
library('dplyr')

roc2microM = 1.887505
sug2microM = .499451
if (exists('adf')){rm(adf)}
if (exists('params')){rm(params)}

# Check if data exists, otherwise create it
dataSourceDir = [yourPROJECTdirectory]

if (!('mydata.csv' %in% list.files())) {
  
  print('Datafile doesnt exist. Creating it')
  targetScen = c('G', 'H', 'I' , 'J', 'K', 'L') # Scenarios to be joined and analyzed
  
  setwd(dataSourceDir)
  scenarios = tibble(read.csv('scenarios.csv'))
  
  for (letter in targetScen){
    scenDsug = scenarios[scenarios$sceName==letter,]$d0SUG
    scenDsug = ifelse(scenDsug=='var','adjusted','fixed') # change for easier factoring
    
    # Read tofAnalyzer
    tdf_tofAn  = tibble(
      read.csv(
        paste0('./simulations/scenario_',letter,'/analyzedTOF.csv'))) |>
      select(isRNB, t90, tIniRNB, durRNB, lowestTOF ) |>
      mutate(sug=scenDsug) 
    
    # Read forRegression
    tdf_forReg = tibble(
      read.csv(
        paste0('./simulations/scenario_',letter,'/forRegression.csv')))
    tdf = cbind(tdf_tofAn, tdf_forReg)
    
    # Read params
    tdf_params  = tibble(
      read.csv(
        paste0('./simulations/scenario_',letter,'/params.csv'))) |>
      select(!c('simId'))
    #tdf_params = tibble(
    #  cbind(tdf_params, tdf_tofAn |> select(isRNB, sug)))
    
    if (!(exists('adf'))){
      adf = tdf
    } else {
      adf = tibble(rbind(adf, tdf))
    }
    
    if (!(exists('params'))){
      params = tdf_params
    } else {
      params = tibble(rbind(params, tdf_params))
    }
  }

  
  setwd('[yourPROJECTdirectory]/statAnalysis')
  
  mydata = tibble(cbind(tibble(adf), params)) |> 
    mutate( #rSUGamt is rate sugammadex dose to remaining roc in A1 and A2
            #calculated as molar ratio, that should be close to 1.
            rSUGamt = (79 * dSUG * sug2microM)/((A1ro + A2ro)), #79kg scenG-L
            A1ro = A1ro / roc2microM, # was imported as microM, needed mg so reader understands
            A2ro = A2ro / roc2microM, # was imported as microM, needed mg so reader understands
	    sug  = factor(sug, levels=c('adjusted','fixed')),
            dROC = factor(dROC),
            dSUG = factor(round(dSUG,2))) 
  
  write.csv(mydata, 'mydata.csv', row.names = FALSE)
} else {
  mydata = tibble(read.csv('mydata.csv')) |> 
	    mutate( 
            sug  = factor(sug, levels=c('adjusted','fixed')),
            dROC = factor(dROC),
            dSUG = factor(round(dSUG,2))) 
}  

#--------------Support functions------------------
qtMe = function(x, probs = c(0.05, 0.5, 0.95)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = c('p5','p50','p95')
  )
}

mySummary = function(df) {
  summarized = df |> summarise(across(where(is.numeric), .fns = 
                          list(median = ~quantile(.,.5),
                               p5 = ~quantile(., 0.05),
                               p95 = ~quantile(., 0.95)))) |>
  tidyr::pivot_longer(everything(), 
                      names_sep='_', 
                      names_to=c('variable', '.value'))
  
  return (summarized)
}
mySummary2 = function(df) {
  summarized = df |> summarise(across(where(is.numeric), .fns = 
                                        list(median = ~quantile(.,.5),
                                             p5 = ~quantile(., 0.05),
                                             p95 = ~quantile(., 0.95)))) |>
    tidyr::pivot_longer(everything())
  
  return (summarized)
}

