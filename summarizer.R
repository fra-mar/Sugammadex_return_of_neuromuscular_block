# Goal1: Summarize results regarding TOF and pkpd parameters, for each scenario
rm(list=ls())
library(dplyr)
library(ggplot2)
#library(patchwork)

# Supporting functions
quantileMe = function(sr, srName){
  toReturn = quantile(sr, c(0.05,.5,.95))
  names(toReturn) = paste0(srName,'_',c('p5','p50','p95'))
  return (toReturn)}
  
whichEvent = function(spntRecover, isRNB, isDIP, didRecover){
    if (spntRecover==TRUE){event = 1}
    else if (isRNB==TRUE){event=2}
    else if (isDIP==TRUE){event=3}
    else if (didRecover==FALSE){event=4}
    else {event=0} 
  return (event)
}

wilcoxMe = function(df,eventType) {
  ar = c()
  for (rn in colnames(df)[1:8]){
    tempdf = df[,c(rn,eventType)]
    colnames(tempdf) = c('Value','EventType')
    wt = wilcox.test( Value ~ EventType, data=tempdf)
    ar = c(ar, c(wt$p.value)) }

  toReturn = data.frame(array(ar,dim=c(1,length(ar))))
  colnames(toReturn) = paste0('P_',colnames(df)[1:8])
  return (toReturn)
}
# Functions for human readable data
npMe = function(n, p){
  toRet = paste0(as.character(round(n,3)),
                 ' (',
                 as.character(round(p*100,1)),
                 '%)')
  return (toRet)
}

qMe = function(p5,p50,p95){
  if (any(is.na(c(p5,p50,p95)))) {
    toRet = '-'
  }
  else{
    toRet = paste0(as.character(round(p50,1)),
                   ' (',
                   as.character(round(p5,1)),
                   ',',
                   as.character(round(p95,1)),
                   ')')
  }
  return(toRet)
}
# For these params  variability was reported.
targetPars = c("simId", "CLro" ,"V1ro", "V2ro", "CLsu", 
               "ke0", "EC50", "Hill", "ks")
# To visuzalize in last table scenarios with no RNB events
scenWithNoEvents = c()

scenDirs = list.dirs('simulations', recursive = FALSE)
outputDir = 'outputs'
dir.create(outputDir)

#_________________Main Loop___________________________________
for (scenName in scenDirs) {
  print(paste0("Working with ",scenName))
  tdata_raw = tibble(read.csv(paste0(scenName,'/analyzedTOF.csv')))
  nSub = length(tdata_raw$sim.id)
  
  # Who recovered sponatnously
  spntRec = tibble(
    nSpntRec = dim(tdata_raw |> filter(spntRecover==TRUE))[1],
    pSpntRec = nSpntRec / nSub)
  
  
  # exclude thse who didnt and actually never recovered
  tdata = tdata_raw |> 
    filter((spntRecover==FALSE) & (didRecover==TRUE))
  
  t90_q = quantileMe(tdata$t90,"t90") 
  
  RNB = tdata_raw |> filter((isRNB==TRUE) & (isDip==FALSE))
  npRNB = tibble(
    nRNB = dim(RNB)[1],
    pRNB = nRNB / dim(tdata)[1])
  
  initRNB_q = quantileMe(RNB$tIniRNB, "initRNB")
  durRNB_q = quantileMe(RNB$durRNB, "durRNB")
  lowestRNB_q = quantileMe(RNB$lowestTOF, "lowestRNB")
  
  DIP = tdata_raw |> filter((isDip==TRUE) & (spntRecover==FALSE))
  npDIP = tibble(
    nDIP = dim(DIP)[1],
    pDIP = nDIP / dim(tdata)[1] )
  
  initDIP_q = quantileMe(DIP$tDip, "initDIP")
  #lowestDIP_q = quantileMe(DIP$lowestTOF, "lowestDIP")
  
  scenLetter = substr(scenName,nchar(scenName),nchar(scenName))
  sumup = tibble(scen=scenLetter,
                     cbind(spntRec, t(t90_q),
                           npRNB, t(initRNB_q), t(durRNB_q), t(lowestRNB_q),
                           npDIP, t(initDIP_q)))
  if (!exists("scenSumup")) {
    scenSumup = sumup }
  
  else {
    scenSumup = rbind(scenSumup, sumup)
  }
  
  #----------------PKPD params-----------------------------------
  #Build a table with percentiles 0,50,95 for parameter and scenario
  params = tibble(read.csv(paste0(scenName,'/params.csv'))) |>
    select(all_of(targetPars))
  # Sort events in 0(normal) 1(spntRec) 2(RNB) 3(dip) 4(neverRec)
  idEvents = tdata_raw |> 
    rowwise() |> 
    mutate(whichEv = whichEvent(spntRecover,isRNB,isDip,didRecover)) |>
    select(whichEv)
  
  # Add as new column to params
  params = params |>
    mutate(whichEv = idEvents$whichEv)
  pars_grp = params|> group_by(whichEv) |> 
    summarise(across(targetPars[2:9],
                     list(p5=~quantile(.,.05, na.rm=TRUE),
                          p50=~quantile(.,.5, na.rm=TRUE),
                          p95=~quantile(.,.95, na.rm=TRUE)
                          )))
  letterCol = tibble(scen=rep(scenLetter,dim(pars_grp)[1]))
  pars_grp = tibble(cbind(letterCol, pars_grp))
  
  if (!exists("parsSumup")) {
    parsSumup = pars_grp }
  
  else {
    parsSumup = rbind(parsSumup, pars_grp)
  }
  
  # compare pkpd between RNB and noRNB 
  # Run wilcoxon tests
  params_sel = params |> filter(whichEv %in% c(0,2))
  
  
  if (any(params$whichEv==2)){
    w=wilcoxMe(params_sel |> select(colnames(params_sel)[2:10]), 
             'whichEv')
    is2 = pars_grp[pars_grp$whichEv==2,3:26]
    is0 = pars_grp[pars_grp$whichEv==0,3:26]
    dpars = tibble((is2 - is0)/is0) |>
      select (colnames(pars_grp)[grep('_p50',colnames(pars_grp))])
    diffpars = tibble(cbind(scenLetter,dpars,w)) 
    
    if (!exists("parsCompare")) {
      parsCompare = diffpars }
    else {
      parsCompare = rbind(parsCompare, diffpars)
    }
  }
  else{
    scenWithNoEvents = c(scenWithNoEvents, c(scenLetter))
  }
  
  #---------------Plot pkpdPars------------------------------
 
  params_sel = params_sel |> 
    mutate(whichEv = factor(whichEv,
                            levels=c(0,2),
                            labels=c('noRNB','isRNB')))
  
  pivoted = tidyr::pivot_longer(params_sel ,
                                cols=names(params_sel)[2:9])
  pivoted = pivoted |>
    mutate(nameF = factor(name, 
                          levels=name, 
                          labels=name)) 
   
  gg_recur = ggplot(pivoted, aes(x=value, fill=whichEv)) + 
    geom_density(alpha=.5) +
    facet_wrap(.~nameF, scales='free', nrow=2)+ ggtitle(scenName) +
    labs(x='', y = 'Density') +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank()) 
  ggsave(filename=paste0(scenName,'/pkpd_plot_',scenLetter,'.png'),
         plot=gg_recur,
         width = 12, height = 6,
         units = "in",
         dpi=1000)
  

 
}
write.csv(scenSumup, 	paste0(outputDir,'/TOFsumup.csv'), 	
          row.names = FALSE)
write.csv(parsSumup, 	paste0(outputDir,'/PARSsumup.csv'), 	
          row.names = FALSE)
write.csv(parsCompare, 	paste0(outputDir,'/PARScompare.csv'), 	
          row.names = FALSE)

#_________________End of Main Loop_____________________

# Create a human readable TOF data
TOFprettyPrinted = scenSumup |>
  rowwise() |>
  mutate(spntRecovered = npMe( nSpntRec,  pSpntRec),
         t90           = qMe( t90_p5, t90_p50, t90_p95),
         nRNB          = npMe( nRNB, pRNB),
         initRNB       = qMe( initRNB_p5, initRNB_p50, initRNB_p95),
         durRNB        = qMe( durRNB_p5, durRNB_p50, durRNB_p95),
         lowestRNB     = qMe( lowestRNB_p5, lowestRNB_p50, lowestRNB_p95),
         nDIP          = npMe( nDIP,  pDIP),
         initDIP      = qMe( initDIP_p5, initDIP_p50, initDIP_p95)) |>
  select(c(scen, spntRecovered,
           t90,
           nRNB,
           initRNB,
           durRNB,
           lowestRNB,
           nDIP,
           initDIP))

write.csv(TOFprettyPrinted, paste0(outputDir,'/TOFsumup_h.csv'), 
          row.names = FALSE)



