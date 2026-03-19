library(ggplot2)
library(patchwork)
library(dplyr)

scenDirs = list.dirs( 'simulations', recursive = FALSE)
nrep=200
for (scenDir in scenDirs)  {
  
  print(paste0("Working on ",scenDir))
  
  time = read.csv(paste0(scenDir,'/timeline.csv'))

  rnbData = tibble(read.csv(paste0(scenDir,'/','analyzedTOF.csv')))
  scenLetter = substr(scenDir,nchar(scenDir),nchar(scenDir))

  isEvent = rnbData |> 
    filter( (isRNB == TRUE) | (isDip == TRUE)) |> 
    select(sim.id)
  
  # n of subjects to be shown not exceed number of subjects, max 200
  maxRep = dim(rnbData)[1] - dim(isEvent)[1]
  nRep   = ifelse(maxRep<=200, maxRep, 200)
  
  notEvent = sample(rnbData$sim.id[!(rnbData$sim.id %in% isEvent$sim.id)], nRep)
  
  ids2plot = c(isEvent$sim.id, notEvent)

   # Plot
  tempData = tibble(cbind(time,
                          read.csv(paste0(scenDir,'/C1ro.csv')))) |>
    select(time,all_of(ids2plot))
  tempData = tidyr::pivot_longer(tempData, 
                                 cols=colnames(tempData)[colnames(tempData) != 
                                                           'time'] )

  g_c1ro <- ggplot(data=tempData) + 
    geom_line(aes(x = time, y = log10(value), group = name), 
              color='black', alpha=.3, linewidth=.3) + 
    #scale_y_log10(guide = "axis_logticks") +
    ylab('log10(C1ro), microM') + xlab('Time, minutes')
  
  tempData = tibble(cbind(time,
                          read.csv(paste0(scenDir,'/C2ro.csv')))) |>
    select(time,all_of(ids2plot))
  tempData = tidyr::pivot_longer(tempData, 
                                 cols=colnames(tempData)[colnames(tempData) != 
                                                           'time'] )
  g_c2ro <- ggplot(data=tempData) + 
    geom_line(aes(x = time, y = log10(value), group = name), 
              color='black', alpha=.3, linewidth=.3) + 
    #scale_y_log10(guide = "axis_logticks") +
    ylab('log10(C2ro), microM') + xlab('Time, minutes')
  
  tempData = tibble(cbind(time,
                          read.csv(paste0(scenDir,'/Ceff_ro.csv')))) |>
    select(time,all_of(ids2plot))
  tempData = tidyr::pivot_longer(tempData, 
                                 cols=colnames(tempData)[colnames(tempData) != 
                                                           'time'] )
  g_ceff <- ggplot(data=tempData) + 
    geom_line(aes(x = time, y = log10(value), group = name), 
              color='black', alpha=.3, linewidth=.3) + 
    #scale_y_log10(guide = "axis_logticks") +
    ylab('log10(Ceff), microM') + xlab('Time, minutes')
  
  tempData = tibble(cbind(time,
                          read.csv(paste0(scenDir,'/TOFr.csv')))) |>
    select(time,all_of(ids2plot))
  tempData = tidyr::pivot_longer(tempData, 
                                 cols=colnames(tempData)[colnames(tempData) != 
                                                           'time'] )
  g_tof <- ggplot(data=tempData) + 
    geom_line(aes(x = time, y = value, group = name), 
              color='black', alpha=.3, linewidth=.3) + 
    geom_hline(yintercept=90, linetype='dashed', color='darkgreen',linewidth=.5) +
    ylab('TOFr, %') + xlab('Time, minutes')
  
  toPlot = (g_c1ro +  g_c2ro +g_ceff) / g_tof
  ggsave(filename=paste0(scenDir,'/TOFplot_',scenLetter,'.png'), 
         plot=toPlot,
	 width = 10,
	 height= 7,
	 units= "in",
	 dpi=1000)
  
  rm(rnbData, tempData, toPlot, g_c1ro, g_ceff, g_tof)

}
