library(ggplot2)
library(reshape2)
library(tools)
library(wesanderson)
library(plyr)
library(gridExtra)
library(ggpubr)
library(magrittr)
library(sjPlot)
library(stargazer)
library(grid)

# VARIABLES SPECIFIED BY USER ---------------------------------------------

peakThreshold1 = 0.070          # Threshold for classifying trial as CR
peakThreshold2 = 0.050           # Onset Threshold, or trials already classified as CRs
peakThreshold1_startle = 0.05   # Threshold for classifying trial as Startle
baseline_noise = 0.025
delay = 20                      # moves the analysis window to delay ms after onset of CS before looking for CRs
analysisWindow = 200            # time in ms after startle tone onset that is analyzed for startle
halfW = 0                       # Half-width used by smoothing kernel in analysis
halfW_plotting = 3                      # Half-width used by smoothing kernel in plotting
file_substring_end = 15         # Length of file name (minus the '.txt' part)
begin_paired = c(10,50)       # The next four lines are vectors specifying data used in onset/peak bootstrapping
end_paired = c(60,100)
begin_probe = c(1,20)
end_probe = c(25,55)
numberOfIterations = 5000       # Number of iterations used in bootstrapping onset/peak values
min_test = 2
sliding_window_size = 40

color1 = wes_palette("FantasticFox1")
color2 = wes_palette("Royal1")
color3 = wes_palette("BottleRocket2")
color4 = wes_palette("Zissou1")
color5 = wes_palette("Chevalier1")

baseline_color = color5[4]
startle_color = color1[3]


# Optional Variables if running without shiny APP -----------------------

# dataFolder = '2019-03-SEFL'
# doWaterfall = 1
# cue_color = color1[1]
# baselineDur = 200
# trialDur = 1200
# onset_histogram = 'black'
# airpuff_color = '#D3DDDC' # Gray
# peak_histogram = 'black'
# group1_color = 'green'
# group2_color = 'black'
# include_session_waterfall = 1
# postscript_flag = 0

# CREATE PATHS AND FOLDERS  ----------------------------------------------

path_eyelid_behavior = '/Users/Cheasequah/Desktop/Eyelid_Behavior/'
path_waterfall_raw = paste(path_eyelid_behavior, 'waterfall_raw/', dataFolder, sep = '')
path_eyelid_analysis = paste(path_eyelid_behavior, 'eyelid_analysis/', dataFolder, sep = '')
path_session_dataframes = paste(path_eyelid_analysis, '/session_dataframes/', sep = '')
path_session_dataframes_startle = paste(path_eyelid_analysis, '/session_dataframes_startle/', sep = '')
path_waterfall_postcript = paste(path_eyelid_analysis, '/waterfall_poscript', sep = '')
path_waterfall_png = paste(path_eyelid_analysis, '/waterfall_png', sep = '')
path_waterfall_startle_postcript = paste(path_eyelid_analysis, '/waterfall_startle_postcript', sep='')
path_waterfall_startle_png = paste(path_eyelid_analysis, '/waterfall_startle_png', sep='')
path_trajectories = paste(path_eyelid_analysis, '/trajectories', sep = '')
path_trajectories_startle = paste(path_eyelid_analysis, '/trajectories_startle', sep = '')
path_summary_csv = paste(path_eyelid_analysis, '/summary_csv/', sep = '')
path_plots = paste(path_eyelid_analysis, '/plots/', sep = '')

setwd(path_waterfall_raw)
fileList <- list.files()
if (!file.exists(path_eyelid_analysis)){
  dir.create(path_eyelid_analysis) }
if (!file.exists(path_session_dataframes)){
  dir.create(path_session_dataframes) }
if (!file.exists(path_session_dataframes_startle) && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 ){
  dir.create(path_session_dataframes_startle) }
if (!file.exists(path_waterfall_postcript) && include_session_waterfall == 1 && postscript_flag == 1){
  dir.create(path_waterfall_postcript) }
if (!file.exists(path_waterfall_png) && include_session_waterfall == 1){
  dir.create(path_waterfall_png) }
if (!file.exists(path_waterfall_startle_postcript) && include_session_waterfall == 1 && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 && postscript_flag == 1){
  dir.create(path_waterfall_startle_postcript) }
if (!file.exists(path_waterfall_startle_png) && include_session_waterfall == 1 && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 ){
  dir.create(path_waterfall_startle_png) }
if (!file.exists(path_trajectories)  && include_session_waterfall == 1){
  dir.create(path_trajectories) }
if (!file.exists(path_trajectories_startle)  && include_session_waterfall == 1 && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 ){
  dir.create(path_trajectories_startle) }
if (!file.exists(path_summary_csv)){
  dir.create(path_summary_csv) }
if (!file.exists(path_plots)){
  dir.create(path_plots) }

# DEFINE FUNCTIONS --------------------------------------------------------
readTimeIntervals = function(f) {
  lines = readLines(f)
  waveStarts = sort(grep('^WAVES/', lines))
  minStart = min(waveStarts)
  waveEnds = sort(grep('^END', lines))[1:length(waveStarts)]
  waves = list()
  for (t in 1:length(waveStarts)) {
    waves[[t]] = do.call(rbind, strsplit(
      sub('^\\s+', '', lines[(waveStarts[t]+2):(waveEnds[t]-1)]),
      '\\s+'
    ))
    mode(waves[[t]]) = 'numeric'
    waves[[t]] = apply(
      X = waves[[t]],
      MARGIN = 1,
      FUN = function(c) {rgb(c[1], c[2], c[3], maxColorValue=255)}
    )
  }
  return(waves)
}

readColorMapping = function(f) {
  lines = readLines(f)
  lines = lines[grepl('ModifyGraph zColor', lines)]
  lines = unlist(strsplit(lines, split='},zColor', fixed=TRUE))
  colorMap = as.integer(gsub('^.*\\{colorwave([1234]),.*', '\\1', lines))
  names(colorMap) = gsub('^.*\\((.*?)\\).*', '\\1', lines)
  return(colorMap)
}

readWaterFall = function(f) {
  lines = readLines(f)
  trajectories = list()
  waveStarts = sort(grep('^WAVES\\t', lines)) 
  minStart = min(waveStarts)
  waveEnds = sort(grep('^END', lines))
  waveEnds = waveEnds[waveEnds > minStart]
  waveNames = gsub('^WAVES\\t', '', grep('^WAVES\\t', lines, value=TRUE))
  trajectories = list()
  for (t in 1:length(waveStarts)) {
    trajectories[[ waveNames[[t]] ]] =
      as.numeric(gsub('\\s+', '', lines[(waveStarts[t]+2):(waveEnds[t]-1)]))
  }
  return(trajectories)
}

baselines = function(x, maxt=baselineDur) {
  colMeans(x[x$t <= maxt, setdiff(colnames(x), 't')])
}

subtractBaselines = function(x, maxt=baselineDur) {
  origOrder = colnames(x)
  x = x[ , c(setdiff(colnames(x), 't'), 't')]
  x = sweep(x, 2, c(baselines(x, maxt), 0), `-`)
  return(x[ , origOrder])
}

colSds = function(x, na.rm=FALSE) {
  n = nrow(x)
  return(sqrt((n/(n-1)) * (
    colMeans(x*x, na.rm=na.rm) - colMeans(x, na.rm=na.rm)^2)))
}

basenoise = function(x, maxt=200) {
  ## 0.025-0.05 appears to be a reasonable range for filter value...
  colSds(x[x$t <= maxt, setdiff(colnames(x), 't')])
}

smoothTrajectory = function(traj, halfWidth = halfW) {
  kernapply(traj, k=kernel('daniell', halfWidth))
}

allPeaks = function(traj, halfWidth = halfW, minPeak = 0.15) {
  if (is.matrix(traj) || is.data.frame(traj)) {
    out = lapply(
      X = data.frame(traj[ , setdiff(colnames(traj), 't')]),
      FUN = allPeaks,
      halfWidth = halfWidth
    )
    out$t = NA
    return(out[colnames(traj)])
  }
  if (length(halfWidth) > 0) {
    traj = smoothTrajectory(traj, halfWidth)
  }
  d1 = diff(traj)
  d1pos = ifelse(d1>0, 1, 0)
  peaks = which(diff(d1pos) == -1)
  if (length(halfWidth) > 0) {
    peaks = peaks + halfWidth
  }
  peaks = peaks + 1
  if (length(minPeak) > 0) {
    peaks = peaks[traj[peaks] >= minPeak]
  }
  return(peaks)
}

topPeaks = function(traj, n=1, halfWidth = halfW) {
  if (is.matrix(traj) || is.data.frame(traj)) {
    out = data.frame(lapply(
      X = data.frame(traj[ , setdiff(colnames(traj), 't'), drop=FALSE]),
      FUN = topPeaks,
      n = n,
      halfWidth = halfWidth
    ))
    out[['t']] = NA
    return(out[colnames(traj)])
  }
  peaks = allPeaks(traj, halfWidth=halfWidth)
  peakHeights = sort(structure(traj[peaks], names=peaks), decreasing=TRUE)
  if (length(peakHeights) > 0) {
    out = sort(as.integer(names(peakHeights[1:n])))
  } else {
    out = rep(as.numeric(NA), n)
  }
  if (length(halfWidth) > 0) {
    unsmoothedPeaks = allPeaks(traj, halfWidth=NULL)
    for (i in 1:n) {
      if (!is.na(out[i])) {
        out[i] = unsmoothedPeaks[which.min(abs(unsmoothedPeaks-out[i]))]
      }
    }
  }
  return(out)
}

peakOnset1 = function(traj, peak, threshold1, threshold2) {
  if (all(traj < threshold1)) {return(NA)}
  if (is.na(peak)) {return(NA)}
  belowThreshold = which(traj[1:(peak-1)] < threshold2)
  return(max(belowThreshold) + 1)
}

maxAmp = function(traj, threshold=NULL) {
  if (is.vector(traj)) {
    out = max(traj)
  } else if (length(traj$t) == length(unique(traj$t))) {
    out = apply(traj[ , colnames(traj) != 't'], 2, max)
  } else {
    out = traj %>%
      group_by(trajectory) %>%
      summarize(max_amp = max(value))
    out = structure(out$max_amp, names=as.character(out$trajectory))
  }
  if (length(threshold) > 0) {out[out < threshold] = NA}
  return(out)	
}

integrateArea = function(traj, start, end, fps, threshold) {
  start = start*(fps/1000)
  end = end*(fps/1000)
  region= end-start
  samples = 1:round(region,0)
  out=0
  if (is.vector(traj)) {
    for(i in samples){
      if(traj[start+samples[i]]>threshold){
        out = out+traj[start+samples[i]]
      }
    }
  } else {return(NULL)}
  return(out)	
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# CREATE DATAFRAMES TO COMPILE ACROSS ALL SUBJECTS ------------------------

setwd(path_waterfall_raw)
fileList <- list.files()

# Counts Unique Animals in Experiment
x=c(1:(length(fileList)-1))
for (i in 1:(length(fileList)-1)) {
  x[i] = (substring(fileList[i],1,4) == substring(fileList[i+1],1,4))
}
numAnimals = sum(x==0)+1
numSessions = length(unique(substring(fileList,6,8)))

sessionNames = as.character(unique(substring(fileList,6,8)))
animalNames = unique(substring(fileList,1,4))

for(cueisi in unique(substring(fileList,10,13))){
  # Rename rows and columns in dataframe for paired trials
  summaryDF = data.frame(matrix(0, numAnimals, numSessions + 2))
  names(summaryDF)[1:2] = c('Animal', 'Group')
  summaryDF$Animal = animalNames
  summaryDF$Group = c(substring(animalNames,1,1))
  for (i in 3:(numSessions + 2)) {
    names(summaryDF)[i] <- sessionNames[i-2]
  }
  setwd(path_summary_csv)
  write.csv(summaryDF, paste(cueisi,'-Summary.csv', sep = '')) 
  setwd(path_waterfall_raw)
}

if(sum(as.numeric(substring(fileList,file_substring_end,file_substring_end)))>0){
  # Rename rows and columns in dataframe for startle trials
  startleDF = data.frame(matrix(0, numAnimals, numSessions + 4))
  names(startleDF)[1:2] = c('Animal', 'Group')
  startleDF$Animal = animalNames
  startleDF$Group = c(substring(animalNames,1,1))
  names(startleDF)[3] <- 'numTrials'
  names(startleDF)[4] <- 'countStartle'
  for (i in 5:(numSessions + 4)) {
    names(startleDF)[i] <- sessionNames[i-4]
  }
  setwd(path_summary_csv)
  write.csv(summaryDF, 'xxxx-StartleSummary.csv')  
}

# ANALYZE RAW FILES FOR CONDITIONED RESPONSES -----------------------------

setwd(path_waterfall_raw)
fileList <- list.files()

if(doWaterfall == 1) {
  for (j in 1:length(fileList)) {
    wf = readWaterFall(fileList[j])
    wf = data.frame(do.call(cbind, wf))
    wf$t = trialDur * (1:nrow(wf)) / nrow(wf) 
    
    # smooth a duplicated set of traces called "wf_plot"
    wf_plot = wf[((halfW_plotting+1):(length(wf[,1])-halfW_plotting)),]
    for(i in 1:(length(wf_plot[1,])-1)){
      wf_plot[,i] = smoothTrajectory(wf[,i],halfW_plotting)
    }
    wfsub_plot = subtractBaselines(wf_plot)
    wfsub_plot[is.na(wfsub_plot)] = 0
    wfsub_plot = wfsub_plot[ , c(basenoise(wfsub_plot) < baseline_noise, TRUE), drop=FALSE]
    
    wfsub = subtractBaselines(wf)
    wfsub[is.na(wfsub)] = 0
    wfsub = wfsub[ , c(basenoise(wfsub) < baseline_noise, TRUE), drop=FALSE]
    fps = length(wfsub[,1])/trialDur*1000
    USOnset = as.numeric(substring(fileList[j],11,13)) + baselineDur
    traceStart = round(((baselineDur + delay)*fps/1000),0) + 1    # moves the analysis window delay ms after CS onset. outputs the frame number where this happens
    traceEnd = round(((trialDur-50)*fps/1000),0)           # moves the end of analysis window to 50 before end of trial
    
    # Locate highest peak in each trajectory
    peakTiming_probe = topPeaks(wfsub[traceStart:(traceEnd),], n=1, halfWidth = halfW) + (traceStart-1)
    peakTiming_probe[is.na(peakTiming_probe)] = 0
    peakTiming_paired = topPeaks(wfsub[traceStart:(round((fps/1000)*(USOnset-50))),], n=1, halfWidth = halfW) + (traceStart-1)
    peakTiming_paired[is.na(peakTiming_paired)] = 0
    sessionDF = data.frame(Animal = substring(fileList[j],1,4), Group = substring(fileList[j],1,1), Session = substring(fileList[j],6,8),
                           Trial = as.numeric(substring(colnames(peakTiming_probe)[1:(length(peakTiming_probe)-1)],7,9)),
                           Paired = 0, CR = 0, Onset = 0, Peak = 0, UR = 0, UR_Integrate = 0, maxAmp = 0, CR_Area = 0, 
                           ID = substring(fileList[j],file_substring_end,file_substring_end),trialName=colnames(peakTiming_probe)[1:(length(peakTiming_probe)-1)], 
                           Cue = as.character(substring(fileList[j],10,10)), ISI = as.numeric(substring(fileList[j],11,13)))
    
    # Fill in data frame
    for(k in 1:length((sessionDF$Trial))) {
      if(sessionDF$Trial[k]%%5 == 0){
        if(peakTiming_probe[1,k] > 0){
          sessionDF$Peak[k] = topPeaks(wfsub[traceStart:(traceEnd),k], n=1, halfWidth = halfW)*(1000/fps) + delay
          sessionDF$CR[k] = 1 
          sessionDF$Onset[k] = peakOnset1(wfsub[ , k], peakTiming_probe[1, k], peakThreshold1, peakThreshold2)*(1000/fps) - baselineDur
          sessionDF$maxAmp[k] = maxAmp(wfsub[1:round(((fps/1000)*(USOnset-50)),0), k])
          sessionDF$CR_Area[k] = integrateArea(wfsub[,k], baselineDur, trialDur, fps , peakThreshold1)
        }}
      else{
        sessionDF$Paired[k] = 1
        sessionDF$UR[k] = maxAmp(wfsub[round((fps/1000)*(USOnset),0):traceEnd, k])
        sessionDF$UR_Integrate[k] = integrateArea(wfsub[,k], USOnset, trialDur, fps, peakThreshold1)
        if(peakTiming_paired[1, k] > 0){
          sessionDF$Onset[k] = peakOnset1(wfsub[1:round(((fps/1000)*(USOnset-50)),0), k], peakTiming_paired[1, k], 
                                          peakThreshold1, peakThreshold2)*(1000/fps) - baselineDur 
          sessionDF$CR[k] = 1
        }
      }}
    
    CR_Pct = round(sum(sessionDF$CR)/length(sessionDF$CR)*100,1)
    CR_Onst = round(median(sessionDF$Onset[sessionDF$Onset>0]),0)
    
    
    ## WATERFALL PLOTS
    
    if(include_session_waterfall == 1) {
      # Get time intervals from fileList[j] file
      ti = readTimeIntervals(fileList[j])
      names(ti[[1]]) = as.character(wf_plot$t)
      
      # Reshape wfsub for plotting
      wfmelt = wfsub
      wfmelt_plot = wfsub_plot
      wfmelt_plot$index = 1:nrow(wf_plot)
      wfmelt_plot = melt(wfmelt_plot, id.vars=c('t', 'index'), variable.name='trajectory')
      wfmelt_plot$'time color' = ti[[1]][as.character(wfmelt_plot$t)]
      
      # Fan trajectories out for waterfall plotting
      wfmelt_plot$shifted_value = wfmelt_plot$value + (as.integer(wfmelt_plot$trajectory)-1) * 0.1
      wfmelt_plot$shifted_t = wfmelt_plot$t - (as.integer(wfmelt_plot$trajectory)-1) *
        (400 / length(unique(wfmelt_plot$trajectory)))
      
      # Distinguish between colorwave 1/2 and 3/4
      wfmelt_plot$wavetype = readColorMapping(fileList[j])[as.character(wfmelt_plot$trajectory)]
      
      if(wfmelt_plot$wavetype[1] == 2 || wfmelt_plot$wavetype[1] == 1){
        # Label regions for plot
        wfmelt_plot$region = c('#000000'='initial',
                               '#0000FF'='LED',
                               '#00FF00'='TONE',
                               '#808080'='post-US')[wfmelt_plot$'time color']
        for (i in 1:length(sessionDF$trialName)) {
          if(sessionDF$Peak[i]>0){
            pko1 = peakOnset1(wfsub[[as.character(sessionDF$trialName[i])]],
                              sessionDF$Peak[i],
                              peakThreshold1,
                              peakThreshold2)
            if (!is.na(pko1)) {
              wfmelt_plot[wfmelt_plot$trajectory == as.character(sessionDF[i, 'trial']) &
                            wfmelt_plot$'time color' == '#808080' &
                            wfmelt_plot$index <= sessionDF[i, 'peak'] &
                            wfmelt_plot$index >= pko1, 'region'] = 'rising'
            }
          }}
        
        wfmelt_plot$region = factor(wfmelt_plot$region,
                                    levels = c('initial', 'LED','post-US', 'rising'))
        
        wfmelt_plot[wfmelt_plot$wavetype == 2 & wfmelt_plot$region != 'LED', 'region'] = 'initial'
        
        # Full waterfall plot
        full_waterfall = ggplot(wfmelt_plot, aes(x=shifted_t, y=shifted_value*1.5,
                                                 color=region,group=trajectory))
        full_waterfall = full_waterfall + theme_classic()
        full_waterfall = full_waterfall + theme(
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = 14)
        )
        full_waterfall = full_waterfall + geom_line(size=0.4) + theme(text=element_text(face="bold", size= 12))
        full_waterfall = full_waterfall + scale_color_manual(values=c(
          baseline_color, cue_color, airpuff_color, airpuff_color))
        full_waterfall = full_waterfall + xlab('\n                   time (ms)') +
          scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100)))
        full_waterfall = full_waterfall + ggtitle((paste(substring(fileList[j],1,file_substring_end), ':         ', CR_Pct, '% CR\nMedian CR Onset:        ',
                                                         CR_Onst, ' ms'))) 
        # Collapsed waterfall plots
        if(any(wfmelt_plot$wavetype == 2)){
          probeTrials = wfmelt_plot[wfmelt_plot$wavetype == 2,]
          collapsed1 =  ggplot(probeTrials, aes(x=t, y=value,
                                                color=region,group=trajectory))
          collapsed1 = collapsed1 + theme_classic() + theme(legend.position = 'none')
          collapsed1 = collapsed1 + geom_line(size=0.4)
          collapsed1 = collapsed1 + scale_color_manual(values=c(
            baseline_color, cue_color, airpuff_color, airpuff_color))
          collapsed1 = collapsed1 + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else{collapsed1 = ggplot()}
        
        if(any(wfmelt_plot$wavetype == 1)){
          pairedTrials = wfmelt_plot[wfmelt_plot$wavetype == 1,]
          collapsed2 =  ggplot(pairedTrials, aes(x=t, y=value,
                                                 color=region,group=trajectory))
          collapsed2 = collapsed2 + theme_classic() + theme(legend.position = 'none')
          collapsed2 = collapsed2 + geom_line(size=0.4)
          collapsed2 = collapsed2 + scale_color_manual(values=c(
            baseline_color, cue_color, airpuff_color, airpuff_color))
          collapsed2 = collapsed2 + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else {collapsed2 = ggplot()}
        
        
        if(postscript_flag == 1){
          setwd(path_waterfall_postcript)
          setEPS()
          postscript((paste(substring(fileList[j],1,file_substring_end), '-wf.eps',sep = '')), height=8, width=10)
          capture.output(print(grid.arrange(ggarrange(full_waterfall, nrow=1,ncol=1),ggarrange(collapsed1,collapsed2, nrow=1,ncol=2),heights=c(1.75,1))))
          dev.off()
        }
        setwd(path_waterfall_png)
        png(paste(substring(fileList[j],1,file_substring_end), '-wf.png',sep = ''),
            width     = 7.5,
            height    = 6.5,
            units     = "in",
            res       = 200
        )
        capture.output(grid.arrange(ggarrange(full_waterfall, nrow=1,ncol=1),ggarrange(collapsed1,collapsed2, nrow=1,ncol=2),heights=c(1.75,1)))
        dev.off()
        
        setwd(path_trajectories)
        wfmelt_plot$Animal = substring(fileList[j],1,4)
        wfmelt_plot$Geno = substring(fileList[j],1,1)
        wfmelt_plot$session = substring(fileList[j],6,8)
        wfmelt_plot$fullTrajectory = paste(wfmelt_plot$trajectory, wfmelt_plot$session, sep = '')
        for(traj in levels(wfmelt_plot$trajectory)){
          if(any(traj %in% sessionDF$trialName) == TRUE){
            wfmelt_plot$CR[wfmelt_plot$trajectory == traj] = sessionDF$CR[sessionDF$trialName == traj]
            wfmelt_plot$Paired[wfmelt_plot$trajectory == traj] = sessionDF$Paired[sessionDF$trialName == traj]
          } else{
            wfmelt_plot$CR[wfmelt_plot$trajectory == traj] = 0
            wfmelt_plot$Paired[wfmelt_plot$trajectory == traj] = 2
          }
        }
        write.csv(wfmelt_plot,paste(substring(fileList[j],1,file_substring_end),'-trajectories.txt', sep = ''))
      }
      
      if(wfmelt_plot$wavetype[1] == 3 || wfmelt_plot$wavetype[1] == 4){
        # Label regions for plot
        wfmelt_plot$region = c('#000000'='initial',
                               '#00FF00'='TONE',
                               '#808080'='post-US')[wfmelt_plot$'time color']
        for (i in 1:length(sessionDF$trialName)) {
          if(sessionDF$Peak[i]>0){
            pko1 = peakOnset1(wfsub[[as.character(sessionDF$trialName[i])]],
                              sessionDF$Peak[i],
                              peakThreshold1,
                              peakThreshold2)
            if (!is.na(pko1)) {
              wfmelt_plot[wfmelt_plot$trajectory == as.character(sessionDF[i, 'trial']) &
                            wfmelt_plot$'time color' == '#808080' &
                            wfmelt_plot$index <= sessionDF[i, 'peak'] &
                            wfmelt_plot$index >= pko1, 'region'] = 'rising'
            }
          }} # END OF LOOP THAT READS COLORWAVES 1/2
        
        wfmelt_plot$region = factor(wfmelt_plot$region,
                                    levels = c('initial', 'TONE','post-US', 'rising'))
        
        wfmelt_plot[wfmelt_plot$wavetype == 4 & wfmelt_plot$region != 'TONE', 'region'] = 'initial'
        
        # Full waterfall plot
        full_waterfall = ggplot(wfmelt_plot, aes(x=shifted_t, y=shifted_value*1.5,
                                                 color=region,group=trajectory))
        full_waterfall = full_waterfall + theme_classic()
        full_waterfall = full_waterfall + theme(
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = 14)
        )
        full_waterfall = full_waterfall + geom_line(size=0.4) + theme(text=element_text(face="bold", size= 12))
        full_waterfall = full_waterfall + scale_color_manual(values=c(
          baseline_color, cue_color, airpuff_color, airpuff_color))
        full_waterfall = full_waterfall + xlab('\n                   time (ms)') +
          scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100)))
        full_waterfall = full_waterfall + ggtitle((paste(substring(fileList[j],1,file_substring_end), ':         ', CR_Pct, '% CR\nMedian CR Onset:        ',
                                                         CR_Onst, ' ms'))) 
        # Collapsed waterfall plots
        if(any(wfmelt_plot$wavetype == 4)){
          probeTrials = wfmelt_plot[wfmelt_plot$wavetype == 4,]
          collapsed1 =  ggplot(probeTrials, aes(x=t, y=value,
                                                color=region,group=trajectory))
          collapsed1 = collapsed1 + theme_classic() + theme(legend.position = 'none')
          collapsed1 = collapsed1 + geom_line(size=0.4)
          collapsed1 = collapsed1 + scale_color_manual(values=c(
            baseline_color, cue_color, airpuff_color, airpuff_color))
          collapsed1 = collapsed1 + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else{collapsed1 = ggplot()}
        
        if(any(wfmelt_plot$wavetype == 3)){
          pairedTrials = wfmelt_plot[wfmelt_plot$wavetype == 3,]
          collapsed2 =  ggplot(pairedTrials, aes(x=t, y=value,
                                                 color=region,group=trajectory))
          collapsed2 = collapsed2 + theme_classic() + theme(legend.position = 'none')
          collapsed2 = collapsed2 + geom_line(size=0.4)
          collapsed2 = collapsed2 + scale_color_manual(values=c(
            baseline_color, cue_color, airpuff_color, airpuff_color))
          collapsed2 = collapsed2 + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else {collapsed2 = ggplot()}
        
        
        if(postscript_flag == 1){
          setwd(path_waterfall_postcript)
          setEPS()
          postscript((paste(substring(fileList[j],1,file_substring_end), '-wf.eps',sep = '')), height=8, width=10)
          capture.output(print(grid.arrange(ggarrange(full_waterfall, nrow=1,ncol=1),ggarrange(collapsed1,collapsed2, nrow=1,ncol=2),heights=c(1.75,1))))
          dev.off()
        }
        setwd(path_waterfall_png)
        png(paste(substring(fileList[j],1,file_substring_end), '-wf.png',sep = ''),
            width     = 7.5,
            height    = 6.5,
            units     = "in",
            res       = 200
        )
        capture.output(grid.arrange(ggarrange(full_waterfall, nrow=1,ncol=1),ggarrange(collapsed1,collapsed2, nrow=1,ncol=2),heights=c(1.75,1)))
        dev.off()
        
        setwd(path_trajectories)
        wfmelt_plot$Animal = substring(fileList[j],1,4)
        wfmelt_plot$Geno = substring(fileList[j],1,1)
        wfmelt_plot$session = substring(fileList[j],6,8)
        wfmelt_plot$fullTrajectory = paste(wfmelt_plot$trajectory, wfmelt_plot$session, sep = '')
        for(traj in levels(wfmelt_plot$trajectory)){
          if(any(traj %in% sessionDF$trialName) == TRUE){
            wfmelt_plot$CR[wfmelt_plot$trajectory == traj] = sessionDF$CR[sessionDF$trialName == traj]
            wfmelt_plot$Paired[wfmelt_plot$trajectory == traj] = sessionDF$Paired[sessionDF$trialName == traj]
          } else{
            wfmelt_plot$CR[wfmelt_plot$trajectory == traj] = 0
            wfmelt_plot$Paired[wfmelt_plot$trajectory == traj] = 2
          }
        }
        write.csv(wfmelt_plot,paste(substring(fileList[j],1,file_substring_end),'-trajectories.txt', sep = ''))
      } # END OF LOOP THAT READS COLORWAVES 3/4

    # Save dataframe, trajectory for each raw file
    
    setwd(path_session_dataframes)
    write.csv(sessionDF, paste(substring(fileList[j],1,file_substring_end),'-DF.csv', sep = ''))
    
    setwd(path_waterfall_raw)
    if(j%%20 == 0){
      pct = round(((j)/length(fileList))*100,0)
      print(paste('waterfalling is ', pct, ' % Complete'))
    }
  }}} # End brackets that reads in individual raw data files and doWaterfall if statement

# FILL DATAFRAME FOR EYELID CONDITIONING & CR TIMING HISTOGRAMS ----------------

setwd(path_summary_csv)
csvList = list.files()
hist_onset_list = list()
hist_peak_list = list()

for(cueisi in unique(substring(csvList,1,4))){
  if(cueisi != 'xxxx'){
    file = grep(cueisi, csvList)
    summaryDF = read.csv(csvList[file])
    for (animal in unique(substring(fileList,1,4))) {
      setwd(path_session_dataframes)
      fileList = list.files()
      temp = read.csv(fileList[1])
      temp = temp[0,]
      for(i in 1:length(fileList)){
        if(substring(fileList[i],1,4) == animal){
          csv = read.csv(fileList[i]) 
          if(csv$Cue[1] == TRUE){
            csv$Cue = 'T'
          }
          # If the temp session has startle probes, exclude them
          if(substring(fileList[i],file_substring_end,file_substring_end) > 0){
            temp = temp[temp$Paired == 1,]
          }
          temp = rbind(temp,csv)
        }
      }
      
      temp$X = 1:length(temp[,1])
      temp = temp[temp$Cue == substring(cueisi,1,1),]
      temp = temp[temp$ISI == substring(cueisi,2,4),]
      
      temp_probe = temp[temp$Paired == 0,]
      temp_probe = temp_probe[temp_probe$ID == 0,]
      
      # Calculate learning curve for each animal before removing probe trials
      
      if(length(temp[,1]) > 0){
        for(Session in unique(temp$Session)){
          tempSession = temp[temp$Session == Session,]
          summaryDF[[Session]][summaryDF$Animal == animal] = round(100*sum(tempSession$CR)/length(tempSession$CR), 1)
        }
        
        #Calculate bootstrapped UR for first 2 days of training:
        if (length(temp[,1]) > 0){
          X = substring(temp$Session,1,1)[1]
          tempSession1 = temp[temp$Session == paste(X, "01", sep=''),]
          tempSession2 = temp[temp$Session == paste(X, "02", sep=''),]
          tempSession3 = temp[temp$Session == paste(X, "03", sep=''),]
            
          tempSession = rbind(tempSession1, tempSession2,tempSession3)
          #tempSession = temp
          bootstrappedCoefficients_UR <- numeric(numberOfIterations)
          bootstrappedCoefficients_UR_Integrate <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            UR = tempSession$UR
            UR_Integrate = tempSession$UR_Integrate
            bootstrappedDataset_UR <- sample(UR, length(UR), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_UR[n] <- mean(bootstrappedDataset_UR)
            bootstrappedDataset_UR_Integrate <- sample(UR_Integrate, length(UR_Integrate), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_UR_Integrate[n] <- mean(bootstrappedDataset_UR_Integrate)
          }
          UR_median = round(as.numeric(quantile(bootstrappedCoefficients_UR, .500)),3) # median
          UR_Integrate_median = round(as.numeric(quantile(bootstrappedCoefficients_UR_Integrate, .500)),3) # median
          summaryDF$UR[summaryDF$Animal==animal] = UR_median
          summaryDF$UR_Integrate[summaryDF$Animal==animal] = UR_Integrate_median
        } else{
          summaryDF$UR[summaryDF$Animal==animal] = 0
          summaryDF$UR_Integrate[summaryDF$Animal==animal] = 0
          }
        
        # Find number of trials it takes each animal to reach some pre-selected % CR
        CR_vector = temp$CR
        
        sliding = 0
        for(i in 1:(length(CR_vector)-sliding_window_size)){
          sliding[i] = sum(CR_vector[i:(i+(sliding_window_size-1))])
        }
        
        trial25 = which(sliding >= round(sliding_window_size*.25,0))[1]
        if(is.na(trial25) == 1){
          trial25 = 0
          summaryDF$trial25[summaryDF$Animal== animal] = 0
        } else{
          summaryDF$trial25[summaryDF$Animal== animal] = trial25 + sliding_window_size-1
        }
        
        trial50 = which(sliding >= round(sliding_window_size*.5,0))[1]
        if(is.na(trial50) == 1){
          trial50 = 0
          summaryDF$trial50[summaryDF$Animal== animal] = 0
        } else{
          summaryDF$trial50[summaryDF$Animal== animal] = trial50 + sliding_window_size-1
        }
        
        trial75 = which(sliding >= round(sliding_window_size*.75,0))[1]
        if(is.na(trial75) == 1){
          trial75 = 0
          summaryDF$trial75[summaryDF$Animal== animal] = 0
        } else{
          summaryDF$trial75[summaryDF$Animal== animal] = trial75 + sliding_window_size-1
        }
        
        
        summaryDF$totalCR[summaryDF$Animal== animal] = (sum(CR_vector)/length(CR_vector))
        
        # Analyze paired trials for %CR, onset
        summaryDF$numTrials[summaryDF$Animal== animal] = length(temp[,2])
        temp_CR = temp[temp$Onset>0,]
        summaryDF$countCRs[summaryDF$Animal== animal] = sum(temp_CR$CR)
        if (length(temp_CR$Onset)>=end_paired[1]){
          temp1 = temp_CR[begin_paired[1]:end_paired[1],]
          bootstrappedCoefficients <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            onset = temp1$Onset
            bootstrappedDataset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
            bootstrappedCoefficients[n] <- mean(bootstrappedDataset)
          }
          onset_median = round(as.numeric(quantile(bootstrappedCoefficients, .500)),1) # median
          summaryDF$Onset1[summaryDF$Animal==animal] = onset_median
        } else{summaryDF$Onset1[summaryDF$Animal==animal] = 0}
        
        # Do second round if enough data
        if (length(temp_CR$Onset)>=end_paired[2]){
          temp1 = temp_CR[begin_paired[2]:end_paired[2],]
          bootstrappedCoefficients <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            onset = temp1$Onset
            bootstrappedDataset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
            bootstrappedCoefficients[n] <- mean(bootstrappedDataset)
          }
          onset_median = round(as.numeric(quantile(bootstrappedCoefficients, .500)),1) # median
          summaryDF$Onset2[summaryDF$Animal==animal] = onset_median
        } else{summaryDF$Onset2[summaryDF$Animal==animal] = 0}
      } else{temp_CR = NULL}
      
      # Now calculate stuff for probe trials
      if(length(temp_probe[,1]) > 0){
        summaryDF$numTrials_probe[summaryDF$Animal==animal] = length(temp_probe[,2])
        temp_probeCR = temp_probe[temp_probe$Peak>0,]
        summaryDF$countCRs_probe[summaryDF$Animal==animal] = sum(temp_probeCR$CR)
        
        if (length(temp_probeCR$Peak) >= end_probe[1]) {
          temp_probeCR1 = temp_probeCR[begin_probe[1]:end_probe[1],]
          temp_probeCR1 = temp_probeCR
          bootstrappedCoefficients_peak <- numeric(numberOfIterations)
          bootstrappedCoefficients_onset <- numeric(numberOfIterations)
          bootstrappedCoefficients_amp <- numeric(numberOfIterations)
          bootstrappedCoefficients_CR_Area <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            peak = temp_probeCR1$Peak
            onset = temp_probeCR1$Onset
            ampProbe1 = temp_probeCR1$maxAmp
            CR_Area = temp_probeCR1$CR_Area
            bootstrappedDataset_peak <- sample(peak, length(peak), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_peak[n] <- mean(bootstrappedDataset_peak)
            bootstrappedDataset_onset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_onset[n] <- mean(bootstrappedDataset_onset)
            bootstrappedDataset_ampProbe1 <- sample(ampProbe1, length(ampProbe1), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_amp[n] <- mean(bootstrappedDataset_ampProbe1)
            bootstrappedDataset_CR_Area <- sample(CR_Area, length(CR_Area), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_CR_Area[n] <- mean(bootstrappedDataset_CR_Area)
          }
          peak_median = round(as.numeric(quantile(bootstrappedCoefficients_peak, .500)),1) # median
          onset_median = round(as.numeric(quantile(bootstrappedCoefficients_onset, .500)),1) # median
          ampProbe1_median = round(as.numeric(quantile(bootstrappedCoefficients_amp, .500)),3) # amplitude
          CR_Area_median = round(as.numeric(quantile(bootstrappedCoefficients_CR_Area, .500)),3) # Area Under CR
          summaryDF$Peak1[summaryDF$Animal==animal] = peak_median
          summaryDF$onsetProbe1[summaryDF$Animal==animal] = onset_median
          summaryDF$ampProbe1[summaryDF$Animal==animal] = ampProbe1_median
          summaryDF$CR_Area1[summaryDF$Animal==animal] = CR_Area_median
        } else{
          summaryDF$Peak1[summaryDF$Animal==animal] = 0
          summaryDF$onsetProbe1[summaryDF$Animal==animal] = 0
          summaryDF$ampProbe1[summaryDF$Animal==animal] = 0
          summaryDF$CR_Area1[summaryDF$Animal==animal] = 0
        }
        
        # Do second round if enough data
        if (length(temp_probeCR$Peak) >= end_probe[2]) {
          temp_probeCR2 = temp_probeCR[begin_probe[2]:end_probe[2],]
          bootstrappedCoefficients_peak <- numeric(numberOfIterations)
          bootstrappedCoefficients_onset <- numeric(numberOfIterations)
          bootstrappedCoefficients_amp <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            peak = temp_probeCR2$Peak
            onset = temp_probeCR2$Onset
            ampProbe2 = temp_probeCR2$maxAmp
            bootstrappedDataset_peak <- sample(peak, length(peak), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_peak[n] <- mean(bootstrappedDataset_peak)
            bootstrappedDataset_onset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_onset[n] <- mean(bootstrappedDataset_onset)
            bootstrappedDataset_ampProbe2 <- sample(ampProbe2, length(ampProbe2), replace = TRUE, prob=NULL)
            bootstrappedCoefficients_amp[n] <- mean(bootstrappedDataset_ampProbe2)
          }
          peak_median = round(as.numeric(quantile(bootstrappedCoefficients_peak, .500)),1) # median
          onset_median = round(as.numeric(quantile(bootstrappedCoefficients_onset, .500)),1) # median
          ampProbe2_median = round(as.numeric(quantile(bootstrappedCoefficients_amp, .500)),3) # median
          summaryDF$Peak2[summaryDF$Animal==animal] = peak_median
          summaryDF$onsetProbe2[summaryDF$Animal==animal] = onset_median
          summaryDF$ampProbe2[summaryDF$Animal==animal] = ampProbe2_median
        } else{
          summaryDF$Peak2[summaryDF$Animal==animal] = 0
          summaryDF$onsetProbe2[summaryDF$Animal==animal] = 0
          summaryDF$amp2[summaryDF$Animal==animal] = 0
          }
      } else{temp_probeCR = NULL}
      
      # Make histogram of peak (from probe) and onset (from paired and probe) for histogram list 
      if(length(temp_probeCR$Onset) >= 25){
        histogram_onset = ggplot(temp_CR, aes(x = Onset, color = group)) + 
          geom_histogram(binwidth = 40, color = onset_histogram, fill = onset_histogram, na.rm = TRUE) + 
          theme_classic() + ggtitle(paste(animal, cueisi)) + xlim(0,700) + geom_vline(xintercept = (as.numeric(substring(cueisi,2,4))))
        hist_onset_list[[animal]] = histogram_onset
      } else {hist_onset_list[[animal]] = nullGrob()}
      
      if(length(temp_probeCR$Peak) >= 25){
        histogram_peak = ggplot(temp_probeCR, aes(x = Peak, color = group)) + 
          geom_histogram(binwidth = 40, color = peak_histogram, fill = peak_histogram, na.rm = TRUE) + 
          theme_classic() + ggtitle(paste(animal, cueisi)) + xlim(0,700) + geom_vline(xintercept = (as.numeric(substring(cueisi,2,4))))
        hist_peak_list[[animal]] = histogram_peak
      } else {hist_peak_list[[animal]] = nullGrob()}
    } # End of for loop that reads in each animal
    
    setwd(path_plots)
    if(length(hist_peak_list) > 0){
      pdf((paste(cueisi,'-PeakHistogram.pdf',sep = '')), height=(1.5*(ceiling(length(hist_peak_list))/4)), width=7.5)
      capture.output(print(multiplot(plotlist=hist_peak_list, cols = 4)))
      dev.off()
    }
    if(length(hist_onset_list) > 0 ){
      pdf((paste(cueisi,'-OnsetHistogram.pdf',sep = '')), height=(1.5*(ceiling(length(hist_peak_list))/4)), width=7.5)
      capture.output(print(multiplot(plotlist=hist_onset_list, cols = 4)))
      dev.off()
    }
    
    summaryDF[is.na(summaryDF)] = 0
    setwd(path_summary_csv)
    write.csv(summaryDF, csvList[file])
  } # End of if cueisi != 'xxxx' statement
} # End of loop that reads in each Summary csv file

# ANALYZE RAW FILES FOR STARTLE  ------------------------------------------

if(doStartle == 1){
  allPeaks = function(traj, halfWidth = halfW, minPeak=0.10) {
    if (is.matrix(traj) || is.data.frame(traj)) {
      out = lapply(
        X = data.frame(traj[ , setdiff(colnames(traj), 't')]),
        FUN = allPeaks,
        halfWidth = halfWidth
      )
      out$t = NA
      return(out[colnames(traj)])
    }
    if (length(halfWidth) > 0) {
      traj = smoothTrajectory(traj, halfWidth)
    }
    d1 = diff(traj)
    d1pos = ifelse(d1>0, 1, 0)
    peaks = which(diff(d1pos) == -1)
    if (length(halfWidth) > 0) {
      peaks = peaks + halfWidth
    }
    peaks = peaks + 1
    if (length(minPeak) > 0) {
      peaks = peaks[traj[peaks] >= minPeak]
    }
    return(peaks)
  }
  
  setwd(path_waterfall_raw)
  fileList <- list.files()
  
  if(any(substring(fileList,file_substring_end,file_substring_end)>0)){
    for (j in 1:length(fileList)) {
      if(substring(fileList[j],file_substring_end,file_substring_end) == 4){
        wf_startle = readWaterFall(fileList[j])
        wf_startle = data.frame(do.call(cbind, wf_startle))
        wf_startle$t = trialDur * (1:nrow(wf_startle)) / nrow(wf_startle) 
        fps = length(wf_startle[,1])/trialDur*1000
        traceStart = round((baselineDur*fps/1000),0)   # moves the analysis window to CS onset. outputs the frame number where this happens
        traceEnd = round(((baselineDur + analysisWindow)*fps/1000),0)  # moves the end of analysis window to analysisWindow ms after startle tone
        wf_startle = wf_startle[1:traceEnd,]
        wfsub_startle = subtractBaselines(wf_startle)
        
        # Get rid of all trials that aren't probe trials (ending in either a 0 or a 5)
        temp1 = wfsub_startle[grep('0$', colnames(wfsub_startle))]
        temp2 = wfsub_startle[grep('5$', colnames(wfsub_startle))]
        wfsub_startle = cbind(temp1,temp2,'t' = wfsub_startle$t)
        wfsub_startle[is.na(wfsub_startle)] = 0
        wfsub_startle = wfsub_startle[ , c(basenoise(wfsub_startle) < baseline_noise, TRUE), drop=FALSE]
        
        # Locate highest peak in each trajectory
        peakTiming_probe = topPeaks(wfsub_startle[traceStart:(traceEnd),], n=1, halfWidth = halfW) + (traceStart-1)
        peakTiming_probe[is.na(peakTiming_probe)] = 0
        sessionDF = data.frame(Animal = substring(fileList[j],1,4), Group = substring(fileList[j],1,1), Session = substring(fileList[j],6,8),
                               Trial = as.numeric(substring(colnames(peakTiming_probe)[1:(length(peakTiming_probe)-1)],7,9)),
                               Paired = 0, Startle = 0, Onset = 0, Peak = 0, ID = substring(fileList[j],file_substring_end,file_substring_end),
                               trialName=colnames(peakTiming_probe)[1:(length(peakTiming_probe)-1)], Cue = as.character(substring(fileList[j],10,10)), ISI = as.numeric(substring(fileList[j],11,13)))
        
        # Fill in data frame
        for(k in 1:length((sessionDF$Trial))) {
          if(peakTiming_probe[1,k] > 0){
            sessionDF$Peak[k] = topPeaks(wfsub_startle[traceStart:(traceEnd),k], n=1, halfWidth = halfW)*(1000/fps)
            sessionDF$Startle[k] = 1 
            sessionDF$Onset[k] = peakOnset1(wfsub_startle[ , k], peakTiming_probe[1, k], peakThreshold1_startle, peakThreshold2)*(1000/fps) - baselineDur
          } else{sessionDF$Startle = 0}
        }
        
        ## STARTLE WATERFALL PLOTS 
        
        if(include_session_waterfall == 1) {
          # Get time intervals from fileList[j] file
          ti = readTimeIntervals(fileList[j])
          ti[[1]] = ti[[1]][1:traceEnd]
          names(ti[[1]]) = as.character(wf_startle$t)[1:traceEnd]
          
          # Reshape wfsub_startle for plotting
          wfmelt_startle = wfsub_startle
          wfmelt_startle$index = 1:nrow(wf_startle)
          wfmelt_startle = wfmelt_startle[(100*fps/1000):traceEnd,]
          wfmelt_startle = melt(wfmelt_startle, id.vars=c('t', 'index'), variable.name='trajectory')
          wfmelt_startle$'time color' = ti[[1]][as.character(wfmelt_startle$t)]
          
          # Fan trajectories out for waterfall plotting
          wfmelt_startle$shifted_value = wfmelt_startle$value + (as.integer(wfmelt_startle$trajectory)-1) * 0.1
          wfmelt_startle$shifted_t = wfmelt_startle$t - (as.integer(wfmelt_startle$trajectory)-1) *
            (400 / length(unique(wfmelt_startle$trajectory)))
          
          # Label regions for plot
          wfmelt_startle$region = c('#000000'='initial',
                                    '#0000FF'='CS',
                                    '#808080'='post-US')[wfmelt_startle$'time color']
          for (i in 1:length(sessionDF$trialName)) {
            if(sessionDF$Peak[i]>0){
              pko1 = peakOnset1(wfsub_startle[[as.character(sessionDF$trialName[i])]],
                                sessionDF$Peak[i],
                                peakThreshold1,
                                peakThreshold2)
              if (!is.na(pko1)) {
                wfmelt_startle[wfmelt_startle$trajectory == as.character(sessionDF[i, 'trial']) &
                                 wfmelt_startle$'time color' == '#808080' &
                                 wfmelt_startle$index <= sessionDF[i, 'peak'] &
                                 wfmelt_startle$index >= pko1, 'region'] = 'rising'
              }
            }}
          
          wfmelt_startle$region = factor(wfmelt_startle$region,
                                         levels = c('initial', 'CS', 'post-US', 'rising'))
          wfmelt_startle$wavetype = readColorMapping(fileList[j])[as.character(wfmelt_startle$trajectory)]
          wfmelt_startle[wfmelt_startle$wavetype == 2 & wfmelt_startle$region != 'CS', 'region'] = 'initial'
          
          # Full waterfall plot
          full_waterfall = ggplot(wfmelt_startle, aes(x=shifted_t, y=shifted_value*1.5,
                                                      color=region,group=trajectory))
          full_waterfall = full_waterfall + theme_classic()
          full_waterfall = full_waterfall + theme(
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none',
            plot.title = element_text(size = 14)
          )
          full_waterfall = full_waterfall + geom_line(size=0.45) + theme(text=element_text(face="bold", size= 12))
          full_waterfall = full_waterfall + scale_color_manual(values=c(
            baseline_color, startle_color))
          full_waterfall = full_waterfall + xlab('\n                   time (ms)') +
            scale_x_continuous(breaks= seq(from = baselineDur, to = traceEnd, by = 25), labels = as.character(seq(from = 0, to = (traceEnd-baselineDur), by = 25)))
          full_waterfall = full_waterfall + ggtitle((paste(substring(fileList[j],1,file_substring_end), ':                       ', 'startlePct', '% CR\nMedian CR Onset:            ',
                                                           'startleOnset', ' ms'))) 
          # Collapsed startle waterfall plots
          probeTrials = wfmelt_startle[wfmelt_startle$wavetype == 2,]
          collapsed1 =  ggplot(probeTrials, aes(x=t, y=value,
                                                color=region,group=trajectory))
          collapsed1 = collapsed1 + theme_classic() + theme(legend.position = 'none')
          collapsed1 = collapsed1 + geom_line(size=0.4)
          collapsed1 = collapsed1 + scale_color_manual(values=c(
            baseline_color, startle_color))
          collapsed1 = collapsed1 + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = traceEnd, by = 25), labels = as.character(seq(from = 0, to = (traceEnd-baselineDur), by = 25))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
          
          subplot1 = ggarrange(collapsed1,ggplot(), nrow = 1,ncol = 2)
          if(postscript_flag == 1){
            setwd(path_waterfall_startle_postcript)
            setEPS()
            postscript((paste(substring(fileList[j],1,file_substring_end), '-waterfallStartle.eps',sep = '')), height = 7, width = 6)
            capture.output(print(ggarrange(full_waterfall, subplot1, nrow = 2, ncol = 1, heights=c(2,1))))
            dev.off()
          }
          setwd(path_waterfall_startle_png)
          png(paste(substring(fileList[j],1,file_substring_end), '-waterfallStartle.png', sep = ''),
              width     = 5,
              height    = 4,
              units     = "in",
              res       = 200
          )
          capture.output(print(ggarrange(full_waterfall, subplot1, nrow = 2, ncol = 1, heights=c(2,1))))
          dev.off()
          
          setwd(path_trajectories_startle)
          wfmelt_startle$Animal = substring(fileList[j],1,4)
          wfmelt_startle$Geno = substring(fileList[j],1,1)
          wfmelt_startle$session = substring(fileList[j],6,8)
          wfmelt_startle$fullTrajectory = paste(wfmelt_startle$trajectory, wfmelt_startle$session, sep = '')
          for(traj in levels(wfmelt_startle$trajectory)){
            
            wfmelt_startle$CR[wfmelt_startle$trajectory == traj] = sessionDF$CR[sessionDF$trialName == traj]
          }
          write.csv(wfmelt_startle,paste(substring(fileList[j],1,file_substring_end),'-trajectoriesStartle.txt', sep = ''))
        }
        
        # Save dataframe & trajectory CSV and postscript of raw startle traces
        setwd(path_session_dataframes_startle)
        write.csv(sessionDF, paste(substring(fileList[j],1,file_substring_end),'-startleDF.csv', sep = ''))
        
        setwd(path_waterfall_raw)
        if(j%%20 == 0){
          pct = round(((j)/length(fileList))*100,0)
          print(paste('startle waterfalling is ', pct, ' % Complete'))
        }
      }
    } # End bracket for loop that reads in the individual raw waterfall files for startle  
  } # End bracket for end statement that reads in files if there are any startle files
} # End bracket for doStartle flag


# FILL DATAFRAME FOR STARTLE ---------------------------------------------------

if(any((substring(fileList,file_substring_end,file_substring_end) > 0))){
  
  setwd(path_session_dataframes_startle)
  fileList = list.files()
  
  for (animal in unique(substring(fileList,1,4))) {
    temp = read.csv(fileList[1])
    temp = temp[0,]
    for(i in 1:length(fileList)){
      if(substring(fileList[i],1,4) == animal){
        csv = read.csv(fileList[i])
        if(csv$Cue[1] == TRUE){
          csv$Cue = 'T'
        }
        temp = rbind(temp,csv)
      }
    }
    
    temp$X = 1:length(temp[,1])
    temp = temp[temp$Paired == 0,]
    temp = temp[temp$ID == 4,]
    
    for(Session in unique(temp$Session)){
      tempSession = temp[temp$Session == Session,]
      startleDF[[Session]][startleDF$Animal == animal] = round(100*sum(tempSession$Startle)/length(tempSession$Startle), 1)
    }
    
    startleDF$numTrials[startleDF$Animal==animal] = length(temp[,2])
    startleDF$countStartle[startleDF$Animal==animal] = length(which(temp$Peak>0))
    
    # Save Summary CSV for Startle Trials
    setwd(path_summary_csv)
    write.csv(startleDF, paste('xxxx-StartleSummary.csv', sep = '')) 
    setwd(path_session_dataframes_startle)
  }
}


# PLOT AVERAGE TRAJECTORIES ACROSS ANIMALS --------------------------------

# Specify fps, just in case the waterfall analysis isn't performed
setwd(path_waterfall_raw)
fileList <- list.files()
wf = readWaterFall(fileList[1])
wf = data.frame(do.call(cbind, wf))
fps = length(wf[,1])/trialDur*1000

setwd(path_trajectories)
fileList <- list.files()

if(doTrajectories == 1){
  
  setwd(path_trajectories)
  fileList_all <- list.files()
  
  for(cueisi in unique(substring(fileList_all,10,13))){
    fileList = fileList_all[substring(fileList_all,10,13) == cueisi]
    animalList = unique(substring(fileList,1,4))
    n = 0
    responseAverage_Paired = data.frame(Animal = 'animal', Geno = 'geno', t = 0, responseAmplitude = 0, Paired = 1)
    responseAverage_Probe = data.frame(Animal = 'animal', Geno = 'geno', t = 0, responseAmplitude = 0, Paired = 0)
    
    for(i in 1:length(animalList)){
      dfAnimal = data.frame()
      for (k in 1:length(fileList)){
        if(substring(fileList[k],1,4) == animalList[i]){
          temp = read.csv(fileList[k])
          dfAnimal = rbind(dfAnimal,temp)
        }
      }
      dfAnimal = dfAnimal[dfAnimal$t != 'NA',]
      dfAnimal = dfAnimal[dfAnimal$CR > 0,]
      if (length(dfAnimal$value) >= fps*trialDur*(trajMin)/1000){
        if(restrict == 1){
          if (length(dfAnimal$value) >= fps*trialDur*((trajRestrict[1]+1)+trajRestrict[2])/1000){
            dfAnimal = dfAnimal[(fps*trialDur*(trajRestrict[1]+1)/1000):(fps*trialDur*((trajRestrict[1]+1)+trajRestrict[2])/1000),]
          } else {dfAnimal = 0}
          
        } 
        if(length(dfAnimal) > 1){
          dfAnimal_Paired = dfAnimal[dfAnimal$Paired == 1,]
          dfAnimal_Probe = dfAnimal[dfAnimal$Paired == 0,]
          tempPaired = data.frame()
          tempProbe = data.frame()
          for(t in unique(dfAnimal$t)){
            tpaired = data.frame(Animal = dfAnimal_Paired$Animal[i], Geno = dfAnimal_Paired$Geno[i], 't' = t, responseAmplitude = mean(dfAnimal_Paired$value[dfAnimal_Paired$t == t]), Paired = 1)
            #tpaired$responseAmplitude = smoothTrajectory(tpaired$responseAmplitude, halfWidth = 4)
            tprobe = data.frame(Animal = dfAnimal_Probe$Animal[i], Geno = dfAnimal_Probe$Geno[i], 't' = t, responseAmplitude = mean(dfAnimal_Probe$value[dfAnimal_Probe$t == t]), Paired = 0)
           # tprobe$responseAmplitude = smoothTrajectory(tprobe$responseAmplitude, halfWidth = 4)
            tempPaired = rbind(tempPaired, tpaired)
            tempProbe = rbind(tempProbe, tprobe)
          }
          responseAverage_Paired = rbind(responseAverage_Paired, tempPaired)
          responseAverage_Probe = rbind(responseAverage_Probe, tempProbe)
        }}
      n = n+1
      print(paste(n, ' out of ', length(unique(animalList)), ' animal trajectory averages complete'))
    }
    
    responseAverage_Paired = responseAverage_Paired[2:length(responseAverage_Paired$t),]
    responseAverage_Probe = responseAverage_Probe[2:length(responseAverage_Probe$t),]
    
    avgPaired_plot = ggplot(responseAverage_Paired, aes(x=t, y=responseAmplitude, group = Animal, color = Geno)) + geom_line(alpha = .5) + scale_x_continuous(breaks = c(200,400,600,800,1000,1200), 
                    labels=c('0', '200', '400', '600', '800', '1000')) + theme_classic() + ylim(-0.15,1) + scale_color_manual(values=c(group1_color,group2_color, 'pink')) + theme(legend.position="none")
    
    avgProbe_plot = ggplot(responseAverage_Probe, aes(x=t, y=responseAmplitude, group = Animal, color = Geno)) + geom_line(alpha = .5) + scale_x_continuous(breaks = c(200,400,600,800,1000,1200),
                    labels=c('0', '200', '400', '600', '800', '1000')) + theme_classic() + ylim(-0.15,1) + scale_color_manual(values=c(group1_color,group2_color, 'pink'))
    
    trajPlot = ggarrange(avgPaired_plot,avgProbe_plot, nrow = 1,ncol = 2)
    
    setwd(path_plots)
    pdf(paste(cueisi,'-TrajPlots.pdf', sep =''), height=4, width=10)
    capture.output(print(trajPlot))
    dev.off()
    
    setwd(path_trajectories)
  }
}

# PLOTS WITH STATISTICAL TESTS ---------------------------------------

setwd(path_summary_csv)
csvList = list.files()

for(cueisi in unique(substring(csvList,1,4))){
  file = grep(cueisi, csvList)
  summaryDF = read.csv(csvList[file])
  
  # If summary csv is anything but the StartleSummary.csv
  if(cueisi != 'xxxx'){
    # t-tests for onset and peak data, and wilcox test for CR data. 
    Onset1 = summaryDF[summaryDF$Onset1>0,]
    Onset2 = summaryDF[summaryDF$Onset2>0,]
    trial25 = summaryDF[summaryDF$trial25 > 0,]
    trial50 = summaryDF[summaryDF$trial50 > 0,]
    totalCR = summaryDF
    Peak1 = summaryDF[summaryDF$Peak1>0,]
    Peak2 = summaryDF[summaryDF$Peak2>0,]
    onsetProbe1 = summaryDF[summaryDF$onsetProbe1>0,]
    onsetProbe2 = summaryDF[summaryDF$onsetProbe2>0,]
    Group = unique(summaryDF$Group)
    
    summaryDFm = melt(summaryDF)
    names(summaryDFm)[3] = 'TrainingSession'
    names(summaryDFm)[4] = 'ConditionedResponse'
    #IMPROVE
    summaryDFm1 = summaryDFm[substring(summaryDFm$TrainingSession,1,1) == 'A',]
    summaryDFm2 = summaryDFm[substring(summaryDFm$TrainingSession,1,1) == 'D',]
    summaryDFm = rbind(summaryDFm1, summaryDFm2)
    # Create Group Average Dataframe
    groupAvg = summaryDFm
    groupAvg = groupAvg[0,2:4]
    
    for(group in Group){
      temp1 = summaryDFm[summaryDFm$Group == group,]
      for(session in unique(summaryDFm$TrainingSession)){
        temp2 = temp1[temp1$TrainingSession == session,]
        averageCR = mean(temp2$ConditionedResponse)
        groupAvg = rbind(groupAvg, data.frame('Group' = group, 'TrainingSession' = session, 'ConditionedResponse' = averageCR))
      }
    }
    
    # Make sure each has more than min values so that they can be compared, otherwise do not caluclate test
    group1 = Onset1$Onset1[Onset1$Group == Group[1]]
    group2 = Onset1$Onset1[Onset1$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      op1_ttest = t.test(group1,group2)
      op1 = ggplot(Onset1, aes(x = Group, y = Onset1, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE)  + 
        geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
        ylim(75,400) + ggtitle(paste('t.test, p = ', round(as.numeric(op1_ttest[3]),3)))
    }else{
      op1_ttest = 0
      op1 = ggplot()
    }
    
    group1 = Onset2$Onset2[Onset2$Group == Group[1]]
    group2 = Onset2$Onset2[Onset2$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      op2_ttest = t.test(group1,group2)
      op2 = ggplot(Onset2, aes(x = Group, y = Onset2, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE) + 
        geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
        ylim(75,400) + ggtitle(paste('t.test, p = ', round(as.numeric(op2_ttest[3]),3)))
    }else{
      op2_ttest = 0
      op2 = ggplot()
    }
    
    group1 = Peak1$Peak1[Peak1$Group == Group[1]]
    group2 = Peak1$Peak1[Peak1$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      pk1_ttest = t.test(group1,group2)
      pk1 = ggplot(Peak1, aes(x = Group, y = Peak1, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE)  + 
        geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
        ylim(50,750) + ggtitle(paste('t.test, p = ', round(as.numeric(pk1_ttest[3]),3)))
    }else{
      pk1_ttest = 0
      pk1 = ggplot()
    }
    
    group1 = Peak2$Peak2[Peak2$Group == Group[1]]
    group2 = Peak2$Peak2[Peak2$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      pk2_ttest = t.test(group1,group2)
      pk2 = ggplot(Peak2, aes(x = Group, y = Peak2, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE) + 
        geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
        ylim(50,750) + ggtitle(paste('t.test, p = ', round(as.numeric(pk2_ttest[3]),3)))
    }else{
      pk2_ttest = 0
      pk2 = ggplot()
    }
    
    group1 = onsetProbe1$onsetProbe1[onsetProbe1$Group == Group[1]]
    group2 = onsetProbe1$onsetProbe1[onsetProbe1$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      po1_ttest = t.test(group1,group2)
      po1 = ggplot(onsetProbe1, aes(x = Group, y = onsetProbe1, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE)  + 
        geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
        ylim(50,750) + ggtitle(paste('t.test, p = ', round(as.numeric(po1_ttest[3]),3)))
    }else{
      po1_ttest = 0
      po1 = ggplot()
    }
    
    group1 = onsetProbe2$onsetProbe2[onsetProbe2$Group == Group[1]]
    group2 = onsetProbe2$onsetProbe2[onsetProbe2$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      po2_ttest = t.test(group1,group2)
      po2 = ggplot(onsetProbe2, aes(x = Group, y = onsetProbe2, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE) + 
        geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
        ylim(50,750) + ggtitle(paste('t.test, p = ', round(as.numeric(po2_ttest[3]),3)))
    }else{
      po2_ttest = 0
      po2 = ggplot()
    }
    
    group1 = trial25$trial25[trial25$Group == Group[1]]
    group2 = trial25$trial25[trial25$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      trial25_wilcox = wilcox.test(group1,group2)
      t25 = ggplot(trial25, aes(x = Group, y = trial25, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE)  + 
        geom_point(color = 'black', shape=1, size=2, na.rm = TRUE) + ylim(0,1250)
        ggtitle(paste('rank-sum, p =', round(as.numeric(trial25_wilcox[3]),4))) + theme(plot.title = element_text(hjust = -.75))
    }else{
      trial25_wilcox = 0
      t25 = ggplot()
    }
    
    group1 = totalCR$totalCR[totalCR$Group == Group[1]]
    group2 = totalCR$totalCR[totalCR$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      totalCR_wilcox = wilcox.test(group1,group2)
      totalCR_plot = ggplot(totalCR, aes(x = Group, y = totalCR*100, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5, na.rm = TRUE)  + 
        geom_point(color = 'black', shape=1, size=2, na.rm = TRUE) + ylim(0,100) +
        ggtitle(paste('rank-sum, p =', round(as.numeric(totalCR_wilcox[3]),5))) + theme(plot.title = element_text(hjust = -.75))
    }else{
      trial75_wilcox = 0
      totalCR_plot = ggplot()
    }
    
    cr_plot = ggplot(summaryDFm, aes(x = TrainingSession, y= ConditionedResponse, group = Animal, color = Group))  + theme_classic() +
      geom_line(size=.5, alpha= .3) + theme(legend.position="none") + scale_color_manual(values = c(group1_color,group2_color)) +
      geom_line(data = groupAvg, aes(x=TrainingSession, y=ConditionedResponse, group = Group, color= Group), size = 3) + ylim(0,100)
    
    summaryDF$Categorical.Analysis25 = ifelse(summaryDF$trial25 > 0, 'learner', 'non-learner')
    summaryDF$Categorical.Analysis50 = ifelse(summaryDF$trial50 > 0, 'learner', 'non-learner')
    
    subplot1 = ggarrange(t25,totalCR_plot, nrow = 2,ncol = 1)
    subplot2 = ggarrange(cr_plot, subplot1, nrow = 1, ncol = 2, widths = c(2,1))
    subplot3 = ggarrange(op1,op2,po1,pk1, nrow = 2, ncol = 2)
    
    setwd(path_plots)
    pdf((paste(cueisi,'-AnalysisPlots.pdf',sep = '')), height=7, width=7)
    capture.output(print(ggarrange(subplot2,subplot3, nrow = 2, ncol = 1)))
    dev.off()
    
    print(sjt.xtab(summaryDF$Group, summaryDF$Categorical.Analysis25, title = paste(dataFolder, cueisi), file = paste(cueisi, '-Contingency-25.doc', sep='')))
    dev.off()
    
    print(sjt.xtab(summaryDF$Group, summaryDF$Categorical.Analysis50, title = paste(dataFolder, cueisi), file = paste(cueisi, '-Contingency-50.doc', sep='')))
    dev.off()

    setwd(path_summary_csv)
  } else { # This else statement goes with the if csv file is != 'xxxx'
    for(i in 1:length(summaryDF[,1])){
      summaryDF$totalStartle[i] = round((summaryDF$countStartle[i])/(summaryDF$numTrials[i])*100,1)
    }
    Group = unique(summaryDF$Group)
    
    summaryDFm = melt(summaryDF)
    names(summaryDFm)[3] = 'TrainingSession'
    names(summaryDFm)[4] = 'StartlePercent'
    summaryDFm1 = summaryDFm[substring(summaryDFm$TrainingSession,1,1) == 'A',]
    summaryDFm2 = summaryDFm[substring(summaryDFm$TrainingSession,1,1) == 'D',]
    summaryDFm = rbind(summaryDFm1, summaryDFm2)
    # Create Group Average Dataframe
    groupAvg = summaryDFm
    groupAvg = groupAvg[0,2:4]
    
    for(group in Group){
      temp1 = summaryDFm[summaryDFm$Group == group,]
      for(session in unique(summaryDFm$TrainingSession)){
        temp2 = temp1[temp1$TrainingSession == session,]
        averageCR = mean(temp2$StartlePercent)
        groupAvg = rbind(groupAvg, data.frame('Group' = group, 'TrainingSession' = session, 'StartlePercent' = averageCR))
      }
    }
    
    group1 = summaryDF$totalStartle[summaryDF$Group == Group[1]]
    group2 = summaryDF$totalStartle[summaryDF$Group == Group[2]]
    if(length(group1) >= min_test && length(group2) >= min_test){
      totalStartle_wilcox = wilcox.test(group1,group2)
      totalStartle = ggplot(summaryDF, aes(x = Group, y = totalStartle, group = Group)) + theme_classic() +
        geom_boxplot(fill= c(group1_color,group2_color), width = .75, alpha=.5)  + geom_point(color = 'black', shape=1, size=2) +
        ggtitle(paste('rank-sum, p =', round(as.numeric(totalStartle_wilcox[3]),4))) + theme(plot.title = element_text(hjust = 0))
    }else{
      totalStartle_wilcox = 0
      totalStartle = ggplot()
    }
    
    startle_plot = ggplot(summaryDFm, aes(x = TrainingSession, y= StartlePercent, group = Animal, color = Group))  + theme_classic() +
      geom_line(size=1, alpha= .5) + theme(legend.position="none") + scale_color_manual(values = c(group1_color,group2_color)) +
      geom_line(data = groupAvg, aes(x=TrainingSession, y=StartlePercent, group = Group, color= Group), size = 3)
    
    setwd(path_plots)
    pdf((paste(cueisi,'-StartlePlots.pdf',sep = '')), height=3, width=7)
    capture.output(print(ggarrange(startle_plot,totalStartle, nrow = 1, ncol = 2, widths = c(2,1))))
    dev.off()
    setwd(path_summary_csv)
  }
} # End loop that reads in each cueisi summary csv

