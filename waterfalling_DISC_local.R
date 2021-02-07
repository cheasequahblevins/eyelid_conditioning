#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(tools)
library(plyr)
library(gridExtra)
library(ggpubr)
library(magrittr)
library(sjPlot)
library(stargazer)
library(grid)
library(snakecase)
library(RGraphics)
library(coin)
options(digits = 4)

# VARIABLES SPECIFIED BY USER ---------------------------------------------

dataFolder = 'AgedPS1x-Mates-AC-collapsed'
baselineDur = 200
trialDur = 1200

analyzeTraces = 1
doTrajectories = 1
doTrajectories_aligned = 1
doIndividual_trajPlots = 1
doIndividual_trajPlots_aligned = 1
doTimingPlots = 0
doSessionWaterfall = 0
doPostscriptWaterfall = 0
Discrim_Session = 0
peakThreshold1 = 0.065             # Threshold for classifying trial as CR
peakThreshold2 = 0.040              # Onset Threshold, or trials already classified as CRs
baseline_noise = 0.025
offsetStart = 0                    # Moves analysis window to offsetStart ms after onset of CS before looking for CRs
offsetEnd = 100                    # Moves analysis window to offsetEnd ms before trace end
usOffset = 25                      # Moves analysis window usOffset ms before airpuff 
halfW_analysis = 0                 # Half-width used by smoothing kernel in analysis
halfW_plotting = 0                 # Half-width used by smoothing kernel in plotting
file_substring_end = 15            # Length of file name (minus the '.txt' part)

# Define and Restrict Analysis of CRs
restrict = 1
trajMin = 80
trajRestrict = c(31,80)
trajRestrict_Probe = c(6,16)
min_test = 3
sliding_window_size = 50           # Size of sliding window for trials to learning onset analysis
numberOfIterations = 1000       # Number of iterations used in bootstrapping onset/peak values
lowerCI = 0.1
upperCI = 0.9
traj_truc = 500                     # How many ms to collect from each CR-aligned-to-Onset
traj_offset = 100                      # How many ms to plot before CR-aligned onset
deltaSlope_ms = 50                 # How many ms from CR Onset to measure amplitude again
deltaSlope_cutoff = 0.25

# Startle Analysis
doStartle = 0
analysisWindow_startle = 200
peakThreshold_startle = .05

# Colors of waterfall and analysis plots

LED_color = 'blue'
TONE_color = 'pink'

baseline_color = 	'#A0522D'   # Sienna Brown
airpuff_color = '#D3DDDC'     # Gray
cue_color = 'blue'
startle_color = 'blue'

onset_histogram = 'black'
peak_histogram = 'black'

group1_color = 'red'
group2_color = 'black'
group3_color = 'green'


# CREATE PATHS AND FOLDERS  

path_eyelid_behavior = '/Users/Cheasequah/Desktop/eyelid_behavior/'
path_waterfall_raw = paste(path_eyelid_behavior, 'waterfall_raw/', dataFolder, sep = '')
path_eyelid_analysis = paste(path_eyelid_behavior, 'eyelid_analysis/', dataFolder, sep = '')
path_session_dataframes = paste(path_eyelid_analysis, '/session_dataframes/', sep = '')
path_animal_dataframes = paste(path_eyelid_analysis, '/animal_dataframes/', sep = '')
path_timing_plots = paste(path_eyelid_analysis, '/timing_plots/', sep = '')
path_session_dataframes_startle = paste(path_eyelid_analysis, '/session_dataframes_startle/', sep = '')
path_waterfall_postcript = paste(path_eyelid_analysis, '/waterfall_poscript', sep = '')
path_waterfall_png = paste(path_eyelid_analysis, '/waterfall_png', sep = '')
path_waterfall_startle_postcript = paste(path_eyelid_analysis, '/waterfall_startle_postcript', sep='')
path_waterfall_startle_png = paste(path_eyelid_analysis, '/waterfall_startle_png', sep='')
path_trajectories = paste(path_eyelid_analysis, '/trajectories', sep = '')
path_trajectories_aligned = paste(path_eyelid_analysis, '/trajectories_aligned', sep = '')
path_trajectories_startle = paste(path_eyelid_analysis, '/trajectories_startle', sep = '')
path_summary_csv = paste(path_eyelid_analysis, '/summary_csv/', sep = '')
path_plots = paste(path_eyelid_analysis, '/plots/', sep = '')
path_traj_plots = paste(path_eyelid_analysis, '/traj_plots/', sep = '')

setwd(path_waterfall_raw)
fileList <- list.files()
if (!file.exists(path_eyelid_analysis)){
  dir.create(path_eyelid_analysis) }
if (!file.exists(path_session_dataframes)){
  dir.create(path_session_dataframes) }
if (!file.exists(path_animal_dataframes)){
  dir.create(path_animal_dataframes) }
if (!file.exists(path_timing_plots) && doTimingPlots == 1){
  dir.create(path_timing_plots) }
if (!file.exists(path_session_dataframes_startle) && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 ){
  dir.create(path_session_dataframes_startle) }
if (!file.exists(path_waterfall_postcript) && doSessionWaterfall == 1 && doPostscriptWaterfall == 1){
  dir.create(path_waterfall_postcript) }
if (!file.exists(path_waterfall_png) && doSessionWaterfall == 1){
  dir.create(path_waterfall_png) }
if (!file.exists(path_waterfall_startle_postcript) && doSessionWaterfall == 1 && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 && doPostscriptWaterfall == 1){
  dir.create(path_waterfall_startle_postcript) }
if (!file.exists(path_waterfall_startle_png) && doSessionWaterfall == 1 && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 ){
  dir.create(path_waterfall_startle_png) }
if (!file.exists(path_trajectories)){
  dir.create(path_trajectories) }
if (!file.exists(path_trajectories_aligned)){
  dir.create(path_trajectories_aligned) }
if (!file.exists(path_trajectories_startle) && sum(as.numeric(substring(fileList,file_substring_end,file_substring_end))) > 0 ){
  dir.create(path_trajectories_startle) }
if (!file.exists(path_summary_csv)){
  dir.create(path_summary_csv) }
if (!file.exists(path_plots)){
  dir.create(path_plots) }
if (!file.exists(path_traj_plots) && doTrajectories == 1){
  dir.create(path_traj_plots) }

# DEFINE FUNCTIONS --------------------------------------------------------
readTimeIntervals = function(f) {
  lines = readLines(f)
  waveStarts = sort(grep('^WAVES/', lines))
  minStart = min(waveStarts)
  waveEnds = sort(grep('^END', lines))[1:length(waveStarts)]
  waves = list()
  waveTypes = unique(name_rows(melt(readColorMapping(f)))$value)
  i=0
  for (t in waveTypes) {
    i = i+1
    waves[[t]] = do.call(rbind, strsplit(
      sub('^\\s+', '', lines[(waveStarts[i]+2):(waveEnds[i]-1)]),
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
  colSds(x[x$t <= maxt, setdiff(colnames(x), 't')])
}

smoothTrajectory = function(traj, halfWidth = halfW_analysis) {
  kernapply(traj, k=kernel('daniell', halfWidth))
}

allPeaks = function(traj, halfWidth = halfW_analysis, minPeak = 0.15) {
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

topPeaks = function(traj, n=1, halfWidth = halfW_analysis) {
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

peakOnset = function(traj, peak, threshold1, threshold2) {
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
  start = start
  end = end
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

# ANALYZE RAW TEXT FILES FOR EACH ANIMAL CR ANALYSIS  --------------------------

setwd(path_waterfall_raw)
fileList <- list.files()

if(analyzeTraces == 1) {
  for (j in 1:length(fileList)) {
    wf = readWaterFall(fileList[j])
    wf = data.frame(do.call(cbind, wf))
    wf$t = trialDur * (1:nrow(wf)) / nrow(wf) # adds a column to the dataframe of the actual time (ms) of the trial
    fps = length(wf[,1])/trialDur*1000
    traceStart = round(((baselineDur + offsetStart)*fps/1000),0) + 1    # moves the analysis window offsetStart ms after CS onset. outputs the frame number where this happens
    traceEnd = round(((trialDur-offsetEnd)*fps/1000),0)           # moves the end of analysis window to offsetEnd before end of trial
    # Subtract baseline from waterfall(wf) object
    wfsub = round(subtractBaselines(wf),4)
    wfsub[is.na(wfsub)] = 0
    wfsub = wfsub[ , c(basenoise(wfsub) < baseline_noise, TRUE), drop=FALSE]
    
    # Get time intervals from fileList[j] file
    ti = readTimeIntervals(fileList[j])
    #re-arrange time intervals so that 1 always equals blue paired, 2 blue probe, 3 green paired, 4 green probe
    length_ti = 0
    for(i in 1:length(ti)){
      if(!is.null(ti[[i]])){
        length_ti = c(length_ti,i)
      }}
    length_ti = length_ti[length_ti > 0]
    ti_hold = ti
    modifiedColormap = data.frame(original = 1:length(ti), modified = 1:length(ti))
    for(i in length_ti){
      if(identical(unique(ti[[i]]), c("#000000", "#0000FF", "#808080"))){ # blue paired
        ti_hold[[1]] = ti[[i]]
        modifiedColormap$modified[i] = 1
      } else if(identical(unique(ti[[i]]),c("#000000", "#0000FF"))){ # blue probe
        ti_hold[[2]] = ti[[i]]
        modifiedColormap$modified[i] = 2
      } else if(identical(unique(ti[[i]]),c("#000000", "#00FF00", "#808080"))){ # green paired
        ti_hold[[3]] = ti[[i]]
        modifiedColormap$modified[i] = 3
      } else if(identical(unique(ti[[i]]),c("#000000", "#00FF00"))){ # green probe
        ti_hold[[4]] = ti[[i]]
        modifiedColormap$modified[i] = 4
      }}
    ti = ti_hold
    
    # pull trialNames so that sessionDF can be initiated at the correct size
    # Only pull trial names that survive baseline noise exclusion
    trialName = names(wfsub)[names(wfsub)!='t']
    # get all colorwaves (1,2,3,4) and then keep only those that correspond to kept trials
    allColorwaves = name_rows(melt(readColorMapping(fileList[j])))
    keptColorwaves = allColorwaves$value[allColorwaves$.rownames %in% trialName]

    # Create sessionDF object and fill in basic information obtained from raw file name
    sessionDF = data.frame(trialName = trialName, trialNum = as.numeric(substring(trialName,7,9)), Animal = substring(fileList[j],1,4), 
                           Group = substring(fileList[j],1,1), Session = substring(fileList[j],6,8), Session_UR = paste(substring(fileList[j],6,8), '_UR', sep=''),
                           Paired = 0, ISI = 0, ID = substring(fileList[j],file_substring_end,file_substring_end), wavetype = keptColorwaves,
                           CR = 0, Onset = 0, Peak = 0, Slope = 0, UR = 0, UR_Int = 0, maxAmp = 0, ampInt=0)
   
    # Add "Cue1" or "Cue2" to sessionDF
    for(i in 1:length(sessionDF$wavetype)){
      if(sessionDF$wavetype[i] == 1 || sessionDF$wavetype[i] == 2){
        sessionDF$Cue[i] = 'Cue1'
      } else if(sessionDF$wavetype[i] == 3 || sessionDF$wavetype[i] == 4){
        sessionDF$Cue[i] = 'Cue2'
      }}

    # Get ISI from full timeIntervals/colorwave list object (ti), fill in ISI & Paired
    for(i in 1:length(ti)){
      # Only look at the timeIntervals/colorwaves that actually exist in session
      if(!is.null(ti[[i]])){
        uniqueColors = unique(ti[[i]])
        if(length(uniqueColors) > 2){ # paired trial
          wave = ti[[i]]
          sessionDF$Paired[sessionDF$wavetype == i] = 1
          sessionDF$ISI[sessionDF$wavetype == i] = as.numeric(wfsub$t[which(wave == uniqueColors[3])[1]]) + 1000/fps
        } else{
          sessionDF$Paired[sessionDF$wavetype == i] = 0
          sessionDF$ISI[sessionDF$wavetype == i] = 0
        }}}
    
    # Locate peak and onset timing in each trialName (paired and probe)
    sessionDF = sessionDF[!is.na(sessionDF$trialName),]
    for(k in sessionDF$trialName){
      i = which(sessionDF$trialName == k)
      if(sessionDF$Paired[i]==1){ # Paired Trial
        if(!is.na(topPeaks(wfsub[traceStart:((sessionDF$ISI[i]-usOffset)/(1000/fps)),k], n=1, halfWidth = halfW_analysis))){
          sessionDF$Peak[i] = (topPeaks(wfsub[traceStart:((sessionDF$ISI[i]-usOffset)*(fps/1000)),k], n=1, halfWidth = halfW_analysis))*(1000/fps)
          sessionDF$Onset[i] = peakOnset(wfsub[traceStart:((sessionDF$ISI[i]-usOffset)*(fps/1000)),k], sessionDF$Peak[i]*(fps/1000),
                                         peakThreshold1, peakThreshold2)*(1000/fps)
          sessionDF$Slope[i] = wfsub[(sessionDF$Onset[i] + deltaSlope_ms + baselineDur)*(fps/1000),k]
          sessionDF$maxAmp[i] = maxAmp(wfsub[traceStart:((sessionDF$ISI[i]-usOffset)*(fps/1000)),i])
          sessionDF$UR[i] = maxAmp(wfsub[(sessionDF$ISI[i]/(1000/fps)):traceEnd,i])
          sessionDF$UR_Int[i] = integrateArea(wfsub[,i], sessionDF$ISI[i]/(1000/fps), traceEnd, fps, peakThreshold1)
          sessionDF$CR[i] = 1
        } else{ # Not CR, but we still want UR information
          sessionDF$UR[i] = maxAmp(wfsub[(sessionDF$ISI[i]*(fps/1000)):traceEnd,i])
          sessionDF$UR_Int[i] = integrateArea(wfsub[,i], sessionDF$ISI[i]*(fps/1000), traceEnd, fps, peakThreshold1)
        }
      } else if(sessionDF$Paired[i] == 0){ # Probe Trial 
        if(!is.na(topPeaks(wfsub[traceStart:(traceEnd), k], n=1, halfWidth = halfW_analysis))){
        sessionDF$Peak[i] = topPeaks(wfsub[traceStart:(traceEnd),k], n=1, halfWidth = halfW_analysis)*(1000/fps)
        sessionDF$Onset[i] = peakOnset(wfsub[traceStart:(traceEnd),k], sessionDF$Peak[i]*(fps/1000),
                                    peakThreshold1, peakThreshold2)*(1000/fps)
        sessionDF$Slope[i] = wfsub[(sessionDF$Onset[i] + deltaSlope_ms + baselineDur)*(fps/1000),k]
        sessionDF$ampInt[i] = integrateArea(wfsub[,i], traceStart, traceEnd, fps, peakThreshold1)
        sessionDF$maxAmp[i] = maxAmp(wfsub[(traceStart:traceEnd),i])
        sessionDF$CR[i] = 1
        }
      } 
        else{ # Unknown trial type
        sessionDF$trialName[i] = 'unknown'
        }
    } # End for(k in sessionDF$trialName) statement

    # calculate a %CR and mean response onset to print on waterfall plots
    # First, initialize variables to 'n/a' so incorrect information is never printed on plot
    CROnst_Cue1 = 'n/a'
    CROnst_Cue2 = 'n/a'
    CRPct_Cue1 = 0
    CRPct_Cue2 = 0
    for(cue in unique(sessionDF$Cue)){
      assign(paste('CRPct_', cue, sep=''), round(sum(sessionDF$CR[sessionDF$Cue == cue])/length(sessionDF$CR[sessionDF$Cue == cue])*100,0))
      onsets = sessionDF$Onset[sessionDF$Cue == cue]
      onsets = onsets[onsets>0]
      if(length(onsets)>4){
        assign(paste('CROnst_', cue, sep=''), round(mean(onsets),0))
      } else{
        assign(paste('CROnst_', cue, sep=''), 'n/a')
        }}

    ## MAKE WATERFALL PLOTS

      # smooth a waterfall object for plotting
    if(halfW_plotting>0){
      wfsub_trunc = wfsub[((halfW_plotting+1):(length(wfsub[,1])-halfW_plotting)),]
      wfsub_plot = wfsub_trunc
      for(i in 1:(length(wfsub_trunc[1,])-1)){
        wfsub_plot[,i] = smoothTrajectory(wfsub[,i],halfW_plotting)
      }
    } else{wfsub_plot = wfsub}

      # Reshape wfsub for plotting
      wfsub_plot$index = 1:nrow(wfsub_plot)
      wfsub_plot = melt(wfsub_plot, id.vars=c('t', 'index'), variable.name='trialName')
      
      # Distinguish between colorwave 1/2 and 3/4
        for(k in sessionDF$trialName){
          wfsub_plot$wavetype[wfsub_plot$trialName == k] = sessionDF$wavetype[sessionDF$trialName == k]
        }

      # Get Time Intervals
      for(i in sessionDF$trialName){
        if(sessionDF$wavetype[sessionDF$trialName == i] == 1){
          wfsub_plot$'time color'[wfsub_plot$trialName == i] = ti[[1]][((halfW_plotting+1):(length(wfsub[,1])-halfW_plotting))]
        }
        else if(sessionDF$wavetype[sessionDF$trialName == i] == 2){ 
          wfsub_plot$'time color'[wfsub_plot$trialName == i] = ti[[2]][((halfW_plotting+1):(length(wfsub[,1])-halfW_plotting))]
        }
        else if(sessionDF$wavetype[sessionDF$trialName == i] == 3){
          wfsub_plot$'time color'[wfsub_plot$trialName == i] = ti[[3]][((halfW_plotting+1):(length(wfsub[,1])-halfW_plotting))]
        }
        else if(sessionDF$wavetype[sessionDF$trialName == i] == 4){
          wfsub_plot$'time color'[wfsub_plot$trialName == i] = ti[[4]][((halfW_plotting+1):(length(wfsub[,1])-halfW_plotting))]
        }}
      
      # Fan trajectories out for waterfall plotting
      wfsub_plot$shifted_value = round((wfsub_plot$value + (as.integer(wfsub_plot$trialName)-1) * 0.1),3)
      wfsub_plot$shifted_t = wfsub_plot$t - (as.integer(wfsub_plot$trialName)-1) *
        (400 / length(unique(wfsub_plot$trialName)))
      
      # Label regions for plot
      wfsub_plot$region = c('#000000'='initial',
                             '#0000FF'='LED',
                             '#00FF00'='TONE',
                             '#808080'='post-US')[wfsub_plot$'time color']
      for (i in 1:length(sessionDF$trialName)) {
        if(sessionDF$Peak[i]>0){
          pko1 = peakOnset(wfsub[[as.character(sessionDF$trialName[i])]],
                            sessionDF$Peak[i],
                            peakThreshold1,
                            peakThreshold2)
          if (!is.na(pko1)) {
            wfsub_plot[wfsub_plot$trialName == as.character(sessionDF[i, 'trial']) &
                          wfsub_plot$'time color' == '#808080' &
                          wfsub_plot$index <= sessionDF[i, 'peak'] &
                          wfsub_plot$index >= pko1, 'region'] = 'rising'
          }}}
      
      wfsub_plot$region = factor(wfsub_plot$region,
                                  levels = c('initial', 'LED','TONE', 'post-US', 'rising'))
      
      wfsub_plot[wfsub_plot$wavetype == 2 & wfsub_plot$region != 'LED', 'region'] = 'initial'
      wfsub_plot[wfsub_plot$wavetype == 4 & wfsub_plot$region != 'TONE', 'region'] = 'initial'
      
      if((doPostscriptWaterfall + doSessionWaterfall) > 0){
        # Full waterfall plot
        full_waterfall = ggplot(wfsub_plot, aes(x=shifted_t, y=shifted_value*1.5,
                                                color=region,group=trialName))
        full_waterfall = full_waterfall + theme_classic()
        full_waterfall = full_waterfall + theme(
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = 14)
        )
        full_waterfall = full_waterfall + geom_line(size=0.5) + theme(text=element_text(face="bold", size= 12))
        full_waterfall = full_waterfall + scale_color_manual(values=c(
          baseline_color, LED_color, TONE_color,airpuff_color))
        full_waterfall = full_waterfall + xlab('\n                   time (ms)') +
          scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100)))
        full_waterfall = full_waterfall + ggtitle((paste(substring(fileList[j],1,file_substring_end), '\n\n', '%CR Cue 1: ',CRPct_Cue1, '\nOnset Cue1: ',
                                                         CROnst_Cue1, ' ms\n', '%CR Cue 2: ', CRPct_Cue2, '\nOnset Cue2: ', CROnst_Cue2, ' ms' ))) 
        # Collapsed waterfall plots
        if(any(wfsub_plot$wavetype == 2)){
          probeTrials_LED = wfsub_plot[wfsub_plot$wavetype == 2,]
          Probe_LED =  ggplot(probeTrials_LED, aes(x=t, y=value,
                                                   color=region,group=trialName))
          Probe_LED = Probe_LED + theme_classic() + theme(legend.position = 'none')
          Probe_LED = Probe_LED + geom_line(size=0.4)
          Probe_LED = Probe_LED + scale_color_manual(values=c(
            baseline_color, LED_color ,airpuff_color))
          Probe_LED = Probe_LED + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else{Probe_LED = ggplot()}
        
        if(any(wfsub_plot$wavetype == 4)){
          probeTrials_TONE = wfsub_plot[wfsub_plot$wavetype == 4,]
          Probe_TONE =  ggplot(probeTrials_TONE, aes(x=t, y=value,
                                                     color=region,group=trialName))
          Probe_TONE = Probe_TONE + theme_classic() + theme(legend.position = 'none')
          Probe_TONE = Probe_TONE + geom_line(size=0.4)
          Probe_TONE = Probe_TONE + scale_color_manual(values=c(
            baseline_color, TONE_color,airpuff_color))
          Probe_TONE = Probe_TONE + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else{Probe_TONE = ggplot()}
        
        if(any(wfsub_plot$wavetype == 1)){
          pairedTrials_LED = wfsub_plot[wfsub_plot$wavetype == 1,]
          Paired_LED =  ggplot(pairedTrials_LED, aes(x=t, y=value,
                                                     color=region,group=trialName))
          Paired_LED = Paired_LED + theme_classic() + theme(legend.position = 'none')
          Paired_LED = Paired_LED + geom_line(size=0.4)
          Paired_LED = Paired_LED + scale_color_manual(values=c(
            baseline_color, LED_color,airpuff_color))
          Paired_LED = Paired_LED + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else {Paired_LED = ggplot()}
        
        if(any(wfsub_plot$wavetype == 3)){
          pairedTrials_TONE = wfsub_plot[wfsub_plot$wavetype == 3,]
          Paired_TONE =  ggplot(pairedTrials_TONE, aes(x=t, y=value,
                                                       color=region,group=trialName))
          Paired_TONE = Paired_TONE + theme_classic() + theme(legend.position = 'none')
          Paired_TONE = Paired_TONE + geom_line(size=0.4)
          Paired_TONE = Paired_TONE + scale_color_manual(values=c(
            baseline_color, TONE_color,airpuff_color))
          Paired_TONE = Paired_TONE + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
            theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
            scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
            scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100'))
        } else {Paired_TONE = ggplot()}
        
        
        if(doPostscriptWaterfall == 1){
          setwd(path_waterfall_postcript)
          setEPS()
          postscript((paste(substring(fileList[j],1,file_substring_end), '-wf.eps',sep = '')), height=13, width=10)
          capture.output(print(ggarrange(ggarrange(full_waterfall, nrow=1,ncol=1), ggarrange(Probe_LED,Paired_LED,Probe_TONE,Paired_TONE, nrow=2,ncol=2),nrow=2,ncol=1,heights=c(1.25,1))))
          dev.off()
        }
        if(doSessionWaterfall == 1){
          setwd(path_waterfall_png)
          png(paste(substring(fileList[j],1,file_substring_end), '-wf.png',sep = ''),
              width     = 10,
              height    = 13,
              units     = "in",
              res       = 150
          )
          capture.output(print(ggarrange(ggarrange(full_waterfall, nrow=1,ncol=1), ggarrange(Probe_LED,Paired_LED,Probe_TONE,Paired_TONE, nrow=2,ncol=2),nrow=2,ncol=1,heights=c(1.25,1))))
          dev.off()
        }
      }
      
       setwd(path_trajectories)
       wfsub_plot$Animal = substring(fileList[j],1,4)
       wfsub_plot$Geno = substring(fileList[j],1,1)
       wfsub_plot$session = substring(fileList[j],6,8)
       wfsub_plot$fullTrajectory = paste(wfsub_plot$trialName, wfsub_plot$session, sep = '')
       for(traj in levels(wfsub_plot$trialName)){
         if(any(traj %in% sessionDF$trialName) == TRUE){
           wfsub_plot$CR[wfsub_plot$trialName == traj] = sessionDF$CR[sessionDF$trialName == traj]
           wfsub_plot$Paired[wfsub_plot$trialName == traj] = sessionDF$Paired[sessionDF$trialName == traj]
         } else{
           wfsub_plot$CR[wfsub_plot$trialName == traj] = 0
           wfsub_plot$Paired[wfsub_plot$trialName == traj] = 5
         }
       }
       write.csv(wfsub_plot,paste(substring(fileList[j],1,file_substring_end),'-trajectories.txt', sep = ''))
       
       ## Go into each trace and snip out a 500 ms window surrounding response onset, so that they can be plotted aligned to onset.
       CR_df = sessionDF[sessionDF$Onset > 0,]
       wfsub_plot_onsetAligned = wfsub_plot[1,]
       wfsub_plot_onsetAligned$t_truc = 0
       wfsub_plot_onsetAligned = wfsub_plot_onsetAligned[0,]

       for(traj in sessionDF$trialName){
         trajDF = wfsub_plot[wfsub_plot$trialName == traj,]
         onset_t = CR_df$Onset[CR_df$trialName == traj] + baselineDur
         trajDF = trajDF[trajDF$t >= (onset_t - traj_offset),]
         if(length(trajDF$t) > traj_truc*(fps/1000)){
           trajDF = trajDF[1:(traj_truc*(fps/1000)),]
           trajDF$t_truc = seq((1000/fps),traj_truc,(1000/fps))
           wfsub_plot_onsetAligned = rbind(wfsub_plot_onsetAligned, trajDF)
         }
       }
       
       setwd(path_trajectories_aligned)
       write.csv(wfsub_plot_onsetAligned,paste(substring(fileList[j],1,file_substring_end),'-Aligned_trajectories.txt', sep = ''))
       setwd(path_trajectories)
      
    # Save dataframe, trialName for each raw file
     setwd(path_session_dataframes)
     write.csv(sessionDF, paste(substring(fileList[j],1,file_substring_end),'-DF.csv', sep = ''))

    setwd(path_waterfall_raw)
    if(j == 1){
      pct = round(((j)/length(fileList))*100,0)
      print(paste('           begin analysis: ', pct, '% Complete'))
    } else if(j%%80 == 0){
      pct = round(((j)/length(fileList))*100,0)
      print(paste('your wildest dreams await: ', pct, '% Complete'))
    } else if(j%%20 == 0){
      pct = round(((j)/length(fileList))*100,0)
      print(paste('     analyzing experiment: ', pct, '% Complete'))
    } else if(j%%5 == 0){
      pct = round(((j)/length(fileList))*100,0)
      print(paste('     ..................... ', pct, '% Complete'))
    }
  } # End bracket that reads in individual raw data files  
}


# CREATE & FILL IN DATAFRAMES FOR WHOLE EXPERIMENT ANALYSIS ------------------------

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

for(cue in unique(sessionDF$Cue)){
  for(cueisi in unique(substring(fileList,10,13))){
    # Rename rows and columns in dataframe for paired trials
    df = data.frame(matrix(0, numAnimals, numSessions + 3))
    names(df)[1:3] = c('Animal', 'Group', 'Cue')
    df$Animal = animalNames
    df$Group = c(substring(animalNames,1,1))
    df$Cue = cue
    # Make entries for Session CR rates per animal
    for (i in 4:(numSessions + 3)) {
      names(df)[i] <- sessionNames[i-3]
    }
    # # Make entries for Session UR-Integrate values per animal
    # for (i in (numSessions + 4):(2*numSessions + 3)) {
    #   names(df)[i] <- paste(sessionNames[(i-numSessions-3)],'_UR', sep='')
    # }
    setwd(path_summary_csv)
    write.csv(df, paste(cueisi,'-', cue, '-', substring(dataFolder, 1, nchar(dataFolder)), '.csv', sep = '')) 
    setwd(path_waterfall_raw)
  }}

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
  write.csv(startleDF, 'StartleSummary.csv')  
}


# FILL DATAFRAME FOR WHOLE EXPERIMENT ANALYSIS 

setwd(path_summary_csv)
csvList = list.files()

for(cue in unique(sessionDF$Cue)){
for(cueisi in unique(substring(csvList,1,4))){
  if(cueisi != 'Star'){
    file = intersect(grep(cueisi, csvList), grep(cue, csvList)) # Finds unique file name that is for cueisi & cue
    df = read.csv(csvList[file]) # then opens the partially filled-in csv
    for (animal in unique(substring(fileList,1,4))) {
      setwd(path_session_dataframes)
      fileList = list.files()
      temp = read.csv(fileList[1])
      temp = temp[0,]
      for(i in 1:length(fileList)){
        if(substring(fileList[i],1,4) == animal){
          csv = read.csv(fileList[i]) 
          csv = csv[csv$Cue == cue,] # Only keep trials with the relevant cue
          
          # If the temp session has startle probes, exclude them
          if(substring(fileList[i],file_substring_end,file_substring_end) > 0){
            csv = csv[csv$Paired == 1,]
          }
          temp = rbind(temp,csv)
        }}
      
      # Do only for non-discrimination sessions
      # if(!Discrim_Session){
      #   temp = temp[temp$ISI == substring(cueisi,2,4),]
      # }
      
      # temp = temp[temp$ID == 0,] # gets rid of startle sessions. 
      
      # save total CSV
      setwd(path_animal_dataframes)
      write.csv(temp, paste(animal, '-allSessions.csv', sep=''))
      
      setwd(path_session_dataframes)
      temp_probe = temp[temp$Paired == 0,]
      temp_paired = temp[temp$Paired == 1,]
      
      # Calculate learning curve for each animal
      
      if(length(temp[,1]) > 0){
        
        CR_vector = temp$CR
        temp_CR = temp[temp$Onset>0,]
        df$totalCR[df$Animal == animal] = (sum(CR_vector)/length(CR_vector))
        df$totalUR[df$Animal == animal] = median(temp_paired$UR_Int[temp_paired$Animal == animal])
        df$numTrials[df$Animal== animal] = length(temp[,2])
        df$countCRs[df$Animal== animal] = sum(temp_CR$CR)
        
        # Fill in %CR per day
        for(Session in unique(temp$Session)){
          tempSession = temp[temp$Session == Session,]
          df[[Session]][df$Animal == animal] = round(100*sum(tempSession$CR)/length(tempSession$CR), 1)
        }
        
         X = substring(temp$Session,1,1)[1]
         tempSession1 = temp[temp$Session == paste(X, "01", sep=''),]
         tempSession2 = temp[temp$Session == paste(X, "02", sep=''),]
         tempSession = rbind(tempSession1, tempSession2)
            
         if(length(temp[,1]) > 0){
            #Calculate bootstrapped UR for first 2 days of training:
            # bootCoeff_UR <- numeric(numberOfIterations)
            # bootCoeff_UR_Int <- numeric(numberOfIterations)
            # for (n in 1:numberOfIterations){
            #   UR = tempSession$UR
            #   UR_Int = tempSession$UR_Int
            #   bootDataset_UR <- sample(UR, length(UR), replace = TRUE, prob=NULL)
            #   bootCoeff_UR[n] <- median(bootDataset_UR)
            #   bootDataset_UR_Int <- sample(UR_Int, length(UR_Int), replace = TRUE, prob=NULL)
            #   bootCoeff_UR_Int[n] <- median(bootDataset_UR_Int)
            # }
            # UR_boot = round(as.numeric(quantile(bootCoeff_UR, .500)),3) # median
            # UR_median = median(UR)
            # UR_Int_boot = round(as.numeric(quantile(bootCoeff_UR_Int, .500)),3) # median
            # UR_Int_median = median(UR_Int)
            #df$UR_boot[df$Animal==animal] = UR_boot
            df$UR_median[df$Animal==animal] = median(tempSession$UR)
            #df$UR_Int_boot[df$Animal==animal] = UR_Int_boot
            df$UR_Int_median[df$Animal==animal] = median(tempSession$UR_Int)
            df$SlopeCutoff[df$Animal==animal] = length(temp_CR$Slope[temp_CR$Slope > deltaSlope_cutoff])
          } else{
            #df$UR_boot[df$Animal==animal] = 0
            df$UR_median[df$Animal==animal] = 0
            #df$UR_Int_boot[df$Animal==animal] = 0
            df$UR_Int_median[df$Animal==animal] = 0
            df$SlopeCutoff[df$Animal==animal] = 0
          }
          
          # Find number of trials it takes each animal to reach 25% CR
          sliding = 0
          if(length(CR_vector) > sliding_window_size*2){
          for(i in 1:((length(CR_vector)-sliding_window_size))){
            sliding[i] = sum(CR_vector[i:(i+(sliding_window_size-1))])
          }
          trial25 = which(sliding >= round(sliding_window_size*.25,0))[1]
          if(is.na(trial25) == 1){
            trial25 = 0
            df$trial25[df$Animal== animal] = 0
          } else{
            df$trial25[df$Animal== animal] = trial25 + sliding_window_size-1
          }} else{df$trial25[df$Animal== animal] = 0}
  
          # Analyze paired trials for onset, amp at 100 ms after onset
         if(restrict == 1){
            if (length(temp_CR$Onset) >= trajRestrict[2]){
              temp1 = temp_CR[trajRestrict[1]:trajRestrict[2],]
              bootCoeff <- numeric(numberOfIterations)
              for (n in 1:numberOfIterations){
                onset = temp1$Onset
                bootDataset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
                bootCoeff[n] <- median(bootDataset)
              }
              onset_lower = round(as.numeric(quantile(bootCoeff, lowerCI)),1) # lowerCI 
              onset_upper = round(as.numeric(quantile(bootCoeff, upperCI)),1) # upperCI
              
              df$Onset1_median[df$Animal==animal] = median(onset)
              df$Onset1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff, .500)),1) # median
              df$Onset1_ConfInt[df$Animal==animal] = onset_upper - onset_lower
              df$Onset1_Var[df$Animal==animal] = round(var(onset),0)
              
              # Now do slope 
              bootCoeff <- numeric(numberOfIterations)
              for (n in 1:numberOfIterations){
                slope = temp1$Slope
                bootDataset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
                bootCoeff[n] <- median(bootDataset)
              }
              slope_lower = round(as.numeric(quantile(bootCoeff, lowerCI)),1) # lowerCI 
              slope_upper = round(as.numeric(quantile(bootCoeff, upperCI)),1) # upperCI
              
              df$Slope1_median[df$Animal==animal] = median(slope)
              df$Slope1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff, .500)),1) # median
              df$Slope1_ConfInt[df$Animal==animal] = slope_upper - slope_lower
              df$Slope1_Var[df$Animal==animal] = round(var(slope),0)
              
            } 
           else{
              df$Onset1_median[df$Animal==animal] = 0
              df$Onset1_boot[df$Animal==animal] = 0
              df$Onset1_ConfInt[df$Animal==animal] = 0
              df$Onset1_Var[df$Animal==animal] = 0
              df$Slope1_median[df$Animal==animal] = 0
              df$Slope1_boot[df$Animal==animal] = 0
              df$Slope1_ConfInt[df$Animal==animal] = 0
              df$Slope1_Var[df$Animal==animal] = 0
                }
      } # End of if(restrict == 1) statement
      else if(restrict == 0){ # restict = 0, do all CRs
        if(length(temp_CR[,1]) > trajMin){
          temp1 = temp_CR[1:length(temp_CR[,1]),]
          bootCoeff <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            onset = temp1$Onset
            bootDataset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
            bootCoeff[n] <- median(bootDataset)
          }
          onset_lower = round(as.numeric(quantile(bootCoeff, lowerCI)),1) # lowerCI 
          onset_upper = round(as.numeric(quantile(bootCoeff, upperCI)),1) # upperCI
          
          df$Onset1_median[df$Animal==animal] = median(onset)
          df$Onset1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff, .500)),1) # median
          df$Onset1_ConfInt[df$Animal==animal] = onset_upper - onset_lower
          df$Onset1_Var[df$Animal==animal] = round(var(onset),0)
          
          # Now do slope 
          bootCoeff <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            slope = temp1$Slope
            bootDataset <- sample(onset, length(onset), replace = TRUE, prob=NULL)
            bootCoeff[n] <- median(bootDataset)
          }
          slope_lower = round(as.numeric(quantile(bootCoeff, lowerCI)),1) # lowerCI 
          slope_upper = round(as.numeric(quantile(bootCoeff, upperCI)),1) # upperCI
          
          df$Slope1_median[df$Animal==animal] = median(slope)
          df$Slope1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff, .500)),1) # median
          df$Slope1_ConfInt[df$Animal==animal] = slope_upper - slope_lower
          df$Slope1_Var[df$Animal==animal] = round(var(slope),0)
        } # end  if(length(temp_CR[,1]) > trajMin)
        else{ 
          df$Onset1_median[df$Animal==animal] = 0
          df$Onset1_boot[df$Animal==animal] = 0
          df$Onset1_ConfInt[df$Animal==animal] = 0
          df$Onset1_Var[df$Animal==animal] = 0
          df$Slope1_median[df$Animal==animal] = 0
          df$Slope1_boot[df$Animal==animal] = 0
          df$Slope1_ConfInt[df$Animal==animal] = 0
          df$Slope1_Var[df$Animal==animal] = 0
        }
      } # End else if(restrict == 0){
        
        # take median of all responses
        if(length(temp_CR$Onset >= trajMin)){
          df$allOnset[df$Animal==animal] = median(temp_CR$Onset)
          df$allSlope[df$Animal==animal] = median(temp_CR$Slope)
        }

        # Plot all onsets and slopes across trials as well as their histograms
        if(length(temp_CR) > 1){ 
          temp_CR$trialNumFull = as.numeric(substring(temp_CR$Session,2,3))*temp_CR$trialNum
          OnsetxTrial_point = ggplot(temp_CR, aes(x = trialNumFull, y = Onset, group = Animal, color = Animal)) + 
            geom_point(alpha = 0.5) + theme_classic() + xlim(0,600) + ylim(0,600)
          OnsetxTrial_histogram = ggplot(temp_CR, aes(x = Onset, group = Animal, color = Animal)) + 
            geom_histogram(binwidth = 10, alpha = 0.5) + theme_classic() + xlim(0,600) + ylim(0,55)
          SlopexTrial_point = ggplot(temp_CR, aes(x = trialNumFull, y = Slope, group = Animal, color = Animal)) + 
            geom_point(alpha = 0.5) + theme_classic() + xlim(0,600) + ylim(0,1)
          SlopexTrial_histogram = ggplot(temp_CR, aes(x = Slope, group = Animal, color = Animal)) + 
            geom_histogram(binwidth = .02, alpha = 0.5) + theme_classic() + xlim(0,1) + ylim(0,55)
          
          if(doTimingPlots == 1){
            setwd(path_timing_plots)
            pdf(paste(temp_CR$Animal, '-', cueisi, '-TimingxTrial.pdf',sep = ''), height=6, width=10)
            capture.output(print(ggarrange(OnsetxTrial_point,OnsetxTrial_histogram,SlopexTrial_point,SlopexTrial_histogram,nrow = 2, ncol = 2)))
            dev.off()
            setwd(path_summary_csv)
          }
        }
        
      # Now calculate stuff for probe trials
    
     if(length(temp_probe[,1]) > 0){
         df$numTrials_probe[df$Animal==animal] = length(temp_probe[,2])
         temp_probeCR = temp_probe[temp_probe$Peak>0,]
         df$countCRs_probe[df$Animal==animal] = sum(temp_probeCR$CR)  
         
         if(restrict == 1){
          if (length(temp_probeCR$Peak) >= trajRestrict_Probe[2] && length(temp_CR$Onset) >= trajRestrict[2]){
            temp_probeCR1 = temp_probeCR[trajRestrict_Probe[1]:trajRestrict_Probe[2],]
            bootCoeff_peak <- numeric(numberOfIterations)
            bootCoeff_onset <- numeric(numberOfIterations)
            bootCoeff_amp <- numeric(numberOfIterations)
            bootCoeff_ampInt <- numeric(numberOfIterations)
            for (n in 1:numberOfIterations){
              peak = temp_probeCR1$Peak
              ampProbe = temp_probeCR1$maxAmp
              ampInt = temp_probeCR1$ampInt
              bootDataset_peak <- sample(peak, length(peak), replace = TRUE, prob=NULL)
              bootCoeff_peak[n] <- median(bootDataset_peak)
              bootDataset_ampProbe <- sample(ampProbe, length(ampProbe), replace = TRUE, prob=NULL)
              bootCoeff_amp[n] <- median(bootDataset_ampProbe)
              bootDataset_ampInt <- sample(ampInt, length(ampInt), replace = TRUE, prob=NULL)
              bootCoeff_ampInt[n] <- median(bootDataset_ampInt)
            }
            
            peak_lower = round(as.numeric(quantile(bootCoeff_peak, lowerCI)),1) # lowerCI
            peak_upper = round(as.numeric(quantile(bootCoeff_peak, upperCI)),1) # upperCI
            
            amp_lower = round(as.numeric(quantile(bootCoeff_amp, lowerCI)),3) # lowerCI
            amp_upper = round(as.numeric(quantile(bootCoeff_amp, upperCI)),3) # upperCI
            
            ampInt_lower = round(as.numeric(quantile(bootCoeff_ampInt, lowerCI)),3) # lowerCI
            ampInt_upper = round(as.numeric(quantile(bootCoeff_ampInt, upperCI)),3) # upperCI
            
            df$Peak1_median[df$Animal==animal] = median(peak)
            df$Peak1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff_peak, .500)),1) # median
            df$PeakCI[df$Animal==animal] = peak_upper - peak_lower
            df$PeakVar[df$Animal==animal] = round(var(peak),0)
            
            df$ampProbe1_median[df$Animal==animal] = median(ampProbe)
            df$ampProbe1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff_amp, .500)),3) # median
            df$ampProbeCI[df$Animal==animal] = amp_upper - amp_lower
            df$ampProbeVar[df$Animal==animal] = round(var(ampProbe),3)
            
            df$ampInt1_median[df$Animal==animal] = median(ampInt)
            df$ampInt1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff_ampInt, .500)),3) # Area Under CR
            df$ampIntCI[df$Animal==animal] = ampInt_upper - ampInt_lower 
            df$ampIntVar[df$Animal==animal] = round(var(ampInt),1)
            
          } else{
            df$Peak1_median[df$Animal==animal] = 0
            df$Peak1_boot[df$Animal==animal] = 0
            df$PeakCI[df$Animal==animal] = 0
            df$PeakVar[df$Animal==animal] = 0
            df$ampProbe1_median[df$Animal==animal] = 0
            df$ampProbe1_boot[df$Animal==animal] = 0
            df$ampProbeCI[df$Animal==animal] = 0 
            df$ampProbeVar[df$Animal==animal] = 0
            df$ampInt1_median[df$Animal==animal] = 0
            df$ampInt1_boot[df$Animal==animal] = 0
            df$ampIntCI[df$Animal==animal] = 0
            df$ampIntVar[df$Animal==animal] = 0
          }
        } # end if(restrict == 1) statement
         else if(restrict == 0){ #  Restrict = 0, not restricted, analyze all CRs
           if (length(temp_CR$Onset) > trajMin){
             temp_probeCR1 = temp_probeCR[1:length(temp_probeCR[,1]),]
             temp_probeCR1 = temp_probeCR
             bootCoeff_peak <- numeric(numberOfIterations)
             bootCoeff_onset <- numeric(numberOfIterations)
             bootCoeff_amp <- numeric(numberOfIterations)
             bootCoeff_ampInt <- numeric(numberOfIterations)
             for (n in 1:numberOfIterations){
               peak = temp_probeCR1$Peak
               ampProbe = temp_probeCR1$maxAmp
               ampInt = temp_probeCR1$ampInt
               bootDataset_peak <- sample(peak, length(peak), replace = TRUE, prob=NULL)
               bootCoeff_peak[n] <- median(bootDataset_peak)
               bootDataset_ampProbe <- sample(ampProbe, length(ampProbe), replace = TRUE, prob=NULL)
               bootCoeff_amp[n] <- median(bootDataset_ampProbe)
               bootDataset_ampInt <- sample(ampInt, length(ampInt), replace = TRUE, prob=NULL)
               bootCoeff_ampInt[n] <- median(bootDataset_ampInt)
             }
           
             peak_lower = round(as.numeric(quantile(bootCoeff_peak, lowerCI)),1) # lowerCI
             peak_upper = round(as.numeric(quantile(bootCoeff_peak, upperCI)),1) # upperCI
             amp_lower = round(as.numeric(quantile(bootCoeff_amp, lowerCI)),3) # lowerCI
             amp_upper = round(as.numeric(quantile(bootCoeff_amp, upperCI)),3) # upperCI
             ampInt_lower = round(as.numeric(quantile(bootCoeff_ampInt, lowerCI)),3) # lowerCI
             ampInt_upper = round(as.numeric(quantile(bootCoeff_ampInt, upperCI)),3) # upperCI
             
             df$Peak1_median[df$Animal==animal] = median(peak)
             df$Peak1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff_peak, .500)),1) # median
             df$PeakCI[df$Animal==animal] = peak_upper - peak_lower
             df$PeakVar[df$Animal==animal] = round(var(peak),0)
             
             df$ampProbe1_median[df$Animal==animal] = median(ampProbe)
             df$ampProbe1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff_amp, .500)),3) # median
             df$ampProbeCI[df$Animal==animal] = amp_upper - amp_lower
             df$ampProbeVar[df$Animal==animal] = round(var(ampProbe),3)
             
             df$ampInt1_median[df$Animal==animal] = median(ampInt)
             df$ampInt1_boot[df$Animal==animal] = round(as.numeric(quantile(bootCoeff_ampInt, .500)),3) # Area Under CR
             df$ampIntCI[df$Animal==animal] = ampInt_upper - ampInt_lower 
             df$ampIntVar[df$Animal==animal] = round(var(ampInt),1)
             
           } else{ # there aren't trajMin CRs
             df$Peak1_median[df$Animal==animal] = 0
             df$Peak1_boot[df$Animal==animal] = 0
             df$PeakCI[df$Animal==animal] = 0
             df$PeakVar[df$Animal==animal] = 0
             df$ampProbe1_median[df$Animal==animal] = 0
             df$ampProbe1_boot[df$Animal==animal] = 0
             df$ampProbeCI[df$Animal==animal] = 0 
             df$ampProbeVar[df$Animal==animal] = 0
             df$ampInt1_median[df$Animal==animal] = 0
             df$ampInt1_boot[df$Animal==animal] = 0
             df$ampIntCI[df$Animal==animal] = 0
             df$ampIntVar[df$Animal==animal] = 0
           } 
         } # end  else if(restrict == 0){
      } # end if(length(temp_probe[,1]) > 0) statement
   else{temp_probeCR = NULL} 
        
       # Now calculate %CR on last two days for all animals (then manually exclude based on other criteria if desired)
        last2_SessionNames = sessionNames[(length(sessionNames)-1):length(sessionNames)]
        last2_dat = temp[0,]
        
        if(length(temp$Onset[temp$Onset>0]) >= trajMin){ # can change this number to restrict analysis
            for(i in 1:length(last2_SessionNames)){
              hold = temp[temp$Session == last2_SessionNames[i],]
              last2_dat = rbind(last2_dat, hold)
            }
              df$final2[df$Animal == animal] = round(sum(last2_dat$CR)/length(last2_dat$CR),3)
        } else{df$final2[df$Animal == animal] = 0}
        
    } # Ends if temp[,1] is greater than 0 statement
     else{temp_CR = NULL} # End if(length(temp[,1]) > 0) Statement
    } # End of for loop that reads in each animal
    
    df[is.na(df)] = 0
    setwd(path_summary_csv)
    write.csv(df, csvList[file])
  } # End of if cueisi != 'Star' statement
} # End of loop that reads in each Summary csv file
} # End loop that does a separate analysis for each Cue

# ANALYZE RAW FILES FOR STARTLE  ------------------------------------------

if(doStartle == 1){
  allPeaks = function(traj, halfWidth = halfW_analysis, minPeak= 0.05) {
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
        traceStart = round((baselineDur*(fps/1000)),0)   # moves the analysis window to CS onset. outputs the frame number where this happens
        traceEnd = round(((baselineDur + analysisWindow_startle)*(fps/1000)),0)  # moves the end of analysis window to analysisWindow_startle ms after startle tone
        wf_startle = wf_startle[1:traceEnd,]
        wfsub_startle = subtractBaselines(wf_startle)
        
        # Get rid of all trials that aren't probe trials (ending in either a 0 or a 5)
        temp1 = wfsub_startle[grep('0$', colnames(wfsub_startle))]
        temp2 = wfsub_startle[grep('5$', colnames(wfsub_startle))]
        wfsub_startle = cbind(temp1,temp2,'t' = wfsub_startle$t)
        wfsub_startle[is.na(wfsub_startle)] = 0
        wfsub_startle = wfsub_startle[ , c(basenoise(wfsub_startle) < baseline_noise, TRUE), drop=FALSE]
        
        # Locate highest peak in each trajectory
        peakTiming_probe = topPeaks(wfsub_startle[traceStart:(traceEnd),], n=1, halfWidth = halfW_analysis)
        peakTiming_probe[is.na(peakTiming_probe)] = 0
        
        if(length(peakTiming_probe) > 0){
          sessionDF = data.frame(Animal = substring(fileList[j],1,4), Group = substring(fileList[j],1,1), Session = substring(fileList[j],6,8),
                                 Trial = as.numeric(substring(colnames(peakTiming_probe)[1:(length(peakTiming_probe)-1)],7,9)),
                                 Paired = 0, Startle = 0, Onset = 0, Peak = 0, maxAmp = 0, ampInt = 0, ID = substring(fileList[j],file_substring_end,file_substring_end),
                                 trialName=colnames(peakTiming_probe)[1:(length(peakTiming_probe)-1)], Cue = as.character(substring(fileList[j],10,10)), ISI = as.numeric(substring(fileList[j],11,13)))
          # Fill in data frame
          for(k in 1:length((sessionDF$Trial))) {
            if(peakTiming_probe[1,k] > 0){ # If there is a startle 
              sessionDF$Peak[k] = topPeaks(wfsub_startle[traceStart:(traceEnd),k], n=1, halfWidth = halfW_analysis)*(1000/fps)
              sessionDF$Onset[k] = peakOnset(wfsub_startle[traceStart:(traceEnd),k], sessionDF$Peak[k]*(fps/1000), 
                                             peakThreshold_startle, peakThreshold2)*(1000/fps)
              sessionDF$Startle[k] = 1 
              sessionDF$maxAmp[k] = maxAmp(wfsub_startle[traceStart:(traceEnd),k])
              sessionDF$ampInt[k] = integrateArea(wfsub[,k], traceStart, traceEnd, fps, peakThreshold_startle)
            } else{}
          }
        }
        
        ## STARTLE WATERFALL PLOTS 
        
        if(doSessionWaterfall == 1) {
          # Get time intervals from fileList[j] file
          ti = readTimeIntervals(fileList[j])
          ti[[1]] = ti[[1]][1:traceEnd]
          names(ti[[1]]) = as.character(wf_startle$t)[1:traceEnd]
          
          # Reshape wfsub_startle for plotting
          wfmelt_startle = wfsub_startle
          wfmelt_startle$index = 1:nrow(wf_startle)
          wfmelt_startle = wfmelt_startle[(100*(fps/1000)):traceEnd,]
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
              pko1 = peakOnset(wfsub_startle[[as.character(sessionDF$trialName[i])]],
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
          if(doPostscriptWaterfall == 1){
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
    startleDF$maxAmp[startleDF$Animal==animal] = round(median(temp$maxAmp),3)
    
    # Save Summary CSV for Startle Trials
    setwd(path_summary_csv)
    write.csv(startleDF, paste('StartleSummary.csv', sep = '')) 
    setwd(path_session_dataframes_startle)
  }
}

# PLOTS WITH STATISTICAL TESTS ---------------------------------------

setwd(path_summary_csv)
csvList = list.files()

# Read in all CSV's and pull out unique cues 
cues = NULL

 for(i in 1:length(csvList)){
   df = read.csv(csvList[i])
   cues = c(cues, as.character(unique(df$Cue)))
 }

for(cue in unique(cues)){
  for(cueisi in unique(substring(csvList,1,4))){
    # If summary csv is anything but the StartleSummary.csv
     if(cueisi != 'Star'){
      file = intersect(grep(cueisi, csvList), grep(cue, csvList)) # Finds unique file name that is for cueisi & cue
      df = read.csv(csvList[file])
      
      # K-S test for all comparisons. 
      Onset1 = df[df$Onset1_median>0,]
      maxAmp1 = Onset1[Onset1$ampProbe1_median>0,]
      ampInt1 = Onset1[Onset1$ampInt1_median>0,]
      totalUR = df
      totalCR = df
      slopeDF = df[df$Slope1_median>0,]
      trial25 = df[df$trial25 > 0,]
      Peak1 = df[df$Peak1_median>0,]
      Group = unique(df$Group)
      final2 = df[df$final2>0,]
      
      if(length(Group) == 2){
        group_colors = c(group1_color, group2_color)
      } else if(length(Group) == 3){
        group_colors = c(group1_color, group2_color, group3_color)
      }
      
      # melts dataframe & deletes unncessary info so that only %CR data across sessions remains
      # This dataframe will be used to plot learning curve (median %CR over training)
      df_learningCurve = melt(df, id=c('Animal', 'Group', 'Cue'))
      names(df_learningCurve)[4] = 'TrainingSession'
      names(df_learningCurve)[5] = 'ConditionedResponse'
      toDelete = NULL
      for(i in 1:length(df_learningCurve$TrainingSession)){
        if(nchar(as.character(df_learningCurve$TrainingSession[i])) > 3 || nchar(as.character(df_learningCurve$TrainingSession[i])) < 3){
          toDelete = c(toDelete,i)
        }}
      df_learningCurve = df_learningCurve[-c(toDelete),]
      df_learningCurve = df_learningCurve[substring(df_learningCurve$TrainingSession,2,2)!= '.',] # Gets rid of annoying X.1 column artifact

      # Create Group Average Dataframe
      groupMedian = df_learningCurve
      groupMedian = groupMedian[0,2:4]
      for(group in Group){
        temp1 = df_learningCurve[df_learningCurve$Group == group,]
        for(session in unique(df_learningCurve$TrainingSession)){
          temp2 = temp1[temp1$TrainingSession == session,]
          medianCR = median(temp2$ConditionedResponse)
          groupMedian = rbind(groupMedian, data.frame('Group' = group, 'TrainingSession' = session, 'ConditionedResponse' = medianCR))
        }
      }
      
      # Perform statistical tests on data and plot
      
      if(length(df_learningCurve[,1]) > 0){
        cr_plot = ggplot(df_learningCurve, aes(x = TrainingSession, y= ConditionedResponse, group = Animal, color = Group))  + theme_classic() +
          geom_line(size=1, alpha= .25) + theme(legend.position="none") + scale_color_manual(values = group_colors) + ylim(0,100) +
          geom_line(data = groupMedian, aes(x=TrainingSession, y=ConditionedResponse, group = Group, color= Group), size = 2)
      } else{cr_plot = ggplot()}
      
      
      group1 = df$totalCR[df$Group == Group[1]]
      group2 = df$totalCR[df$Group == Group[2]]
      if(length(Group) > 2){
        group3 = df$totalCR[df$Group == Group[3]]
      } else{ group3 = 1:10}
      if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
        totalCR_test = kruskal_test(totalCR ~ Group, data = df, distrubution = 'exact')
        totalCR_boxplot = ggplot(df, aes(x = Group, y = totalCR, group = Group)) + theme_classic() +
          ylab('Total % CR') +
          geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
          geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(0, 0.5) + 
          ggtitle(paste('ks-test, p = ', round(pvalue(totalCR_test),3)))
      }else{
        totalCR_test = 'n/a'
        totalCR_boxplot = ggplot()
      }
      
      group1 = df$totalUR[df$Group == Group[1]]
      group2 = df$totalUR[df$Group == Group[2]]
      if(length(Group) > 2){
        group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
      } else{ group3 = 1:10}
      if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
        ur_test = kruskal_test(totalUR ~ Group, data = df, distrubution = 'exact')
        totalUR = ggplot(df, aes(x = Group, y = totalUR, group = Group)) + theme_classic() + ylab("UR Magnitude \n (all URs)") + 
          geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + ylim(25,175) + 
          geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
          ggtitle(paste('ks-test, p = ', round(pvalue(ur_test),4)))
      }else{
        ur_test = 'n/a'
        totalUR = ggplot()
      }
      
      group1 = trial25$trial25[trial25$Group == Group[1]]
      group2 = trial25$trial25[trial25$Group == Group[2]]
      if(length(Group) > 2){
        group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
      } else{ group3 = 1:10}
      if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
        trial25_test = kruskal_test(trial25 ~ Group, data = trial25, distrubution = 'exact')
        t25 = ggplot(trial25, aes(x = Group, y = trial25, group = Group)) + theme_classic() + ylim(0,1200) + 
          geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
          geom_point(color = 'black', shape=1, size=2, na.rm = TRUE) + ylab("Trials to 25% CR") + 
          ggtitle(paste('      ks-test, p =', round(pvalue(trial25_test),3))) + theme(plot.title = element_text(hjust = -.75))
      }else{
        trial25_test = 'n/a'
        t25 = ggplot()
      }
      
      group1 = final2$final2[final2$Group == Group[1]]
      group2 = final2$final2[final2$Group == Group[2]]
      if(length(Group) > 2){
        group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
      } else{ group3 = 1:10}
      if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
        cr_test = kruskal_test(final2 ~ Group, data = final2, distrubution = 'exact')
        last2_boxplot = ggplot(final2, aes(x = Group, y = final2, group = Group)) + theme_classic() +
          ylab('%CR Plateau \n Final 2 Days') + 
          geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
          geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
          ylim(0,1) + ggtitle(paste('ks-test, p = ', round(pvalue(cr_test),4)))
      }else{
        cr_test = 'n/a'
        last2_boxplot = ggplot()
      }
      
      # Plot differently depending on whether the data is restricted or not
      
      if(restrict == 1){
        group1 = slopeDF$Slope1_median[slopeDF$Group == Group[1]]
        group2 = slopeDF$Slope1_median[slopeDF$Group == Group[2]]
        if(length(Group) > 2){
          group3 = slopeDF$Slope1_median[slopeDF$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          slope_test = kruskal_test(Slope1_median ~ Group, data = slopeDF, distrubution = 'exact')
          slope_boxplot = ggplot(slopeDF, aes(x = Group, y = Slope1_median, group = Group)) + theme_classic() +
            ylab(paste('amp at ', deltaSlope_ms, ' ms\n', 'trials ',trajRestrict[1], '-', trajRestrict[2], sep='')) +
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
            ggtitle(paste('ks-test, p = ', round(pvalue(slope_test),3)))
        }else{
          slope_test = 'n/a'
          slope_boxplot = ggplot()
        }
        
        group1 = Onset1$Onset1_median[Onset1$Group == Group[1]]
        group2 = Onset1$Onset1_median[Onset1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          op1_test = kruskal_test(Onset1_median ~ Group, data = Onset1, distrubution = 'exact')
          op1 = ggplot(Onset1, aes(x = Group, y = Onset1_median, group = Group)) + theme_classic() +
            ylab(paste("CR Onset (ms)\n", 'trials ',trajRestrict[1], '-', trajRestrict[2], sep='')) +
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(75,400) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(op1_test),3)))
        }else{
          op1_test = 'n/a'
          op1 = ggplot()
        }
        
        group1 = ampInt1$ampInt1_median[ampInt1$Group == Group[1]]
        group2 = ampInt1$ampInt1_median[ampInt1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          ampInt1_test = kruskal_test(ampInt1_median ~ Group, data = ampInt1, distrubution = 'exact')
          ampInt1_boxplot = ggplot(ampInt1, aes(x = Group, y = ampInt1_median, group = Group)) + theme_classic() +
            ylab(paste("CR Magnitude\n", 'probe ',trajRestrict_Probe[1], '-', trajRestrict_Probe[2], sep='')) +
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE) + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(0,100) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(ampInt1_test),3)))
        }else{
          ampInt1_test = 'n/a'
          ampInt1_boxplot = ggplot()
        }
        
        group1 = maxAmp1$ampProbe1_median[maxAmp1$Group == Group[1]]
        group2 = maxAmp1$ampProbe1_median[maxAmp1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          maxAmp1_test = kruskal_test(ampProbe1_median ~ Group, data = maxAmp1, distrubution = 'exact')
          maxAmp1_boxplot = ggplot(maxAmp1, aes(x = Group, y = ampProbe1_median, group = Group)) + theme_classic() +
            ylab(paste("CR Amplitude\n", 'probe ',trajRestrict_Probe[1], '-', trajRestrict_Probe[2], sep='')) +
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE) + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(0,1) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(maxAmp1_test),3)))
        }else{
          maxAmp1_test = 'n/a'
          maxAmp1_boxplot = ggplot()
        }
        
        group1 = Peak1$Peak1_median[Peak1$Group == Group[1]]
        group2 = Peak1$Peak1_median[Peak1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          pk1_test = kruskal_test(Peak1_median ~ Group, data = Peak1, distrubution = 'exact')
          pk1 = ggplot(Peak1, aes(x = Group, y = Peak1_median, group = Group)) + theme_classic() + 
            ylab(paste("Peak Timing (ms)\n", 'probe ', trajRestrict_Probe[1], '-', trajRestrict_Probe[2], sep='')) +
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(250,(500)) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(pk1_test),3)))
        }else{
          pk1_test = 'n/a'
          pk1 = ggplot()
        }
        
        if(length(ampInt1[,1]>0)){
          if(length(Onset1[,1]>0)){
            ampInt1_scatterplot = ggplot(ampInt1, aes(x = Onset1_median, y= ampInt1_median, group = Group, color = Group)) + theme_classic() +
              geom_point(size = .75, alpha = .75) + xlim(75,350) + ylim(0,100) + ylab(paste("CR Magnitude \n", 'probe ',trajRestrict_Probe[1], '-', trajRestrict_Probe[2], sep='')) + 
              xlab(paste("CR Onset (ms) \n", 'trials ', trajRestrict[1], '-', trajRestrict[2], sep='')) + 
              theme(legend.position="none") + scale_color_manual(values = group_colors)
          }}else{ampInt1_scatterplot = ggplot()}
        
        if(length(maxAmp1[,1]>0)){
          if(length(Onset1[,1]>0)){
            maxAmp1_scatterplot = ggplot(maxAmp1, aes(x = Onset1_median, y= ampProbe1_median, group = Group, color = Group)) + theme_classic() +
              geom_point(size = .75, alpha = .75) + xlim(75,350) + ylim(0,1) + ylab(paste("CR Amplitude \n",'probe ', trajRestrict_Probe[1],'-', trajRestrict_Probe[2], sep='')) + 
              xlab(paste("CR Onset (ms) \n", 'trials ', trajRestrict[1], '-',trajRestrict[2], sep='')) + 
              theme(legend.position="none") + scale_color_manual(values = group_colors)
          }}else{maxAmp1_scatterplot = ggplot()}
        
        if(length(slopeDF[,1]>0)){
          if(length(Onset1[,1]>0)){
            slopeXonset_scatterplot = ggplot(Onset1, aes(x = Onset1_median, y= Slope, group = Group, color = Group)) + theme_classic() +
              geom_point(size = .75, alpha = .75) + xlim(75,350) + ylim(0,1) + ylab(paste("Amp 100ms after onset \n", sep='')) + 
              xlab(paste("CR Onset (ms) \n", 'trials ', trajRestrict[1], '-', trajRestrict[2], sep='')) + 
              theme(legend.position="none") + scale_color_manual(values = group_colors)
          }}else{slopeXonset_scatterplot = ggplot()}
      } # end if(restrict == 1) 
      
      
      
      else if(restrict == 0){
        
        group1 = slopeDF$Slope1_median[slopeDF$Group == Group[1]]
        group2 = slopeDF$Slope1_median[slopeDF$Group == Group[2]]
        if(length(Group) > 2){
          group3 = slopeDF$Slope1_median[slopeDF$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          slope_test = kruskal_test(Slope1_median ~ Group, data = slopeDF, distrubution = 'exact')
          slope_boxplot = ggplot(slopeDF, aes(x = Group, y = Slope1_median, group = Group)) + theme_classic() +
            ylab(paste('amp at ', deltaSlope_ms, ' ms\n', 'all trials')) +
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) +
            ggtitle(paste('ks-test, p = ', round(pvalue(slope_test),3)))
        }else{
          slope_test = 'n/a'
          slope_boxplot = ggplot()
        }
        
        group1 = Onset1$Onset1_median[Onset1$Group == Group[1]]
        group2 = Onset1$Onset1_median[Onset1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          op1_test = kruskal_test(Onset1_median ~ Group, data = Onset1, distrubution = 'exact')
          op1 = ggplot(Onset1, aes(x = Group, y = Onset1_median, group = Group)) + theme_classic() +
            ylab(paste("CR Onset (ms)\n", 'all trials')) + 
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(75,400) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(op1_test),3)))
        }else{
          op1_test = 'n/a'
          op1 = ggplot()
        }
        
        group1 = ampInt1$ampInt1_median[ampInt1$Group == Group[1]]
        group2 = ampInt1$ampInt1_median[ampInt1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          ampInt1_test = kruskal_test(ampInt1_median ~ Group, data = ampInt1, distrubution = 'exact')
          ampInt1_boxplot = ggplot(ampInt1, aes(x = Group, y = ampInt1_median, group = Group)) + theme_classic() +
            ylab(paste("CR Magnitude\n", 'all probe trials')) + 
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE) + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(0,100) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(ampInt1_test),3)))
        }else{
          ampInt1_test = 'n/a'
          ampInt1_boxplot = ggplot()
        }
        
        group1 = maxAmp1$ampProbe1_median[maxAmp1$Group == Group[1]]
        group2 = maxAmp1$ampProbe1_median[maxAmp1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          maxAmp1_test = kruskal_test(ampProbe1_median ~ Group, data = maxAmp1, distrubution = 'exact')
          maxAmp1_boxplot = ggplot(maxAmp1, aes(x = Group, y = ampProbe1_median, group = Group)) + theme_classic() +
            ylab(paste("CR Amplitude\n", 'all probe trials')) + 
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE) + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(0,1) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(maxAmp1_test),3)))
        }else{
          maxAmp1_test = 'n/a'
          maxAmp1_boxplot = ggplot()
        }
        
        group1 = Peak1$Peak1_median[Peak1$Group == Group[1]]
        group2 = Peak1$Peak1_median[Peak1$Group == Group[2]]
        if(length(Group) > 2){
          group3 = Onset1$Onset1_median[Onset1$Group == Group[3]]
        } else{ group3 = 1:10}
        if(length(group1) >= min_test && length(group2) >= min_test && length(group3) >= min_test){
          pk1_test = kruskal_test(Peak1_median ~ Group, data = Peak1, distrubution = 'exact')
          pk1 = ggplot(Peak1, aes(x = Group, y = Peak1_median, group = Group)) + theme_classic() + 
            ylab(paste("Peak Timing (ms)\n", 'all probe trials')) + 
            geom_boxplot(fill= group_colors, width = .75, alpha=.5, na.rm = TRUE)  + 
            geom_point(color = 'black', shape=1, size=1.5, na.rm = TRUE) + ylim(250,(500)) + 
            ggtitle(paste('ks-test, p = ', round(pvalue(pk1_test),3)))
        }else{
          pk1_test = 'n/a'
          pk1 = ggplot()
        }
        
        if(length(ampInt1[,1]>0)){
          if(length(Onset1[,1]>0)){
            ampInt1_scatterplot = ggplot(ampInt1, aes(x = Onset1_median, y= ampInt1_median, group = Group, color = Group)) + theme_classic() +
              geom_point(size = .75, alpha = .75) + xlim(75,400) + ylim(0,100) + ylab(paste("CR Magnitude \n", 'all probe trials')) + 
              xlab(paste("CR Onset (ms) \n", 'all trials')) + 
              theme(legend.position="none") + scale_color_manual(values = group_colors)
          }}else{ampInt1_scatterplot = ggplot()}
        
        if(length(maxAmp1[,1]>0)){
          if(length(Onset1[,1]>0)){
            maxAmp1_scatterplot = ggplot(maxAmp1, aes(x = Onset1_median, y= ampProbe1_median, group = Group, color = Group)) + theme_classic() +
              geom_point(size = .75, alpha = .75) + xlim(0,400) + ylim(0,1) + ylab(paste("CR Amplitude \n",'all probe trials')) + 
              xlab(paste("CR Onset (ms) \n", 'all trials')) + 
              theme(legend.position="none") + scale_color_manual(values = group_colors)
          }}else{maxAmp1_scatterplot = ggplot()}
        
        if(length(slopeDF[,1]>0)){
          if(length(Onset1[,1]>0)){
            slopeXonset_scatterplot = ggplot(Onset1, aes(x = Onset1_median, y= Slope, group = Group, color = Group)) + theme_classic() +
              geom_point(size = .75, alpha = .75) + xlim(75,400) + ylim(0,1) + ylab(paste("Amp 100ms after onset \n", sep='')) + 
              xlab(paste("CR Onset (ms) \n", 'all trials')) + 
              theme(legend.position="none") + scale_color_manual(values = group_colors)
          }}else{slopeXonset_scatterplot = ggplot()}
        
      } # end else if(restrict == 0)
      
      df$Categorical.Analysis = ifelse(df$trial25 > 0, 'learner', 'non-learner')
      subplot1 = ggarrange(t25, last2_boxplot, nrow = 2,ncol = 1)
      subplot2 = ggarrange(cr_plot, subplot1, nrow = 1, ncol = 2, widths = c(2,1))
      subplot3 = ggarrange(totalCR_boxplot, op1, pk1, ampInt1_boxplot, maxAmp1_boxplot, totalUR, slope_boxplot, ampInt1_scatterplot, maxAmp1_scatterplot,
                           nrow = 3, ncol = 3)
      
      setwd(path_plots)
      pdf((paste(cueisi,'-', cue, '-', substring(dataFolder,1,nchar(dataFolder)),'-AnalysisPlots.pdf',sep = '')), height=10.5, width=8)
      capture.output(print(ggarrange(subplot2,subplot3, nrow = 2, ncol = 1)))
      dev.off()
      
      if(length(which(df$Categorical.Analysis == 'learner')) > 0 && length(which(df$Categorical.Analysis == 'non-learner')) > 0){
        print(sjt.xtab(df$Group, df$Categorical.Analysis, title = paste(dataFolder, cueisi), file = paste(cueisi, '-', cue, '-Contingency.doc', sep='')))
      }
      setwd(path_summary_csv)
    } #  Ends if(cueisi != 'Star')
else { # This else statement goes with the if csv file is != 'Star'
      for(i in 1:length(df[,1])){
        df$totalStartle[i] = round((df$countStartle[i])/(df$numTrials[i])*100,1)
      }
      Group = unique(df$Group)

      for(group in Group){
        temp1 = dfm[dfm$Group == group,]
        for(session in unique(dfm$TrainingSession)){
          temp2 = temp1[temp1$TrainingSession == session,]
          medianCR = median(temp2$StartlePercent)
          groupMedian = rbind(groupMedian, data.frame('Group' = group, 'TrainingSession' = session, 'StartlePercent' = medianCR))
        }
      }

      group1 = df$totalStartle[df$Group == Group[1]]
      group2 = df$totalStartle[df$Group == Group[2]]
      if(length(group1) >= min_test && length(group2) >= min_test){
        totalStartle_wilcox = wilcox.exact(group1,group2)
        totalStartle = ggplot(df, aes(x = Group, y = totalStartle, group = Group)) + theme_classic() +
          geom_boxplot(fill= group_colors, width = .75, alpha=.5)  + geom_point(color = 'black', shape=1, size=2) +
          ggtitle(paste('ks-test, p =', round(as.numeric(totalStartle_wilcox[3]),4))) + theme(plot.title = element_text(hjust = 0))
      }else{
        totalStartle_wilcox = 0
        totalStartle = ggplot()
      }

      startle_plot = ggplot(dfm, aes(x = TrainingSession, y= StartlePercent, group = Animal, color = Group))  + theme_classic() +
        geom_line(size=1, alpha= .5) + theme(legend.position="none") + scale_color_manual(values = group_colors) +
        geom_line(data = groupMedian, aes(x=TrainingSession, y=StartlePercent, group = Group, color= Group), size = 3)

      setwd(path_plots)
      pdf((paste(cueisi,'-StartlePlots.pdf',sep = '')), height=3, width=7)
      capture.output(print(ggarrange(startle_plot,totalStartle, nrow = 1, ncol = 2, widths = c(2,1))))
      dev.off()
      setwd(path_summary_csv)
    }
   } # End loop that considers cueisi summary csv
 } # End loop that considers for(cue in Cues)

# PLOT MEDIAN TRAJECTORIES ACROSS ANIMALS --------------------------------

if(doTrajectories == 1){
  
  # Specify fps, just in case the waterfall analysis isn't performed
  setwd(path_waterfall_raw)
  fileList <- list.files()
  wf = readWaterFall(fileList[1])
  wf = data.frame(do.call(cbind, wf))
  fps = length(wf[,1])/trialDur*1000
  
  setwd(path_trajectories)
  fileList_all <- list.files()
  
  for(cueisi in unique(substring(fileList_all,10,13))){
    fileList = fileList_all[substring(fileList_all,10,13) == cueisi]
    animalList = unique(substring(fileList,1,4))
    nAnimal = 0
    responseMedian_Paired = data.frame(Animal = 'animal', Geno = 'geno', t = 0, responseAmplitude = 0, Paired = 1)
    responseMedian_Probe = data.frame(Animal = 'animal', Geno = 'geno', t = 0, responseAmplitude = 0, Paired = 0)
    
    for(animal in animalList){
      dfAnimal = data.frame()
      for (k in 1:length(fileList)){
        if(substring(fileList[k],1,4) == animal){
          temp = read.csv(fileList[k])
          dfAnimal = rbind(dfAnimal,temp)
        }
      }
      dfAnimal = dfAnimal[dfAnimal$t != 'NA',]
      dfAnimal = dfAnimal[dfAnimal$CR > 0,]
      
      # exclude trials from dfAnimal based on how the analysis is restricted
      if(restrict == 1){
        if (length(dfAnimal$value) >= fps*trialDur*(trajRestrict[2]+1)/1000){
          dfAnimal = dfAnimal[(fps*(trialDur/1000)*trajRestrict[1]+1):(fps*(trialDur/1000)*trajRestrict[2]),]
        } else {dfAnimal = 0}
      } else if (restrict == 0){
        if(length(dfAnimal$value) >= fps*(trialDur/1000)*trajMin){
          dfAnimal = dfAnimal
        } else {dfAnimal = 0}
      } else {dfAnimal = 0}
      
      if(length(dfAnimal) > 1){ # If there are any trials in dfAnimal
        dfAnimal_Paired = dfAnimal[dfAnimal$Paired == 1,]
        dfAnimal_Probe = dfAnimal[dfAnimal$Paired == 0,]
        tempPaired = data.frame()
        tempProbe = data.frame()
        for(t in unique(dfAnimal$t)){
          tpaired = data.frame(Animal = animal, Geno = dfAnimal_Paired$Geno[dfAnimal_Paired$Animal == animal][1], 't' = t, 
                               responseAmplitude = round(median(dfAnimal_Paired$value[dfAnimal_Paired$t == t]), 3), Paired = 1)
          tprobe = data.frame(Animal = animal, Geno = dfAnimal_Probe$Geno[dfAnimal_Probe$Animal == animal][1], 't' = t, 
                              responseAmplitude = round(median(dfAnimal_Probe$value[dfAnimal_Probe$t == t]), 3), Paired = 0)
          tempPaired = rbind(tempPaired, tpaired)
          tempProbe = rbind(tempProbe, tprobe)
        }
        responseMedian_Paired = rbind(responseMedian_Paired, tempPaired)
        responseMedian_Probe = rbind(responseMedian_Probe, tempProbe)
        
      if(doIndividual_trajPlots == 1){
      # Make individual animal median + confidence interval plots 
        Paired = data.frame(t = unique(dfAnimal_Paired$t), median = 0, upper = 0, lower = 0)
        for(t in unique(dfAnimal_Paired$t)){
          timepoint = dfAnimal_Paired$value[dfAnimal_Paired$t == t]
          bootCoeff <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
            bootCoeff[n] <- median(bootDataset)
          }
          Paired$t[Paired$t == t] = t
          Paired$median[Paired$t == t] = round(as.numeric(quantile(bootCoeff, .500)),3) # median
          Paired$upper[Paired$t == t] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
          Paired$lower[Paired$t == t] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
        }
        Paired = Paired[2:length(Paired[,1]),]
        
        Probe = data.frame(t = unique(dfAnimal_Probe$t), median = 0, upper = 0, lower = 0)
        for(t in unique(dfAnimal_Probe$t)){
          timepoint = dfAnimal_Probe$value[dfAnimal_Probe$t == t]
          bootCoeff <- numeric(numberOfIterations)
          for (n in 1:numberOfIterations){
            bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
            bootCoeff[n] <- median(bootDataset)
          }
          Probe$t[Probe$t == t] = t
          Probe$median[Probe$t == t] = round(as.numeric(quantile(bootCoeff, .500)),3) # median
          Probe$upper[Probe$t == t] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
          Probe$lower[Probe$t == t] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
        }
        Probe = Probe[2:length(Probe[,1]),]
        
        indPaired_plot = ggplot(Paired, aes(x=t, y=median)) + geom_line(size = .5, alpha = 1) +
          geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25, size = .1) + ggtitle(paste(length(unique(dfAnimal_Paired$fullTrajectory)), 'trials')) +
          scale_x_continuous(breaks = c(200,400,600,800,1000,1200),labels=c('0', '200', '400', '600', '800', '1000')) +
          theme_classic() + ylim(-0.05,1) + scale_color_manual(values = 'black') + xlab(paste('time (ms)', length(unique(dfAnimal_Paired$fullTrajectory)), ' trials'))
        
        indProbe_plot = ggplot(Probe, aes(x=t, y=median)) + geom_line(size = .5, alpha = 1) + 
          geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25, size = .1) + ggtitle(paste(length(unique(dfAnimal_Probe$fullTrajectory)), 'trials')) +
          scale_x_continuous(breaks = c(200,400,600,800,1000,1200),labels=c('0', '200', '400', '600', '800', '1000')) +
          theme_classic() + ylim(-0.05,0.75) + scale_color_manual(values = 'black')  + xlab(paste('time (ms)', length(unique(dfAnimal_Probe$fullTrajectory)), ' trials'))
        
        trajPlot = ggarrange(indPaired_plot, indProbe_plot, nrow = 1,ncol = 2)
        
        setwd(path_traj_plots)
        pdf(paste(animal, '-', cueisi,'-TrajPlots.pdf', sep =''), height=4, width=6)
        capture.output(print(trajPlot))
        dev.off()
        setwd(path_trajectories)
      } # end doIndividual_trajPlots
      }
      nAnimal = nAnimal+1
      print(paste(nAnimal, ' out of ', length(unique(animalList)), ' animal trajectory medians complete'))
    } # end for(i in 1:length(animalList)) loop that reads in each animal
    
    responseMedian_Paired = responseMedian_Paired[2:length(responseMedian_Paired$t),]
    responseMedian_Probe = responseMedian_Probe[2:length(responseMedian_Probe$t),]
    
    # Loop to add an average-of-average trajectory and upper and lower confidence intervals 
    summaryPaired = data.frame(t = 0, group = 'group', mean = 0, upper = 0, lower = 0)
    
    for(group in unique(responseMedian_Paired$Geno)){
      groupDF = responseMedian_Paired[responseMedian_Paired$Geno == group,]
      Paired = data.frame(t = unique(responseMedian_Paired$t), group = group, mean = 0, upper = 0, lower = 0)
      for(t in unique(responseMedian_Paired$t)){
        timepoint = groupDF$responseAmplitude[groupDF$t == t]
        bootCoeff <- numeric(numberOfIterations)
        for (n in 1:numberOfIterations){
          bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
          bootCoeff[n] <- mean(bootDataset)
        }
        Paired$t[Paired$t == t] = t
        Paired$mean[Paired$t == t] = round(as.numeric(quantile(bootCoeff, .500)),3) # median
        Paired$upper[Paired$t == t] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
        Paired$lower[Paired$t == t] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
      }
      Paired = Paired[2:length(Paired[,1]),]
      summaryPaired = rbind(summaryPaired, Paired)
    }
    summaryPaired = summaryPaired[2:length(summaryPaired[,1]),]
    
    
    # Now do Probe
    summaryProbe = data.frame(t = 0, group = 'group', mean = 0, upper = 0, lower = 0)
    
    for(group in unique(responseMedian_Probe$Geno)){
      groupDF = responseMedian_Probe[responseMedian_Probe$Geno == group,]
      Probe = data.frame(t = unique(responseMedian_Probe$t), group = group, mean = 0, upper = 0, lower = 0)
      for(t in unique(responseMedian_Probe$t)){
        timepoint = groupDF$responseAmplitude[groupDF$t == t]
        bootCoeff <- numeric(numberOfIterations)
        for (n in 1:numberOfIterations){
          bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
          bootCoeff[n] <- mean(bootDataset)
        }
        
        Probe$t[Probe$t == t] = t
        Probe$mean[Probe$t == t] = round(as.numeric(quantile(bootCoeff, .500)),3) # median
        Probe$upper[Probe$t == t] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
        Probe$lower[Probe$t == t] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
      }
      Probe = Probe[2:length(Probe[,1]),]
      summaryProbe = rbind(summaryProbe, Probe)
    }
    summaryProbe = summaryProbe[2:length(summaryProbe[,1]),]
    
    Group = unique(summaryProbe$group)
    
    if(length(Group) == 2){
      group_colors = c(group1_color, group2_color)
    } else if(length(Group) == 3){
      group_colors = c(group1_color, group2_color, group3_color)
    }
    
    avgPaired_plot = ggplot(summaryPaired, aes(x=t, y=mean, group = group, color = group)) + geom_line(size = .5, alpha = 1) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = .25, size = .1) + 
      scale_x_continuous(breaks = c(200,400,600,800,1000,1200),labels=c('0', '200', '400', '600', '800', '1000')) +
      theme_classic() + ylim(-0.05,1) + scale_color_manual(values = group_colors) 
    
    avgProbe_plot = ggplot(summaryProbe, aes(x=t, y=mean, group = group, color = group)) + geom_line(size = .5, alpha = 1) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = .25, size = .1) + 
      scale_x_continuous(breaks = c(200,400,600,800,1000,1200),labels=c('0', '200', '400', '600', '800', '1000')) +
      theme_classic() + ylim(-0.05,0.75) + scale_color_manual(values = group_colors) 
    
    Paired_plot = ggplot(responseMedian_Paired, aes(x=t, y=responseAmplitude, group = Animal, color = Geno)) + geom_line(alpha = .75) + 
      scale_x_continuous(breaks = c(200,400,600,800,1000,1200),labels=c('0', '200', '400', '600', '800', '1000')) +
      theme_classic() + ylim(-0.05,1) + scale_color_manual(values=group_colors) + theme(legend.position="none")
    
    Probe_plot = ggplot(responseMedian_Probe, aes(x=t, y=responseAmplitude, group = Animal, color = Geno)) + geom_line(alpha = .75) + 
      scale_x_continuous(breaks = c(200,400,600,800,1000,1200), labels=c('0', '200', '400', '600', '800', '1000')) + theme_classic() + 
      ylim(-0.05,1) + scale_color_manual(values=group_colors)
    
    trajPlot = ggarrange(Paired_plot,Probe_plot, avgPaired_plot, avgProbe_plot, nrow = 2,ncol = 2)
    
    setwd(path_plots)
    pdf(paste(dataFolder, '-', cueisi,'-TrajPlots.pdf', sep =''), height=8, width=10)
    capture.output(print(trajPlot))
    dev.off()
    
    # setwd(path_summary_csv)
    # write.csv(responseAverage_Paired, paste(dataFolder,cueisi, 'avgTraj_Paired.csv', sep = ''))
    # write.csv(responseAverage_Probe, paste(dataFolder,cueisi, 'avgTraj_Probe.csv', sep = ''))
    
    setwd(path_trajectories)
  } # end for each cueisi in experiment statement
} # end if(doTrajectories == 1) 


# PLOT MEDIAN TRAJECTORIES ACROSS ANIMALS ALIGNED TO CR ONSET ---------------------------------

if(doTrajectories_aligned){
  setwd(path_trajectories_aligned)
  fileList_all <- list.files()
  
  for(cueisi in unique(substring(fileList_all,10,13))){
    
    fileList = fileList_all[substring(fileList_all,10,13) == cueisi]
    animalList = unique(substring(fileList,1,4))
    nAnimal = 0
    responseMedian_Paired = data.frame(Animal = 'animal', Geno = 'geno', t_truc = 0, responseAmplitude = 0, Paired = 1)
    responseMedian_Probe = data.frame(Animal = 'animal', Geno = 'geno', t_truc = 0, responseAmplitude = 0, Paired = 0)
    
    for(animal in animalList){
      dfAnimal = data.frame()
      for (k in 1:length(fileList)){
        if(substring(fileList[k],1,4) == animal){
          temp = read.csv(fileList[k])
          dfAnimal = rbind(dfAnimal,temp)
        }
      }
      dfAnimal = dfAnimal[dfAnimal$t_truc != 'NA',]
      dfAnimal = dfAnimal[dfAnimal$CR > 0,]
      
      # exclude trials from dfAnimal based on how the analysis is restricted
      if(restrict == 1){
        if (length(dfAnimal$value) >= fps*(traj_truc/1000)*(trajRestrict[2])){ 
          dfAnimal = dfAnimal[(fps*(traj_truc/1000)*trajRestrict[1]+1):(fps*(traj_truc/1000)*trajRestrict[2]),]
        } else {dfAnimal = 0}
      } else if (restrict == 0){
        if(length(dfAnimal$value) >= fps*(traj_truc/1000)*trajMin){
          dfAnimal = dfAnimal
        } else {dfAnimal = 0}
      } else {dfAnimal = 0}
      
      if(length(dfAnimal) > 1){
        dfAnimal_Paired = dfAnimal[dfAnimal$Paired == 1,]
        dfAnimal_Probe = dfAnimal[dfAnimal$Paired == 0,]
        tempPaired = data.frame()
        tempProbe = data.frame()
        for(t_truc in unique(dfAnimal$t_truc)){
          tpaired = data.frame(Animal = animal, Geno = dfAnimal_Paired$Geno[dfAnimal_Paired$Animal == animal][1], 't_truc' = t_truc, 
                               responseAmplitude = round(median(dfAnimal_Paired$value[dfAnimal_Paired$t_truc == t_truc]), 3), Paired = 1)
          tprobe = data.frame(Animal = animal, Geno = dfAnimal_Probe$Geno[dfAnimal_Probe$Animal == animal][1], 't_truc' = t_truc, 
                              responseAmplitude = round(median(dfAnimal_Probe$value[dfAnimal_Probe$t_truc == t_truc]),3), Paired = 0)
          tempPaired = rbind(tempPaired, tpaired)
          tempProbe = rbind(tempProbe, tprobe)
        }
        responseMedian_Paired = rbind(responseMedian_Paired, tempPaired)
        responseMedian_Probe = rbind(responseMedian_Probe, tempProbe)
        
        if(doIndividual_trajPlots_aligned == 1){
          # Make individual animal median + confidence interval plots 
          Paired = data.frame(t = unique(dfAnimal_Paired$t_truc), median = 0, upper = 0, lower = 0)
          for(t in unique(dfAnimal_Paired$t_truc)){
            timepoint = dfAnimal_Paired$value[dfAnimal_Paired$t_truc == t]
            bootCoeff <- numeric(numberOfIterations)
            for (n in 1:numberOfIterations){
              bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
              bootCoeff[n] <- median(bootDataset)
            }
            Paired$t[Paired$t == t] = t
            Paired$median[Paired$t == t] = round(as.numeric(quantile(bootCoeff, .500)),3) # median
            Paired$upper[Paired$t == t] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
            Paired$lower[Paired$t == t] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
          }
          Paired = Paired[2:length(Paired[,1]),]
          
          Probe = data.frame(t = unique(dfAnimal_Probe$t_truc), median = 0, upper = 0, lower = 0)
          for(t in unique(dfAnimal_Probe$t_truc)){
            timepoint = dfAnimal_Probe$value[dfAnimal_Probe$t_truc == t]
            bootCoeff <- numeric(numberOfIterations)
            for (n in 1:numberOfIterations){
              bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
              bootCoeff[n] <- median(bootDataset)
            }
            
            Probe$t[Probe$t == t] = t
            Probe$median[Probe$t == t] = round(as.numeric(quantile(bootCoeff, .500)),3) # median
            Probe$upper[Probe$t == t] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
            Probe$lower[Probe$t == t] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
          }
          Probe = Probe[2:length(Probe[,1]),]
          
          indPaired_plot = ggplot(Paired, aes(x=t, y=median)) + geom_line(size = .5, alpha = 1) +
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25, size = .1) + ggtitle(paste(length(unique(dfAnimal_Paired$fullTrajectory)), 'trials')) +
            scale_x_continuous(breaks= seq(from = 0, to = traj_truc, by = 100), labels = as.character(seq(from = -traj_offset, to = (traj_truc-traj_offset), by = 100))) +
            theme_classic() + ylim(-0.05,1) + scale_color_manual(values = 'black') + xlab('time (ms)') 
          
          indProbe_plot = ggplot(Probe, aes(x=t, y=median)) + geom_line(size = .5, alpha = 1) + 
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25, size = .1) + ggtitle(paste(length(unique(dfAnimal_Probe$fullTrajectory)), 'trials')) +
            scale_x_continuous(breaks= seq(from = 0, to = traj_truc, by = 100), labels = as.character(seq(from = -traj_offset, to = (traj_truc-traj_offset), by = 100))) +
            theme_classic() + ylim(-0.05,0.75) + scale_color_manual(values = 'black') + xlab('time (ms)')
          
          trajPlot = ggarrange(indPaired_plot, indProbe_plot, nrow = 1,ncol = 2)
          
          setwd(path_traj_plots)
          pdf(paste(animal, '-', cueisi,'-TrajPlots_aligned.pdf', sep =''), height=4, width=6)
          capture.output(print(trajPlot))
          dev.off()
          setwd(path_trajectories_aligned)
        } # end doIndividual_trajPlots_aligned
      } # end if animaldf > 1
      nAnimal = nAnimal+1
      print(paste(nAnimal, ' out of ', length(unique(animalList)), ' CR-Aligned animal trajectory medians complete'))
    } # end loop read in animals
    
    responseMedian_Paired = responseMedian_Paired[2:length(responseMedian_Paired$t_truc),]
    responseMedian_Probe = responseMedian_Probe[2:length(responseMedian_Probe$t_truc),]
    
    # Loop to add an average-of-average trajectory and upper and lower confidence intervals 
    # Do Paired First 
    
    summaryPaired = data.frame(t_truc = 0, group = 'group', mean = 0, upper = 0, lower = 0)
    
    for(group in unique(responseMedian_Paired$Geno)){
      groupDF = responseMedian_Paired[responseMedian_Paired$Geno == group,]
      Paired = data.frame(t_truc = unique(responseMedian_Paired$t_truc), group = group, mean = 0, upper = 0, lower = 0)
      for(t_truc in unique(responseMedian_Paired$t_truc)){
        timepoint = groupDF$responseAmplitude[groupDF$t_truc == t_truc]
        bootCoeff <- numeric(numberOfIterations)
        for (n in 1:numberOfIterations){
          bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
          bootCoeff[n] <- mean(bootDataset)
        }
        Paired$t_truc[Paired$t_truc == t_truc] = t_truc
        Paired$mean[Paired$t_truc == t_truc] = round(as.numeric(quantile(bootCoeff, .500)),3) # mean
        Paired$upper[Paired$t_truc == t_truc] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
        Paired$lower[Paired$t_truc == t_truc] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
      }
      Paired = Paired[2:length(Paired[,1]),]
      summaryPaired = rbind(summaryPaired, Paired)
    }
    
    summaryPaired = summaryPaired[2:length(summaryPaired[,1]),]
    
    
    # Now do Probe
    summaryProbe = data.frame(t_truc = 0, group = 'group', mean = 0, upper = 0, lower = 0)
    
    for(group in unique(responseMedian_Probe$Geno)){
      groupDF = responseMedian_Probe[responseMedian_Probe$Geno == group,]
      Probe = data.frame(t_truc = unique(responseMedian_Probe$t_truc), group = group, mean = 0, upper = 0, lower = 0)
      for(t_truc in unique(responseMedian_Probe$t_truc)){
        timepoint = groupDF$responseAmplitude[groupDF$t_truc == t_truc]
        bootCoeff <- numeric(numberOfIterations)
        for (n in 1:numberOfIterations){
          bootDataset <- sample(timepoint, length(timepoint), replace = TRUE, prob=NULL)
          bootCoeff[n] <- mean(bootDataset)
        }
        Probe$t_truc[Probe$t_truc == t_truc] = t_truc
        Probe$mean[Probe$t_truc == t_truc] = round(as.numeric(quantile(bootCoeff, .500)),3) # mean
        Probe$upper[Probe$t_truc == t_truc] = round(as.numeric(quantile(bootCoeff, upperCI)),3) # upper
        Probe$lower[Probe$t_truc == t_truc] = round(as.numeric(quantile(bootCoeff, lowerCI)),3) # lower
      }
      Probe = Probe[2:length(Probe[,1]),]
      summaryProbe = rbind(summaryProbe, Probe)
    }
    
    summaryProbe = summaryProbe[2:length(summaryProbe[,1]),]
    
    avgPaired_plot = ggplot(summaryPaired, aes(x=t_truc, y=mean, group = group, color = group)) + geom_line(size = .5, alpha = 1) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = .25, size = .1) + 
      scale_x_continuous(breaks= seq(from = 0, to = traj_truc, by = 100), labels = as.character(seq(from = -traj_offset, to = (traj_truc-traj_offset), by = 100))) +
      theme_classic() + ylim(-0.05,1) + scale_color_manual(values = group_colors) 
    
    avgProbe_plot = ggplot(summaryProbe, aes(x=t_truc, y=mean, group = group, color = group)) + geom_line(size = .5, alpha = 1) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = .25, size = .1) + 
      scale_x_continuous(breaks= seq(from = 0, to = traj_truc, by = 100), labels = as.character(seq(from = -traj_offset, to = (traj_truc-traj_offset), by = 100))) +
      theme_classic() + ylim(-0.05,0.75) + scale_color_manual(values = group_colors) 
    
    Paired_plot = ggplot(responseMedian_Paired, aes(x=t_truc, y=responseAmplitude, group = Animal, color = Geno)) + geom_line(alpha = .75) + 
      scale_x_continuous(breaks= seq(from = 0, to = traj_truc, by = 100), labels = as.character(seq(from = -traj_offset, to = (traj_truc-traj_offset), by = 100))) +
      theme_classic() + ylim(-0.05,1) + scale_color_manual(values=group_colors) + theme(legend.position="none")
    
    Probe_plot = ggplot(responseMedian_Probe, aes(x=t_truc, y=responseAmplitude, group = Animal, color = Geno)) + geom_line(alpha = .75) + 
      scale_x_continuous(breaks= seq(from = 0, to = traj_truc, by = 100), labels = as.character(seq(from = -traj_offset, to = (traj_truc-traj_offset), by = 100))) +
      theme_classic() + ylim(-0.05,1) + scale_color_manual(values=group_colors) + theme(legend.position="none")
    
    trajPlot = ggarrange(Paired_plot,Probe_plot, avgPaired_plot, avgProbe_plot, nrow = 2,ncol = 2)
    
    setwd(path_plots)
    pdf(paste(dataFolder, '-', cueisi,'-TrajPlots-Aligned.pdf', sep =''), height=8, width=10)
    capture.output(print(trajPlot))
    dev.off()
    
    setwd(path_trajectories_aligned)
  } # end for each cueisi in experiment statement
} # End doTrajectories_aligned statement