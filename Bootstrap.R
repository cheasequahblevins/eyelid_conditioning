sampleZ = 0
for(i in 1:100000){
  sampleZ[i] = mean(sample(IsolatDisc_CR, 8, replace = TRUE))
}

hist(sampleZ,100)
quantile(sampleZ,c(.01,.05,.95,.99))
ecdf(sampleZ)(.1963)
ecdf(sampleZ)(.08622)

## Perform PCA per-animal  ------------
library(ggfortify)
library(pca3d)

setwd('/Users/Cheasequah/Desktop/Eyelid_Behavior/eyelid_analysis/2019-03-SEFL/trajectories/')
fileList = list.files()
animalList = unique(substring(fileList,1,4))

for(j in animalList){
  traj= data.frame()
  for(k in 1:length(fileList)){
    if(substring(fileList[k],1,4) == j){
      hold = read.csv(fileList[k])
      hold = hold[hold$Paired == 0,]
      hold = hold[hold$CR == 1,]
      hold = hold[c('value', 'fullTrajectory', 'Geno')]
      #hold = hold[2:length(hold[,1]),2:length(hold[1,])]
      traj = rbind(traj,hold)
    }}
      uniqueFT = unique(traj$fullTrajectory)
      uniqueFT = uniqueFT[!is.na(uniqueFT)]
      traj_PCA = data.frame()
      for(i in uniqueFT){
        fullTrajectory = traj$fullTrajectory[traj$fullTrajectory == i][1]
        hold1 = traj[traj$fullTrajectory == i,]
        hold2 = hold1[c('value')]
        hold2 = as.data.frame(t(hold2$value))
        hold2 = cbind(fullTrajectory,hold2)
        hold2 = hold2[,1:595]
        traj_PCA = rbind(traj_PCA,hold2)
      }
      if(length(traj_PCA[,1])>4){
        setwd('/Users/Cheasequah/Desktop/PCA/')
        pdf(paste(j,'-PCA.pdf', sep=''))
        ggo = autoplot(prcomp(traj_PCA[,2:595]), data = traj_PCA, x = 1, y =2, size = 1, alpha  = .25, colour = 'black',
                       xlim = c(-1,1), ylim = c(-1,1),loadings = TRUE, loadings.colour = 'blue') + theme_classic()
        print(ggo)
        dev.off()
        setwd('/Users/Cheasequah/Desktop/Eyelid_Behavior/eyelid_analysis/2019-03-SEFL/trajectories/')
      }}


project.b = predict(prcomp(traj_PCA[,2:595]), b)
traj.pca <- prcomp(traj_PCA[,2:595], center = TRUE,scale. = TRUE)


## Take traj and separate out steep traces ------------
library(ggplot2)
library(ggfortify)
library(pca3d)
library(rgl)

baselineDur = 200
trialDur = 1200
baseline_color = 'black'
cue_color = 'red'
airpuff_color = 'gray'
setwd('/Users/Cheasequah/Desktop/Eyelid_Behavior/eyelid_analysis/2019-03-SEFL/trajectories/')
fileList = list.files()
animalList = unique(substring(fileList,1,4))

for(j in animalList){
  traj = data.frame()
  for(k in 1:length(fileList)){
    if(substring(fileList[k],1,4) == j){
      hold = read.csv(fileList[k])
      hold = hold[hold$Paired == 0,]
      hold = hold[hold$CR == 1,]
      traj = rbind(traj, hold)
    }}
  uniqueFT = unique(traj$fullTrajectory)
  uniqueFT = uniqueFT[!is.na(uniqueFT)]
  
reject = traj[0,]
accept = traj[0,]
for(i in uniqueFT){
  sample = traj$value[traj$fullTrajectory == i]
  diff = diff(x = sample[1:584],lag = 15, differences = 1) # Leave off last 20 ms
  if(max(rev(sort(diff))[1]) >.40){
    reject = rbind(reject, traj[traj$fullTrajectory == i,])
  } else{
    accept = rbind(accept, traj[traj$fullTrajectory == i,])
  }}
  if(length(accept[,1])>0){
    accept_plot =  ggplot(accept, aes(x=t, y=value,
                                      color=region,group=fullTrajectory))
    accept_plot = accept_plot + theme_classic() + theme(legend.position = 'none')
    accept_plot = accept_plot + geom_line(size=0.4)
    accept_plot = accept_plot + scale_color_manual(values=c(
      baseline_color, cue_color, airpuff_color, airpuff_color))
    accept_plot = accept_plot + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
      theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
      scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
      scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100')) +
      ggtitle(paste(length(unique(accept$fullTrajectory)), 'accepted'))
  } else{accept_plot = ggplot()}
  if(length(reject[,1])>0){
    reject_plot =  ggplot(reject, aes(x=t, y=value,
                                      color=region,group=fullTrajectory))
    reject_plot = reject_plot + theme_classic() + theme(legend.position = 'none')
    reject_plot = reject_plot + geom_line(size=0.4)
    reject_plot = reject_plot + scale_color_manual(values=c(
      baseline_color, cue_color, airpuff_color, airpuff_color))
    reject_plot = reject_plot + xlab('\n             time (ms)') + ylab('Percent Eyelid Closure') +
      theme(text=element_text(face="bold", size=12)) + coord_cartesian(ylim=c(0,1)) +
      scale_x_continuous(breaks= seq(from = baselineDur, to = trialDur, by = 100), labels = as.character(seq(from = 0, to = (trialDur-baselineDur), by = 100))) +
      scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c('0', '25', '50', '75', '100')) + 
      ggtitle(paste(length(unique(reject$fullTrajectory)), 'rejected'))
  } else{reject_plot = ggplot()}
    
    setwd('/Users/Cheasequah/Desktop/PCA/')
    pdf(paste(j,'-trajectories.pdf', sep=''), height = 2, width = 7)
    print(ggarrange(accept_plot,reject_plot, nrow=1,ncol=2))
    dev.off()
    setwd('/Users/Cheasequah/Desktop/Eyelid_Behavior/eyelid_analysis/2019-03-SEFL/trajectories/')
    
    accept$Include = 'Accept'
    reject$Include = 'Reject'
    traj = rbind(accept,reject)
    traj = traj[c('value', 'fullTrajectory', 'Include')]
    uniqueFT = unique(traj$fullTrajectory)
    uniqueFT = uniqueFT[!is.na(uniqueFT)]
    traj_PCA = data.frame()
    for(i in uniqueFT){
      Include = traj$Include[traj$fullTrajectory == i][1]
      hold1 = traj[traj$fullTrajectory == i,]
      hold2 = hold1[c('value')]
      hold2 = as.data.frame(t(hold2$value))
      hold2 = cbind(Include,hold2)
      hold2 = hold2[,1:595]
      traj_PCA = rbind(traj_PCA,hold2)
    }
    if(length(traj_PCA[,1])>4){
      setwd('/Users/Cheasequah/Desktop/PCA/')
      pdf(paste(j,'-PCA.pdf', sep=''))
      ggo = autoplot(prcomp(traj_PCA[,150:500]), data = traj_PCA, x = 1, y =2, size = 2, alpha  = .5, colour = 'Include',
                     xlim = c(-1,1), ylim = c(-1,1)) + theme_classic()
      print(ggo)
      dev.off()
      setwd('/Users/Cheasequah/Desktop/Eyelid_Behavior/eyelid_analysis/2019-03-SEFL/trajectories/')
    }
}

pc = princomp(traj_PCA[,150:175])

plot3d(pc$scores[,1:3])