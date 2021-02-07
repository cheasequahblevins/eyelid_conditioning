sampleZ = 0
for(i in 1:100000){
  sampleZ[i] = mean(sample(C, 3))
}

hist(sampleZ,100)
quantile(sampleZ,c(.01,.05,.95,.99))
ecdf(sampleZ)(mean(CI))