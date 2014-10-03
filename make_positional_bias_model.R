## use estimated positional bias profiles for fragment selection

### RNA fragmentation model:
rnaf = read.csv('rnaf_pos_bias.csv', header=FALSE, colClasses=c('numeric', 'numeric'))
names(rnaf) = c('pospct', 'prob')
rnaf$pospct = seq(0.01, 1, by=0.01) #x-axis should be calibrated from 0-99%, but it was hard to click exactly :) 
rnaf$prob[rnaf$prob < 0] = 0 #minor adjustment to y-axis calibration
# make sure probabilities sum to 1:
rnaf$prob = rnaf$prob/sum(rnaf$prob) #it was super close anyway, but let's be careful.

save(rnaf, file='data/rnaf.rda', compress='xz')

### cDNA fragmentation model:
cdnaf = read.csv('cdna_pos_bias.csv', header=FALSE, colClasses=c('numeric', 'numeric'))
names(cdnaf) = c('pospct', 'prob')
cdnaf$pospct = seq(0.01, 1, by=0.01)
cdnaf$prob[cdnaf$prob < 0] = 0
cdnaf$prob = cdnaf$prob/sum(cdnaf$prob)

save(cdnaf, file='data/cdnaf.rda', compress='xz')


