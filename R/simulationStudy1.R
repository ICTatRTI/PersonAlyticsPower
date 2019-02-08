# inputs for 'simulation study 1' rare disease ICT R21

library(tidyr)
library(ggplot2)

# settings ####

# sample sizes
sampleSizes <- c(
  1,                    # not originally proposed, added to ensure n=1 works from the begining
  seq(10, 40, by = 5),  # as originally prposed
  seq(50, 100, by = 10) # as originally prposed
)

# observations per participant
nObservations <- c(
  seq(20, 50, by = 5),
  seq(60, 100, by = 10)
)

# proportion of observations that comprise the baseline phase
propBaseline <- c(1:5, 10, 15, 20)/100

#
propNobs <- expand.grid(nObservations, propBaseline)
names(propNobs) <- c('nObservations', 'propBaseline')
propNobs$nBaseline <- propNobs[,1]*propNobs[,2]
write.csv(propNobs, file='Sim1_nBaseline.csv')

# degree of variability with and between participants
# not sure what Ty means by this, so let's decompose the total variance
# and see what our options are

# variance of time, assuming an integer sequence starting at 0
timeVariances <- list()
for(i in seq_along(nObservations)) timeVariances[[i]] <- var(0:(nObservations[i]-1))
timeVariances <- data.frame(nObservations=nObservations, timeVariances=unlist(timeVariances))
ggplot(timeVariances, aes(x=nObservations, y=timeVariances)) +
     geom_line() +
     xlab('Number of Observations per Person') +
     ylab('Variance of the Time Variable')

# variance of phase
phaseVariances <- list()
phaseTimeVariances <- list()
for(i in seq_along(nObservations))
{
  times <- 0:(nObservations[i]-1)
  for(j in seq_along(propBaseline))
  {
    nTimes    <- length(times)
    nBaseline <- ceiling( propBaseline[j]*nTimes )
    phase     <- rep(1, nTimes)
    phase[1:nBaseline] <- 0
    phaseVarName <- paste('O', nTimes, 'BL', propBaseline[j], sep='')
    phaseVariances[[phaseVarName]] <- var(phase)
    phaseTimeVariances[[phaseVarName]] <- var(times*phase)
  }
}
phaseVariances <- matrix(unlist(phaseVariances), length(nObservations),
                         length(propBaseline), byrow=TRUE)
phaseVariances <- data.frame(phaseVariances)
names(phaseVariances) <- paste('BL', propBaseline, sep='')
phaseVariances <- data.frame(nObservations, phaseVariances)
phaseVariances <- gather(phaseVariances, 'propBL', 'Variance', BL0.01:BL0.2)
phaseVariances$propBL <- factor(as.numeric(gsub('BL', '', phaseVariances$propBL)))
ggplot(phaseVariances, aes(x=nObservations, y=Variance, group=propBL, col=propBL)) +
  geom_line(position=position_jitter(w=0.02, h=0.01), size = 2) +
  xlab('Number of Observations per Person') +
  ylab('Variance of the Phase Variable')

# variance of the phase X time interaction
phaseTimeVariances <- matrix(unlist(phaseTimeVariances), length(nObservations),
                         length(propBaseline), byrow=TRUE)
phaseTimeVariances <- data.frame(phaseTimeVariances)
names(phaseTimeVariances) <- paste('BL', propBaseline, sep='')
phaseTimeVariances <- data.frame(nObservations, phaseTimeVariances)
phaseTimeVariances <- gather(phaseTimeVariances, 'propBL', 'Variance', BL0.01:BL0.2)
phaseTimeVariances$propBL <- factor(as.numeric(gsub('BL', '', phaseTimeVariances$propBL)))
ggplot(phaseTimeVariances, aes(x=nObservations, y=Variance, group=propBL, col=propBL)) +
  geom_line(position=position_jitter(w=0.02, h=10), size = 2) +
  xlab('Number of Observations per Person') +
  ylab('Variance of Phase X Time')



