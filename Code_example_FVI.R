#### Example for calculating the FVI of C. hippurus (see table 1 for references). Marín-Enríquez et al.Quantifying fishing vulnerability from life-history traits: the Fishing Vulnerability Index ####
#Dr. Emigdio Marin-Enriquez (emarin@cibnor.mx)

#Case 1: no uncertainity estimates
K <- 0.57 #yr-1
Linf <- 122.15 #cm
L50 <- 80 #cm
Tmax <- 5 #yr
Cf <- 1

Cf*(Tmax/K)*(L50/Linf)

#Case 2: with uncertainity estimate, with no correlation between parameters (or no available correlation information)
source("FVI_MC.R") #call the function
c.hip.unc <- FVI_MC(K=K,KLow = 0.52,KUp = 0.62,TmaxUp = 5.5, TmaxLow = 4.3,L50=L50, Tmax=Tmax, L50Low = 70, L50Up = 92, Linf = Linf, LinfLow = 110, LinfUp = 135) #please note that the uncertainity values are not real, just established for demonstrational purposes.
names(c.hip.unc) #explore the results
hist(c.hip.unc$simulations, probability=T ,xlab="FVI", main="FVI of C. hippurus") #simple plot
lines(density(c.hip.unc$simulations), lwd=2)

#case 3: with uncertainity estimate, and with correlation between parameters 
c.hip.corr <- FVI_MC(K=K,KLow = 0.52,KUp = 0.62,TmaxUp = 5.5, TmaxLow = 4.3,L50=L50, Tmax=Tmax, L50Low = 70, L50Up = 92, Linf = Linf, LinfLow = 110, LinfUp = 135, cor_K_Linf = -0.3)
hist(c.hip.corr$simulations, probability=T ,xlab="FVI", main="FVI of C. hippurus with correlation") #simple plot
lines(density(c.hip.unc$simulations), lwd=2)

