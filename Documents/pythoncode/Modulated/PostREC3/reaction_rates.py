##Reaction rate plotter over a ange of different temperatures, this will help us determine the timescale of some of the rates and help determine whether any of the reactants are in equilibrium



##Idea is to see which rates dominate over a certain temperature range

##Rates have to be fast over the whole temperature range in order for them to dominate


from matplotlib import *
from pylab import *
from math import exp, log
import numpy as np


matter_temperature = []
radiation_temperature = []

for i in arange (1,8000):
	matter_temperature.append(i)

for i in arange (2,3000):
	radiation_temperature.append(i)
	
Z_initial = 1

k0 = []
k1 = []
k2 = []
k3 = []
k4 = []
k5 = []
k6 = []
k7 = []
k8 = []
k9 = []
k10 = []
k11 = []
k12 = []
k13 = []
k14 = []
k15 = []
k16 = []
k17 = []
k18 = []
k19 = []
k20 = []
k51 = []
k52 = []




for i in arange (0,(len(matter_temperature))):
 
	k0.append(4.91*(10**-22)*(matter_temperature[i]**4))								##Abell et al 1997									cm3s-1
	k1.append(1.4*(10**-18)*(matter_temperature[i]**0.928)*exp(-matter_temperature[i]/16200))   		##Galli D & Palla 1998								cm3s-1
		
	if (matter_temperature[i] < 300):
		k2.append(1.5*(10**-9))			 							##Galli D & Palla 1998 								cm3s-1
	else:
		k2.append(4.0*(10**-9)*(matter_temperature[i]**-0.17))
		
	k3.append(1.85*(10**-23)*(matter_temperature[i]**1.8))   						##Shapiro + Kong 1987								cm3s-1
		
	k4.append(6.4*(10**-10))			 								##Karpos et al 1979, Galli D & Palla 1998			cm3s-1
	k5.append(5.5*(10**-29)*(1/matter_temperature[i]))								##PSS83
	k6.append(3.0*(10**-10)*exp(-21050/matter_temperature[i])	)  				    ##Galli D & Palla 1998								cm3s-1
	k7.append(1.91*(10**-9)*(matter_temperature[i]**0.136)*exp(-53407.1/matter_temperature[i])) 		##Trevian, C.S & Tennyson J.2002					cm3s-1
	k8.append(1.067*(10**-10)*(matter_temperature[i]**2.012)*exp(-(4.463/matter_temperature[i])*((1+0.2472*matter_temperature[i])**3.512))) ##Dove and Mandy 1986												##Fit to data of MSM96
	k9.append(0.0)                                                   #Janev et al 1987 
	k10.append(2.5634*(10**-9)*(matter_temperature[i]**1.78186))         ##Janev et al 1987
	k11.append(1.40*(10**-7)*((matter_temperature[i]/300)**0.487)*exp(-matter_temperature[i]/29300)) 	##Lepp,S. Stancil P.C., & Dalgarno 2002				cm3s-1
	k12.append(6.9*(10**-9)*(matter_temperature[i]**-0.35))							##Galli D & Palla 1998								cm3s-1
	k13.append(2.0*(10**-7)*(matter_temperature[i]**-0.5))							    ##Galli D & Palla 1998 								cm3s-1
	k14.append(5*(10**-6)*(matter_temperature[i]**-0.5))					##Abell et al 1996 (modelling primordial gas in numerical cosmology)
	k15.append(36.7*(matter_temperature[i]**-2.28)*exp(-47172/matter_temperature[i]))
	k51.append(1.88*(10**-10)*(matter_temperature[i]**-0.64)) 
	
	
for i in arange(0,(len(radiation_temperature))):	
	## Capitelli, M., Coppola, C. M., Diomede, P., & Longo, S. 2007, A&A, 470, 811	
	k16.append(0.144*(radiation_temperature[i]**2.13)*exp(-8650/radiation_temperature[i]))     				##Abell et al 1997									cm3s-1
	k17.append(1.63*(10**7)*exp(-32400/radiation_temperature[i]))  						##Galli D & Palla 1998 								cm3s-1
	k18.append(2.9*(10**2)*(radiation_temperature[i]**1.56)*exp(-178500/radiation_temperature[i]))	 		##Galli D & Palla 1998								cm3s-1
	k19.append(90*(radiation_temperature[i]**1.48)*(exp(-335000/radiation_temperature[i]))) ##Galli D & Palla 1998
	k20.append(1.13*(10**6)*(radiation_temperature[i]**0.369)*exp(-140000/radiation_temperature[i]))			  ##Glover, S.C.0 & Jappsen 2007						cm3s-1
	k52.append((2.41*(10**15)*(radiation_temperature[i]**1.5)*exp(-39472/radiation_temperature[i]))*(8.76*10**-11)*(1+(Z_initial**-0.58)))


from matplotlib import pyplot


##some radiation based process are currently excluded


##Plotting H- destruction mechanisms
figure(1)
plot(matter_temperature,k2,'b')
plot(matter_temperature,k9,'g')
plot(matter_temperature,k10,'m')
plot(matter_temperature,absolute(k11),'c')
plot(matter_temperature,k12,'k')
plot(matter_temperature,k14,'r')
pyplot.yscale('log')
pyplot.xscale('log')

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k2$', r'$k9$', r'$k10$', r'$k11$', r'$k11$', r'$k14$']	
legend_location = 'upper left'
colors = ['b','g','m','c','k','r']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the destruction of $H^-$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')	
show()



##Plotting H- creation mechanisms
figure(2)
plot(matter_temperature,k1,'b')
plot(matter_temperature,k15,'r')
pyplot.yscale('log')
pyplot.xscale('log')

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k1$', r'$k15$']	
legend_location = 'lower left'
colors = ['b','r']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the creation of $H^-$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')	
show()

##Plotting H+ creation mechanisms

figure(3)
plot(matter_temperature,k4,'b')
#plot(matter_temperature,k17,'r')
pyplot.yscale('log')
pyplot.xscale('log')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k4$']	
legend_location = 'upper left'
colors = ['b']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the creation of $H^+$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')	
show()


##Plotting H+ destruction mechanisms

figure(4)
plot(matter_temperature,k3,'b')
plot(matter_temperature,k6,'r')
plot(matter_temperature,k11,'k')
plot(matter_temperature,k12,'c')
plot(matter_temperature,k51,'m')
pyplot.yscale('log')
pyplot.xscale('log')

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k3$', r'$k6$', r'$k11$', r'$k12$', r'$k51$']	
legend_location = 'lower left'
colors = ['b','r','k','c','m']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the destruction of $H^+$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')	
show()

##Plotting H2+ formation mechamisms

figure(5)
plot(matter_temperature,k3,'b')
plot(matter_temperature,k6,'r')
plot(matter_temperature,k12,'k')
#plot(matter_temperature,k18,'c')
pyplot.yscale('log')
pyplot.xscale('log')

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k3$', r'$k6$', r'$k12$']	
legend_location = 'lower left'
colors = ['b','r','k']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the formation of $H_{2}^{+}$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')
	
show()

##Plotting H2+ destruction mechanisms

figure(6)
plot(matter_temperature,k4,'b')
plot(matter_temperature,k13,'r')
plot(matter_temperature,k14,'k')
#plot(matter_temperature,k17,'c')
#plot(matter_temperature,k19,'m')
pyplot.yscale('log')
pyplot.xscale('log')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k4$', r'$k13$', r'$k14$']	
legend_location = 'upper right'
colors = ['b','r','k']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the destruction of $H_{2}^{+}$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')

show()

##Plotting H2 formation mechanisms

figure(7)
plot(matter_temperature,k2,'b')
plot(matter_temperature,k4,'r')
plot(matter_temperature,k5,'k')
plot(matter_temperature,k14,'c')

pyplot.yscale('log')
pyplot.xscale('log')

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k2$', r'$k4$', r'$k5$', r'$k14$']	
legend_location = 'lower left'
colors = ['b','r','k','c']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the formation of $H_2$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')

	
show()

##Plotting H2 destruction mechanisms

figure(8)
plot(matter_temperature,k6,'b')
plot(matter_temperature,k7,'r')
plot(matter_temperature,k8,'k')
plot(matter_temperature,k15,'c')
#plot(matter_temperature,k18,'m')
#plot(matter_temperature,k20,'m')
pyplot.yscale('log')
pyplot.xscale('log')

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})	
labels = [r'$k6$', r'$k7$', r'$k8$', r'$k15$']	
legend_location = 'lower left'
colors = ['b','r','k','c']	
plt.legend(labels, loc=legend_location)	
title(r'Timescales of rates relating to the destruction of $H_2$')	
xlabel('Temperature (K)')	
ylabel(r'Rate of reaction ($cm^{-3}s^{-1}$)')


show()

