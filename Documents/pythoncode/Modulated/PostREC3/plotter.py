####Python plotting program for PostREC, takes values from integration and plots their evolution as a function of redhsift for the user. Hydrogen, Deuterium, Lithium and Helium are all plotted
##as well as the evolution of the matter temperature. Each figure is shown in turn with an option for the user to save or dismiss the file.

##Inputs:
##Z_initial : The chosen start point of the Integration (used for the title of the plot)
##Z_final : The chosen end point of the Integration (used for the title of the plot)
## y : the data values obatined from the integration, size is 25 by X, where X is (Z_initial-Z_final)/step_size

def plotting_program(Z_initial,Z_final,y):
	
	from pylab import *                                                                         #Needed for plotting commands
	from matplotlib import rc                                                                   #Needed for titles of graphs
	from numpy import absolute                                                                  #Needed to calculate absolute values
	from matplotlib import pyplot                                                               #Needed for plotting commands
	
	##Begin plotting of Hydrogen species

	figure()
	
	##Each species is plotted as a fractional abundance of the total number of nuceli. The list of chemical species and thier corresponding array columns are shown for reference purposes.
	##The fractional abudnace is plotted as a function of redshift, which is calculated by taking the CMB temperature, dividing it by today's value and subtracting 1 (Becuase T_0 = 2.725 at Z = 0).
	##The scale is logarithmic for each axis, for the best visual representation. The limits are also set with this in mind. The lines are color coded and labelled as such.
	
	#species_list = ["[H]","[H+]","[H-]","[H2]","[H2+]","[e-]","[hv]","[D]","[D-]","[D+]","[HD]","[HD+]","[He]","[He+]","[He++]","[HeH+]",[Li], [Li+], [Li-], [LiH], [LiH+], [H2D+],  [D2]]	
		          #    0     1      2      3      4       5      6      9    10     11      12      13     14     15       16      17      18     19    20      21     22  ,  23   ,   24
	subplot(221)	
	plot(((y[:,7]/2.725)-1),absolute(y[:,0]/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'r')
	plot(((y[:,7]/2.725)-1),absolute((y[:,1])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'b')
	plot(((y[:,7]/2.725)-1),absolute((y[:,2])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'c')
	plot(((y[:,7]/2.725)-1),absolute((y[:,3])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'m')
	plot(((y[:,7]/2.725)-1),absolute((y[:,4])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'g')
	
	##axes scales and limits	
	pyplot.xscale('log')
	pyplot.yscale('log')
	plt.gca().invert_xaxis()
	ylim(10**-25,10)
	xlim(Z_initial,Z_final+1)
	pyplot.yticks([10**-25,10**-20,10**-15,10**-10,10**-5,10**0])
	
	##Font settings for title and labels
	rc('font',**{'family':'serif','sans-serif':['Helvetica'] ,'size': 20 })
	
	##create a legend and plot it
	labels = [r'$H$', r'$H^+$', r'$H^-$', r'$H_2$', r'$H_2^+$']
	legend_location = 'lower left'
	colors = ['r','b','c','m','g']
	
	plt.legend(labels, loc=legend_location, prop={'size':10})
	
	##The same is repeated for each chemical species, as shown below.
	
	##Plot Deuterium Species
	subplot(222)
	plot(((y[:,7]/2.725)-1),absolute((y[:,9])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'r')
	plot(((y[:,7]/2.725)-1),absolute((y[:,10])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'b')
	plot(((y[:,7]/2.725)-1),absolute((y[:,11])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'c')
	plot(((y[:,7]/2.725)-1),absolute((y[:,12])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'k')
	plot(((y[:,7]/2.725)-1),absolute((y[:,13])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'g')

	pyplot.xscale('log')
	pyplot.yscale('log')
	plt.gca().invert_xaxis()
	ylim(10**-25,10)
	xlim(Z_initial,Z_final+1)
	pyplot.yticks([10**-25,10**-20,10**-15,10**-10,10**-5,10**0])

	rc('font',**{'family':'serif','sans-serif':['Helvetica']})
	
	labels = [r'$D$', r'$D^-$', r'$D^+$', r'$HD$', r'$HD^+$']
	
	legend_location = 'lower left'
	
	colors = ['r','b','c','m','g']
	
	plt.legend(labels, loc=legend_location, prop={'size':10})
	
	##Begin plotting Helium Species

	subplot(223)

	plot(((y[:,7]/2.725)-1),absolute((y[:,14])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'r')
	plot(((y[:,7]/2.725)-1),absolute((y[:,15])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'b')
	plot(((y[:,7]/2.725)-1),absolute((y[:,17])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'k')
	
	pyplot.xscale('log')
	pyplot.yscale('log')
	plt.gca().invert_xaxis()
	ylim(10**-25,10**1)
	xlim(Z_initial,Z_final+1)
	pyplot.yticks([10**-25,10**-20,10**-15,10**-10,10**-5,10**0])

	rc('font',**{'family':'serif','sans-serif':['Helvetica'] ,'size': 20 })
	
	labels = [r'$He$', r'$He^+$', r'$HeH^{+}$']
	
	legend_location = 'lower left'
	
	colors = ['r','b','k']
	
	xlabel('Redshift (Z)', fontsize=24)
	
	plt.gca().xaxis.set_label_coords(1.10, -0.105)
	
	ylabel(r' log$_{10}$(Fractional Abundance)', fontsize=24)
	
	plt.gca().yaxis.set_label_coords(-0.15, 1.105)

	plt.legend(labels, loc=legend_location, prop={'size':10})
	
	##Begin plotting Lithium Species
	subplot(224)

	plot(((y[:,7]/2.725)-1),absolute((y[:,18])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'r')
	plot(((y[:,7]/2.725)-1),absolute((y[:,19])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'b')
	plot(((y[:,7]/2.725)-1),absolute((y[:,20])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'c')
	plot(((y[:,7]/2.725)-1),absolute((y[:,21])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'k')
	plot(((y[:,7]/2.725)-1),absolute((y[:,22])/(y[:,0]+y[:,1]+y[:,2]+(2*y[:,3])+(2*y[:,4])+y[:,9]+y[:,10]+y[:,11]+(2*y[:,12])+(2*y[:,13])+y[:,15]+y[:,16]+(2*y[:,17])+y[:,18]+y[:,19]+y[:,20]+(2*y[:,21])+(2*y[:,22])+(3*y[:,23])+y[:,14]+(2*y[:,24]))),'g')
	
	pyplot.xscale('log')
	pyplot.yscale('log')
	plt.gca().invert_xaxis()
	ylim(10**-25,10**1)
	xlim(Z_initial,Z_final+1)
	pyplot.yticks([10**-25,10**-20,10**-15,10**-10,10**-5,10**0])

	rc('font',**{'family':'serif','sans-serif':['Helvetica'] ,'size': 12 })
	
	labels = [r'$Li$', r'$Li^+$', r'$Li^-$', r'$LiH$', r'$LiH^+$']
	
	legend_location = 'lower left'
	
	colors = ['r','b','c','k','g']
	suptitle(r'Change in Number Density of Chemical Species between ' + str(Z_initial) + '$\leq$ Z $\leq$' + str(Z_final), fontsize=20)
	
	plt.legend(labels, loc=legend_location, prop={'size':10})
	
	show()		
	
	
	
	
	##Begin plotting of matter temperature evolution
	figure()
	plot(((y[:,7]/2.725)-1),absolute(y[:,8]),'r')
	
	pyplot.xscale('log')
	pyplot.yscale('log')
	plt.gca().invert_xaxis()
	xlim(Z_initial,Z_final+1)
	pyplot.yticks([10**-25,10**-20,10**-15,10**-10,10**-5,10**0])

	rc('font',**{'family':'serif','sans-serif':['Helvetica']})
	
	title(r'The change in gas temperature between ' + str(Z_initial) + ' $\leq$ Z $\leq$ '+ str(Z_final))
	xlabel('Redshift (Z)')

	
	ylabel(r'Gas temperature (K)')

	#show() 				##disabled by default, as this plot is not commonly used apart from checking it has the correct form.
		
	exit = raw_input("\nProcess complete. Press any key to exit...")
	
	raise SystemExit