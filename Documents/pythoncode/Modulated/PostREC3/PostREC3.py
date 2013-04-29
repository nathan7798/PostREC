#PostREC3.0.py

##Nathan Bush, April 2013

#A python program designed to process user defined reaction equations, organise them into an array,
#calculate the rate of change of each reactant species, and integrate them over a certain time period.
#Results can be returned numerically and displayed graphically.

##DEPENDANCIES: The following modules and files are needed in the same directory for this function to
##run properly

##Sundials and assimulo need to be installed if running as python files
##Sundials: https://computation.llnl.gov/casc/sundials/main.html
##Assimulo: http://www.jmodelica.org/assimulo
##The programs also need: SCIPY, NUMPY and MATPLOTLIB

#####Python Modules#####

## system_arguments.py                     Evaluates system arguments and processes them for the main program
## new_ageconverter.py                     Calculates the age of the Universe for a given redshift value
## reaction_and_rate_selector.py           Calculates the required reactions and rates that are passed to the solver
## plotter.py                              Plots the results

#####Data files #####

## data.dat                                Has the ionised fraction and matter to radiation temperature ratio stored for a 
##                                         Range of redshift values

#####Pickled Files#####                    Information on pickled files:http://docs.python.org/2/library/pickle.html

##cooling_rates.p                          Pickled values for the cooling rates
##photoionisation.p                        Pickled values for photoionisation cross sections
##previous_use.p                           Information on the previous run of the solver (used to save time at points)
##reaction_equations.p                     Pickled reaction equations ready to solve
##reactionratearray.p                      Pickled array containing the reaction rates to be used

##FUNCTION INPUTS: The following is a list of the function inputs and thier descriptions

##Plot_request:                            Whether or not the results need to be plotted (1 = No, 2 = Yes)

##Z_initial:                               Initial redshift value to integrate between. (No larger than 2000)

##Z_final:                                 End point of the integration (No smaller than 1)

##step_size:                               Step size for the integration (A step size larger than 2 causes instability, lower than 0.1
##                                         means increased run-time. Between 0.1-1 is recommended)

##number_of_particles:                     Number density of particles to use, defualt is calculated using comsological model parameters

##manual_matter_temperature:               Value for manual matter temperature, otherwise it is set using cosmological model parameters

##source_function:                         Source function for radiation flux (optional)

##chemical_model:                          Chemical model to use for the integration. 1 = full model 2 = minimal model 3 = Abel et al 1996
##                                         4 = Galli & Palla 1998

def PostREC(plot_request,Z_initial,Z_final,step_size,number_of_particles, manual_matter_temperature, source_function,chemical_model):
	
	global c,kb,h_p,H0,W0,T_0,kbev,Z,cooling_function,flux,k,T_m, hev	###Global parameters defined for use throughout the function
	                                                                    ###Mianly cosmological parameter values, numerical constants etc.

	c = 2.99792458*(10**8)			#Speed of light     (SI units)
	kb = 1.38065*(10**-23)			#Boltzmann constant (SI units)
	h_p = 6.626068*(10**-34)		#Planck's constant  (SI units)
	H0 = 69.7						#100h
	W0 = 1.0						#Total matter density (flat universe)
	T_0 = 2.725						#Modern day temperature of the CMB
	kbev = 8.6173324*(10**-5)       #Boltzman constant in eV
	hev = 4.135667*(10**-15)        #Plancks constant in eV

	##################################    Program Functions    #################################################
	
	##Residual is the main function that the solver, IDA, uses to solve the system of equations. Each "residual"
	##(res_1, res_2 etc) corresponds to a certain chemical species. The reaction eqautions array contains the 
	##reaction equation for that particular species. Values for the rates are substituted in (k) along with the
	##matter and radiation temperatures. IDA then trys to make the value of the residual as close to 0 as possible
	##by doing so it finds the equilibrium state of the system.

	##list of chemical species and the array elements that each corresponds to
	
	#species_list = ["[H]","[H+]","[H-]","[H2]","[H2+]","[e-]","[hv]","[D]","[D-]","[D+]","[HD]","[HD+]","[He]","[He+]","[He++]","[HeH+]",[Li], [Li+], [Li-], [LiH], [LiH+], [H2D+], [D2]]
			      #    0     1      2      3      4       5      6      9    10     11      12      13     14     15       16      17      18     19    20      21     22  ,  23   ,   24
	
	def residual(t,y,yd):
				
		res_0 = reaction_equations[0](k,Z,y[8],y[7],y) - yd[0]
		
		res_1 = reaction_equations[1](k,Z,y[8],y[7],y) - yd[1]
		
		res_2 = reaction_equations[2](k,Z,y[8],y[7],y) - yd[2]
		
		res_3 = reaction_equations[3](k,Z,y[8],y[7],y) - yd[3]
		
		res_4 = reaction_equations[4](k,Z,y[8],y[7],y) - yd[4]
		
		res_5 = reaction_equations[5](k,Z,y[8],y[7],y) - yd[5]

		photon_density = 0.0 -yd[6]
		
		radiation_temp = 0.0 -yd[7]
		
		matter_temp = 0.0 -yd[8]

		res_9 = reaction_equations[7](k,Z,y[8],y[7],y) - yd[9]
		
		res_10 =reaction_equations[8](k,Z,y[8],y[7],y) -yd[10]
		
		res_11 = reaction_equations[9](k,Z,y[8],y[7],y)-yd[11]
		
		res_12 =  reaction_equations[10](k,Z,y[8],y[7],y) - yd[12]
		
		res_13 =  reaction_equations[11](k,Z,y[8],y[7],y) -yd[13]
		
		res_14 = reaction_equations[12](k,Z,y[8],y[7],y)-yd[14]
		
		res_15 =  reaction_equations[13](k,Z,y[8],y[7],y) - yd[15]
		
		res_16 =  reaction_equations[14](k,Z,y[8],y[7],y)- yd[16]
		
		res_17 = reaction_equations[15](k,Z,y[8],y[7],y) -yd[17]
		
		res_18 = reaction_equations[16](k,Z,y[8],y[7],y) -yd[18]
		
		res_19 = reaction_equations[17](k,Z,y[8],y[7],y) -yd[19]
		
		res_20 = reaction_equations[18](k,Z,y[8],y[7],y) -yd[20]
		
		res_21 =reaction_equations[19](k,Z,y[8],y[7],y) -yd[21]
		
		res_22 =reaction_equations[20](k,Z,y[8],y[7],y)-yd[22]
	
		res_23 =reaction_equations[21](k,Z,y[8],y[7],y)-yd[23]
		
		res_24 =reaction_equations[22](k,Z,y[8],y[7],y)-yd[24]
		
		return array([res_0,res_1,res_2,res_3,res_4,res_5,photon_density,radiation_temp,matter_temp,res_9,res_10,res_11,res_12,res_13,res_14,res_15,res_16,res_17,res_18,res_19,res_20,res_21,res_22,res_23,res_24])	
	
	##deriv provides a jacobian matrix for IDA so that it has apprpriate initial conditions and can solve the equations more efficiently.
	
	def deriv(y,t): 

		res_0 = reaction_equations[0](k,Z,y[8],y[7],y)
						
		res_1 = reaction_equations[1](k,Z,y[8],y[7],y)
		
		res_2 = reaction_equations[2](k,Z,y[8],y[7],y)
		
		res_3 = reaction_equations[3](k,Z,y[8],y[7],y)
		
		res_4 = reaction_equations[4](k,Z,y[8],y[7],y)
		
		res_5 = reaction_equations[5](k,Z,y[8],y[7],y)

		photon_density = 0.0
		
		radiation_temp = 0.0
		
		matter_temp = -cooling_function(y[8],y[7],y[0],y[3],y[12],y[5],y[1],y[4],y[14],y[15],y[16],y[21])

		res_9 = reaction_equations[7](k,Z,y[8],y[7],y) 
		
		res_10 =reaction_equations[8](k,Z,y[8],y[7],y) 
		
		res_11 = reaction_equations[9](k,Z,y[8],y[7],y)
		
		res_12 =  reaction_equations[10](k,Z,y[8],y[7],y) 
		
		res_13 =  reaction_equations[11](k,Z,y[8],y[7],y)
		
		res_14 = reaction_equations[12](k,Z,y[8],y[7],y)
		
		res_15 =  reaction_equations[13](k,Z,y[8],y[7],y)
		
		res_16 =  reaction_equations[14](k,Z,y[8],y[7],y)
		
		res_17 = reaction_equations[15](k,Z,y[8],y[7],y)
		
		res_18 = reaction_equations[16](k,Z,y[8],y[7],y)
		
		res_19 = reaction_equations[17](k,Z,y[8],y[7],y)
		
		res_20 = reaction_equations[18](k,Z,y[8],y[7],y)
		
		res_21 =reaction_equations[19](k,Z,y[8],y[7],y) 
		
		res_22 =reaction_equations[20](k,Z,y[8],y[7],y)
		
		res_23 =reaction_equations[21](k,Z,y[8],y[7],y)
		
		res_24 =reaction_equations[22](k,Z,y[8],y[7],y)
		
		return array([res_0,res_1,res_2,res_3,res_4,res_5,photon_density,radiation_temp,matter_temp,res_9,res_10,res_11,res_12,res_13,res_14,res_15,res_16,res_17,res_18,res_19,res_20,res_21,res_22,res_23,res_24])	
		
	##The error_checking function performs an error test on the results. It does this by comparing the final total number of particles to the initial number.
	##If there is a large discrepancy, there is a bug somewhere within the code or the solution is unstable. The user is then warned.

	def error_checking(error_test,end_line,Z_initial,Z_final):
				
		atoms_lost = error_test*(((1+Z_final)**3)/((1+Z_initial)**3))-(end_line[0]+end_line[1]+end_line[2]+(2*end_line[3])+(2*end_line[4])+end_line[9]+end_line[10]+end_line[11]+(2*end_line[12])+(2*end_line[13])+end_line[15]+end_line[16]+(2*end_line[17])+end_line[18]+end_line[19]+end_line[20]+(2*end_line[21])+(2*end_line[22])+(3*end_line[23])+end_line[14]+(2*end_line[24]))

		return(atoms_lost)
	
	##The initial conditions function calculates the initial conditions ready for integration
	
	def initial_conditions(Z_initial,number_of_particles, manual_matter_temperature):		#This function returns the initial conditions of the chosen epoch ready to be passed to the numerical integrator
		
		no_H = number_of_particles/((2.68*(10**-5))+0.240+(10**-9.9)+1.0)           #Calculate the hydrogen fraction using the fractional abundance of other elements
		
		no_He = 0.240*no_H                                                          #Calculate the fractional abundances of elements using nucelosynthesis results
		no_Deu = 2.68*(10**-5)*no_H
		no_Lithium = (10**-9.9)*no_H
				
		recombination_data = loadtxt('data.dat', delimiter=' ', dtype = (str))		#Loads in an array giving us the recombination fraction and ratio of matter to radiation temperature
		
		with open('data.dat','r') as f:
			rows=[map(float,L.strip().split(' ')) for L in f]  
		
		arr=array(rows)  
		
		free_electron_fraction = arr[8000-Z_initial][1]                             #Reads in the free electron fraction from the loaded data
	
		matter_radiation_temperature_ratio = arr[8000-Z_initial][2]                 #Reads the matter radiation temperature ratio from the loaded data
		
		radiation_temperature = (2.725)*(1+Z_initial)	                            #Calculates the initial radiation temperature using the initial redshift

		if manual_matter_temperature == 0:                                          #Calculates the matter temperature by first checking if a manual
		                                                                            #value was entered, if not then it is calculated from the ratio.
			matter_temperature = matter_radiation_temperature_ratio*radiation_temperature
		else:
			matter_temperature = manual_matter_temperature
				
		no_e = no_H*free_electron_fraction                                          #Calculate number of electrons from free electron fraction

		if free_electron_fraction > 1.0:                                            #Calculate the ionised fraction of elements using ionised fraction.
			no_He_plus = (free_electron_fraction - 1.0)*no_He
			no_He = no_He - no_He_plus
			no_H_plus = no_H
			no_H = 0.0
			no_Deu_plus = no_Deu
			no_Deu = 0.0
		else:
			no_He_plus = 0.0 
			no_H_plus = free_electron_fraction*no_H
			no_H = no_H - no_H_plus
			no_Deu_plus = free_electron_fraction*no_Deu
			no_Deu = no_Deu - no_Deu_plus
		
		photon_density = (number_of_particles*((10**10)/4.5))                       # Calculates the number density of photons using the photon to baryon ratio

		return ([no_H,no_H_plus,0,0,0,no_e,photon_density,radiation_temperature,matter_temperature,no_Deu,0.0,no_Deu_plus,0.0,0.0,no_He,no_He_plus,0.0,0.0,0.0,no_Lithium,0.0,0.0,0.0,0.0,0.0],(no_H+no_H_plus+no_Deu+no_Deu_plus+no_He+no_He_plus+no_Lithium))
	
	##It retuns the abundances of each chemical species, in accordance with the labelling below. It also returns number densities to be used by the error checking function later on.
	
	#species_list = ["[H]","[H+]","[H-]","[H2]","[H2+]","[e-]","[hv]","[D]","[D-]","[D+]","[HD]","[HD+]","[He]","[He+]","[He++]","[HeH+]",[Li], [Li+], [Li-], [LiH], [LiH+], [H2D+], [D2]]					
				  #    0     1      2      3      4       5      6      9    10     11      12      13     14     15       16      17      18     19    20      21     22  ,  23   ,   24
	
	def photo_ionisation_calculator(source_function): ##This function calculates the effect of photoionising flux on the pre-galacti gas
		
		##refer to reaction_and_rte_writer.py for more info and sources for these rates
		
		
		flux =  [0] * 12     ##Prepare list ready to be filled with values

		##Reaction H   + hv -> H+   + e-
		
		##Integrate the cross section over the required frequency range through use of the source function, repeat for other reactions
		flux[0] = quad(lambda v: ((6.30*(10**-18))*(((3.28846551*(10**15))/v)**4.0)*((exp(4-(4*atan((((v/(3.28846551*(10**15)))-1.0)**0.5))/(((v/(3.28846551*(10**15)))-1)**0.5))))/(1-exp((-2*pi)/(((v/(3.28846551*(10**15)))-1)**0.5)))))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (3.28846551*(10**15)), 10**20)[0]
		
		##Reaction He+ + hv -> He++ + e-
		
		flux[1] = quad(lambda v: ((6.30*(10**-18))*(((3.28846551*(10**15))/v)**4.0)*((exp(4-(4*atan((((v/(3.28846551*(10**15)))-1)**0.5))/(((v/(3.28846551*(10**15)))-1)**0.5))))/(1-exp((-2*pi)/(((v/(3.28846551*(10**15)))-1)**0.5)))))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (3.28846551*(10**15)), 10**20)[0]

		##Reaction He +  hv -> He+  + e-
		
		flux[2] = quad(lambda v: (7.42*(10**-18)*(((1.66*(v/(24.6/(4.13566*(10**-15))))**-2.05))-(0.66*((v/(24.6/(4.13566*(10**-15))))**-3.05))))*(4.0*pi*((h_p*v)**-1))*eval(source_function),(5.948253*(10**15)), 10**20)[0]

		##reaction H- +  hv -> H    + e-
		
		flux[3] = quad(lambda v: (7.928*(10**5)*((v-(0.755/(4.13566*(10**-15))))**1.5)*(v**-3))*(4.0*pi*((h_p*v)**-1))*eval(source_function),(1.8255*(10**14)), 10**20)[0]

		##reaction D + hv -> D+ + e-
		
		flux[4] = quad(lambda v: ((6.30*(10**-18))*(((3.28846551*(10**15))/v)**4.0)*((exp(4-(4*atan((((v/(3.28846551*(10**15)))-1.0)**0.5))/(((v/(3.28846551*(10**15)))-1)**0.5))))/(1-exp((-2*pi)/(((v/(3.28846551*(10**15)))-1)**0.5)))))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (3.28846551*(10**15)), 10**20)[0]

		##reaction: H2 +  hv -> H2+  + e- (predisscoaition reaction) there are 3 elements for this element, one for 
		##the effect of each frequency band
		

		flux[5] = quad(lambda v: (6.2*(10**-18)*((4.135667516*(10**-15))*v)-(9.4*(10**-17)))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (3.728539*(10**15)),(3.9896824*(10**15)))[0]
		
		flux[6] = quad(lambda v: (1.4*(10**-18)*((4.135667516*(10**-15))*v)-(1.48*(10**-17)))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (3.9896824*(10**15)), (4.2798411*(10**15)))[0]
		
		flux[7] = quad(lambda v: (2.5*(10**-14)*(((4.135667516*(10**-15))*v)**-2.71))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (4.2798411*(10**15)), 10**20)[0]
		
		##Reaction: H2+ + hv -> H + H+
		
		flux[8] = quad(lambda v: (10**((-1.6547717*(10**6))+(1.866033*(10**5)*log(v))-(7.8986431*(10**3)*(log(v)**2))+(148.73693*(log(v)**3))-(1.0513032*(log(v)**4))))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (6.407671772*(10**14)), 10**20)[0]
	
		##Reaction: H2+ + hv -> 2H+  + e-
		
		flux[9] = quad(lambda v: (10**((-16.926)-((4.528*(10**-2))*((4.135667516*(10**-15))*v))+(2.238*(10**-4)*(((4.135667516*(10**-15))*v)**2))+(4.425*(10**-7)*(((4.135667516*(10**-15))*v)**3))))*(4.0*pi*((h_p*v)**-1))*eval(source_function), (7.253968*(10**15)), (2.176190413*(10**16)))[0]
		
		
		##reaction: H2  + hv -> H2* -> H + H"
		
		source_function_jw = source_function
		
		v = 12.87/hev
		j_w = eval(source_function)
		del v
		
		flux[10] = ((1.1*(10**8)*j_w))
		
		##reaction H2 +  hv -> H2+  + e- (Direct photodissociation into lyman warner continuum bands)
		##different cross sections need to be used fordifferent frequeny bands, these are shown below, along with the limits (in eV).
		##These are all taken from tom Abels 1996 paper "Modelling Primordial gas in Numerical Cosmology"
		
		##Total ionisation rate is the sum of the ionisation from each frequency band.
		
		sigmaL01 = lambda v: ((10**-18)*(10**(15.1289-1.05139*hev*v))) ##14.675 < hv < 16.820
		sigmaL02 = lambda v:((10**-18)*(10**(-31.41 + (1.8042*(10**-2)*((hev*v)**3)) - (4.2339*(10**-5))*((hev*v)**5))))  ##16.820 < hv < 17.6
		sigmaW0  = lambda v:((10**-18)*(10**(13.5311-(0.9182618*hev*v)))) ##14.675 < hv < 17.7
		sigmaL1 = lambda v:((10**-18)*(10**(12.0218406-(0.819429*hev*v)))) ## 14.159 < hv < 15.302
		sigmaL2 = lambda v:((10**-18)*(10**(16.04644-(1.082438*hev*v)))) ##15.302 < hv < 17.2
		sigmaW1 = lambda v:((10**-18)*(10**(12.87367-(0.85088597*hev*v)))) ##14.159 < hv < 17.2

		band1 = quad(lambda v: ((0.25*(sigmaL01(v)+sigmaW0(v))+(0.75*sigmaL1(v)+sigmaW1(v)))*(4.0*pi*((h_p*v)**-1))*eval(source_function)), (3.4236*(10**15)),(3.7000*(10**15)))[0]
		band2 = quad(lambda v: ((0.25*(sigmaL01(v)+sigmaW0(v))+(0.75*sigmaL2(v)+sigmaW1(v)))*(4.0*pi*((h_p*v)**-1))*eval(source_function)), (3.7000*(10**15)),(4.067*(10**15)))[0]
		band3 = quad(lambda v: ((0.25*(sigmaL02(v)+sigmaW0(v))+(0.75*sigmaL2(v)+sigmaW1(v)))*(4.0*pi*((h_p*v)**-1))*eval(source_function)), (4.067*(10**15)), (4.2798*(10**15)))[0]
		
		flux[11] = band1+band2+band3
				
		return flux
		
	
	##The age_converter gives the age of the Universe as any redshift Z in Gyrs. This 
	##is converted into seconds using the formula below.
	
	def range_calculator(Z_initial,Z_final):

		range = (Z_initial-Z_final)*(10**9)*365.25*24*60*60

		return range
	
	def file_unpacking(chemical_model):	                                            ##Performs necessary file unpacking for reaction rates and cooling rates
		
		pickle_file0 = open("previous_use.p","rb")                                  ##Checks to see what the previous chemical model in use was. If it is the same as the current use
		                                                                            ##then the arrays are loaded straight away - nothing needs to be recalculated
		previous_model_used = pickle.load(pickle_file0)

		if (previous_model_used == chemical_model):
			pickle_file1 =open("reactionratearray.p","rb")
			pickle_file4 = open("reaction_equations.p","rb")
		else:                                                                       ##Otherwise, the chemical equations and rates are loaded in accordance with the model chosen.
			reaction_and_rate_selector.reaction_writer(chemical_model)
			pickle_file1 =open("reactionratearray.p","rb")
			pickle_file4 = open("reaction_equations.p","rb")
			pickle.dump(chemical_model,open("previous_use.p","wb"))					     

		array_size_lookup = {1:92  , 2:39 , 3 : 20, 4 : 42}                          ##Dictionary containing the size of the array for each model.
		
		array_size = array_size_lookup[chemical_model]                              ##Returns the size of the array for the chosen model.
		 
		rate_array = pickle.load(pickle_file1)                                      ##Loads the reactions, cooling rates, photoionisation rates and regular rates

		reactions = pickle.load(pickle_file4)

		return (rate_array,array_size,reactions)

	###################################    end of program functions     ##############################################	

	global reaction_equations ##Define a variable ready to be filled with the reaction equations
	
	rate_array, array_size,reactions  = file_unpacking(chemical_model)       ##Unpack the equations
	
	reaction_equations = [0] * len(reactions)      ##Turn the unpacked chemical equations into something the program can use,
	                                               ##In this case each reaction is converted into a lambda function
	for i in range(0,len(reactions)):
		if (reactions[i] == ''):
			reaction_equations[i] = lambda k,Z,T_m,T_r,y : 0.0
		else:
			reaction_equations[i] = eval('lambda k,Z,T_m,T_r,y: ' + reactions[i])

	
	k = [0] * array_size                          ##Turn the unpacked reaction rates into something the program can use,
	                                              ##lambda functions once again.
		
	for i in range(0,array_size):
		k[i] = eval(rate_array[i][9]) 
	
	##Now we need a cooling function for the gas, each term within the function is due to a different cooling/heating process
	
	cooling_function = lambda T_m,T_r,n_H,n_H2,n_HD,n_e,n_Hplus,n_H2plus,n_He,n_Heplus,n_Heplusplus,n_LiH: ((-1.017*(10**-37)*(T_r**4)*(T_m-T_r)*n_e*n_H)  + \
	                                                                                                 (1.27*(10**-21)*(T_m**0.5)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-157809.1/T_m))*n_e*n_H + \
                                                                                                     (9.38*(10**-22)*(T_m**0.5)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-285335.4/T_m))*n_e*n_He + \
                                                                                                     (4.95*(10**-22)*(T_m**0.5)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-631515.0/T_m))*n_e*n_Heplus + \
	                                                                                                 (8.70*(10**-27)*(T_m**0.5)*((T_m*10**-3)**-0.2)*((1+(T_m*10**-6)**0.7)**-1.0))*n_e*n_Hplus + \
																									 (1.55*(10**-26)*(T_m**0.3647))*n_e*n_Heplus + \
																									 (3.48*(10**-26)*(T_m**0.5)*((T_m*10**-3)**-0.2)*((1+(T_m*10**-6)**0.7)**-1.0))*n_e*n_Heplusplus + \
																									 (1.24*(10**-13)*(T_m**-1.5)*exp(-47000.0/T_m)*(1+0.3*exp(-94000.0/T_m)))*n_e*n_Heplus + \
																									 (7.50*(10**-19)*((1+(T_m*10**-5)**0.5)**-1)*exp(-118348.0/T_m))*n_e*n_H +\
																									 (5.54*(10**-17)*(T_m**-0.397)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-473638.0/T_m))*n_e*(n_Heplus) + \
																									 (1.42*(10**-27)*(T_m**0.5)*(1.10+0.34*exp(-((5.50-log10(T_m))**2)/3.0)))*(n_e*(n_Hplus+n_Heplus+4*n_Heplusplus)) + \
																									 (10**(-103.0 + (97.59*log10(T_m))-(48.05*((log10(T_m))**2.0))+(10.80*((log10(T_m))**3))-(0.9032*((log10(T_m))**4))))*n_H2*n_H + \
																									 (1.1*(10**-19)*(T_m**-0.34)*exp(-3025.0/T_m))*n_H2plus*n_e + \
																									 (1.36*(10**-22)*exp(-3152.0/T_m))*n_H2plus*n_H + \
																									 (10**(-36.42+5.95*log10(T_m)-0.526*((log10(T_m))**2)))*n_Hplus*n_H + \
																									 (10**(-42.45906+(21.90083*T_m)-(10.1954*(T_m**2))+(2.19788*(T_m**3))-(0.17286*(T_m**4))))*n_H*n_HD + \
																									 (10**(-31.47+(8.817*(log10(T_m)))-(4.1444*((log10(T_m))**2))+(0.8292*((log10(T_m))**3))-(0.04996*((log10(T_m))**4))))*n_LiH)/(kbev)
	
	
	##Determine whether or not a source function was entered, if it was, then find the relevant flux for the reactions
	
	if source_function != "0.0":
		flux = photo_ionisation_calculator(source_function)
	else:
		flux = [0] * 12
	
	#################################      BEGINNING OF INTEGRATION PROCEDURE      ###########################################################
	
	##Begin integration procedure
	
	Z = Z_initial
			
	yinit, error_test_begin = initial_conditions(Z_initial,number_of_particles,manual_matter_temperature) ##Get the initial conditions, and get the initial values ready 
	                                                                                                      ##for the error testing function to make use of
	yd0init = deriv(yinit,0.0)                                          ##Get the jacobian for IDA to use (faster solving)      
  
	Z_initial_s = new_ageconverter.age_converter(Z_initial)				##Convert Z value into age of universe

	Z_final_s = new_ageconverter.age_converter(Z_initial-step_size)		##Convert Z value into age of universe

	time_range = range_calculator(Z_final_s,Z_initial_s)				##Calculate the actual time in seconds to integrate over by using the step size.

	model = Implicit_Problem(residual,yinit,yd0init,0.0)                ##Set up the problem

	sim = IDA(model)                                                 

	sim.verbosity = 50                                                  ##Set verbosity to low(not needed)
	
	print float(Z)
	
	t,y,yd = sim.simulate(time_range)                                   ##Perform the simulation for one step
	
	eoa = (shape(y)[0])-1 
	
	end_line = y[(shape(y)[0])-1,:]       

	yd0init = yd[eoa,:]                                                  ##Provide the new jacobian from the newly calculated results, ready for the next step

	plot_values = y[eoa,:]

	previous_time_range = time_range
	
	for i in arange((Z_initial-step_size),Z_final,-step_size):      ##Z_final+step_size is used since it then subtracts the step size in the iteration
				
		##Recompute the flux for some reactions as Z has changed
				
		print i
		Z = i
		Z_initial_s = new_ageconverter.age_converter(i)				    ##Convert Z value into age of universe
		Z_final_s = new_ageconverter.age_converter(i-step_size)		    ##Convert Z value into age of universe
		time_range = range_calculator(Z_final_s,Z_initial_s)	        ##Calculate the actual time in seconds to integrate over

		
		
		 ##re-initialise the values ready for integration during the next step
		##each of the number densities is reduced by a factor of step-size**3
		##due to the expansion of the universe in this time
		##the radiation temperature and matter temperature are set according
		##to the redshift, over Z = 200 the two are coupled together, and so
		##are equal, after this the matter temperature evolves adiabtically
		##The cooling function also comes into effect at the lower Z values
		
		if i > 200:                                                    
			yinit = [end_line[0]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[1]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[2]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[3]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[4]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[5]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[6]*(((1.0+i)/(1.0+i+step_size))**4), \
					 (T_0*(1.0+i)), \
					 (T_0*(1.0+i)), \
					 end_line[9]*(((1.0+i)/(1.0+i+step_size))**3),  \
					 end_line[10]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[11]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[12]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[13]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[14]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[15]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[16]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[17]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[18]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[19]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[20]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[21]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[22]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[23]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[24]*(((1.0+i)/(1.0+i+step_size))**3)] 
		else:
			yinit = [end_line[0]*(((1.0+i)/(1.0+i+step_size))**3), \
			         end_line[1]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[2]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[3]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[4]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[5]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[6]*(((1.0+i)/(1.0+i+step_size))**4), \
					 (T_0*(1+i)), \
					 ((T_0*((1.0+i)**2))/201.0)-cooling_function(end_line[8],end_line[7],end_line[0],end_line[3],end_line[12],end_line[5],end_line[1],end_line[4],end_line[14],end_line[15],end_line[16],end_line[21]),\
					 end_line[9]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[10]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[11]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[12]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[13]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[14]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[15]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[16]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[17]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[18]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[19]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[20]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[21]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[22]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[23]*(((1.0+i)/(1.0+i+step_size))**3), \
					 end_line[24]*(((1.0+i)/(1.0+i+step_size))**3)] 

		yd0init = yd[eoa,:]                                                  ##Provide the new jacobian from the newly calculated results, ready for the next step
	
		model = Implicit_Problem(residual,yinit,yd0init,previous_time_range) ##Simulate again, using new step
		sim = IDA(model)
		sim.verbosity = 50
			
		t,y,yd = sim.simulate(time_range+previous_time_range)
	
		eoa = (shape(y)[0])-1                                      

		end_line = y[eoa,:]                                                  ##store plot values, alter time range and get ready for the next step.
		previous_time_range = time_range + previous_time_range
		plot_values = vstack((plot_values,y[eoa,:]))
		
	print("\nError test results:\n")
	print(error_checking(error_test_begin,end_line,Z_initial,i))		

	final_calculated_values = array(plot_values[(shape(plot_values)[0])-1][:])          ##Finds the final calculated abundacnes, ready to return to the user
	print("\n Final calculated abundances:\n")
	print end_line
	
	if plot_request == 2:																##If the user wanted to plot the results, this series of commands is run
		plotter.plotting_program(Z_initial,Z_final,plot_values)
			
	return (final_calculated_values) ##Return the final calculated values to the user

	
##Below is the start process for the function(its at the bottom of the code since the functions it makes use of need to be defined
##before it can be run). It checks if this function is being run directly from the command line (__name__ == '__main__'). If it is
##then it checks for input and also prompts the user for input if needed. It imports all the required modules, processes all of the
##inputs and then passes them to the main function (POSTREC).
	
if __name__ == '__main__':

	from numpy import array, shape, math, loadtxt, log10, vstack, arange
	from scipy.integrate import quad       
	from pylab import *                
	from numpy import pi as pi			 
	import system_arguments
	import new_ageconverter
	import reaction_and_rate_selector
	import plotter
	import sys
	from assimulo.solvers.sundials import IDA 
	from assimulo.problem import Implicit_Problem
	from math import exp, log10, fabs, atan, log
	import pickle
	
	input_test, integration_arguments = system_arguments.check_system_arguments()																				   
		
	if input_test == False:
		
		import argparse                     #Needed to check the system arguments
		
		##Below are a series of functions that check whether or not the input is sensible (prevents complicated errors
		##and possible file corruption later in the process)
		
		def check_negative(value):
			ivalue = int(value)
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			return value
			
		def check_chemical_model(value):
			value = int(value)
			if value != 1 and value != 2 and value != 3 and value != 4:
				raise argparse.ArgumentTypeError("%s is an invalid chemical model" % value)
			return value
			
		def check_integration_range(value):
			ivalue = int(value)
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			if value > 2000:
				raise argparse.ArgumentTypeError("%s is an invalid value, select between 2000 and 1" % value)
			return value
		
		def check_step_size(value):
			ivalue = int(value)
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			if value > 2 :
				print("\nStep size accepted, however be warned step sizes larger than 2 can lead to instability")
			if value < 0.01: 
				print("\nStep size accepted, however be warned small step sizes cause long completion times")
			return value
		
		def check_matter_temp_choice(value):
			if value != 1 and value != 2:
				raise argparse.ArgumentTypeError("%s is an invalid choice" % value)
			return value
		
		def check_matter_temp(value):
			ivalue = int(value)
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			if value > 10**6:
				print("Accepted, but be warned large values can cause instability")
			return value
		
		##Below are a series of prompts for the user, these are only activated if no input parameters are provided
		
		print("\n\nWelcome to to PostREC v.3.0, this program simulates the post recombination chemical evolution of the Universe.\n\nFor help running the program run 'PostREC2.0.py -h' " )
			
		chemical_model = input("\nSelect your chemical model:\n\nFull Model = 1\nMinimal Model = 2\nAbel et al 1996 = 3\nGalli & Palla 1998 = 4\n\n")
		
		check_chemical_model(chemical_model)
		
		Z_initial = input("\nPlease enter the startpoint of the integration(Z): \n")
		
		check_integration_range(Z_initial)
		
		Z_final = input("\nPlease enter the endpoint of the integration(Z): \n")
		
		check_integration_range(Z_final)

		step_size = input("\nEnter your integration step size:\n")
		
		check_step_size(step_size)
		
		matter_temp_choice = input("\nWould you like to enter a value for the matter temperature?  (1 = no 2 = yes)\n")
		
		check_matter_temp_choice(matter_temp_choice)
		
		if matter_temp_choice == 2:
			manual_matter_temperature = input("\nEnter the matter temperature you would like to start with:\n")
		else:
			manual_matter_temperature = 0	
		
		check_matter_temp(manual_matter_temperature)
		
		source_function_choice = input("\nWould you like to enter a source function for external flux?  (1 = no 2 = yes)\n")
		
		if source_function_choice == 2:
			source_function = raw_input("\nEnter the source function you would like to use\n")
		
		else:
			source_function = "0.0"
			
		particles_choice = input("\nWould you like to enter a custom number density? (1 = no 2 = yes)\n")
		
		check_matter_temp_choice(particles_choice)
		
		if particles_choice ==2:
			number_of_particles = input("\nEnter the number density (cm-3):\n")
		else:
			number_of_particles = ((1+Z_initial)**3)*(3.66*(10**-3)*4.5*(0.67**-2)*(9.19*(10**-6)*(0.67**2)))			
		
		check_matter_temp(number_of_particles)
		
		plot_request = input("\nPlot the results? (1 = no 2 = yes):\n")
		
		check_matter_temp_choice(plot_request)
		
	else:
		
		chemical_model = integration_arguments[0]     ##Function mode runs automatically with the minimal model
		
		Z_initial = integration_arguments[1]          ##Proces the inputs
		Z_final = integration_arguments[2]
		step_size = integration_arguments[3]
			
		if integration_arguments[4] == 'D':           ##D stands for default
			manual_matter_temperature = 0
		else:
			manual_matter_temperature = integration_arguments[3]
			
		if integration_arguments[5] == 'D':
			number_of_particles = number_of_particles = ((1+Z_initial)**3)*(3.66*(10**-3)*4.5*(0.67**-2)*(9.19*(10**-6)*(0.67**2)))
		else: 
			number_of_particles = eval(integration_arguments[4])
			
		if integration_arguments[6] == 'D':
			source_function = "0.0"

		else:
			source_function = integration_arguments[6]
			
		plot_request = integration_arguments[7]
				
	##The input arguments have been processed, now run the function
	PostREC(plot_request,Z_initial,Z_final,step_size,number_of_particles, manual_matter_temperature, source_function,chemical_model)
	
	exit = raw_input("\nProcess complete. Press any key to exit\n")
	
	raise SystemExit