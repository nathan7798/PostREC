#Python Module to check for input system arguments, it takes care of all of the input values and processes them ready for use by
#POSTREC. It also takes care of any "nonsense values" ie values that are not feasible to use, ormay have been entered in error.
#It also provides help messages for how to use the program in "function form" where it can be run directly from the terminal/
#command line

##Information on argparse: http://docs.python.org/3.4/library/argparse.html

import sys							#Needed to check the input commands
import argparse                     #Needed to check the system arguments

def check_system_arguments():		#All arguments are put under a single function so that it can be called from the main program.
	
	sys_array = [True, False , 0, False, False]
	
	if (len(sys.argv) > 1):         ##If some input arguments have been used, check them
		
		def check_negative(value):  ##This function checks if a negative value has been entered
			
			ivalue = int(value)     ##converts the input into a numerical value (was a string)

			if value == "D" :       ##default is accepted
				return value
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			return value
			
		def check_integration_range(value):  ##Integration range is a max of 2000 to 0, this
			ivalue = int(value)          ##function checks the range
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			if ivalue > 2000:
				raise argparse.ArgumentTypeError("%s is an invalid value, select between 2000 and 1" % value)
			return float(value)
		
		def check_step_size(value):     ##Gives warning about using small step sizes
			value = eval(value)
			if value < 0.01: 
				print("\nStep size accepted, however be warned small step sizes cause long completion times")
			return float(value)
		
		def check_matter_temp_choice(value):  
			value = int(value)
			if value != 1 and value != 2:
				raise argparse.ArgumentTypeError("%s is an invalid choice" % value)
			return value
			
		def check_chemical_model(value):
			value = int(value)
			if value != 1 and value != 2 and value != 3 and value != 4:
				raise argparse.ArgumentTypeError("%s is an invalid chemical model" % value)
			return value			
		
		def check_matter_temp(value):
			if value == "D" :
				return value
			ivalue = int(value)
			if ivalue < 0:
				raise argparse.ArgumentTypeError("%s is an invalid positive integer" % value)
			if value > 10**6:
				print("Accepted, but be warned large values can cause instability")
			return value
		
		
		##Below the possible arguments to the function are represented as classes, each with different attributes.
		##Each variable has a "help" parameter, which can give the user insight into how to use the function
		
		parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, description='Inputs ready for integration, note to use a default setting for a particular variable, type "D"')

		parser.add_argument('chemical_model',metavar='chemical_model',type=check_chemical_model,
							 help ='Chemical Model to use. 1 = full model 2=minimal model (faster) 3 = Abel et al 1996 4 = Galli & Palla 1998')		
	
		parser.add_argument('Start_point',metavar='Z_initial',type=check_integration_range,
							 help ='Start point of integration')

		parser.add_argument('End_point', metavar='Z_final', type=check_integration_range,
						   help='End point of the integration')

		parser.add_argument('Step_size', metavar='Step_size', type=float,
						   help='Step size used within the integration')				   
   
		parser.add_argument('Matter_temperature', metavar='Matter_temperature', type=check_matter_temp,
						   help='Manual matter temperature at start of integration')				   
						   
		parser.add_argument('Particle_density', metavar='Particle_density', type=check_matter_temp,
						   help='Particle density used in integration (cm3)')	
		
		parser.add_argument('Plotting', metavar='Plotting', type=int, choices = {1,2},
						   help='Plot the results? 1 = no 2 = yes')
						   
		parser.add_argument('Source_function', metavar='Source_function', type=str,
						   help='Radiation source function to be used in the integration')

		parser.add_argument('Startflux', help='Startpoint of the flux spectrum (must be positive)')

		parser.add_argument('Endflux', help='endpoint of the flux spectrum (must be positive)')

		args = parser.parse_args()  ##creates all possible arguments defined above
		
		##If arguments are present on the command line, they are processed and passed back to PostREC
				
		return [True,[args.chemical_model,args.Start_point,args.End_point,args.Step_size, args.Matter_temperature, args.Particle_density ,args.Source_function,args.Plotting]]

	else:
		##If not, the boolean "False" value prompts PostREC into entering verbose mode, where the user is prompted for input.																		
		return [False,0]