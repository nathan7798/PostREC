###This is an external function which calculates the reaction equations for each chemical species depending on the "model" chosen by the User. 
##It then loads the appropriate reaction rates to use for the model, based on either what previous authors used or through selection according
##to the project report. The input is "chemical model" - a number between 1 and 4 which is tagged to a certain model to use. The key is:
## 1 The Full model, as detailed in the project report.
## 2 The minimal model, as detailed in the project report.
## 3.The model used by Abel et al (1996) in thier paper "Modelling primordial gas in numerical cosmology"
## 4.The model used by Galli & Palla (1998) in thier paper "The chemistry of the early Universe"

##The function is essentially a series of arrays, the correct array is chosen according to the chemical model chosen in "PostREC.py". 
##The reaction equations are then calculated, and stored in a pickle file. The rates are then also stored in a pickle file.

##Array information:

##Column######Information
#  1            This is the reaction number, this is not used explicitly anywhere, it is essentially just a reference number to keep track of
#               of which reactions are used and where
#  2            The reaction in question in text form. Not used anywhere explicitly, but useful for identification purposes and checking for 
#               possible errors
#  3-8          These columns contain the chemical species involved for the reaction shown by column two. Only one chemical species can be assigned per column. 
#               Columns 3-5 conatin the reactants for the reaction, columns 6-8 contain the prducts of the reaction. They are each enclosed in square brackets for 
#               clarity.
#  9            The rate reference for the reaction. The rates are contained with an array (rate_array), each element in the array corresponds to a different rate.
#               The rates themselves are lambda functions (for portability) hence the appropriate inputs are also shown in this column. For example, y[8] is the 
#               matter temperature, which many rates make use of. Others make use of the radiation temperature y[7], some make use of the redshift Z. Note that also
#               some rates use the temperature in units of eV, in which case a conversion is made before substitution into the function. (multiply by the boltzmann constant
#               that is in units of electron volts (kbev).
#  10           The actual rate for the reaction, in the form of a lambda function.
#  11           The source of the rate.
#  12           Range of applicability of the rate, 0-8000K is standard.

from numpy import array, shape		##Import for array functionality									
import pickle                       ##Needed to store files

##Begin function

def reaction_writer(chemical_model):
	
	##Perform comparison to use the array appropriate to the selected chemical model
	
	if (chemical_model ==1) :

		reaction_matrix=array([
							   [ 1,"H   + e- ->  H-  + hv"   , "[H]"  , "[e-]" , ""      , "[H-]"  , "[hv]" , ""     , "k[0](y[8])", "lambda T_m: 1.4*(10**-18)*(T_m**0.928)*exp(-T_m/16200.0)"                                          , "Galli D & Palla 1998"                                         , "0-8000"], #chemical equations we are interested in.											
							   [ 2,"H-  + H  ->  H2  + e-"    ,"[H-]" , "[H]"  , ""      , "[H2]"  , "[e-]" , ""     , "k[1](y[8])", "lambda T_m: 4.0*(10**-9)*(T_m**-0.17) if T_m > 300 else (1.5*10**-9)" 					         , "Galli D & Palla 1998"                                         , "0-8000"],
							   [ 3,"H+  + H  ->  H2+ + hv"    ,"[H+]" , "[H]"  , ""      , "[H2+]" , "[hv]" , ""     , "k[2](y[8])", "lambda T_m: 10**(-19.38-(1.523*log10(T_m))+(1.118*((log10(T_m))**2))-(0.1269*((log10(T_m))**3)))" , "Glover & Abel 2008"                                           , "0-8000"], 
							   [ 4,"H2+ + H  ->  H2  + H+"    ,"[H2+]", "[H]"  , ""      , "[H2]"  , "[H+]" , ""     , "k[3](y[8])", "lambda T_m: 6.4*(10**-10)"					                                                       , "Galli D & Palla 1998"                                         , "0-8000"], 
							   [ 6,"H2  + H+ ->  H2+ + H "    ,"[H2]" , "[H+]" , ""      , "[H2+]" , "[H]"  , ""     , "k[4](y[8])", "lambda T_m: 3.0*(10**-10)*exp(-21050.0/T_m)"				                                       , "Galli D & Palla 1998"                                         , "0-8000"], 
							   [ 7,"H2  + e- ->  2H  + e-"    ,"[H2]" , "[e-]" , ""      , "[H]"   , "[H]"  , "[e-]" , "k[5](y[8])", "lambda T_m: 1.91*(10**-9)*(T_m**0.136)*exp(-53407.1/T_m)"                                          , "Trevian, C.S & Tennyson J.2002"                               , "0-8000"], 
							   [ 8,"H2  + H  ->  2H  + H "    ,"[H2]" , "[H]"  , ""      , "[H]"   , "[H]"  , "[H]"  , "k[6](y[8])", "lambda T_m: 6.67*(10**-12)*(T_m**0.5)*exp(-(1.0+(63593.0/T_m)))"                                   , "Dove and Mandy 1986"                                          , "0-8000"], 
							   [ 9,"H-  + e- ->  H   + 2e-"   ,"[H-]" , "[e-]" , ""      , "[H]"   , "[e-]" , "[e-]" , "k[7](y[8]*kbev)", "lambda T_m: exp((-18.01849334)+(2.3608522*log(T_m))-(0.28274430*((log(T_m))**2))+(1.62331664*(10**-2)*(log(T_m))**3)-(3.36501203*(10**-2)*((log(T_m))**4))+(1.17832978*(10**-2)*((log(T_m))**5))-(1.65619470*(10**-3)*((log(T_m))**6))+(1.06827520*(10**-4)*((log(T_m))**7))-(2.63128581*(10**-6)*((log(T_m))**8)))", "Glover & Abel 2008", "0-8000"],
							   [10,"H-  + H  ->  2H  + e-"    ,"[H-]" , "[H]"  , ""      , "[H]"   , "[H]"  , "[e-]" ,"k[8](y[8]*kbev)","lambda T_m: 2.5634*(10**-9)*(T_m**1.78186) if T_m < 0.1 else exp((-20.37260896) + (1.13944933*log(T_m)) - (0.14210135*((log(T_m))**2)) + (8.4644554*(10**-3)*(log(T_m))**3) - (1.4327641*(10**-3)*((log(T_m))**4)) + (2.0122503*(10**-4)*((log(T_m))**5)) + (8.6639632*(10**-5)*((log(T_m))**6)) - (2.5850097*(10**-5)*((log(T_m))**7)) + (2.4555012*(10**-6)*((log(T_m))**8)) - (8.0683825*(10**-8)*((log(T_m))**9))) " , "Glover & Abel 2008"  , "0-8000"], 
							   [11,"H-  + H+ ->  H   + H "    ,"[H-]" , "[H+]" , ""      , "[H]"   , "[H]"  , ""     ,"k[9](y[8])", "lambda T_m: 1.40*(10**-7)*((T_m/300.0)**-0.487)*exp(T_m/29300.0)"                              , "Lepp,S. Stancil P.C., & Dalgarno 2002"                        , "0-8000"], 
							   [12,"H-  + H+ ->  H2+ + e-"    ,"[H-]" , "[H+]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[10](y[8])", "lambda T_m: 6.9*(10**-9)*(T_m**-0.35) if T_m < 8000 else 9.6*(10**-7)*(T_m**-0.9)"                                                         , "Glover & Abel 2008"                                         , "0-8000"], 
							   [13,"H2+ + e- ->  H   + H "    ,"[H2+]", "[e-]" , ""      , "[H]"   , "[H]"  , ""     ,"k[11](y[8])", "lambda T_m: 1.0*(10**-8) if T_m < 617 else 1.32*(10**-6)*(T_m**-0.76)"                        , "Glover & Abel 2008" 							               , "0-8000"], 
							   [14,"H2+ + H- ->  H   + H2"    ,"[H2+]", "[H-]" , ""      , "[H]"   , "[H2]" , ""     ,"k[12](y[8])", "lambda T_m: 1.4*(10**-7)*((T_m/300.0)**-0.5)"                                                  , "Glover & Abel 2008"                                             , "0-8000"],
							   [15,"H2  + e- ->  H   + H-"    ,"[H2]" , "[e-]" , ""      , "[H]"   , "[H-]" , ""     ,"k[13](y[8])", "lambda T_m: 36.7*(T_m**-2.28)*exp(-47172.0/T_m)"                                                 , "Capitelli, M., Coppola, C. M., Diomede, P., & Longo, S. 2007" , "0-8000"],
							   [16,"H-  + hv ->  H   + e-"    ,"[H-]" , "[hv]" , ""      , "[H]"   , "[e-]" , ""     ,"k[14](y[7])", "lambda T_r: 0.144*(T_r**2.13)*exp(-8650.0/T_r)"                                                  , "Abell et al 1997"                                             , "0-8000"], 
							   [17,"H2+ + hv ->  H   + H+"    ,"[H2+]", "[hv]" , ""      , "[H]"   , "[H+]" , ""     ,"k[15](y[7])", "lambda T_r: 2.0*10*(T_r**1.59)*exp(-82000.0/T_r)"                                              , "Stancil et al 1994"										   , "0-8000"], 
							   [18,"H2  + hv ->  H2+ + e-"    ,"[H2]" , "[hv]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[16](y[7])", "lambda T_r: 2.9*(10**2)*(T_r**1.56)*exp(-178500.0/T_r)"                                          , "Galli D & Palla 1998"                                         , "0-8000"], 
							   [19,"H2+ + hv ->  2H+ + e-"    ,"[H2+]", "[hv]" , ""      , "[H+]"  , "[H+]" , "[e-]" ,"k[17](y[7])", "lambda T_r: 90*(T_r**1.48)*(exp(-335000.0/T_r))"                                                 , "Galli D & Palla 1998"                                         , "0-8000"],
							   [20,"H2  + hv ->  H   + H "    ,"[H2]" , "[hv]" , ""      , "[H]"   , "[H]"  , ""     ,"k[18](y[7])", "lambda T_r: 1.13*(10**6)*(T_r**0.369)*exp(-140000.0/T_r)"                                        , "Glover, S.C.0 & Jappsen 2007"                                 , "0-8000"], 
							   [21,"D-  + hv ->  D   + e-"    ,"[D-]" , "[hv]" , ""      , "[D]"   , "[e-]" , ""     ,"k[19](y[7])", "lambda T_r: 1.1*(10**-1)*(T_r**2.13)*exp(-8823.0/T_r)"                                           , "Galli D & Palla 1998"                                            , "0-8000"],		
							   [22,"HD+ + hv ->  D   + H+"    ,"[HD+]", "[hv]" , ""      , "[D]"   , "[H+]" , ""     ,"k[20](y[7])", "lambda T_r: 0.5*1.63*(10**-7)*exp(-32400.0/T_r)"                                                  , "Galli D & Palla 1998"                                            , "0-8000"],
							   [23,"HD+ + hv ->  H   + D+"    ,"[HD+]", "[hv]" , ""      , "[H]"   , "[D+]" , ""     ,"k[21](y[7])", "lambda T_r: 0.5*1.63*(10**-7)*exp(-32400.0/T_r)"                                                  , "Galli D & Palla 1998"                                            , "0-8000"],
							   [24,"HD+ + hv ->  H+  + D+ + e","[HD+]", "[hv]" , ""      , "[H+]"  , "[D+]" , "[e-]" ,"k[22](y[7])", "lambda T_r: 90.0*(T_r**1.48)*exp(-335000.0/T_r)"                                                   , "Galli D & Palla 1998"                                           , "0-8000"],  	
							   [25,"HD  + hv ->  HD+ + e-"    ,"[HD]" , "[hv]" , ""      , "[HD+]" , "[e-]" , ""     ,"k[23](y[7])" , "lambda T_r: 2.9*100*(T_r**1.56)*exp(-178500.0/T_r)"                                              , "Galli D & Palla 1998"                                            , "0-8000"],
							   [26,"H+  +  e- -> H  +  hv"    ,"[D+]" , "[e-]" , ""      , "[D]"   , "[hv]" , ""     ,"k[24](Z)", "lambda Z: (8.76*(10**-11))*((1+Z)**-0.58)"                                                                    , "Galli & Palla 1998"         , "0-8000"],
							   [27,"D   + H+ ->  D+  + H"     ,"[D]"  , "[H+]" , ""      , "[D+]"  , "[H]"  , ""     ,"k[25](y[8])", "lambda T_m: 2.0*(10**-10)*(T_m**0.402)*exp(-37.1/T_m) - 3.31*(10**-17)*(T_m**1.48)"              , "Savin 2002"                                                   , "0-8000"],
							   [28,"D+  + H  ->  D   + H+"    ,"[D+]" , "[H]"  , ""      , "[D]"   , "[H+]" , ""     ,"k[26](y[8])", "lambda T_m: 2.06*(10**-10)*(T_m**0.396)*exp(-33.0/T_m) + 2.03*(10**-9)*(T_m**-0.332)"            , "Savin 2002"                                                   , "0-8000"],
							   [29,"D   + H  ->  HD  + hv"    ,"[D]"  , "[H]"  , ""      , "[HD]"  , "[hv]" , ""     ,"k[27](y[8])", "lambda T_m: (10**-25)*(2.80202-(6.63697*log(T_m))+(4.75619*(log(T_m)**2))-(1.39325*(log(T_m)**3))+(0.178259*(log(T_m)**4))-(0.00817097*(log(T_m)**5)))"     , "Glover & Abel 2008" , "0-8000"],
							   [30,"D   + H2 ->  HD  + H"     ,"[D]"  , "[H2]" , ""      , "[HD]"  , "[H]"  , ""     ,"k[28](y[8])", "lambda T_m: 9.0*(10**-11)*exp(-3876.0/T_m) if T_m < 200 else 1.69*(10**-10)*exp((-4680/T_m) +(198800.0/(T_m**2)))"                                                      , "Galli & Palla 1998"                                           , "0-8000"],
							   [31,"HD+ + H  ->  HD  + H+"    ,"[HD+]", "[H]"  , ""      , "[HD]"  , "[H+]" , ""     ,"k[29](y[8])", "lambda T_m: 6.4*(10**-10)"                                                                     , "SLD98"                                                        , "0-8000"],
							   [32,"D+  + H2 ->  HD  + H+"    , "[D+]" , "[H2]" , ""      , "[HD]"  , "[H+]" , ""     ,"k[30](y[8])", "lambda T_m: 1.0*(10**-9)*(0.417+0.846*log10(T_m)-0.137*(log10(T_m)**2))"                           , "GP02"                                                         , "0-8000"],
							   [33,"HD  + H  ->  D   + H2"    , "[HD]" , "[H]"  , ""      , "[D]"   , "[H2]" , ""     ,"k[31](y[8])" , "lambda T_m: 3.2*(10**-11)*exp(-3624.0/T_m) if T_m < 200 else 5.25*(10**-11)*exp((-4430.0/T_m)+(173900.0/(T_m**2)))"                                                      , "Galli & Palla 1998"                                           , "0-8000"],
							   [34,"HD  + H+ ->  D+  + H2"    ,"[HD]" , "[H+]" , ""      , "[D+]"  , "[H2]" , ""     ,"k[32](y[8])", "lambda T_m: 1.1*(10**-9)*exp(-488.0/T_m)"                                                        , "Galli & Palla 2002"                                           , "0-8000"],
							   [35,"D   + H+ ->  HD+ + hv"    , "[D]"  , "[H+]" , ""      , "[HD+]" , "[hv]" , ""     ,"k[33](y[8])", "lambda T_m: 10**(-19.38-1.523*log(T_m)+1.118*(log(T_m)**2)-0.1269*(log(T_m)**3))"              , "Galli & Palla 1998"                                           , "0-8000"],
							   [36,"D+  + H  ->  HD+ + hv"    ,"[D+]" , "[H]"  , ""      , "[HD+]" , "[hv]" , ""     ,"k[34](y[8])", "lambda T_m: 10**(-19.38-1.523*log(T_m)+1.118*(log(T_m)**2)-0.1269*(log(T_m)**3))"              , "Galli & Palla 1998"                                           , "0-8000"],
							   [37,"HD+ + e- ->  D   + H "    , "[HD+]", "[e-]" , ""      , "[D]"   , "[H]"  , ""     ,"k[35](y[8])", "lambda T_m: 7.2*(10**-8)*(T_m**-0.5)"                                                          , "SLD98"                                                        , "0-8000"],
							   [38,"D   + e- ->  D-  + hv"    , "[D]"  , "[e-]" , ""      , "[D-]"  , "[hv]" , ""     ,"k[36](y[8])", "lambda T_m: 3.0*(10**-16)*(T_m/300.0)**0.95*exp(-T_m/9320.0)"                                      , "SLD98"                                                        , "0-8000"],
							   [39,"D+  + D- ->  D   + D"     , "[D+]" , "[D-]" , ""      , "[D]"   , "[D]"  , ""     ,"k[37](y[8])", "lambda T_m: 1.96*(10**-7)*((T_m/300.0)**-0.487)*exp(-T_m/29300.0)"                                   , "LSD02"                                                        , "0-8000"],
							   [40,"H+  + D- ->  D   + H"     , "[H+]" , "[D-]" , ""      , "[D]"   , "[H]"  , ""     ,"k[38](y[8])", "lambda T_m: 1.61*(10**-7)*((T_m/300.0)**-0.487)*exp(-T_m/29300.0)"                                   , "LSD02"                                                        , "0-8000"],
							   [41,"H-  + D  ->  H   + D-"    , "[H-]" , "[D]"  , ""      , "[H]"   , "[D-]" , ""     ,"k[39](y[8])", "lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41"                                                      , "SLD98"                                                        , "0-8000"],
							   [42,"D-  + H  ->  D   + H-"    , "[D-]" , "[H]"  , ""      , "[D]"   , "[H-]" , ""     ,"k[40](y[8])", "lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41"                                                      , "SLD98"                                                        , "0-8000"],
							   [43,"D-  + H  ->  HD  + e-"    , "[D-]" , "[H]"  , ""      , "[HD]"  , "[e-]" , ""     ,"k[41](y[8])", "lambda T_m: 1.5*(10**-9)*(T_m/300.0)**-0.1"                                                      , "SLD98"                                                        , "0-8000"],
							   [44,"D   + H- ->  HD  + e-"    , "[D]"  , "[H-]" , ""      , "[HD]"  , "[e-]" , ""     ,"k[42](y[8])", "lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41"                                                      , "SLD98"                                                        , "0-8000"],
							   [45,"H-  + D+ ->  D   + H"     , "[H-]" , "[D+]" , ""      , "[D]"   , "[H]"  , ""     ,"k[43](y[8])", "lambda T_m: 1.61*(10**-7)*((T_m/300.0)**-0.487)*exp(T_m/29300.0)"                                  , "LSD02"                                                        , "0-8000"],        
							   [46,"He+ + D  ->  D+  + He"    , "[He+]", "[D]" , ""      , "[D+]"  , "[He]"  , ""    ,"k[44](y[8])" , "lambda T_m: 1.1*(10**-15)*((T_m/300.0)**0.25)"              , "Glover & Abel 2008"                                                 , "0-8000"],
							   [47,"He  + D+ ->  D   + He+"   , "[He]" , "[D+]" , ""      , "[D]"   , "[He+]" , ""     ,"k[45](y[8])", "lambda T_m: 1.85*(10**-9)*(T_m**-0.75)*exp(-127500.0/T_m)" , "Glover & Abel 2008"                                                 , "0-8000"],
							   [48,"HD  + He+->  HD+ + He"    , "[HD]" , "[He+]" , ""     , "[HD+]" , "[He]"  , ""     ,"k[46](y[8])", "lambda T_m: 7.2*(10**-15)"                               , "Glover & Abel 2008"                                                        , "0-8000"],
							   [49,"D2  + He+->  D2+ + He"    , "[D2]" , "[He+]" , ""     , "[D2+]" , "[He]"  , ""     ,"k[47](y[8])", "lambda T_m: 2.5*(10**-14)"                               , "Glover & Abel 2008"                                                        , "0-8000"],
							   [50,"D2  + He+->  He + D+ + D" , "[D2]" , "[He+]" , ""     , "[D+]" , "[He]"  , "[D]"   ,"k[48](y[8])", "lambda T_m: 1.1*(10**-13)*((T_m*(10**-3))**-0.24)"       , "Glover & Abel 2008"                                                        , "0-8000"],
							   [51,"H+  +  e- -> H  +  hv"    , "[H+]" , "[e-]" , ""      , "[H]"   , "[hv]" , ""     ,"k[49](Z)"      , "lambda Z: (8.76*(10**-11))*((1.0+Z)**-0.58)"  , "Galli & Palla 1998"  , "0-8000"], 
							   [52,"H+  +  e- -> H  +  hv"    , "[H]"  , "[hv]" , ""      , "[H+]"  , "[e-]" , ""     ,"k[50](y[7],Z)", "lambda T_r,Z: ((2.41*(10**+15)*(T_r**1.5)*exp(-179472/T_r))*((8.76*(10**-11))*((1.0+Z)**-0.58)))"                         , "Galli & Palla 1998"                                          "0-8000"],
							   [53,"H   + e- ->  H+  + 2e-"   , "[H]"  , "[e-]" , ""      , "[H+]"  , "[e-]" , "[e-]" ,"k[51](y[8])", "lambda T_m: 5.85*(10**-11)*(T_m**0.5)*exp(-157809.1/T_m)"                                       ,"Janev et al 1987"                                             , "0-8000"],
							   [54,"H2  + H2 ->  H2  + H + H" , "[H2]" , "[H2]" , ""      , "[H2]"  , "[H]"  , "[H]"  ,"k[52](y[8])", "lambda T_m: exp(-54657.4/T_m)*(5.996*(10**-30)*(T_m**4.1881))/((1.0+(6.761*((10**-6)*T_m)))**5.6881)","Lepp & Shull 1983"                                        , "0-8000"],  
							   [55,"He++ + e- -> He+ + hv"    , "[He++]","[e-]" , ""      , "[He+]" , "[hv]" , ""     ,"k[53](y[8])", "lambda T_m: 3.36*(10**-10)*(T_m**-0.5)*((T_m/1000.0)**-0.2)*(1.0+(T_m/(10**6))**0.7)**-1"             , "Cen 1992"                                                     , "0-8000"],
							   [56,"He+  + hv -> He++ + e-"   , "[He+]", "[hv]" , ""      , "[He++]", "[e-]" , ""     ,"k[54](y[7])", "lambda T_r: 5.0*(10)*(T_r**1.63)*exp(-590000.0/T_r)"                                             , "Galli & Palla 1998"                                             , "0-8000"],
							   [57,"He+  + e- -> He  + hv"    , "[He+]", "[e-]" , ""      , "[He]"  , "[hv]" , ""     ,"k[55](y[8]*kbev)", "lambda T_m: 3.925*(10**-13)*(T_m**-0.6353)"		                                               , "Abell et al 1997"                                              , "0-8000"],
							   [58,"He   + hv -> He+ + e-"    , "[He]" , "[hv]" , ""      , "[He+]" , "[e-]" , ""     ,"k[56](y[7])", "lambda T_r: 1.0*(10**4)*(T_r**1.23)*exp(-280000.0/T_r)"                                           , "Galli & Palla 1998"                                             , "0-8000"],
							   [59,"He   + H+ -> He+ + H"     ,"[He]" , "[H+]" , ""      , "[He+]" , "[H]"  , ""     ,"k[57](y[8])", "lambda T_m: 1.26*(10**-9)*(T_m**-0.75)*exp(-127500.0/T_m)"                                        , "Glover & Jappsen 2007"                                        , "0-8000"],
							   [60,"He+  + H  -> He  + H+"    , "[He+]", "[H]"  , ""      , "[He]"  , "[H+]" , ""     ,"k[58](y[8])", "lambda T_m: 1.25*(10**-15)*(T_m/300.0)**0.25"                                                    , "Zygelman et al 1989"                                          , "0-8000"],
							   [61,"He   + H+ -> HeH+ + hv"   , "[He]" , "[H+]" , ""      , "[HeH+]", "[hv]" , ""     ,"k[59](y[8])", "lambda T_m: 7.6*(10**-18)*(T_m**-0.5)"                                    , "Galli & Palla 1998"                                           , "0-8000"],
							   [63,"He  + H2+ -> HeH+ + H"    , "[He]" , "[H2+]", ""      , "[HeH+]", "[H]"  , ""     ,"k[60](y[8])", "lambda T_m: 3.0*(10**-10)*exp(-6717.0/T_m)"                                                       , "Galli & Palla 1998"                                           , "0-8000"],
							   [64,"He+  + H  -> HeH+ + hv"   , "[He+]", "[H]"  , ""      , "[HeH+]", "[hv]" , ""     ,"k[61](y[8])", "lambda T_m: 1.6*(10**-14)*(T_m**-0.33) if T_m < 4000 else 1.0*(10**-15)"                                        , "Stancil et al 1998"                                           , "0-8000"],
							   [65,"HeH+ + H  -> He   + H2+"  , "[HeH+]","[H]"  , ""      , "[He]"  , "[H2+]", ""     ,"k[62](y[8])", "lambda T_m: 0.69*(10**-9)*((T_m/300.0)**0.13)*exp(-T_m/33100.0)"                                    , "Linder et al  1995"                                           , "0-8000"],
							   [66,"HeH+ + e- -> He   + H"    , "[HeH+]","[e-]" , ""      , "[He]"  , "[H]"  , ""     ,"k[63](y[8])", "lambda T_m: 3.0*(10**-8)*((T_m/300.0)**-0.47)"                                                    , "Stancil et al 1998"                                           , "0-8000"],
							   [67,"HeH+ + hv -> He   + H+"   , "[HeH+]","[hv]" , ""      , "[He]"  , "[H+]" , ""     ,"k[64](y[7])", "lambda T_r: 220*(T_r**0.9)*exp(-22740.0/T_r)"                                                     , "Jurek et al 1995"                                             , "0-8000"],
							   [68,"HeH+ + hv -> He+  + H"    , "[HeH+]","[hv]" , ""      , "[He+]" , "[H]"  , ""     ,"k[65](y[7])", "lambda T_r: 7.8*(10**3)*(T_r**1.2)*exp(-240000.0/T_r)"                                            , "Galli & Palla 1998"                                           , "0-8000"],
							   [69,"He   + e- -> He+  + 2e-"  , "[He]"  ,"[e-]" , ""      , "[He+]" , "[e-]" , "[e-]" ,"k[66](y[8]*kbev)", "lambda T_m: exp(-44.09864886 + 23.91596563*log(T_m) - 10.7532302*(log(T_m)**2) + 3.05803875*(log(T_m)**3) - 0.56851189*(log(T_m)**4) + 6.79539123*(10**-2)*(log(T_m)**5) - 5.00905610*(10**-3)*(log(T_m)**6)+2.06723616*(10**-4)*(log(T_m)**7) - 3.64916141*(10**-6)*(log(T_m)**8))", "Abel et al 1997", "0-8000"],
							   [70,"He+  + e- -> He++ + 2e-"  , "[He+]" ,"[e-]" , ""      , "[He++]", "[e-]" , "[e-]" ,"k[67](y[8]*kbev)", "lambda T_m: exp(-68.71040990 + 43.93347633*log(T_m) - 18.4806699*(log(T_m)**2) + 4.70162649*(log(T_m)**3) - 0.76924663*(log(T_m)**4) + 8.113042*0.01*(log(T_m)**5) - 5.32402063*(10**-3)*(log(T_m)**6)+1.97570531*(10**-4)*(log(T_m)**7) - 3.16558106*(10**-6)*(log(T_m)**8))", "Abel et al 1997", "0-8000"],   
							   [71,"Li+  + e- -> Li   + hv"   , "[Li+]" ,"[e-]" , ""      , "[Li]"  , "hv"   , ""     ,"k[68](y[8])", "lambda T_m: 1.036*(10**-11)*((((T_m/107.7)**0.5)*((1+((T_m/107.7)**0.5))**0.612)*(1+(T_m/(1.177*10**7))**0.5)**1.388)**-1.0)", "Verner & Ferlan 1996", "0-8000"],
							   [72,"Li   + hv -> Li+  + e-"   , "[Li]"  ,"[hv]" , ""      , "[Li+]" , "[e-]" , ""     ,"k[69](y[7])", "lambda T_r: 1.3*(10**3)*(T_r**1.45)*exp(-60500.0/T_r)"                                                                       , "Galli & Palla 1998"  , "0-8000"],
							   [73,"Li-  + H+ -> Li   + H"    , "[Li-]" ,"[H+]" , ""      , "[Li]"  , "[H]"  , ""     ,"k[70](y[8])", "lambda T_m: (6.3*(10**-6)*(T_m**-0.5)) - (7.6*(10**-9)) + (2.6*(10**-10)*(T_m**0.5)) + (2.7*(10**-14)*T_m)"                  , "Peart & Hayton 1994" , "0-8000"],
							   [74,"Li   + e- -> Li-  + hv"   , "[Li]"  ,"[e-]" , ""      , "[Li-]" , "[hv]" , ""     ,"k[71](y[8])", "lambda T_m: 6.1*(10**-17)*(T_m**0.58)*exp(-T_m/17200.0)"                                                                     , "Ramsbottom et al 1994", "0-8000"],
							   [75,"Li-  + hv -> Li   + e-"   , "[Li-]" ,"[hv]" , ""      , "[Li]"  , "[e-]" , ""     ,"k[72](y[7])" , "lambda T_r: 1.8*(10**2)*(T_r**1.4)*exp(-8100.0/T_r)"                                                                         , "Ramsbottom et al 1994", "0-8000"],
							   [76,"Li   + H  -> LiH  + hv"   , "[Li]"  ,"[H]"  , ""      , "[LiH]" , "[hv]" , ""     ,"k[73](y[8])", "lambda T_m: ((5.6*(10**19)*(T_m**-0.15))+(7.2*(10**15)*(T_m**1.21)))**-1"                                                    , "Dalgarno et al 1996"  , "0-8000"],
							   [77,"LiH  + hv -> Li   + H"    , "[LiH]" ,"[hv]" , ""      , "[Li]"  , "[H]"  , ""     ,"k[74](y[7])", "lambda T_r: 8.3*(10**4)*(T_r**0.3)*exp(-29000.0/T_r)"                                                                        , "Galli & Palla 1998"   , "0-8000"],
							   [78,"Li-  + H  -> LiH  + e-"   , "[Li-]" ,"[H]"  , ""      , "[LiH]" , "[e-]" , ""     ,"k[75](y[8])", "lambda T_m: 4.0*(10**-10)"                                                                                                   , "Stancil et al 1996"   , "0-8000"],
							   [79,"LiH  + H  -> Li   + H2"   , "[LiH]" ,"[H]"  , ""      , "[Li]"  , "[H2]" , ""     ,"k[76](y[8])" , "lambda T_m: 2.0*(10**-11)"                                                                                                   , "Stancil et al 1996"   , "0-8000"],
							   [80,"Li+  + H  -> LiH+ + hv"   , "[Li+]" ,"[H]"  , ""      , "[LiH+]", "[hv]" , ""     ,"k[77](y[8])", "lambda T_m: 10**(-22.4+(0.999*log10(T_m))-(0.351*((log10(T_m))**2)))"                                                        , "Dalgarno et al 1996"  , "0-8000"],
							   [81,"Li   + H+ -> LiH+ + hv"   , "[Li]"  ,"[H+]" , ""      , "[LiH+]", "[hv]" , ""     ,"k[78](y[8])", "lambda T_m: 4.8*(10**-14)*(T_m**-0.49)"                                                                                      , "Dalgarno et al 1996"  , "0-8000"],
							   [82,"LiH+ + e- -> Li   + H"    ,"[LiH+]","[e-]" , ""      , "[Li]"  , "[H]"  , ""     ,"k[79](y[8])" , "lambda T_m: 3.8*(10**-7)*(T_m**-0.47)"                                                                                       , "Stancil et al 1996"   , "0-8000"],
							   [83,"LiH+ + H  -> Li+  + H2"   , "[LiH+]","[H]"  , ""      , "[Li+]" , "[H2]" , ""     ,"k[80](y[8])", "lambda T_m: 3.0*(10**-10)"                                                                                                   , "Stancil et al 1996"   , "0-8000"],
							   [84,"LiH+ + hv -> Li+  + H"    , "[LiH+]","[hv]" , ""      , "[Li+]" , "[H]"  , ""     ,"k[81](y[7])", "lambda T_r: 7.0*(10**2)*exp(-1900.0/T_r) if T_r < 1000 else (3.45*(10**-16)*(T_r**-1.06))"                                   , "Dalgarno et al 1996"  , "0-8000"],	   
							   [85,"H-  + H2+ -> 2H   + H"    ,"[H-]"  ,"[H2+]", ""      , "[H]"   , "[H]"  , "[H]"  ,"k[82](y[8])", "lambda T_m: 1.4*(10**-7)*((T_m/300.0)**-0.5)"                                                                                , "Glover & Abell 2008"  , "0-8000"],
							   [86,"H2  + He+ -> He + H + H+" ,"[H2]"  ,"[He+]", ""      , "[He]"  , "[H]"  , "[H+]" ,"k[83](y[8])", "lambda T_m: 3.7*(10**-14)*exp(35.0/T_m)"                                                                                     , "Glover & Abel 2008"   , "0-8000"],	
							   [87,"H2  + He+ -> H2+  + He"   ,"[H2]"  ,"[He+]", ""      , "[H2+]" , "[He]" , ""     ,"k[84](y[8])", "lambda T_m: 7.2*(10**-15)"                                                                                                   , "Glover & Abel 2008"   , "0-8000"],  
							   [88,"He+ + H-  -> He   + H  "  ,"[He+]" ,"[H-]" , ""      , "[He]" ,  "[H]"  , ""     ,"k[85](y[8])", "lambda T_m: 2.32*(10**-7)*((T_m/300.0)**-0.52)*exp(T_m/22400.0)"                                                             , "Glover & Abel 2008"   , "0-8000"],  
							   [89,"He  + H-  -> He + H + e-" ,"[He]"  ,"[H-]" , ""      , "[He]" ,  "[H]"  , "[e-]" ,"k[86](y[8])", "lambda T_m: 4.1*(10**-17)*(T_m**2)*exp(-19870.0/T_m)"                                                                        , "Glover & Abel 2008"   , "0-8000"],	
							   [90,"H2D+ + e- -> H + H + D"   ,"[H2D+]","[e-]" , ""      , "[H]"  ,  "[H]"  , "[D]"  ,"k[87](y[8])", "lambda T_m: 0.73*1.0*(10**-6)*(T_m**-0.5)"                                  , "Galli & Palla 1998"                                                        , "0-8000"],	
							   [91,"H2D+ + e- -> H2  + D"       , "[H2D+]", "[e-]" , ""      , "[H2]" ,  "[D]"  , ""   ,"k[88](y[8])", "lambda T_m: 0.07*1.0*(10**-6)*(T_m**-0.5)"                                  , "Galli & Palla 1998"                                                        , "0-8000"],		
							   [92,"H2D+ + e- -> HD  + H"       , "[H2D+]", "[e-]" , ""      , "[HD]" ,  "[H]"  , ""   ,"k[89](y[8])", "lambda T_m: 0.20*1.0*(10**-6)*(T_m**-0.5)"                                  , "Galli & Palla 1998"                                                        , "0-8000"],		
							   [93,"HD+  + H2 -> H2D+ + H"      , "[HD+]" , "[H2]" , ""      , "[H2D+]", "[H]"  , ""   ,"k[90](y[8])", "lambda T_m: 2.0*(10**-9)"                                                   , "Galli & Palla 1998"                                                        , "0-8000"],		
							   [94,"H+  +  e- -> H  +  hv"      , "[D]"   , "[hv]" , ""      , "[D+]" ,  "[e-]" , ""   ,"k[91](y[7],Z)", "lambda T_r, Z: 2.41*(10**15)*(T_r**1.5)*exp(-179472.0/T_r)*(8.76*(10**-11))*((1+Z)**-0.58)"          ,  "Galli & Palla 1998"        , "0-8000"],
							   [100, "H   + hv -> H+   + e-"    , "[H]"   ,"" , ""      , "[H+]" ,  "[e-]" , ""     ,"flux[0]"],
							   [101, "He+ + hv -> He++ + e-"    , "[He+]" ,"" , ""      , "[He++]", "[e-]" , ""     ,"flux[1]"],
							   [102, "He +  hv -> He+  + e-"    , "[He]"  ,"" , ""      , "[He+]",  "[e-]" , ""     ,"flux[2]"],
							   [103, "H- +  hv -> H    + e-"    , "[H-]"  ,"" , ""      , "[H]"  ,  "[e-]" , ""     ,"flux[3]"],
							   [104, "D +  hv -> D+  + e-"      ,  "[D]"  ,"" , ""      , "[D+]",  "[e-]" , ""     ,"flux[4]"],
							   [105, "H2 +  hv -> H2+  + e-"    , "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[5]"],
							   [106, "H2 +  hv -> H2+  + e-"    , "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[6]"],
							   [107, "H2 +  hv -> H2+  + e-"    , "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[7]"],
							   [108, "H2+ + hv -> H    + H+"    , "[H2+]" ,"" , ""      , "[H]"  ,  "[H+]" , ""     ,"flux[8]"],
							   [109, "H2+ + hv -> 2H+  + e-"    , "[H2+]" ,"" , ""      , "[H+]" ,  "[H+]" , "[e-]" ,"flux[9]"],
							   [110, "H2  + hv -> H2* -> H + H" , "[H2]"  ,"" , ""      , "[H]"  ,  "[H]"  , ""     ,"flux[10]"],
							   [111, "H2  + hv -> H + H"        , "[H2]"  ,"" , ""      , "[H]"  ,  "[H]"  , ""     ,"flux[11]"]])
							
							
	elif (chemical_model ==  2):			   
		
		reaction_matrix=array([
							   [ 1,"H   + e- ->  H-  + hv"   , "[H]"  , "[e-]" , ""      , "[H-]"  , "[hv]" , ""     , "k[0](y[8])", "lambda T_m: 1.4*(10**-18)*(T_m**0.928)*exp(-T_m/16200.0)"                                          , "Galli D & Palla 1998"                                         , "0-8000"], #chemical equations we are interested in.											
							   [ 2,"H-  + H  ->  H2  + e-"    ,"[H-]" , "[H]"  , ""      , "[H2]"  , "[e-]" , ""     , "k[1](y[8])", "lambda T_m: 4.0*(10**-9)*(T_m**-0.17) if T_m > 300 else (1.5*10**-9)" 					         , "Galli D & Palla 1998"                                         , "0-8000"],
							   [ 3,"H+  + H  ->  H2+ + hv"    ,"[H+]" , "[H]"  , ""      , "[H2+]" , "[hv]" , ""     , "k[2](y[8])", "lambda T_m: 10**(-19.38-(1.523*log10(T_m))+(1.118*((log10(T_m))**2))-(0.1269*((log10(T_m))**3)))" , "Glover & Abel 2008"                                           , "0-8000"], 
							   [ 4,"H2+ + H  ->  H2  + H+"    ,"[H2+]", "[H]"  , ""      , "[H2]"  , "[H+]" , ""     , "k[3](y[8])", "lambda T_m: 6.4*(10**-10)"					                                                       , "Galli D & Palla 1998"                                         , "0-8000"], 
							   [ 6,"H2  + H+ ->  H2+ + H "    ,"[H2]" , "[H+]" , ""      , "[H2+]" , "[H]"  , ""     , "k[4](y[8])", "lambda T_m: 3.0*(10**-10)*exp(-21050.0/T_m)"				                                       , "Galli D & Palla 1998"                                         , "0-8000"], 
							   [ 7,"H2  + e- ->  2H  + e-"    ,"[H2]" , "[e-]" , ""      , "[H]"   , "[H]"  , "[e-]" , "k[5](y[8])", "lambda T_m: 1.91*(10**-9)*(T_m**0.136)*exp(-53407.1/T_m)"                                          , "Trevian, C.S & Tennyson J.2002"                               , "0-8000"], 
							   [ 9,"H-  + e- ->  H   + 2e-"   ,"[H-]" , "[e-]" , ""      , "[H]"   , "[e-]" , "[e-]" , "k[6](y[8]*kbev)", "lambda T_m: exp((-18.01849334)+(2.3608522*log(T_m))-(0.28274430*((log(T_m))**2))+(1.62331664*(10**-2)*(log(T_m))**3)-(3.36501203*(10**-2)*((log(T_m))**4))+(1.17832978*(10**-2)*((log(T_m))**5))-(1.65619470*(10**-3)*((log(T_m))**6))+(1.06827520*(10**-4)*((log(T_m))**7))-(2.63128581*(10**-6)*((log(T_m))**8)))", "Glover & Abel 2008", "0-8000"],
							   [10,"H-  + H  ->  2H  + e-"    ,"[H-]" , "[H]"  , ""      , "[H]"   , "[H]"  , "[e-]" ,"k[7](y[8]*kbev)","lambda T_m: 2.5634*(10**-9)*(T_m**1.78186) if T_m < 0.1 else exp((-20.37260896) + (1.13944933*log(T_m)) - (0.14210135*((log(T_m))**2)) + (8.4644554*(10**-3)*(log(T_m))**3) - (1.4327641*(10**-3)*((log(T_m))**4)) + (2.0122503*(10**-4)*((log(T_m))**5)) + (8.6639632*(10**-5)*((log(T_m))**6)) - (2.5850097*(10**-5)*((log(T_m))**7)) + (2.4555012*(10**-6)*((log(T_m))**8)) - (8.0683825*(10**-8)*((log(T_m))**9))) " , "Glover & Abel 2008"  , "0-8000"], 
							   [11,"H-  + H+ ->  H   + H "    ,"[H-]" , "[H+]" , ""      , "[H]"   , "[H]"  , ""     ,"k[8](y[8])", "lambda T_m: 1.40*(10**-7)*((T_m/300.0)**-0.487)*exp(T_m/29300.0)"                              , "Lepp,S. Stancil P.C., & Dalgarno 2002"                        , "0-8000"], 
							   [12,"H-  + H+ ->  H2+ + e-"    ,"[H-]" , "[H+]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[9](y[8])", "lambda T_m: 6.9*(10**-9)*(T_m**-0.35) if T_m < 8000 else 9.6*(10**-7)*(T_m**-0.9)"                                                         , "Glover & Abel 2008"                                         , "0-8000"], 
							   [13,"H2+ + e- ->  H   + H "    ,"[H2+]", "[e-]" , ""      , "[H]"   , "[H]"  , ""     ,"k[10](y[8])", "lambda T_m: 1.0*(10**-8) if T_m < 617 else 1.32*(10**-6)*(T_m**-0.76)"                        , "Glover & Abel 2008" 							               , "0-8000"], 
							   [14,"H2+ + H- ->  H   + H2"    ,"[H2+]", "[H-]" , ""      , "[H]"   , "[H2]" , ""     ,"k[11](y[8])", "lambda T_m: 1.4*(10**-7)*((T_m/300.0)**-0.5)"                                                  , "Glover & Abel 2008"                                             , "0-8000"],
							   [15,"H2  + e- ->  H   + H-"    ,"[H2]" , "[e-]" , ""      , "[H]"   , "[H-]" , ""     ,"k[12](y[8])", "lambda T_m: 36.7*(T_m**-2.28)*exp(-47172.0/T_m)"                                                 , "Capitelli, M., Coppola, C. M., Diomede, P., & Longo, S. 2007" , "0-8000"],
							   [16,"H-  + hv ->  H   + e-"    ,"[H-]" , "[hv]" , ""      , "[H]"   , "[e-]" , ""     ,"k[13](y[7])", "lambda T_r: 0.144*(T_r**2.13)*exp(-8650.0/T_r)"                                                  , "Abell et al 1997"                                             , "0-8000"], 
							   [17,"H2+ + hv ->  H   + H+"    ,"[H2+]", "[hv]" , ""      , "[H]"   , "[H+]" , ""     ,"k[14](y[7])", "lambda T_r: 2.0*10*(T_r**1.59)*exp(-82000.0/T_r)"                                              , "Stancil et al 1994"										   , "0-8000"], 
							   [18,"H2  + hv ->  H2+ + e-"    ,"[H2]" , "[hv]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[15](y[7])", "lambda T_r: 2.9*(10**2)*(T_r**1.56)*exp(-178500.0/T_r)"                                          , "Galli D & Palla 1998"                                         , "0-8000"], 
							   [20,"H2  + hv ->  H   + H "    ,"[H2]" , "[hv]" , ""      , "[H]"   , "[H]"  , ""     ,"k[16](y[7])", "lambda T_r: 1.13*(10**6)*(T_r**0.369)*exp(-140000.0/T_r)"                                        , "Glover, S.C.0 & Jappsen 2007"                                 , "0-8000"], 
							   [21,"D-  + hv ->  D   + e-"    ,"[D-]" , "[hv]" , ""      , "[D]"   , "[e-]" , ""     ,"k[17](y[7])", "lambda T_r: 1.1*(10**-1)*(T_r**2.13)*exp(-8823.0/T_r)"                                           , "Galli D & Palla 1998"                                            , "0-8000"],		
							   [25,"HD  + hv ->  HD+ + e-"    ,"[HD]" , "[hv]" , ""      , "[HD+]" , "[e-]" , ""     ,"k[18](y[7])" , "lambda T_r: 2.9*100*(T_r**1.56)*exp(-178500.0/T_r)"                                              , "Galli D & Palla 1998"                                            , "0-8000"],
							   [26,"H+  +  e- -> H  +  hv"    ,"[D+]" , "[e-]" , ""      , "[D]"   , "[hv]" , ""     ,"k[19](Z)", "lambda Z: (8.76*(10**-11))*((1+Z)**-0.58)"                                                                    , "Galli & Palla 1998"         , "0-8000"],
							   [27,"D   + H+ ->  D+  + H"     ,"[D]"  , "[H+]" , ""      , "[D+]"  , "[H]"  , ""     ,"k[20](y[8])", "lambda T_m: 2.0*(10**-10)*(T_m**0.402)*exp(-37.1/T_m) - 3.31*(10**-17)*(T_m**1.48)"              , "Savin 2002"                                                   , "0-8000"],
							   [28,"D+  + H  ->  D   + H+"    ,"[D+]" , "[H]"  , ""      , "[D]"   , "[H+]" , ""     ,"k[21](y[8])", "lambda T_m: 2.06*(10**-10)*(T_m**0.396)*exp(-33.0/T_m) + 2.03*(10**-9)*(T_m**-0.332)"            , "Savin 2002"                                                   , "0-8000"],
							   [29,"D   + H  ->  HD  + hv"    ,"[D]"  , "[H]"  , ""      , "[HD]"  , "[hv]" , ""     ,"k[22](y[8])", "lambda T_m: (10**-25)*(2.80202-(6.63697*log(T_m))+(4.75619*(log(T_m)**2))-(1.39325*(log(T_m)**3))+(0.178259*(log(T_m)**4))-(0.00817097*(log(T_m)**5)))"     , "Glover & Abel 2008" , "0-8000"],
							   [30,"D   + H2 ->  HD  + H"     ,"[D]"  , "[H2]" , ""      , "[HD]"  , "[H]"  , ""     ,"k[23](y[8])", "lambda T_m: 9.0*(10**-11)*exp(-3876.0/T_m) if T_m < 200 else 1.69*(10**-10)*exp((-4680/T_m) +(198800.0/(T_m**2)))"                                                      , "Galli & Palla 1998"                                           , "0-8000"],
							   [31,"HD+ + H  ->  HD  + H+"    ,"[HD+]", "[H]"  , ""      , "[HD]"  , "[H+]" , ""     ,"k[24](y[8])", "lambda T_m: 6.4*(10**-10)"                                                                     , "SLD98"                                                        , "0-8000"],
							   [32,"D+  + H2 ->  HD  + H+"    , "[D+]" , "[H2]" , ""      , "[HD]"  , "[H+]" , ""     ,"k[25](y[8])", "lambda T_m: 1.0*(10**-9)*(0.417+0.846*log10(T_m)-0.137*(log10(T_m)**2))"                           , "GP02"                                                         , "0-8000"],
							   [33,"HD  + H  ->  D   + H2"    , "[HD]" , "[H]"  , ""      , "[D]"   , "[H2]" , ""     ,"k[26](y[8])" , "lambda T_m: 3.2*(10**-11)*exp(-3624.0/T_m) if T_m < 200 else 5.25*(10**-11)*exp((-4430.0/T_m)+(173900.0/(T_m**2)))"                                                      , "Galli & Palla 1998"                                           , "0-8000"],
							   [34,"HD  + H+ ->  D+  + H2"    ,"[HD]" , "[H+]" , ""      , "[D+]"  , "[H2]" , ""     ,"k[27](y[8])", "lambda T_m: 1.1*(10**-9)*exp(-488.0/T_m)"                                                        , "Galli & Palla 2002"                                           , "0-8000"],
							   [35,"D   + H+ ->  HD+ + hv"    , "[D]"  , "[H+]" , ""      , "[HD+]" , "[hv]" , ""     ,"k[28](y[8])", "lambda T_m: 10**(-19.38-1.523*log(T_m)+1.118*(log(T_m)**2)-0.1269*(log(T_m)**3))"              , "Galli & Palla 1998"                                           , "0-8000"],
							   [36,"D+  + H  ->  HD+ + hv"    ,"[D+]" , "[H]"  , ""      , "[HD+]" , "[hv]" , ""     ,"k[29](y[8])", "lambda T_m: 10**(-19.38-1.523*log(T_m)+1.118*(log(T_m)**2)-0.1269*(log(T_m)**3))"              , "Galli & Palla 1998"                                           , "0-8000"],
							   [38,"D   + e- ->  D-  + hv"    , "[D]"  , "[e-]" , ""      , "[D-]"  , "[hv]" , ""     ,"k[30](y[8])", "lambda T_m: 3.0*(10**-16)*(T_m/300.0)**0.95*exp(-T_m/9320.0)"                                      , "SLD98"                                                        , "0-8000"],
							   [41,"H-  + D  ->  H   + D-"    , "[H-]" , "[D]"  , ""      , "[H]"   , "[D-]" , ""     ,"k[31](y[8])", "lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41"                                                      , "SLD98"                                                        , "0-8000"],
							   [42,"D-  + H  ->  D   + H-"    , "[D-]" , "[H]"  , ""      , "[D]"   , "[H-]" , ""     ,"k[32](y[8])", "lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41"                                                      , "SLD98"                                                        , "0-8000"],
							   [43,"D-  + H  ->  HD  + e-"    , "[D-]" , "[H]"  , ""      , "[HD]"  , "[e-]" , ""     ,"k[33](y[8])", "lambda T_m: 1.5*(10**-9)*(T_m/300.0)**-0.1"                                                      , "SLD98"                                                        , "0-8000"],
							   [44,"D   + H- ->  HD  + e-"    , "[D]"  , "[H-]" , ""      , "[HD]"  , "[e-]" , ""     ,"k[34](y[8])", "lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41"                                                      , "SLD98"                                                        , "0-8000"],
							   [51,"H+  +  e- -> H  +  hv"    , "[H+]" , "[e-]" , ""      , "[H]"   , "[hv]" , ""     ,"k[35](Z)"      , "lambda Z: (8.76*(10**-11))*((1.0+Z)**-0.58)"  , "Galli & Palla 1998"  , "0-8000"], 
							   [52,"H+  +  e- -> H  +  hv"    , "[H]"  , "[hv]" , ""      , "[H+]"  , "[e-]" , ""     ,"k[36](y[7],Z)", "lambda T_r,Z: ((2.41*(10**+15)*(T_r**1.5)*exp(-179472/T_r))*((8.76*(10**-11))*((1.0+Z)**-0.58)))"                         , "Galli & Palla 1998"                                          "0-8000"],
							   [53,"H   + e- ->  H+  + 2e-"   , "[H]"  , "[e-]" , ""      , "[H+]"  , "[e-]" , "[e-]" ,"k[37](y[8])", "lambda T_m: 5.85*(10**-11)*(T_m**0.5)*exp(-157809.1/T_m)"                                       ,"Janev et al 1987"                                             , "0-8000"],
							   [94,"H+  +  e- -> H  +  hv"      , "[D]"   , "[hv]" , ""  , "[D+]" ,  "[e-]" , ""     , "k[38](y[7],Z)", "lambda T_r, Z: 2.41*(10**15)*(T_r**1.5)*exp(-179472.0/T_r)*(8.76*(10**-11))*((1+Z)**-0.58)"          ,  "Galli & Palla 1998"        , "0-8000"],
							   [100, "H   + hv -> H+   + e-"    , "[H]"   ,"" , ""      , "[H+]" ,  "[e-]" , ""     ,"flux[0]"],
							   [101, "He+ + hv -> He++ + e-"    , "[He+]" ,"" , ""      , "[He++]", "[e-]" , ""     ,"flux[1]"],
							   [102, "He +  hv -> He+  + e-"    , "[He]"  ,"" , ""      , "[He+]",  "[e-]" , ""     ,"flux[2]"],
							   [103, "H- +  hv -> H    + e-"    , "[H-]"  ,"" , ""      , "[H]"  ,  "[e-]" , ""     ,"flux[3]"],
							   [104, "D +  hv -> D+  + e-"      ,  "[D]"  ,"" , ""      , "[D+]",  "[e-]" , ""      ,"flux[4]"],
							   [105, "H2 +  hv -> H2+  + e-"    , "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[5]"],
							   [106, "H2 +  hv -> H2+  + e-"    , "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[6]"],
							   [107, "H2 +  hv -> H2+  + e-"    , "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[7]"],
							   [108, "H2+ + hv -> H    + H+"    , "[H2+]" ,"" , ""      , "[H]"  ,  "[H+]" , ""     ,"flux[8]"],
							   [109, "H2+ + hv -> 2H+  + e-"    , "[H2+]" ,"" , ""      , "[H+]" ,  "[H+]" , "[e-]" ,"flux[9]"],
							   [110, "H2  + hv -> H2* -> H + H" , "[H2]"  ,"" , ""      , "[H]"  ,  "[H]"  , ""     ,"flux[10]"],
							   [111, "H2  + hv -> H + H"        , "[H2]"  ,"" , ""      , "[H]"  ,  "[H]"  , ""     ,"flux[11]"]])
							   
	elif (chemical_model ==  3): 
		reaction_matrix=array([
							   [53,"H   + e- ->  H+  + 2e-"   , "[H]"  , "[e-]" , ""      , "[H+]"  , "[e-]" , "[e-]" ,"k[0](y[8]*kbev)", "lambda T_m : exp((-32.71396786) + (13.536556*log(T_m)) - (5.73932875*((log(T_m))**2)) + (1.56315498*(log(T_m))**3) - (0.2877056*((log(T_m))**4)) + (3.48255977*(10**-2)*((log(T_m))**5)) - (2.63197617*(10**-3)*((log(T_m))**6)) + (1.11954395*(10**-4)*((log(T_m))**7)) - (2.03914985*(10**-6)*((log(T_m))**8)) )","Abel et al 1996", "0-8000"],
							   [51,"H+  + e- ->  H   + hv"    , "[H+]" , "[e-]" , ""      , "[H]"   , "[hv]" , ""     ,"k[1]((y[8]*kbev))", "lambda T_m: exp((-28.6130338) - (0.72411256*log(T_m)) - (2.02604473*(10**-2)*((log(T_m))**2)) - (2.38086188*(10**-3)*(log(T_m))**3) - (3.21260521*(10**-4)*((log(T_m))**4)) - (1.42150291*(10**-5)*((log(T_m))**5)) + (4.98910892*(10**-6)*((log(T_m))**6)) + (5.75561414*(10**-7)*((log(T_m))**7)) - (1.85676704*(10**-8)*((log(T_m))**8)) - (3.07113524*(10**-9)*((log(T_m))**9)))"   , "Abel et al 1996"                                           , "0-8000"],
							   [69,"He   + e- -> He+  + 2e-"  , "[He]"  ,"[e-]" , ""      , "[He+]" , "[e-]" , "[e-]" ,"k[2](y[8]*kbev)", "lambda T_m: exp(-44.09864886 + 23.91596563*log(T_m) - 10.7532302*(log(T_m)**2) + 3.05803875*(log(T_m)**3) - 0.56851189*(log(T_m)**4) + 6.79539123*(10**-2)*(log(T_m)**5) - 5.00905610*(10**-3)*(log(T_m)**6)+2.06723616*(10**-4)*(log(T_m)**7) - 3.64916141*(10**-6)*(log(T_m)**8))", "Abel et al 1996", "0-8000"],
							   [57,"He+  + e- -> He  + hv"    , "[He+]", "[e-]" , ""      , "[He]"  , "[hv]" , ""     ,"k[3](y[8]*kbev)", "lambda T_m: 3.925*(10**-13)*(T_m**-0.6353)"		                                               , "Abel et al 1996"                                              , "0-8000"],
							   [70,"He+  + e- -> He++ + 2e-"  , "[He+]" ,"[e-]" , ""      , "[He++]", "[e-]" , "[e-]" ,"k[4](y[8]*kbev)", "lambda T_m: exp(-68.71040990 + 43.93347633*log(T_m) - 18.4806699*(log(T_m)**2) + 4.70162649*(log(T_m)**3) - 0.76924663*(log(T_m)**4) + 8.113042*0.01*(log(T_m)**5) - 5.32402063*(10**-3)*(log(T_m)**6)+1.97570531*(10**-4)*(log(T_m)**7) - 3.16558106*(10**-6)*(log(T_m)**8))", "Abel et al 1996", "0-8000"],   
							   [55,"He++ + e- -> He+ + hv"    , "[He++]","[e-]" , ""      , "[He+]" , "[hv]" , ""     ,"k[5](y[8])", "lambda T_m: 3.36*(10**-10)*(T_m**-0.5)*((T_m/1000.0)**-0.2)*(1.0+(T_m/(10**6))**0.7)**-1"             , "Abel et al 1996"                                                     , "0-8000"],
							   [ 1,"H   + e- ->  H-  + hv"    , "[H]"  , "[e-]" , ""      , "[H-]"  , "[hv]" , ""     , "k[6](y[8]*kbev)", "lambda T_m:  1.429*(10**-18)*(T_m**0.7620)*(T_m**(0.1523*log10(T_m)))*(T_m**((-3.274*(10**-2))*((log10(T_m))**2))) if T_m < 0.51702 else 3.802*(10**-17)*(T_m**(0.1998*log10(T_m))*10**(4.0415*(10**-5)*((log10(T_m))**6)-(5.447*(10**-3)*((log10(T_m))**4))))"                                          , "Abel et al 1996"                                         , "0-8000"], 											
							   [ 2,"H-  + H  ->  H2  + e-"    ,"[H-]" , "[H]"  , ""      , "[H2]"  , "[e-]" , ""     , "k[7](y[8]*kbev)", "lambda T_m: 1.428*(10**-9) if T_m < 0.1 else exp((-20.06913897) + (0.22898*log(T_m)) + (3.5998377*(10**-2)*((log(T_m))**2)) - (4.55512*(10**-3)*(log(T_m))**3) - (3.10511544*(10**-4)*((log(T_m))**4)) +(1.0732940*(10**-4)*((log(T_m))**5)) - (8.36671960*(10**-6)*((log(T_m))**6)) + (2.23830623*(10**-7)*((log(T_m))**7)))" 					           , "Abel et al 1996"                                        , "0-8000"],
							   [ 3,"H+  + H  ->  H2+ + hv"    ,"[H+]" , "[H]"  , ""      , "[H2+]" , "[hv]" , ""     , "k[8](y[8]*kbev)", "lambda T_m: (3.833*(10**-16)*(T_m*1.8)) if T_m < 0.577 else 5.81*(10**-16)*((0.20651*T_m)**(-0.2891*(log10(0.20651*T_m))))" , "Abel et al 1996"                                           , "0-8000"], 
							   [ 4,"H2+ + H  ->  H2  + H+"    ,"[H2+]", "[H]"  , ""      , "[H2]"  , "[H+]" , ""     , "k[9](y[8])", "lambda T_m: 6.4*(10**-10)"					                                                       , "Abel et al 1996"                                         , "0-8000"], 
							   [ 6,"H2  + H+ ->  H2+ + H "    ,"[H2]" , "[H+]" , ""      , "[H2+]" , "[H]"  , ""     , "k[10](y[8]*kbev)", "lambda T_m:-24.24914687 + (3.40082444*log(T_m)) - (3.89800396*((log(T_m))**2)) + (2.04558782*(log(T_m))**3) - (0.541618285**((log(T_m))**4)) +(8.41077503*(10**-2)*((log(T_m))**5)) - (7.87902615*(10**-3)*((log(T_m))**6)) + (4.13839842*(10**-4)*((log(T_m))**7)) - (9.36345888*(10**-6)*((log(T_m))**8))"				                                       , "Abel et al 1996"                                       , "0-8000"], 
							   [ 7,"H2  + e- ->  2H  + e-"    ,"[H2]" , "[e-]" , ""      , "[H]"   , "[H]"  , "[e-]" , "k[11](y[8])", "lambda T_m: (5.6*(10**-11)*(T_m**0.5)*exp(-102124.0/T_m))"                                          , "Abel et al 1996"                               , "0-8000"], 
							   [ 8,"H2  + H  ->  2H  + H "    ,"[H2]" , "[H]"  , ""      , "[H]"   , "[H]"  , "[H]"  , "k[12](y[8]*kbev)", "lambda T_m: 1.067*(10**-10)*(T_m**2.012)*exp(-(4.463/T_m)*((1+(0.2472*T_m))**3.512))                                            "                                   , "Abel et al 1996"                                          , "0-8000"], 
							   [ 9,"H-  + e- ->  H   + 2e-"   ,"[H-]" , "[e-]" , ""      , "[H]"   , "[e-]" , "[e-]" , "k[13](y[8]*kbev)", "lambda T_m: exp((-18.01849334)+(2.3608522*log(T_m))-(0.28274430*((log(T_m))**2))+(1.62331664*(10**-2)*(log(T_m))**3)-(3.36501203*(10**-2)*((log(T_m))**4))+(1.17832978*(10**-2)*((log(T_m))**5))-(1.65619470*(10**-3)*((log(T_m))**6))+(1.06827520*(10**-4)*((log(T_m))**7))-(2.63128581*(10**-6)*((log(T_m))**8)))", "Abel et al 1996", "0-8000"],
							   [11,"H-  + H+ ->  H   + H "    ,"[H-]" , "[H+]" , ""      , "[H]"   , "[H]"  , ""     ,"k[14](y[8])", "lambda T_m: 7.0*(10**-8)*((T_m/100.0)**-0.5)"                              , "Abel et al 1996"                        , "0-8000"], 
							   [12,"H-  + H+ ->  H2+ + e-"    ,"[H-]" , "[H+]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[15](y[8]*kbev)", "lambda T_m: 2.291*(10**-10)*(T_m**-0.4) if T_m < 1.719 else (8.4258*(10**-10)*(T_m**-1.4)*exp(-1.301/T_m))"                                                         , "Abel et al 1996"                                        , "0-8000"], 
							   [13,"H2+ + e- ->  H   + H "    ,"[H2+]", "[e-]" , ""      , "[H]"   , "[H]"  , ""     ,"k[16](y[8])", "lambda T_m: 1.0*(10**-8) if T_m < 617 else 1.32*(10**-6)*(T_m**-0.76)"                        , "Abel et al 1996" 							               , "0-8000"], 
							   [14,"H2+ + H- ->  H   + H2"    ,"[H2+]", "[H-]" , ""      , "[H]"   , "[H2]" , ""     ,"k[17](y[8])", "lambda T_m: 5*(10**-7)*((100.0/T_m)**0.5)"                                                  , "Abel et al 1996"                                             , "0-8000"],
							   [13,"H-  + H  ->  H   + H + e-","[H-]", "[H]", ""    , "[H]"  ,  "[H]"  , "[e-]"  ,"k[18](y[8]*kbev)", "lambda T_m: exp((-20.37260896) + (1.13944933*log(T_m)) - (0.14210135*((log(T_m))**2)) + (8.4644554*(10**-3)*(log(T_m))**3) - (1.4327641*(10**-3)*((log(T_m))**4)) + (2.0122503*(10**-4)*((log(T_m))**5)) + (8.6639632*(10**-5)*((log(T_m))**6)) - (2.5850097*(10**-5)*((log(T_m))**7)) + (2.4555012*(10**-6)*((log(T_m))**8)) - (8.0683825*(10**-8)*((log(T_m))**9))) if T_m > 0.1 else 2.5634*(10**-9)*(T_m**1.78186)", "Abel et al 1996", "0-8000"],
							   [52,"H+  +  e- -> H  +  hv","[H]"  , "[hv]" , ""      , "[H+]"  , "[e-]" , ""     ,"k[19](y[7],Z)", "lambda T_r,Z: ((2.41*(10**+15)*(T_r**1.5)*exp(-179472/T_r))*((8.76*(10**-11))*((1.0+Z)**-0.58)))"                         , "Abel et al 1996"                                          "0-8000"],							   
							   [100,"H   + hv -> H+   + e-",  "[H]"   ,"" , ""      , "[H+]" ,  "[e-]" , ""     ,"flux[0]"],
							   [101,"He+ + hv -> He++ + e-",  "[He+]" ,"" , ""      , "[He++]", "[e-]" , ""     ,"flux[1]"],
							   [102,"He +  hv -> He+  + e-",  "[He]"  ,"" , ""      , "[He+]",  "[e-]" , ""     ,"flux[2]"],
							   [103,"H- +  hv -> H    + e-",  "[H-]"  ,"" , ""      , "[H]"  ,  "[e-]" , ""     ,"flux[3]"],
							   [104, "D +  hv -> D+  + e-",   "[D]"  ,"" , ""      , "[D+]",  "[e-]" , ""       ,"flux[4]"],
							   [105, "H2 +  hv -> H2+  + e-",  "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[5]"],
							   [106, "H2 +  hv -> H2+  + e-",  "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[6]"],
							   [107, "H2 +  hv -> H2+  + e-", "[H2]"  ,"" , ""      , "[H2+]",  "[e-]" , ""     ,"flux[7]" ],
							   [108, "H2+ + hv -> H    + H+", "[H2+]" ,"" , ""      , "[H]"  ,  "[H+]" , ""     ,"flux[8]" ],
							   [109,"H2+ + hv -> 2H+  + e-", "[H2+]" ,"" , ""      , "[H+]" ,  "[H+]" , "[e-]" ,"flux[9]"  ],
							   [110,"H2  + hv -> H2* -> H + H" ,"[H2]"  ,"" , ""      , "[H]"  ,  "[H]"  , ""     ,"flux[10]"],
							   [111,"H2  + hv -> H + H"        ,"[H2]"  ,"" , ""      , "[H]"  ,  "[H]"  , ""     ,"flux[11]" ]])
							
	elif (chemical_model ==  4): 	   
		reaction_matrix=array([
							   [51,"H+  +  e- -> H  +  hv","[H+]" , "[e-]" , ""      , "[H]"   , "[hv]" , ""     ,"k[0](Z)"      , "lambda Z: (8.76*(10**-11))*((1.0+Z)**-0.58)"  , "Galli & Palla 1998"  , "0-8000"], 
							   [52,"H+  +  e- -> H  +  hv","[H]"  , "[hv]" , ""      , "[H+]"  , "[e-]" , ""     ,"k[1](y[7],Z)", "lambda T_r,Z: ((2.41*(10**+15)*(T_r**1.5)*exp(-179472/T_r))*((8.76*(10**-11))*((1.0+Z)**-0.58)))"                         , "Galli & Palla 1998"                                          "0-8000"],
							   [ 1,"H+  +  e- -> H  +  hv","[H]"  , "[e-]" , ""      , "[H-]"  , "[hv]" , ""     ,"k[2](y[8])", "lambda T_m: 1.4*(10**-18)*(T_m**0.928)*exp(-T_m/16200.0)"                                                       , "Galli & Palla 1998"                                         , "0-8000"], 										
							   [16,"H+  +  e- -> H  +  hv","[H-]" , "[hv]" , ""      , "[H]"   , "[e-]" , ""     ,"k[3](y[7])", "lambda T_r: 0.11*(T_r**2.13)*exp(-8823.0/T_r)"                                                                  , "Galli & Palla 1998"                                             , "0-8000"], 
							   [ 2,"H+  +  e- -> H  +  hv","[H-]" , "[H]"  , ""      , "[H2]"  , "[e-]" , ""     ,"k[4](y[8])", "lambda T_m: 4.0*(10**-9)*(T_m**-0.17) if T_m > 300 else (1.5*10**-9)" 					                        , "Galli & Palla 1998"                                         , "0-8000"], 
							   [12,"H+  +  e- -> H  +  hv","[H-]" , "[H+]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[5](y[8])", "lambda T_m: 6.9*(10**-9)*(T_m**-0.35) if T_m < 8000 else 9.6*(10**-7)*(T_m**-0.9)"                                                         , "Glover & Abel 2008"                                         , "0-8000"], 					       					       
							   [11,"H+  +  e- -> H  +  hv","[H-]" , "[H+]" , ""      , "[H]"   , "[H]"  , ""     ,"k[6](y[8])", "lambda T_m: (5.7*(10**-6)*(T_m**-0.5))+(6.3*(10**-8))-(9.2*(10**-11)*(T_m**0.5))+(4.4*(10**-13)*T_m)"           , "Galli & Palla 1998"   , "0-8000"], 
							   [ 3,"H+  +  e- -> H  +  hv","[H+]" , "[H]"  , ""      , "[H2+]" , "[hv]" , ""     ,"k[7](y[8])", "lambda T_m: 10**(-19.38-(1.523*log10(T_m))+(1.118*((log10(T_m))**2))-(0.1269*((log10(T_m))**3)))"              , "Galli & Palla 1998"                      , "0-8000"], 
							   [17,"H+  +  e- -> H  +  hv","[H2+]", "[hv]" , ""      , "[H]"   , "[H+]" , ""     ,"k[8](y[7])", "lambda T_r: 2.0*10*(T_r**1.59)*exp(-82000.0/T_r)"                                                               , "Galli & Palla 1998" 			                  , "0-8000"], 
							   [ 4,"H+  +  e- -> H  +  hv","[H2+]", "[H]"  , ""      , "[H2]"  , "[H+]" , ""     ,"k[9](y[8])", "lambda T_m: 6.4*(10**-10)"					                                                                    , "Galli & Palla 1998"                           , "0-8000"], 
							   [13,"H+  +  e- -> H  +  hv","[H2+]", "[e-]" , ""      , "[H]"   , "[H]"  , ""     ,"k[10](y[8])", "lambda T_m: 2.0*(10**-7)*(T_m**-0.5)"                                                                           , "Glover & Abel 2008" 							               , "0-8000"], 
							   [19,"H+  +  e- -> H  +  hv","[H2+]", "[hv]" , ""      , "[H+]"  , "[H+]" , "[e-]" ,"k[11](y[7])", "lambda T_r: 90*(T_r**1.48)*(exp(-335000.0/T_r))"                                                 , "Galli D & Palla 1998"                                         , "0-8000"],
							   [15,"H+  +  e- -> H  +  hv","[H2]" , "[e-]" , ""      , "[H]"   , "[H-]" , ""     ,"k[12](y[8])", "lambda T_m: 2.7*(10**-8)*(T_m**-1.27)*exp(-43000.0/T_m)"                                                 , "Capitelli, M., Coppola, C. M., Diomede, P., & Longo, S. 2007" , "0-8000"],					       
							   [ 7,"H+  +  e- -> H  +  hv","[H2]" , "[e-]" , ""      , "[H]"   , "[H]"  , "[e-]" , "k[13](y[8])", "lambda T_m: 4.4*(10**-10)*(T_m**0.35)*exp(-102000.0/T_m)"                                          , "Trevian, C.S & Tennyson J.2002"                               , "0-8000"], 
							   [18,"H+  +  e- -> H  +  hv","[H2]" , "[hv]" , ""      , "[H2+]" , "[e-]" , ""     ,"k[14](y[7])", "lambda T_r: 2.9*(10**2)*(T_r**1.56)*exp(-178500.0/T_r)"                                          , "Galli D & Palla 1998"                                         , "0-8000"], 						   
							   [ 6,"H+  +  e- -> H  +  hv","[H2]" , "[H+]" , ""      , "[H2+]" , "[H]"  , ""     ,"k[15](y[8])", "lambda T_m: 3.0*(10**-10)*exp(-21050.0/T_m)"				                                                    , "Galli & Palla 1998"                           , "0-8000"], 
							   [26,"H+  +  e- -> H  +  hv","[D+]" , "[e-]" , ""      , "[D]"   , "[hv]" , ""     ,"k[16](Z)", "lambda Z: (8.76*(10**-11))*((1+Z)**-0.58)"                                                                    , "Galli & Palla 1998"         , "0-8000"],
							   [94,"H+  +  e- -> H  +  hv","[D]"   ,"[hv]" , ""      , "[D+]" ,  "[e-]" , ""     ,"k[17](y[7],Z)", "lambda T_r, Z: 2.41*(10**15)*(T_r**1.5)*exp(-179472.0/T_r)*(8.76*(10**-11))*((1+Z)**-0.58)"          ,  "Galli & Palla 1998"        , "0-8000"],
							   [27,"H+  +  e- -> H  +  hv","[D]"  , "[H+]" , ""      , "[D+]"  , "[H]"  , ""     ,"k[18](y[8])", "lambda T_m: 3.7*(10**-10)*(T_m**0.28)*exp(-43.0/T_m)"              , "Galli & Palla 1998"                 , "0-8000"],
							   [28,"H+  +  e- -> H  +  hv","[D+]" , "[H]"  , ""      , "[D]"   , "[H+]" , ""     ,"k[19](y[8])", "lambda T_m: 3.7*(10**-10)*(T_m**0.28)"            , "Galli & Palla 1998"               , "0-8000"],						 
							   [29,"H+  +  e- -> H  +  hv","[D]"  , "[H]"  , ""      , "[HD]"  , "[hv]" , ""     ,"k[20](y[8])", "lambda T_m: (10**-25)*(2.80202-(6.63697*log(T_m))+(4.75619*(log(T_m)**2))-(1.39325*(log(T_m)**3))+(0.178259*(log(T_m)**4))-(0.00817097*(log(T_m)**5)))"  , "Glover & Abel 2008"  , "0-8000"],
							   [30,"H+  +  e- -> H  +  hv","[D]"  , "[H2]" , ""      , "[HD]"  , "[H]"  , ""     ,"k[21](y[8])", "lambda T_m: 9.0*(10**-11)*exp(-3876.0/T_m)"                                                      , "Galli & Palla 1998"                                           , "0-8000"],
							   [31,"H+  +  e- -> H  +  hv","[HD+]", "[H]"  , ""      , "[HD]"  , "[H+]" , ""     ,"k[22](y[8])", "lambda T_m: 6.4*(10**-10)"                                                                     , "SLD98"                                                        , "0-8000"], 
							   [32,"H+  +  e- -> H  +  hv","[D+]" , "[H2]" , ""      , "[HD]"  , "[H+]" , ""     ,"k[23](y[8])", "lambda T_m: 2.1*(10**-9)"                           , "Galli & Palla 1998"                   , "0-8000"],
							   [34,"H+  +  e- -> H  +  hv","[HD]" , "[H+]" , ""      , "[D+]"  , "[H2]" , ""     ,"k[24](y[8])", "lambda T_m: 1.0*(10**-9)*exp(-464.0/T_m)"                                                        , "Galli D & Palla 1998"       , "0-8000"],
							   [71,"H+  +  e- -> H  +  hv","[Li+]" ,"[e-]" , ""      , "[Li]"  , "hv"   , ""     ,"k[25](y[8])", "lambda T_m: 1.036*(10**-11)*((((T_m/107.7)**0.5)*((1+((T_m/107.7)**0.5))**0.612)*(1+(T_m/(1.177*10**7))**0.5)**1.388)**-1.0)", "Galli & Palla 1998", "0-8000"],
							   [72,"H+  +  e- -> H  +  hv","[Li]"  ,"[hv]" , ""      , "[Li+]" , "[e-]" , ""     ,"k[26](y[7])", "lambda T_r: 1.3*(10**3)*(T_r**1.45)*exp(-60500.0/T_r)"                                                                       , "Galli & Palla 1998"  , "0-8000"],
							   [73,"H+  +  e- -> H  +  hv","[Li-]" ,"[H+]" , ""      , "[Li]"  , "[H]"  , ""     ,"k[27](y[8])", "lambda T_m: (6.3*(10**-6)*(T_m**-0.5)) - (7.6*(10**-9)) + (2.6*(10**-10)*(T_m**0.5)) + (2.7*(10**-14)*T_m)"                  , "Galli & Palla 1998" , "0-8000"],
							   [74,"H+  +  e- -> H  +  hv","[Li]"  ,"[e-]" , ""      , "[Li-]" , "[hv]" , ""     ,"k[28](y[8])", "lambda T_m: 6.1*(10**-17)*(T_m**0.58)*exp(-T_m/17200.0)"                                                                     , "Galli & Palla 1998", "0-8000"],
							   [75,"H+  +  e- -> H  +  hv","[Li-]" ,"[hv]" , ""      , "[Li]"  , "[e-]" , ""     ,"k[29](y[7])", "lambda T_r: 1.8*(10**2)*(T_r**1.4)*exp(-8100.0/T_r)"                                                                         , "Galli & Palla 1998", "0-8000"],
							   [76,"H+  +  e- -> H  +  hv","[Li]"  ,"[H]"  , ""      , "[LiH]" , "[hv]" , ""     ,"k[30](y[8])", "lambda T_m: ((5.6*(10**19)*(T_m**-0.15))+(7.2*(10**15)*(T_m**1.21)))**-1"                                                    , "Galli & Palla 1998"  , "0-8000"],
							   [77,"H+  +  e- -> H  +  hv","[LiH]" ,"[hv]" , ""      , "[Li]"  , "[H]"  , ""     ,"k[31](y[7])", "lambda T_r: 8.3*(10**4)*(T_r**0.3)*exp(-29000.0/T_r)"                                                                        , "Galli & Palla 1998"   , "0-8000"],
							   [78,"H+  +  e- -> H  +  hv","[Li-]" ,"[H]"  , ""      , "[LiH]" , "[e-]" , ""     ,"k[32](y[8])", "lambda T_m: 4.0*(10**-10)"                                                                                                   , "Galli & Palla 1998"   , "0-8000"],
							   [79,"H+  +  e- -> H  +  hv","[LiH]" ,"[H]"  , ""      , "[Li]"  , "[H2]" , ""     ,"k[33](y[8])", "lambda T_m: 2.0*(10**-11)"                                                                                                   , "Galli & Palla 1998"   , "0-8000"],
							   [80,"H+  +  e- -> H  +  hv","[Li+]" ,"[H]"  , ""      , "[LiH+]", "[hv]" , ""     ,"k[34](y[8])", "lambda T_m: 10**(-22.4+(0.999*log10(T_m))-(0.351*((log10(T_m))**2)))"                                                        , "Galli & Palla 1998"  , "0-8000"],
							   [81,"H+  +  e- -> H  +  hv","[Li]"  ,"[H+]" , ""      , "[LiH+]", "[hv]" , ""     ,"k[35](y[8])", "lambda T_m: 4.8*(10**-14)*(T_m**-0.49)"                                                                                      , "Galli & Palla 1998"  , "0-8000"],
							   [82,"H+  +  e- -> H  +  hv","[LiH+]","[e-]" , ""      , "[Li]"  , "[H]"  , ""     ,"k[36](y[8])", "lambda T_m: 3.8*(10**-7)*(T_m**-0.47)"                                                                                       , "Galli & Palla 1998"   , "0-8000"],
							   [83,"H+  +  e- -> H  +  hv","[LiH+]","[H]"  , ""      , "[Li+]" , "[H2]" , ""     ,"k[37](y[8])", "lambda T_m: 3.0*(10**-10)"                                                                                                   , "Galli & Palla 1998"   , "0-8000"],
							   [84,"H+  +  e- -> H  +  hv","[LiH+]","[hv]" , ""      , "[Li+]" , "[H]"  , ""     ,"k[38](y[7])", "lambda T_r: 7.0*(10**2)*exp(-1900.0/T_r)"                                                                                    , "Galli & Palla 1998"  , "0-8000"],					      
							   [61,"H+  +  e- -> H  +  hv","[He]" , "[H+]" , ""      , "[HeH+]", "[hv]" , ""     ,"k[39](y[8])", "lambda T_m: 7.6*(10**-18)*(T_m-0.5)"                                                                                         , "Galli & Palla 1998"   , "0-8000"],
							   [65,"H+  +  e- -> H  +  hv","[HeH+]","[H]"  , ""      , "[He]"  , "[H2+]", ""     ,"k[40](y[8])", "lambda T_m: 9.1*(10**-10)"                                                                        , "Galli & Palla 1998"                                           , "0-8000"],
							   [67,"H+  +  e- -> H  +  hv","[HeH+]","[hv]" , ""      , "[He]"  , "[H+]" , ""     ,"k[41](y[7])", "lambda T_r: 0.68*(T_r**1.5)*exp(-22750.0/T_r)"                                                     , "Galli & Palla 1998"                                             , "0-8000"]])
							  

	##Below here is the code that calculates the rates of change for each species from the chosen reaction matrix.
    ##It is performed by running a scan for each chemical species over the appropriate part of the reaction array (columns 3 to 8)
	##If the element being scanned is equal to the chemical species it is looking for, then it adds it to the equation for the
	##rate of change of a species. 
	
	##For example, if [H] is found on the LHS of reaction 1 (and the other reactant is [e-]), then it knows that D[H]/dt = -k1[H][e-].
	##By successively adding terms in this way (with the sign depending on if the species is a product or reactant) it can build up
	##the total equation for the rate of change of any particular species. Each species actually corresponds to an element an array
	##which the numerical integrator uses to find a solution, so it is in fact the array elements that are used in the final equation.
	##It knows what array elements o use by using the "dictionary" below.
	
	rlst ={'[H]': '*y[0]','[H+]': '*y[1]',"[H-]" : '*y[2]',"[H2]": '*y[3]',"[H2+]": '*y[4]',
		   "[e-]": '*y[5]',"[hv]": '*y[6]', "[D]":"*y[9]", "[D-]":"*y[10]", "[D+]": '*y[11]' ,
		   "[HD]": '*y[12]',"[HD+]":"*y[13]", "[He]":"*y[14]", "[He+]":"*y[15]", "[He++]":"*y[16]", "[HeH+]" : "*y[17]" , "[Li]":"*y[18]", "[Li+]":"*y[19]","[Li-]":"*y[20]","[LiH]":"*y[21]","[LiH+]":"*y[22]","[H2D+]":"*y[23]","[D2]":"*y[24]", "" : ""}																		#Dictionary converting species names into array elements ready for integration.
	   

	##List of chemical species to perform the scan for.
	species_list = ["[H]","[H+]","[H-]","[H2]","[H2+]","[e-]","[hv]","[D]","[D-]","[D+]","[HD]","[HD+]","[He]","[He+]","[He++]","[HeH+]","[Li]","[Li+]","[Li-]","[LiH]","[LiH+]","[H2D+]","[D2]"]					#List of the different chemical species that are contained within the reactions.
	
	##Generate an output file for debugging if needed.
	file = open('Reaction_writer_output.txt','w')																								#Open the file ready for writing to.
																						
	file.write("List of reaction equations for integration\n")
	
	##Creates an array to store each generated equation in.
	species_array = [0]*len(species_list)
	
    ##Perform the scan for each chemical species.
	
	for k in range(0, len(species_list)):																					#Iterate over the species list
		reactant = species_list[k]
		reaction = ''																										#A variable that contains the differential equation, initialised and set to be blank
		for j in range(0,(shape(reaction_matrix)[0])):																					
																															#Loop over all elements and rows in the matrix
			for i in range(2,8):
		
				if cmp(reaction_matrix[j][i],reactant) == 0:																#Check to see if the reactant has been found
						
					if   i == 2:
						reaction = reaction + '-' + str(reaction_matrix[j][8]) + str(rlst[reaction_matrix[j][i]]) + str(rlst[reaction_matrix[j][i+1]]) + str(rlst[reaction_matrix[j][i+2]])	#Add to the formula, depending upon its position in the reaction 
					elif i == 3:
						reaction = reaction + '-' + str(reaction_matrix[j][8]) + str(rlst[reaction_matrix[j][i]]) + str(rlst[reaction_matrix[j][i+1]]) + str(rlst[reaction_matrix[j][i-1]])
					elif i == 4:
						reaction = reaction + '-' + str(reaction_matrix[j][8]) + str(rlst[reaction_matrix[j][i]]) + str(rlst[reaction_matrix[j][i-1]]) + str(rlst[reaction_matrix[j][i-2]])
					elif i == 5:
						reaction = reaction + '+' + str(reaction_matrix[j][8]) + str(rlst[reaction_matrix[j][i-3]]) + str(rlst[reaction_matrix[j][i-2]]) + str(rlst[reaction_matrix[j][i-1]])
					elif i == 6:
						reaction = reaction + '+' + str(reaction_matrix[j][8]) + str(rlst[reaction_matrix[j][i-4]]) + str(rlst[reaction_matrix[j][i-3]]) + str(rlst[reaction_matrix[j][i-2]])
					elif i == 7:
						reaction = reaction + '+' + str(reaction_matrix[j][8]) + str(rlst[reaction_matrix[j][i-5]]) + str(rlst[reaction_matrix[j][i-4]]) + str(rlst[reaction_matrix[j][i-3]])
		
		species_array[k] = str(reaction)                                                                                    ##Store the generated equation within the array.

		
		file.write("\nThe equation for the rate of change of species ")													    #Write the reaction just created into the text file
		file.write((species_list)[k])
		file.write(" is:\n\n")
		file.write(str(reaction))
		file.write("\n")
		del reaction																										#Delete the reaction we just created, ready for the next iteration.


	file.close()								
	pickle.dump(species_array, open("reaction_equations.p", "wb"))                                                          #Store the generated equations in a pickle file
																												
	print("Different Model selected to previous use, recreating equations and loading rates...")                            #A useful message to the user of POSTREC if this process is occuring.
	
	
	
	###Cooling rates and photoionisation rates are treated differently, they are taken care of in the main function, but the rates and thier sources are shown here for reference purposes.
	
	cooling_rate_array = [[0, "Compton cooling"	                    , "-1.017*(10**-37)*(T_r**4)*(T_m-T_r)"                                               , "Shapiro & Kang 1987"          , "0-8000"],
						  [1, "Collisional ionisation cooling H"	, "1.27*(10**-21)*(T_m**0.5)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-157809.1/T_m)"        , "Black 1981, Cen 1992"         , "0-8000"],
						  [2, "Collisional ionisation cooling He"   , "9.38*(10**-22)*(T_m**0.5)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-285335.4/T_m)" 		  , "Black 1981, Cen 1992"         , "0-8000"],			 
						  [3, "Collisional ionisation cooling He+"  , "4.95*(10**-22)*(T_m**0.5)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-631515.0/T_m)"		  , "Black 1981, Cen 1992"         , "0-8000"],	
						  [4, "Recombination cooling H+"            , "8.70*(10**-27)*(T_m**0.5)*((T_m*10**-3)**-0.2)*((1+(T_m*10**-6)**0.7)**-1.0)"	  , "Black 1981, Cen 1992"         , "0-8000"],
						  [5, "Recombination cooling He+"           , "1.55*(10**-26)*(T_m**0.3647)"                                                      , "Black 1981, Cen 1992"         , "0-8000"],
						  [6, "Recombination cooling He++"          , "3.48*(10**-26)*(T_m**0.5)*((T_m*10**-3)**-0.2)*((1+(T_m*10**-6)**0.7)**-1.0)"	  , "Black 1981, Cen 1992"         , "0-8000"],
						  [7, "Dielectric recombination to He"      , "1.24*(10**-13)*(T_m**-1.5)*exp(-47000.0/T_m)*(1+0.3*exp(-94000.0/T_m))"            , "Black 1981, Cen 1992"         , "0-8000"],
						  [8, "Collisional excitation of H"         , "7.50*(10**-19)*((1+(T_m*10**-5)**0.5)**-1)*exp(-118348.0/T_m)"                     , "Black 1981, Cen 1992"         , "0-8000"],
						  [9, "Collisional excitation of He+"       , "5.54*(10**-17)*(T_m**-0.397)*((1+(T_m*10**-5)**0.5)**-1.0)*exp(-473638.0/T_m)"	  , "Black 1981, Cen 1992"         , "0-8000"],
						  [10,"Bremsstrahlung cooling"              , "(1.42*(10**-27)*(T_m**0.5)*(1.10+0.34*exp(-((5.50-log10(T_m))**2)/3.0)))"          , "Janev et al 1987"             , "0-8000"],
						  [11,"Cooling due to H2-H"                 , "10**(-103.0 + (97.59*log10(T_m))-(48.05*((log10(T_m))**2.0))+(10.80*((log10(T_m))**3))-(0.9032*((log10(T_m))**4)))", "Galli & Palla 1998", "10-10000"],
						  [12,"Collisional excitation H2+ with e-"  , "1.1*(10**-19)*(T_m**-0.34)*exp(-3025.0/T_m)"                                       , "Glover & Abel 2008"           , "0-2000"],
						  [13,"Collisional excitation H2+ with H"   , "1.36*(10**-22)*exp(-3152.0/T_m)"                                                   , "Glover & Abel 2008"           , "0-1000"],
						  [14,"Collisional excitation of H+ with H" , "10**(-36.42+5.95*log10(T_m)-0.526*((log10(T_m))**2))"                               , "Glover & Abel 2008"           , "0-8000"],
						  [15,"Cooling due to HD-H"                 , "10**(-42.45906+(21.90083*T_m)-(10.1954*(T_m**2))+(2.19788*(T_m**3))-(0.17286*(T_m**4)))", "Glover & Abel 2008"      , "0-8000"],
						  [16,"Cooling due to Li"                   , "0.0"                                                                                  , "unknown"                   , "0-8000"]]

						  
	photoionisation_array = [[0, "H   + hv -> H+   + e-",  "((6.30*(10**-18))*(((3.28846551*(10**15))/v)**4.0)*((exp(4-(4*atan((((v/(3.28846551*(10**15)))-1)**0.5))/(((v/(3.28846551*(10**15)))-1)**0.5))))/(1-exp((-2*pi)/(((v/(3.28846551*(10**15)))-1)**0.5)))))", "Abel et al 1997", "all frequencies"],
							 [1, "He+ + hv -> He++ + e-",  "((6.30*(10**-18))*(((3.28846551*(10**15))/v)**4.0)*((exp(4-(4*atan((((v/(3.28846551*(10**15)))-1)**0.5))/(((v/(3.28846551*(10**15)))-1)**0.5))))/(1-exp((-2*pi)/(((v/(3.28846551*(10**15)))-1)**0.5)))))", "Abel et al 1997", "all frequencies"],
							 [2, "He +  hv -> He+  + e-",  "7.42*(10**-18)*(((1.66*(v/(24.6/(4.13566*(10**-15))))**-2.05))-(0.66*((v/(24.6/(4.13566*(10**-15))))**-3.05)))" , "Abel et al 1997" , "v > 5.94825*10**15"],			 
							 [3, "H- +  hv -> H    + e-",  "7.928*(10**5)*((v-(0.755/(4.13566*(10**-15))))**1.5)*(v**-3)"		  , "Abel et al 1997"         , "v > 1.8255*10**14"],	
							 [4, "D +  hv -> D+  + e-",  "((6.30*(10**-18))*(((3.28846551*(10**15))/v)**4.0)*((exp(4-(4*atan((((v/(3.28846551*(10**15)))-1)**0.5))/(((v/(3.28846551*(10**15)))-1)**0.5))))/(1-exp((-2*pi)/(((v/(3.28846551*(10**15)))-1)**0.5)))))", "Abel et al 1997", "all frequencies"],
							 [5, "H2 +  hv -> H2+  + e-",  "6.2*(10**-18)*((4.135667516*(10**-15))*v)-(9.4*(10**-17))" , "Abel et al 1996" , "15.42 < hv < 16.50"],
							 [6, "H2 +  hv -> H2+  + e-",  "1.4*(10**-18)*((4.135667516*(10**-15))*v)-(1.48*(10**-17))", "Abel et al 1996" , "16.50 < hv < 17.70"],
							 [7, "H2 +  hv -> H2+  + e-",  "2.5*(10**-14)*(((4.135667516*(10**-15))*v)**-2.71)"        , "Abel et al 1996" , "hv > 17.7 eV"],
							 [8, "H2+ + hv -> H    + H+",  "10**((-1.6547717*(10**6))+(1.866033*(10**5)*log(v))-(7.8986431*(10**3)*(log(v)**2))+(148.73693*(log(v)**3))-(1.0513032*(log(v)**4)))"  , "Abel et al 1996", "valid for all frequency"],
							 [9, "H2+ + hv -> 2H+  + e-",  "10**((-16.926)-((4.528*(10**-2))*((4.135667516*(10**-15))*v))+(2.238*(10**-4)*(((4.135667516*(10**-15))*v)**2))+(4.425*(10**-7)*(((4.135667516*(10**-15))*v)**3)))"  , "Abel et al 1996" , "30<hv<90eV"],
							 [10,"H2  + hv -> H2* -> H + H" , "1.1*(10**8)"                                                                   , "Abel et al 1996"             , "neglecting selfshielding"],
							 [11,"H2  + hv -> H + H"        , "0.0"                                                            , "Abel et al 1997", "10-10000"]]
							 
	##Below we store all of the results in external files ready for use in POSTREC.						 
	pickle.dump(reaction_matrix,open("reactionratearray.p","wb"))					     
	
	return 0 ##End function
