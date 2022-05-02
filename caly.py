#!/usr/bin/python3
"""
Calibration of ALEX electrical probes.
In principle, it works by adjusting the uncalibrated current and current derivative signals to an exponential function, 
and using these parameters to obtain the adjutment constants.
In general, to do the fitting of all the signals, it is necessary to start after the first peak if the value, 
because of the presence of the Spark Gap that is not on the ideal RLC circuit.

This script must be executed from a folder that contains all the short circuit shots.

coding: UTF-8
Version 2
"""



# Modules:
import os #File management
import numpy as np #Numerical work in Python. Yeah!!!
import matplotlib.pyplot as plt #Plotting capabilities
from scipy.signal import savgol_filter #To smooth the data...
from scipy.optimize import curve_fit #Data fitting.
from scipy.signal import find_peaks #Peak finding




# Short circuit functions for voltage, current, and its derivative, called dcurrent plus error functions:


#Current function:
def curr_function(t, alpha, omega, t0, CA):
	return CA * np.exp(- alpha * (t - t0)) * np.sin( omega * (t - t0) )
	
#Current error function:
def err_curr_function(t, alpha, omega, t0, CA, err_alpha, err_omega, err_t0, err_CA):
	exponente = np.exp( -alpha * (t -t0) )
	seno = np.sin(omega * (t - t0) )
	
	return ( err_CA * exponente * seno )
	+ ( err_alpha * CA * (t - t0) * exponente * seno )
	+ err_t0 * (  ( CA * alpha * exponente * seno )  + (  CA * omega * exponente * np.cos( omega * (t - t0) ) )  )
	+ ( err_omega * CA * (t - t0) * exponente * np.cos( omega * (t - t0) ) )

#Dot Current function: (Its derivative)
def dcurr_function(t, alpha, omega, t0, B, C):
	return C + B * np.exp(- alpha * (t - t0)) * np.cos( omega * (t - t0) )

#Dot current error function:
def err_dcurr_function(t, alpha, omega, t0, B, C, err_alpha, err_omega, err_t0, err_B, err_C):
	exponente = np.exp(  -alpha * ( t - t0 )  )
	coseno = np.cos(  omega * ( t - t0 )  )
	return err_C + ( err_B * exponente * coseno )
	+ ( err_alpha * B * ( t-t0 ) * exponente * coseno )
	+ ( err_omega * B *  ( t-t0 ) * exponente * np.sin( omega * ( t - t0 ) ) )
	+ B * err_t0 * (   ( alpha * exponente * coseno) + ( omega * exponente * np.sin( omega * ( t-t0 ) )  )   )

#Voltage function:
def volt_function(t, alpha, omega, t0, VA):
	return VA * np.exp(- alpha * (t - t0)) * np.cos( omega * (t - t0) )


#To convert np.arrays into strings:
def to_str(var):
    return str(list(np.reshape(np.asarray(var), (1, np.size(var)))[0]))[1:-1]




#Initializing the calibration vectors and their errors...
K_curr = []
K_curr_err  = [] 
K_rog = [] 
K_rog_err = []
K_L = []
K_L_err = []
K_per = []
K_per_err = []

#Path were folders with the shots information data are:
folderdata = os.path.dirname(os.path.abspath(__file__)) 

folderlist = os.walk(folderdata)
#os.walk returns an structure with the files, directories and files inside the directories of a given path,
#   ordered as triplets, with the order used in the next FOR instruction.

for dirpath,dirnames,_ in folderlist: #Make a loop through all the shot folders(dirnames) inside dirpath
	for folders in dirnames: #This is the structure of an iterator on the objects of dirnames, that I name folders, and are the folders within the main folder
		if "RAW" in folders: #Folders with just channel info have the string "RAW" in their name
			#So, for each shot do the calibration:
			
			#Which are the channels in the folder? Answer:
			channels = os.listdir(os.path.join(dirpath,folders)) #Channels files list
			print("Working with shot in folder: ", folders)
			#Now assign the associated probe to every channel. 
			#It is assumed an order that is fixed between shots and previously known and that NO MORE THAN ONE SCOPE IS USED.
			for channels in channels:
				#Read the channels CSV and addapt them to the magnitudes that they have:
				with open(os.path.join(dirpath,folders,channels), 'r') as CHcsv: #This code line assures that the file is closed after its use...
				#CHcsv is the open file object.
					CH = np.loadtxt(CHcsv) #read the dat into a numpy array. It works like a charm.
				#Transforming time from seconds into microseconds for ALL the channels: (It makes the fitting much easier)
				CH[:,0] = CH[:,0]*1e6
				if 'CH1' in channels:
					dot_current= CH #That's the dot current without adjusting anything					
				elif 'CH2' in channels:
					current= CH #That's the current integrator channel
				elif 'CH3' in channels:
					unscaled_voltage = np.array(CH[:,1])	#Unscaled voltage used for calculating the error of V0.
					voltage = CH #That's the DI04 resistive divider channel.
					voltage[:,1] = voltage[:,1]*6930
			#Smoothing voltage to obtain peaks. ONLY ONE ** COLUMN!!!!!		
			voltage_smooth = savgol_filter(voltage[:,1],31,2) #Smoothing voltage to obtain peaks. ONLY ONE ** COLUMN!!!!!
		
			# Finding the signal period:
			#Index positions of positive peaks. It loks like it works quite well
			#'picos' is an array with the indexes positions and the ',_' is for removing the other things that the function produces as output.
			picos , _= find_peaks(voltage_smooth, height = np.max(voltage_smooth)/5, width=10)
			
			#Creating the plot for checking visually the results (VOLTAGE):
			volt_fig, volt_axes = plt.subplots() #With this I create the figure with all inside ('volt_fig') and the axes with the data ('volt_axes')
			
			volt_axes.plot(voltage_smooth) #Smoothed voltage
			volt_axes.plot(picos, voltage_smooth[picos], "+")
			volt_axes.set_title("Voltage in shot "+folders[1:-4])
			volt_axes.set_xlabel("time (AU, points)")
			volt_axes.set_ylabel("Volts")
			volt_fig.savefig("volt"+folders[:-4]+".pdf")	
			


			#Find the differences in time between all the peaks
			#(but the first one, that is related with the spark gap) and calculate the mean of that. 
			#And you have the period of the signal. IN MICROSECONDS!!!!!!!
			#It does not use the first peak, tooo close to the Spark Gap discharge
			Period = np.mean(  np.abs (np.diff(voltage[picos[1:-1], 0 ])  )   )  
			print("Periods: ", Period)
			K_per.append(Period)
			
			#Its error is the standard deviation:
			err_Period = np.std(np.diff(voltage[picos[1:-1], 0] )  )
			print("Err Period: ", err_Period)
			K_per_err.append(err_Period)
			
			#position from where strat alll the fitting: (After Spark Gap discharge, so circuit will be closer to the ideal)
			pos_start_fit = picos[1]
			
			#Calculating inductance:
			L = ( (Period)**2 * 1e-12) / (2.2 * 1e-6 * 4 * np.pi**2 ) #Inductance in Farads
			K_L.append(L)
			
			err_L = (2 * err_Period * Period * 1e-12) / (2.2 * 1e-6 * 4 * np.pi**2 ) #Inductance error in Farads
			K_L_err.append(err_L)
			
			print("L: ", L, "L en nanoFar:", L *1e9)
			print("Error L: ", err_L)
			print("")
			
			#Fitting the current signal not calibrated (after spark gap on)		
			#p0 is the initial guess parametr list for current funcion: alpha, omega, t0, A.
			#Notice that we do the fitting starting at time of the first peak of the voltage, 
			#related with the Spark Gap activation.
			current_fit_params, current_fit_covariance = curve_fit(curr_function, current[pos_start_fit:-1,0],current[pos_start_fit:-1,1], p0 = [0.13,  (2 * np.pi)/Period,  current[pos_start_fit,0],  0.5] )

			#Creating the plot for checking visually the results (CURRENT):
			curr_fig, curr_axes = plt.subplots() #With this I create the figure with all inside ('volt_fig') and the axes with the data ('volt_axes')
			
			curr_axes.plot( current[:,0], current[:,1] ) #current
			curr_axes.plot( current[pos_start_fit:-1,0], curr_function(current[pos_start_fit:-1,0], *current_fit_params) ) #fitting
			curr_axes.set_title("Current in shot "+folders[1:-4])
			curr_axes.set_xlabel("time (µs)")
			curr_axes.set_ylabel("Current (V)")
			curr_fig.savefig("curr"+folders[:-4]+".pdf")				
			
			#Error in the uncalibrated current parameters:
			sigma_curr_fit_params = np.sqrt(np.diagonal(current_fit_covariance) ) 
			
			#Fitting the dotcurrent signal not calibrated (after spark gap on)
			#p0 is the initial guess parametr list for current funcion: alpha, omega, t0, B, C.			
			dcurr_fit_params, dcurr_fit_covariance = curve_fit(dcurr_function, dot_current[pos_start_fit:-1,0],dot_current[pos_start_fit:-1,1], 
			p0 = [0.13,  (2 * np.pi)/Period,  dot_current[pos_start_fit,0],  0.5,  0] )
			# ~ print("Dot current params: ", dcurr_fit_params)
			#Error in the uncalibrated dot current parameters:
			sigma_dcurr_fit_params = np.sqrt(np.diagonal(dcurr_fit_covariance) ) 			
			
			#Creating the plot for checking visually the results (DOT-CURRENT):
			dcurr_fig, dcurr_axes = plt.subplots() #With this I create the figure with all inside ('volt_fig') and the axes with the data ('volt_axes')
			
			dcurr_axes.plot( dot_current[:,0], dot_current[:,1] ) #dot current
			dcurr_axes.plot( dot_current[:,0], dcurr_function(dot_current[:,0], *dcurr_fit_params) ) #fitting
			dcurr_axes.set_title("dot current in shot "+folders[1:-4])
			dcurr_axes.set_xlabel("time (µs)")
			dcurr_axes.set_ylabel("dot current (V)")
			dcurr_fig.savefig("dot_curr"+folders[:-4]+".pdf")		
			
			#Naming some parameters initially for the voltage, following the current parameters (they should be the same...)
			alpha = current_fit_params[0]
			omega = current_fit_params[1]
			t0 = current_fit_params[2]

			
			#Fitting the circuit voltage to find the maximun voltage, initial charging voltage approx.
			#p0 is the initial guess parametr list for current funcion: alpha, omega, t0, VA.	
			volt_params, volt_fit_covariance = curve_fit(volt_function, voltage[pos_start_fit:-1,0],  voltage_smooth[pos_start_fit:-1],  p0 = [alpha,  omega,  t0,  0.5])			
			#Errors on voltage parameters:
			sigma_volt_fit_params = np.sqrt(np.diagonal(volt_fit_covariance))
			
			#Initial voltage. The charging one, supposedely (Volts)
			V0 = voltage_smooth[pos_start_fit]
			#Error with the error of the probe (first term) and the error of the channel (3%, second term)
			err_V0 = 30 * unscaled_voltage[pos_start_fit] + 6930 * 0.03 * unscaled_voltage[pos_start_fit]
			
			#Finding maximum current in no calibrated units to make the calibration of the RC integrator:
			#Pay attention to the fact that the '-50' can produce errors if  'pos_start_fit' is less than this...
			current_max = curr_function(current[pos_start_fit-50:-1, 0], *current_fit_params) .max()
			#Its index position:
			current_max_position = curr_function(current[pos_start_fit-50:-1, 0], *current_fit_params) .argmax() + pos_start_fit-50
			#Its error:
			err_current_max = err_curr_function(  current[current_max_position,0], *current_fit_params, *sigma_curr_fit_params  )

			#Finding maximum dot current in no calibrated units to make the calibration of the rogowski coil:
			#Pay attention to the fact that the '-50' can produce errors if  'pos_start_fit' is less than this...
			dcurrent_max = dcurr_function(dot_current[pos_start_fit-50:-1, 0], *dcurr_fit_params) .max()
			#Its index position:
			dcurrent_max_position = dcurr_function(dot_current[pos_start_fit-50:-1, 0], *dcurr_fit_params) .argmax() + pos_start_fit-50
			#Its error:
			err_dcurrent_max = err_dcurr_function(  dot_current[dcurrent_max_position,0], *dcurr_fit_params, *sigma_dcurr_fit_params  )
			
			#Making the average between the three versions of alpha, omega and t0:
			average_alpha = ( current_fit_params[0] + dcurr_fit_params[0] + volt_params[0] ) / 3			
			err_average_alpha = (1 / 3) * np.sqrt(  sigma_curr_fit_params[0] + sigma_dcurr_fit_params[0] + sigma_volt_fit_params[0] )
			
			average_omega =  ( current_fit_params[1] + dcurr_fit_params[1] + volt_params[1] ) / 3			
			err_average_omega =  (1 / 3) * np.sqrt(  sigma_curr_fit_params[1] + sigma_dcurr_fit_params[1] + sigma_volt_fit_params[1] )
			
			average_t0 =  ( current_fit_params[2] + dcurr_fit_params[2] + volt_params[2] ) / 3			
			err_average_t0 =  (1 / 3) * np.sqrt(  sigma_curr_fit_params[2] + sigma_dcurr_fit_params[2] + sigma_volt_fit_params[2] )

		
			#Calibration constants of the RC integrator and the Rogowsky with their errors:

			#Current probe calibration and error parameter: (1e6 and 1e12 are because the time constants are in µs and SI uses seconds)
			curr_constant = 	V0 / (  average_omega * 1e6 * L * current_max )
			
			curr_constant_err = ( err_V0/ (  average_omega * 1e6 * L * current_max )  )  
			+  (   (V0 * err_L)  /  (average_omega * 1e6 * L**2 * current_max)   ) 
			+ (   (V0 * err_current_max) / (average_omega * 1e6 * L * current_max**2)   ) 
			+ (   (V0 * err_average_omega) / (average_omega**2 * 1e6 * L * current_max )   )					
			
			print("Current scaling and error: ", curr_constant,  curr_constant_err )
			
			#Appending data to the list of calibrations:
			K_curr.append(curr_constant)
			K_curr_err.append(curr_constant_err)
			
			#Rogowsky probe calibration and error parameters:	
			rog_constant  = V0 / ( L * dcurrent_max )			
			rog_constant_err = err_V0 / (L * dcurrent_max )  +  ( V0 * err_L ) /( L**2 * dcurrent_max )  +  ( V0 * curr_constant_err ) / ( L * curr_constant**2 )
			
			#Appending data to the list of calibrations:
			K_rog.append(rog_constant)
			K_rog_err.append(rog_constant_err)
	
			print("Rogowski scaling and error: ", rog_constant,  rog_constant_err )
			print("")
			print("")



#Making weithegthed average and saving data into a TXT file after passing through all the RAW folders:

#For current constant:
K_current = np.sum( K_curr/np.power(K_curr_err,2) )  /  np.sum(1/np.power(K_curr_err,2))
K_current_err = 1/np.sqrt( np.sum(1/np.power(K_curr_err,2)) )

#For Rogowski coil:
K_rogowski = np.sum( K_rog/np.power(K_rog_err,2) )  /  np.sum( 1 / np.power(K_rog_err,2) )
K_rogowski_err = 1  /  np.sqrt(  np.sum( 1 / np.power(K_rog_err,2) )   )

#For the inductance:
K_inductance = np.sum( K_L/np.power(K_L_err,2) )  /  np.sum( 1 / np.power(K_L_err,2) )
K_inductance_err = 1  /  np.sqrt(  np.sum( 1 / np.power(K_L_err,2) )   )

#For the period:
K_period = np.sum( K_per/np.power(K_per_err,2) )  /  np.sum( 1 / np.power(K_per_err,2) )
K_period_err = 1  /  np.sqrt(  np.sum( 1 / np.power(K_per_err,2) )   )

#Saving all the data:
with open("ALEX-Rogowsky-Current-cali.txt","w") as save_file:
	save_file.write("Calibration values for Rogowski and Current channels(S.I. Units)\n\n")
	
	save_file.write("Rogowski constant: (A/s / V)\n")
	save_file.write("{:0.4e}".format(K_rogowski)+ " +- "+"{:0.4e}".format(K_rogowski_err)+"\n")
	
	save_file.write("Current constant: (A / V)\n")
	save_file.write("{:0.4e}".format(K_current)+" +- "+"{:0.4e}".format(K_current_err)+"\n\n")
	
	save_file.write("Inductance (nF)\n")
	save_file.write("{:0.4e}".format(K_inductance*1e9)+" +- "+"{:0.4e}".format(K_inductance_err*1e9)+"\n")
	
	save_file.write("Period (us)\n")
	save_file.write("{:0.4e}".format(K_period)+" +- "+"{:0.4e}".format(K_period_err)+"\n\n\n\n")	
	
	save_file.write(" \n K_curr \n"+to_str(K_curr)+"\n")
	save_file.write("\n K_curr_err \n"+to_str(K_curr_err)+"\n\n")			
	save_file.write("\n K_Rog \n"+to_str(K_rog))
	save_file.write("\n K_Rog_err \n"+to_str(K_rog_err))			

print("Data saved and stored in file `ALEX-Rogowsky-Current-cali.txt'")

#And tha...tha...that's all folks!!
