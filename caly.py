#!/usr/bin/python3
#python 3 env... But only works when called from terminal???
#####################################################################
#
#Calibration of ALEX probes.
#####################################################################
#coding: latin-1
#Version 1.2


######################################
# Importing all the necessary modules.
######################################
import os #File management
import numpy as np #Numerical work in Python. Yeah!!!
import peakutils #Peak finder utilities...From https://bitbucket.org/lucashnegri/peakutils
#PeakUtils is stored in the folder "peakutils" in the folder of this program...
from scipy.signal import savgol_filter #To smooth the data...
from scipy.optimize import curve_fit #Data fitting.



#########################################
# These functions are for a short circuit
#########################################

###########################
# Current function:
###########################
def curr(t, alpha, omega, t0, CA):
	return CA * np.exp(- alpha * (t - t0)) * np.sin( omega * (t - t0) )

###########################
# Dot Current function:
###########################
def dcurr(t, alpha, omega, t0, B, C):
	return C + B * np.exp(- alpha * (t - t0)) * np.cos( omega * (t - t0) )
		
###########################
# Voltage function:
###########################
def volt(t, alpha, omega, t0, VA):
	return VA * np.exp(- alpha * (t - t0)) * np.cos( omega * (t - t0) )

###################################
# To convert np.arrays into strings:
###################################
def to_str(var):
    return str(list(np.reshape(np.asarray(var), (1, np.size(var)))[0]))[1:-1]



#############################################
#Starting program with some intialization:
#############################################

#This dictionary has the structure of the channels in the scope.
#It can be extracted from the HTML shot file, but it takes a lot of time.
chan_list = {'CH1':'Rog','CH2':'2Res','CH3':'3Res','CH4':'Cur'}

#Initializing the calibration vectors and their errors...
K_curr = []
K_curr_err  = [] 
K_Rog = [] 
K_Rog_err = []



###############################################################
# Main boby. Iterative over the RAW folders.
###############################################################

folderdata = os.path.dirname(os.path.abspath(__file__)) #Folder with the folders with the shots information data. Script must executed from it.

folderlist = os.walk(folderdata)
#os.walk returns an structure with the files, directories and fles inside the directories of a given path,
#   ordered as triplets, with the order used in the next FOR instruction.

for dirpath,dirnames,filenames in folderlist: #Make a loop through all the elements
	for folders in dirnames: #This is the structure of an iterator on the objects of dirnames, that I name folders
		if "RAW" in folders: #Folders with just channel info have the string "RAW" in them
			channels = os.listdir(os.path.join(dirpath,folders)) #Channels files list
			#I will assign the associated probe to every channel. 
			#It is assumed an order, chan_list dictionary, and that NO MORE THAN ONE SCOPE IS USED.
			for channels in channels:
				#print 'File: ',os.path.join(folders,channels) #Python 2.7 File being processed
				print('File: ',os.path.join(folders,channels)) #Python 3 File being processed
				#Read the channel CSV:
				with open(os.path.join(dirpath,folders,channels), 'r') as CHcsv: #This code line assures that the file is closed after its use...
				#CHcsv is the open file object.
					CH = np.loadtxt(CHcsv) #read the dat into a numpy array. It works like a charm.
				#Transforming time from seconds into microseconds for ALL the channels:
				CH[:,0] = CH[:,0]*1e6
				if 'CH1' in channels:
					CH1 = CH #That's the dot current without adjusting anything
					print(chan_list.get('CH1'))
					print('\n')
				if 'CH2' in channels:
					print(chan_list.get('CH2'))
					print('\n')	
					Vol_2Res = CH #That's the 2Resistive divider channel.
					Vol_2Res[:,1] = Vol_2Res[:,1]*1359
				if 'CH3' in channels:
					print(chan_list.get('CH3'))	
					print('\n')					
					Vol_3Res = CH #That's the 3Resistive divider channel.
					Vol_3Res[:,1] = Vol_3Res[:,1]*2400
				if 'CH4' in channels:
					CH4 = CH #That's the current without adjusting anything
					print(chan_list.get('CH4'))	
					print('\n')			
			Vol = Vol_3Res #Creating voltage through the wire
			Vol[:,1] = Vol_2Res[:,1] - Vol_3Res[:,1] #Voltage through the wire
			Vol_smooth = savgol_filter(Vol[:,1],31,2) #Smoothing voltage to obtain peaks.
			
			
			###############################################################
			# Finding the signal period
			###############################################################
			indexes = peakutils.indexes(Vol_smooth, thres = 1/1000, min_dist = 2500) #Index positions of peaks. It loks like it works quite well(In some cases)...
			
			Period = np.mean(np.diff(Vol[indexes[1:],0])) #Find the differences in time between all the peaks 
			#(but the first one, that is related with the spark gap) and calculate the mean of that. 
			#And you have the period of the signal.IN MICROSECONDS!!!!!!!
			Err_Period = np.std(np.diff(Vol[indexes[1:],0]))
			
			
			###############################################################
			#Calculating inductance:
			###############################################################
			L = ( (Period)**2 * 1e-12) / (2.2 * 1e-6 * 4 * np.pi**2 ) #Inductance in Farads
			Err_L = ( (Err_Period)**2 * 1e-12) / (2.2 * 1e-6 * 4 * np.pi**2 ) #Inductance error in Farads


			###############################################################
			# Fitting the current signal not calibrated (after spark gap on)
			###############################################################
			curr_params = curve_fit(curr,CH4[indexes[1]:,0],CH4[indexes[1]:,1], p0 = [0.13, (2 * np.pi)/Period, Vol[indexes[0],0], 0.25])
			#p0 is the initial guess parametr list for current funcion: alpha, omega, t0, A.
			#Notice that we do the fitting starting at time of the first peak of the voltage, 
			#related with the sparl Gap activation.
			sigma_curr_params = np.sqrt(np.diagonal(curr_params[1])) #Error in the uncalibrated current parameters			
			#Assinging readable names to the current function parameters:
			alpha = curr_params[0][0]
			omega = curr_params[0][1]
			omega_err = sigma_curr_params[1]
			t0 = curr_params[0][2]
			A = curr_params[0][3]
			
			
			###############################################################
			# Fitting the dotcurrent signal not calibrated (after spark gap on)
			###############################################################
			dcurr_params = curve_fit(dcurr,CH1[indexes[1]:,0],CH1[indexes[1]:,1], p0 = [0.13, (2 * np.pi)/Period, Vol[indexes[0],0], 0.25, 7e-2])
			sigma_dcurr_params = np.sqrt(np.diagonal(dcurr_params[1])) #Error in uncalibrated dcurr parameters
			#p0 is the initial guess parametr list for current funcion: alpha, omega, t0, A.
			#Notice that we do the fitting starting at time of the first peak of the voltage, 
			#related with the sparl Gap activation.

			
			###############################################################
			# Fitting the circuit voltage to find the maximun voltage, 
			#	initial charging voltage approx.
			# Notice that errors are calculated like the substraction
			#   of function(values + sigmas) - function.
			###############################################################
			volt_params = curve_fit(volt, Vol[indexes[1]:,0],Vol[indexes[1]:,1], p0 = [alpha, omega, t0, A])
			sigma_volt_params = np.sqrt(np.diagonal(volt_params[1])) #Error in the voltage parameters
						
			V0 = np.max(volt(Vol[:,0],*volt_params[0])) #Initial voltage.The charging one (Volts)
			V0_up = np.max(volt(Vol[:,0],*(volt_params[0]+sigma_volt_params))) #Calculating the errors...
			V0_err = V0_up - V0 #Voltage error
			
			
			###############################################################
			# Finding maximum current in no calibrated units, 
			#	to make the calibration of the RC filtered Rogowsky signal.
			###############################################################	
			curr_t1 = np.max(curr(CH4[indexes[0]:,0],*curr_params[0])) #Maximum current is the maximum of the fitted current from the begining of the channel data
			curr_t1_up = np.max(curr(CH4[indexes[0]:,0],*(curr_params[0]+sigma_curr_params))) #Calculating the errors...
			curr_t1_err = curr_t1_up - curr_t1 #Uncalibrated current error			
			
			
			###############################################################
			# Finding maximum dotcurrent in no calibrated units, 
			#	to make the calibration of the Rogowsky signal.
			###############################################################	
			dotcurr_t1 = np.max(dcurr(CH1[indexes[0]:,0],*dcurr_params[0])) #The same strategy with the dot current.
			dotcurr_t1_up = np.max(dcurr(CH1[indexes[0]:,0],*(dcurr_params[0]+sigma_dcurr_params)))
			dotcurr_t1_err = dotcurr_t1_up - dotcurr_t1 #Uncalibrated dot current error


			###############################################################
			# Calibration constants of the RC integrator and the Rogowsky:
			#	With error calculation, too
			###############################################################
			#Current probe calibration and error parameters::
			cali_curr = 	V0 / (  omega * 1e6 * L * curr_t1  )
			cali_curr_err = ( V0_err/ (  omega * 1e6 * L * curr_t1  ))  +  ( (V0 * Err_L) / (omega * 1e6 * L**2 * curr_t1) ) + ( (V0 * curr_t1_err )/(omega * 1e6 * L * curr_t1**2) ) + ( (V0 * omega_err)/(omega**2 * 1e12 * L * curr_t1 ) )		
			#Appending parameters if there is nothing extrange with the error:
			if np.isfinite(cali_curr_err):
				K_curr.append(cali_curr) #Calibration of current probe
				K_curr_err.append(cali_curr_err ) 
				
			#Rogowsky probe calibration and error parameters:
			cali_rog = V0 / ( L * dotcurr_t1)
			cali_rog_err =  (V0_err /( L * dotcurr_t1)) + ((V0 * Err_L)/(L**2 * dotcurr_t1)) + ( (V0 * dotcurr_t1_err)/(L * dotcurr_t1**2) )
			#Appending data if there is nothing extrange with the error:
			if np.isfinite(cali_rog_err):
				K_Rog.append(cali_rog) #Calibration of dot current probe, and error next line:
				K_Rog_err.append(cali_rog_err)



##############################################################
#Making weithegthed averaga and saving data into a TXT file:
##############################################################
			
#Weitgthed average of the Current calibration:
K_current = np.sum( K_curr/np.power(K_curr_err,2) ) / np.sum(1/np.power(K_curr_err,2))
K_current_err = 1/np.sqrt( np.sum(1/np.power(K_curr_err,2)) )

#Weitghted average of the Rogowsky calibration:
K_Rogos = np.sum( K_Rog/np.power(K_Rog_err,2) ) / np.sum(1/np.power(K_Rog_err,2))
K_Rogos_err = 1/np.sqrt( np.sum(1/np.power(K_Rog_err,2)) )

#Saving data:
with open("ALEX-Rogowsky-Current-cali.txt","w") as save_file:
	save_file.write("Calibration values for Rogowsky and Current channels(S.I. Units)\n\n")
	save_file.write("Rogowsky constant: (A/s / V)\n")
	save_file.write("{:0.4e}".format(K_Rogos)+ " +- "+"{:0.4e}".format(K_Rogos_err)+"\n\n\n")
	save_file.write("Current constant: (A / V)\n")
	save_file.write("{:0.4e}".format(K_current)+" +- "+"{:0.4e}".format(K_current_err)+"\n")
	save_file.write(" \n K_curr \n"+to_str(K_curr)+"\n")
	save_file.write("\n K_curr_err \n"+to_str(K_curr_err)+"\n\n")			
	save_file.write("\n K_Rog \n"+to_str(K_Rog))
	save_file.write("\n K_Rog_err \n"+to_str(K_Rog_err))			

print("Data saved and stored in file `ALEX-Rogowsky-Current-cali.txt'")


#And tha...tha...that's all folks!!
