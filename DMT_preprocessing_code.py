#!/usr/bin/env python

## Alright people, prepare yourself.
## I am going to:

## Populate my mseed files with event and station information
## Write them to SaC files
## Time shift the SAC files (I KNOW RIGHT)
## Remove the instrument response, detrend, remove mean, filter etc...
## ROTATE THE TRACES
## Then leave them all in the "processed" directory
## Good stuff

#### -----This needs to be run from the "raw" directory as given in the ObspyDMT file----- ####


## import useful packages:

import obspy
import os
from glob import *

## To take in arguments from the terminal
import argparse

# To get relative path
from pathlib import Path

## to only read certain lines in file:
from itertools import islice

##distance baz calculator
from utilities_Array import deg_km_az_baz

# To dictate headers
from obspy.core.util import AttribDict

## Web service packages
from obspy.clients.iris import Client

# TauP
from obspy.taup import TauPyModel

# my model of choice is prem :D
model = TauPyModel(model="prem")

# Make client - just for giggles
client = Client()

import shutil
import numpy

import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description='Preprocessing script for data retrieved from obspy DMT')

parser.add_argument("-fmin","--minimuum_frequency", help="Enter the lowest frequency to be bandpass filtered.", type=float, required=False, action="store", default = 0.05)

parser.add_argument("-fmax","--maximuum_frequency", help="Enter the highest frequency to be bandpass filtered.", type=float, required=False, action="store", default = 1.0)

parser.add_argument("-p","--phase_list", help="Enter the phases you want the travel time predictions for.", type=str, nargs='+', required=True, action="store", default=[])

parser.add_argument("-tap","--taper_percentage", help="Enter the percentage of the trace at each end you want to taper as a decimal.", type=float, required=False, action="store", default=0.05)

parser.add_argument("-rot","--rotate", help="Enter whether you want the traces to be rotated or not.", required=False, action="store_true", default=True)

parser.add_argument("-split","--splitting", help="Enter whether you want the traces stored in a file and be compativle with SHEBA.", type=bool, required=False, action="store", default=False)

parser.add_argument("-detrend","--detrend_T_F", help="Enter whether you want to detrend the trace or not.", type=bool, required=False, action="store", default=True)

parser.add_argument("-v", "--verbose", help="Increase verbosity of output, just --verbose is enough to turn on.",action="store_true")

args = parser.parse_args()

min_freq = args.minimuum_frequency
max_freq = args.maximuum_frequency
phases = args.phase_list[:]
taper_percent = args.taper_percentage
Q_rotate = args.rotate
Q_detrend = args.detrend_T_F
Q_splitting = args.splitting
major=True
print(min_freq, max_freq,phases,taper_percent,Q_rotate,Q_detrend)

### Find the files with the station and event information in them

## Event dir name:
ev_relpath = "EVENTS-INFO/catalog.txt"

## full path of current directory
pwd = os.getcwd()

## path to directory with the Event directory
base_path = str(Path(pwd).parents[1]) + "/"

ev_path = os.path.join(base_path, ev_relpath)


## path to directory with station XML files :D

## gets the dierctory outside the current one:
base_dir = os.path.split(pwd)[0]
station_dir = os.path.join(base_dir,"resp/")
#print(station_dir)

## get event id
ev_name = os.path.basename(base_dir)
#print(ev_name)

## make a list for the stations
station_temp = []

### BEGIN THE LOOP



for M_file in glob('*BH*'):
	print(M_file)
	bits = M_file.split('.')
	station_name=bits[1]
	network_name=bits[0]

### stream the file
	stream = obspy.read(M_file)
## extract the trace
	trace = stream[0]
## copy the trace for no real reason.
	s_time = trace.stats.starttime
	e_time = trace.stats.endtime
	trace2= trace.trim(starttime=s_time)
## open the xml file as an inventory
	inv = obspy.read_inventory("%sSTXML.%s" %(station_dir, M_file))


	network=inv[0]
	station=network[0]
	channel=station[0]

## For each of the channels the station has, if the code (BHZ BHE BHN) is the same
# take the information
# this was the bane of my existance for an hour
	for chnl in station:

		if chnl.code==trace.stats.channel:
			channel_name = chnl.code
			response=channel.response
			dip= chnl.dip
			inclination = chnl.dip + 90

	## take the coordinates of the station
	stla=station.latitude
	stlo=station.longitude
	stel=station.elevation
	station_temp.append(station.code)

	## MORE SANITY CHECKS
	with open(ev_path) as f:
	## only read lines from line 6
	## This will remove all of the nonesense at the start of the file

		for line in islice(f,5,None):
			if not line.startswith('-') and ev_name in line:
				number,event_id,datetime,latitude,longitude,depth,magnitude,magnitude_type,author,flynn_region,mrr,mtt,mpp,mrt,mrp,mtp,stf_func,stf_duration,t1,t2 = line.split(',')
				evla = latitude
				evlo = longitude
				evdp = depth
				mag = magnitude
				event_time = obspy.UTCDateTime(datetime)


	#distances = client.distaz(stalat=stla, stalon=stlo, evtlat=evla, evtlon=evlo)
	dist_deg,dist_km,az,baz = deg_km_az_baz(lat1=float(evla),lon1=float(evlo),lat2=float(stla),lon2=float(stlo))

	#### Want to find the origin time relative to the trace period
	### i.e. we want to find the difference between the event and
	## the start of the trace
	O = obspy.UTCDateTime(event_time) - obspy.UTCDateTime(trace.stats.starttime)
	### I think O is equal to zero...

	### Going to need to loop over all of the phases ###
	## I don't think this is the most efficient way of doing things but I'm ##
	## going to make a list of lists, where i have [phase name, travel time, number (to be put on t) ##
	# good luck my friend #

	ph_time_list = []

#### Populate information into SAC part ####

	trace2.stats.sac = AttribDict({'baz': float(baz), 'az': float(az), 'kstnm': str(station),'cmpinc': inclination, 'gcarc': str(dist_deg),
							'dist': str(dist_km), 'evla': evla, 'evlo': evlo, 'evdp': evdp, 'stla': stla, 'stlo': stlo, 'stel': stel,'o': O,
							'stdp': dip, 'iztype':9, 'mag': mag})


	for t,phase in enumerate(phases):
		time_prediction =  model.get_travel_times(source_depth_in_km=float(evdp), distance_in_degree=float(dist_deg), receiver_depth_in_km=0.0, phase_list=[phase])

		## if the phase is present add the origin time to it
		if time_prediction:
			time_pred_total = O + time_prediction[0].time
		else:
			time_pred_total = None

		print(phase, time_pred_total)
		### add these values to a list ###
		ph_time_list.append([phase, time_pred_total])

		tn = str('t'+str(t+1))
		ktn = str('k'+tn)

		print(tn,ktn)

		# i.e. for SKKS major arc.
		if major == True:
			if len(time_prediction) > 1:
				# only use the first major arc arrival
				time_pred_total_2 = O + time_prediction[1].time
				phase_label2 = "%s_Major" %phase
				ph_time_list.append([phase_label2, time_pred_total_2])



	### for each loop se tthe attributes that I want ###
	# set predicted travel time
	for t,pt_list in enumerate(ph_time_list):
		time_pred = pt_list[1]
		phase_label = pt_list[0]
		tn = str('t'+str(t+1))
		ktn = str('k'+tn)
		setattr(trace2.stats.sac, tn, time_pred)
		setattr(trace2.stats.sac,ktn, phase_label)



	if trace2.stats.channel == "BHZ":
		trace2.stats.sac.cmpaz = 0.0

	print(trace2.stats)

#### process the data, filter etc...

	trace_rr = trace2.copy()

	#question=input("Velocity(v) or Displacement(d)?")

	#trace_rr.remove_response(inventory=inv, output='VEL')

	trace_rr.remove_response(inventory=inv, output='DISP')

	trace_detrend = trace_rr.copy()
	if Q_detrend:
	#### detrend/demean stuff
		trace_detrend.detrend(type='demean')

	#### SO, SO, SO, SO. I'll now filter the trace with the response removed! - Cheers Josh
	#### taper it
	trace_taper = trace_detrend.copy()

	# This gets the same results as the SAC function taper.
	trace_taper.taper(max_percentage=taper_percent) ## Defaults to both sides and 'hann' type.

	trace_filt = trace_taper.copy()

	#trace_filt.filter('bandpass', freqmin=0.05, freqmax=0.3, corners=4, zerophase=True)

	trace_filt.filter('bandpass', freqmin=min_freq, freqmax=max_freq, corners=4, zerophase=True)

	trace_filt.write("%s.%s.SAC" %(station_name,channel_name), format="SAC")

	print("------------Next file -------------")


## Time shift the SAC file
import obspy.io.sac.sactrace

for sac_file in glob('*SAC'):

	sac=obspy.io.sac.sactrace.SACTrace.read(sac_file)
	FILE_NAME = os.path.splitext(sac_file)
	print(FILE_NAME[0])
	sac.reftime= event_time

	print(sac.t1, sac.b, sac.o)

	tr2 = sac.to_obspy_trace()
	print(tr2.stats)
	tr2.write("%s.SAC" %FILE_NAME[0], format="SAC")




print("YAY - It worked, I now have SAC files... WITH POPULATED INFORMATION")

print("I hope this is living up to your expectations")
print("If not, you can leave and write your own script")
print("Then send me it.")
print("It'll be better and more efficient than mine...")

## make a list with only no repeated station name

stations = list(set(station_temp))

print(stations)

Q_rotate = True
#### ROTATE! ####
print("rotate:", Q_rotate)
if Q_rotate == True:

	## Now we rotate ## if the user wants to ##

	## use the lists made at the start of the code

	for station in stations:
	## print station, because, why not?
		print(station)
		## make temporary list for the sac files
		temp_list = []

		## loop over all the sac files, because efficiency is for losers
		for sac_file in glob('*SAC'):
		## print the name, because, see above comment about efficiency and losers
			## take the BHN and BHE files and shove them in the temporary list
			if "%s." %station in sac_file and 'BHN' in sac_file:
				temp_list.append(sac_file)
			if "%s." %station in sac_file and 'BHE' in sac_file:
				temp_list.append(sac_file)
			#print(temp_list)

	## make a stream of those two lists only
		stream = obspy.read(temp_list[0])
		stream += obspy.read(temp_list[1])

		## find back azimuth and inclination
		baz = stream[0].stats.sac.baz
		inc = stream[0].stats.sac.cmpinc

	### check if they are the same length, if not, trim the longer one by the difference in data points

		npts1 = stream[0].stats.npts
		npts2 = stream[1].stats.npts
		samp_rate = stream[0].stats.delta

		s_time1 = stream[0].stats.starttime
		s_time2 = stream[1].stats.starttime

		e_time1 = stream[0].stats.endtime
		e_time2 = stream[1].stats.endtime

		# print(stream)

		## correct for different trace start times

		if s_time1 != s_time2:
			if s_time1 > s_time2:
				stream[0].trim(starttime=s_time1)
				stream[1].trim(starttime=s_time1)
			if s_time1 < s_time2:
				stream[0].trim(starttime=s_time2)
				stream[1].trim(starttime=s_time2)
		if e_time1 != e_time2:
			if e_time1 > e_time2:
				stream[0].trim(endtime=e_time2)
				stream[1].trim(endtime=e_time2)
			if e_time1 < e_time2:
				stream[0].trim(endtime=e_time1)
				stream[1].trim(endtime=e_time1)


		stream.rotate(method = 'NE->RT', back_azimuth = baz, inclination = inc)


	### Get the traces from the stream ###
		trace1 = stream[0]
		trace2 = stream[1]

		s_time3 = trace1.stats.starttime
		s_time4 = trace2.stats.starttime

		e_time3 = trace1.stats.endtime
		e_time4 = trace2.stats.endtime

		## correct for different trace start times

		if s_time3 != s_time4:
			if s_time3 > s_time4:
				trace1.trim(starttime=s_time3)
				trace2.trim(starttime=s_time3)
			if s_time3 < s_time4:
				trace1.trim(starttime=s_time4)
				trace2.trim(starttime=s_time4)
		if e_time3 != e_time4:
			if e_time3 > e_time4:
				trace1.trim(endtime=e_time4)
				trace2.trim(endtime=e_time4)
			if e_time1 < e_time2:
				trace1.trim(endtime=e_time3)
				trace2.trim(endtime=e_time3)

		print(stream)

		ch1 = trace1.stats.channel
		ch2 = trace2.stats.channel

		print(ch1,ch2)

		if ch1 == "BHR":
			cmpaz_1 = float(trace1.stats.sac.baz) + 180
		elif ch1 == "BHT":
			cmpaz_1 = float(trace1.stats.sac.baz) + 270
		else:
			pass

		if ch2 == "BHR":
			cmpaz_2 = float(trace1.stats.sac.baz) + 180
		elif ch2 == "BHT":
			cmpaz_2 = float(trace1.stats.sac.baz) + 270
		else:
			pass

		### write component information into the file
		trace1.stats.sac['cmpaz'] = float(cmpaz_1)
		trace2.stats.sac['cmpaz'] = float(cmpaz_2)

		print("CMPAZ: ", trace1.stats.sac.cmpaz)
	###### get info from the name of the first trace #####
		station_name,channel_name,SAC = temp_list[0].split(".")

	## write the dude ##

		trace1.write("%s.%s.SAC" %(station_name,str(trace1.stats.channel)), format = "SAC")
		trace2.write("%s.%s.SAC" %(station_name,str(trace2.stats.channel)), format = "SAC")

	## Store the BHN and BHE in case I want them for whatever reason.
	# make directory
	H_E_Dir="../processed/BHN_BHE_processed"
	if not os.path.exists(H_E_Dir):
	    os.makedirs(H_E_Dir)
	# move files
	for old_file in glob('*BHN*.SAC'):
		if os.path.exists(str(H_E_Dir) + "/" + str(old_file)):
			os.remove(str(H_E_Dir) + "/" + str(old_file))
		shutil.move(old_file,H_E_Dir)

	for old_file in glob('*BHE*.SAC'):
		if os.path.exists(str(H_E_Dir) + "/" + str(old_file)):
			os.remove(str(H_E_Dir) + "/" + str(old_file))
			shutil.move(old_file,H_E_Dir)


## should be left with only BHR/BHT/BHZ files

### We are, in fact, not done...
## need to move the files to the "processed" directory

if Q_splitting:
### Going to do my splitting work within a seperate directory
	Splitting_Dir = "../Splitting"
	if not os.path.exists(Splitting_Dir):
	    os.makedirs(Splitting_Dir)
else:
	pass

for pro_file in glob('*SAC'):
	pro_file_First = pro_file.split(".")[0]
	pro_file_Second = pro_file.split(".")[1]
	pro_file_no_SAC = pro_file_First + "." + pro_file_Second
	if os.path.exists("../processed/%s" %pro_file):
		os.remove("../processed/%s" %pro_file)
	shutil.copy(pro_file, "../processed")

	if Q_splitting:
		if os.path.exists("%s/%s" %(Splitting_Dir,pro_file)):
			os.remove("%s/%s" %(Splitting_Dir,pro_file))
		shutil.move(pro_file, Splitting_Dir)
		os.rename("%s/%s" %(Splitting_Dir, pro_file), "%s/%s" %(Splitting_Dir,pro_file_no_SAC))
	else:
		pass
