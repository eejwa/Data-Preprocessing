#!/usr/bin/env python

# I'm going to try and develop and develop an algorithm (not very complex one mind) to automate the sorting process...

import obspy
import obspy.signal
import numpy as np
import matplotlib.pyplot as plt
from glob import *
from scipy.fftpack import hilbert
import os
import shutil
from glob import *

### For finding the peaks in the arrays
from scipy import signal

import argparse

parser = argparse.ArgumentParser(description='Automatically sort SAC files based on their signal-noise ratio')

#parser.add_argument("-t","--time_header", help="Enter the time header for the phase you are interested in (e.g. t1).", type=str, nargs=1, required=True, action="store")

parser.add_argument("-d","--dir_wildcard", help="Enter a string to define the directories you want to loop over (e.g. 20*).", type=str, required=True, action="store")

parser.add_argument("-f","--file_wildcard", help="Enter a string to define the files you want to loop over (e.g. *BHR*SAC).", type=str, required=True, action="store")

parser.add_argument("-p","--phase", help="Enter The name of the phase you are analysing (e.g. SKS).", type=str, required=True, action="store")

parser.add_argument("-v", "--verbose", help="Increase verbosity of output, just --verbose is enough to turn on.",action="store_true")

parser.add_argument("-s", "--window_start", help="time in seconds before the predicted arrival for the window of interest. default of -40.",action="store", default=-40, type=float)

parser.add_argument("-e", "--window_end", help="time in seconds after the predicted arrival for the window of interest. default of 30.",action="store", default=30, type=float)

args = parser.parse_args()

#SAC_time_header = args.time_header[0]
dir_def = args.dir_wildcard
file_def = args.file_wildcard
phase = args.phase
win_st = args.window_start
win_end = args.window_end

print(win_st)
print(win_end)


# Open file and write header:
record_file = open("Sorting_Record_%s.txt" %phase, 'w')
record_file.write("Event_Name | Keep | Maybe | Trash \n")

## State name of folder with processed data files in it. s
processed_folder = "processed/"

## full path of current directory
pwd = os.getcwd()



for fold in glob(dir_def):
	relpath = str(fold) + "/" + "processed/"
	os.chdir(relpath)
	cdr = os.getcwd()
	print(cdr)

	## make directories
	if not os.path.exists("Keep_Auto_%s" %phase):
		os.makedirs("Keep_Auto_%s" %phase)

	if not os.path.exists("Trash_Auto_%s" %phase):
		os.makedirs("Trash_Auto_%s" %phase)

	if not os.path.exists("Hmmm_Maybe_Auto_%s" %phase):
		os.makedirs("Hmmm_Maybe_Auto_%s" %phase)

	if not os.path.exists("Confused_Auto_%s" %phase):
		os.makedirs("Confused_Auto_%s" %phase)

	### Read in all files - they should be SAC files.
	st = obspy.read(file_def)

	## make counters
	keep_count =  0
	trash_count = 0
	maybe_count = 0
	confused_count = 0


	for i,tr in enumerate(st):
		### normalising seems like a good idea
		print("The next trace has station: ",tr.stats.station)
		try:
			tr = tr.normalize()
		except ValueError:
			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
				shutil.copy2(sac_file, "Trash_Auto_%s" %phase)
				trash_count += 1
			## Continue should go back to the top of the for loop
			continue

		## Get the time header associated with the phase you want:
		# Need to loop over the file headers and see which has the right phase
		for a in range(1,9,1):
			label="kt"+str(a)
			print(label)
			try:
				phase_label = getattr(tr.stats.sac, label)
			except:
				pass
			if phase_label == phase:
				time_header="t"+str(a)
				break

		## Get the predicted time for the phase you want.
		pred = getattr(tr.stats.sac,time_header) - tr.stats.sac.b

		print(pred)
		s_time = pred + win_st
		e_time = pred + win_end
		print(s_time, e_time)
		window_length = e_time - s_time
		print(window_length)
		## get the number of data points in the window
		pts_in_wind = window_length*tr.stats.sampling_rate
		print(pts_in_wind)
		# define start and end data point
		# Seems to throw up an error when there is no data... typical...
		try:
			start_sample = int((s_time) * tr.stats.sampling_rate + .5)
			print(start_sample)
			print(len(tr.data))
			end_sample = int((tr.stats.endtime - (tr.stats.starttime + e_time)) * tr.stats.sampling_rate + .5)
		except ValueError:
			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
				shutil.copy2(sac_file, "Trash_Auto_%s" %phase)
				trash_count+=1
			## Continue should go back to the top of the for loop
			continue


		## get the data points for the window onlyshould
		dat = tr.data[int(start_sample):(int(start_sample) + int(pts_in_wind))]


		### OKOKOK, we now have the data, now I'll try a few methods to automate the picking process (hopefully it'll be similar to the hand picked sorting).

		### Method 1 - signal-noise of absolute amplitudes
#		try:## Get the average amplitude
#			avg_amp = np.mean(abs(dat))
#		except ValueError:
#			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
#				shutil.copy2(sac_file, "Trash_Auto_%s" %phase)
#				trash_count+=1
#			continue
		## Get the max amplitude
		try:
			max_amp = np.max(abs(dat))
#			sig_noise_ratio = max_amp/avg_amp
		except ValueError:
#			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
#				shutil.copy2(sac_file, "Trash_Auto_%s" %phase)
#				trash_count+=1
			continue



#		print(avg_amp, max_amp, sig_noise_ratio)

		#### Method 2 - Same as above except the mean is with the peaks of the waves, not the whole trace.

		## reset the peaks list
		peaks = []
		peaks_loc = signal.argrelmax(abs(dat))
		for loc in peaks_loc:
			peaks.append(abs(dat[loc]))
		peaks = np.array(peaks)
		peak_avg = np.mean(peaks)


		peak_ratio = max_amp/peak_avg

		print(peak_avg, max_amp,peak_ratio)

		### Method 3 - seismic envelope area or amplitude

#		hilb_dat = hilbert(dat)
#		env_dat = (dat ** 2 + hilb_dat ** 2) ** 0.5
#		avg_amp_env = np.mean(abs(env_dat))

		## Get the max amplitude
#		max_amp_env = np.max(abs(env_dat))

	### Let's try and add a few if statements to sort the traces. This is reasonably generous...

	### First get rid of any traces which have an epicentral distance of less than 80.

	#	fig = plt.figure()
	#	fig.add_subplot(311)
	#	fig.add_subplot(211)
	#	plt.plot(dat)
	#	plt.axhline(peak_avg, color='b', label="Mean")
	#	plt.axhline(max_amp, color='r', label="Max")
	#	plt.title(tr.stats.station)
	#	plt.legend(loc='best')

		dist = tr.stats.sac.gcarc



		#if dist < 80:
		#	print("%s should be in trash" %tr.stats.station)
		#	trash_count += 1
		#	for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
		#		shutil.copy2(sac_file, "Trash_Auto_%s" %phase)


		if peak_ratio >3:
			print("%s in Keep_Auto_%s" %(tr.stats.station, phase))
			keep_count += 1
			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
				shutil.copy2(sac_file, "Keep_Auto_%s" %phase)

		elif peak_ratio >2.5 and peak_ratio < 3:
			print("%s in Hmmm_Maybe_Auto_%s" %(tr.stats.station, phase))
			maybe_count += 1
			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
				shutil.copy2(sac_file, "Hmmm_Maybe_Auto_%s" %phase)

		elif  peak_ratio < 2.5:
			print("%s in Trash_Auto_%s" %(tr.stats.station, phase))
			trash_count +=1
			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
				shutil.copy2(sac_file, "Trash_Auto_%s" %phase)

		else:
			print("%s is weird" %tr.stats.station)
			confused_count += 1
			for sac_file in glob('*%s*%s' %(tr.stats.station, file_def)):
				shutil.copy2(sac_file, "Confused_Auto_%s" %phase)

	#	fig = plt.figure()
	#	fig.add_subplot(311)
	#	fig.add_subplot(211)



	#	plt.plot(dat)
	#	plt.axhline(avg_amp, color='b', label="Mean")
	#	plt.axhline(max_amp, color='r', label="Max")
	#	plt.axhline(peak_avg, color='g', label="Peak_Mean")
	#	plt.title(tr.stats.station)
	#	plt.legend(loc='upper right')

	#	fig.add_subplot(212)
	#	plt.plot(env_dat)
	#	plt.axhline(avg_amp_env, color='b', label="Mean")
	#	plt.axhline(max_amp_env, color='r', label="Max")
	#	plt.title(tr.stats.station)
	#	plt.legend(loc='upper right')
	#	plt.show()

	print("Keep %f | Maybe %f | Trash %f | Confused %f" %(keep_count, maybe_count,trash_count, confused_count))


	### Write the results to a file!
	record_file.write("%s | %s | %s | %s | %s \n" %(fold, keep_count, maybe_count, trash_count, confused_count))


	os.chdir(pwd)



record_file.close()
