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




### Define times before and after predicted arrival:
## Or could do it between the predicted times of SKS and SKKS
win_st = -40
win_end = 30

## make directories
if not os.path.exists("Keep_Auto"):
	os.makedirs("Keep_Auto")

if not os.path.exists("Trash_Auto"):
	os.makedirs("Trash_Auto")

if not os.path.exists("Hmmm_Maybe_Auto"):
	os.makedirs("Hmmm_Maybe_Auto")

if not os.path.exists("Confused_Auto"):
	os.makedirs("Confused_Auto")





## make a file to record the sorting process for each



### Read in all BHR SAC files

st = obspy.read('*B?R*SAC')

keep_count =  0
trash_count = 0
maybe_count = 0
confused_count = 0


for i,tr in enumerate(st):

	### normalising seems like a good idea

	try:
		tr = tr.normalize()
	except ValueError:
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Trash_Auto")
		## Continue should go back to the top of the for loop
		continue
	## Get the predicted SKS time.

	SKS_pred = tr.stats.sac.t1 - tr.stats.sac.b
	SKKS_pred = tr.stats.sac.t2 - tr.stats.sac.b
	print(SKS_pred)
	s_time = SKS_pred + win_st
	e_time = SKS_pred + win_end
	print(s_time, e_time)
	window_length = e_time - s_time
	print(window_length)
	## get the number of data points in the window
	pts_in_wind = window_length*tr.stats.sampling_rate
	print(pts_in_wind)
	# define start and end data point
	# Seems to throw up an error when there is no data...typical...
	try:
		start_sample = int((s_time) * tr.stats.sampling_rate + .5)
		print(start_sample)
		print(len(tr.data))
		end_sample = int((tr.stats.endtime - (tr.stats.starttime + e_time)) * tr.stats.sampling_rate + .5)
	except ValueError:
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Trash_Auto")
		## Continue should go back to the top of the for loop
		continue


	## get the data points for the window only
	dat = tr.data[int(start_sample):(int(start_sample) + int(pts_in_wind))]
	print(tr.data)


	### OKOKOK, we now have the data, now I'll try a few methods to automate the picking process (hopefully it'll be similar to the hand picked sorting).

	### Method 1 - signal-noise of absolute amplitudes


	try:## Get the average amplitude
		avg_amp = np.mean(abs(dat))
	except ValueError:
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Trash_Auto")
		continue
	## Get the max amplitude
	try:
		max_amp = np.max(abs(dat))
		sig_noise_ratio = max_amp/avg_amp
	except ValueError:
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Trash_Auto")
		continue



	print(avg_amp, max_amp, sig_noise_ratio)

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

	hilb_dat = hilbert(dat)
	env_dat = (dat ** 2 + hilb_dat ** 2) ** 0.5
	avg_amp_env = np.mean(abs(env_dat))

	## Get the max amplitude
	max_amp_env = np.max(abs(env_dat))

### Let's try and add a few if statements to sort the traces. This is reasonably generous...

### First get rid of any traces which have an epicentral distance of less than 80.

#	fig = plt.figure()
#	fig.add_subplot(311)
#	fig.add_subplot(211)
#	plt.axhline(0, linestyle=':', label="Zero Reference")
#	plt.plot(abs(dat), color="black")
#	plt.axhline(peak_avg, color='b', label="Mean of Peaks of Absolute Amplitude")
#	plt.axhline(max_amp, color='r', label="Max Amplitude")
#	plt.title("Station: %s Signal-Noise Ratio: %s" %(tr.stats.station, peak_ratio))
#	plt.legend(loc='best')
#	plt.xlabel("Time (s)")
#	plt.ylabel("Displacement (m)")
#	plt.show()

	dist = tr.stats.sac.gcarc



	if dist < 80:
		print("%s should be in trash" %tr.stats.station)
		trash_count += 1
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Trash_Auto")


	elif peak_ratio >3:
		print("%s in Keep_Auto" %tr.stats.station)
		keep_count += 1
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Keep_Auto")

	elif peak_ratio >2.5 and peak_ratio < 3:
		print("%s in Hmmm_Maybe_Auto" %tr.stats.station)
		maybe_count += 1
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Hmmm_Maybe_Auto")

	elif  peak_ratio < 2.5:
		print("%s in Trash_Auto" %tr.stats.station)
		trash_count +=1
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Trash_Auto")

	else:
		print("%s is weird" %tr.stats.station)
		confused_count += 1
		for sac_file in glob('*%s*B?R*' %tr.stats.station):
			shutil.move(sac_file, "Confused_Auto")
	fig = plt.figure()
	fig.add_subplot(311)
	fig.add_subplot(211)



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
