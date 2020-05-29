#!/usr/bin/env python

usage = """Code to filter traces in the directory given a frequency band and file wildcard.

[-fl][-fh][-f][-t] where:

-fl = miniumum frequency value (e.g. 0.05)
-fh = maxiumum frequency value (e.g. 1.0)
-f = filename wildcard (e.g. '*SAC')
-t = type of filtering (e.g. bandpass)
"""

import obspy
import numpy as np
import argparse

# get the arguments from the terminal
parser = argparse.ArgumentParser(description='Preprocessing script for data retrieved from obspy DMT')

parser.add_argument("-f","--file_wildcard", help="Enter the file to be normalised (e.g. *BHR*)", type=str, required=True, action="store", default = "*SAC")

parser.add_argument("-fl","--lower_frequency", help="Enter the lower frequency you want the analysis to be conducted over", type=float, required=True, action="store", default = "0.1")

parser.add_argument("-fh","--upper_frequency", help="Enter the upper frequency you want the analysis to be conducted over", type=float, required=True, action="store", default = "0.4")

parser.add_argument("-t","--filter_type", help="Enter the type of filtering you want to do", type=str, required=True, action="store", default = "bandpass")


args = parser.parse_args()

file_names = args.file_wildcard
flow=args.lower_frequency
fhigh=args.upper_frequency
type=args.filter_type

st = obspy.read(file_names)

# filter the stream

st_filtered = st.filter(type, freqmin=flow, freqmax=fhigh)

for i,tr in enumerate(st_filtered):
    # get information about the trace and rename it
    network=tr.stats.network
    station=tr.stats.station

    tr.write("%s_%s_filtered.SAC" %(network,station), format="SAC")
