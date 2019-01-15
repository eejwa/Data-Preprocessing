#!/usr/bin/bash

### ObspyDMT is an easy to use tool for getting data in whatever format you want. I recommend it to anyone and all my codes are set up using the directory layout. 

### Github page:
# https://github.com/kasra-hosseini/obspyDMT

### the two below are examples that I actually used, they care copied from the log file created by Obspy_DMT!

### get a donut of events around lon:27 lat:-28 between distances 35 and 75 degrees

obspyDMT --datapath DMT_Data_P_S_REVIEWED --min_mag 5.5 --max_mag 7.5 --event_circle 27/-28/35/75 --min_date 1997-01-01 --max_date 2000-01-01 --event_catalog ISC --isc_catalog REVIEWED --data_source IRIS --net XA --cha BH* --preset 0 --offset 1800 --req_parallel --req_np 10 

### get events within a defined rectangle
obspyDMT --datapath ./DMT_data --min_mag 5.5 --max_mag 7.5 --event_rect -130./-60./-60./5 --min_date 1997-01-01 --max_date 2000-01-01 --event_catalog ISC --isc_catalog REVIEWED --data_source IRIS --net XA --loc * --cha BH* --preset 0 --offset 1800 --req_parallel --req_np 10 


