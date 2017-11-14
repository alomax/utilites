
# =========================================================================================================================
# options for applications
# =========================================================================================================================

# =========================================================================================================================
# conditional compiling of different velocity models / travel-time grids / phase types, etc.
#
# global scale AK135 model, time grid and phases generated with iaspei-tau FORTRAN (ttimes_ak135_0-800_10*)
TTIME_DEFS=-D TTIMES_AK135_GLOBAL
export TTIMES_NAME=
#
# global scale AK135 model, time grid and phases generated with edu.sc.seis.TauP Java (ttimes_ak135_0-800_10_TauP*)
#TTIME_DEFS=-D TTIMES_AK135_GLOBAL_TAUP
#export TTIMES_NAME=_AK135_GLOBAL_TAUP
#
# regional scale AK135 model, time grid and phases (ttimes_ak135_regional*)
#TTIME_DEFS=-D TTIMES_AK135_REGIONAL -D TTIMES_REGIONAL
#export TTIMES_NAME=_AK135_REGIONAL
#
# regional scale marsite1 model, time grid and phases (ttimes_marsite1_regional*)
# Estelle CROS e-mail 20130419
#TTIME_DEFS=-D TTIMES_MARSITE1_REGIONAL -D TTIMES_REGIONAL
#export TTIMES_NAME=_MARSITE1_REGIONAL
# AJL adapted with Vp/Vs=1.73 from: 2008__Becel_et__Moho_crustal_architecture_and_deep_deformation_under_the_North_Marmara_Trough_from_the_SEISMARMARA_Leg_1_offshore–onshore_reflection–refraction_survey__Tectonopyhsics
#TTIME_DEFS=-D TTIMES_MARSITE2_REGIONAL -D TTIMES_REGIONAL
#export TTIMES_NAME=_MARSITE2_REGIONAL
#
# regional scale 1D Velocity Model Used by the SCSN (Hadley and Kanamori, 1977; K Hutton et. al, BSSA, 2010)
# HK Model with Vp/Vs=1.73
# Depth to Top of Layer (km), P-Velocity (km/sec)
# 0.0 5.5
# 5.5 6.3
# 16.0 6.7
# 32.0 7.8
#TTIME_DEFS=-D TTIMES_SOCALSN_HK_2010 -D TTIMES_REGIONAL
#export TTIMES_NAME=_SOCALSN_HK_2010
#
# Longmenshan (2008/05/12 Wenchuan, China)
# Zheng Y, Ma H S, Lü J, et al. Source mechanism of strong aftershocks (Ms≥5.6) of the 2008/05/12 Wenchuan earthquake and the implication for seismotectonics. Sci China Ser D-Earth Sci, 2009, 52(6): 739-753, doi: 10.1007/s11430-009-0074-3
#TTIME_DEFS=-D TTIMES_LONGMENSHAN_REGIONAL -D TTIMES_REGIONAL
#export TTIMES_NAME=_LONGMENSHAN_REGIONAL
#
#
export EXE_NAME=$(TTIMES_NAME)
# =========================================================================================================================

# =========================================================================================================================
# AJL_TEST : development / test functions - WARNING: not fully implemented, debugged, tested or memory checked!
# You can comment out the following line to disable all development / test functions
#
DEV_DEFS= $(DEV_DEFS_0) $(DEV_DEFS_1) $(DEV_DEFS_2) $(DEV_DEFS_3) $(DEV_DEFS_4) $(DEV_DEFS_5) $(DEV_DEFS_6) $(DEV_DEFS_7) $(DEV_DEFS_8) $(DEV_DEFS_9)
#
#DEV_DEFS_0= -D COMPILE_OLD_VERSION # do not use!
#
#export EXE_NAME=$(TTIMES_NAME)_ALPHA   # X 20150415
#DEV_DEFS_0= -D ALPHA_VERSION   # X 20150415
#
# 20141203 AJL - Uncomment the following line to enable up-weighting of picks with high S/N for location
#
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_5= -D USE_LOCATION_UPWEIGHT_HIGH_SN_PICKS   # X 20150415
#
# 20150108 AJL - Uncomment the following line to enable use of station corrections
#
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_7= -D USE_STATION_CORRECTIONS   # X 20150415
#
# 20150401 AJL - Uncomment the following line to﻿use all picks for location, irrespective of HF S/N.
#
DEV_DEFS_8= -D USE_USE_ALL_PICKS_FOR_LOCATION   # X 20150415
#
#
# 20140122 AJL - Uncomment the following line to enable check of P amplitude attenuation in association/location
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_1= -D USE_AMPLITUDE_ATTENUATION
#
# 20140624 AJL - Uncomment the following line to check for event persistence - keep previous location results for events with insufficient potential new defining picks
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_2= -D USE_EVENT_PERSISTENCE
#
# 20140801 AJL - Uncomment the following line to﻿correct the computed Mwp and Mwpd based on epicentral distance.
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_4= -D USE_MWP_DISTANCE_CORRECTION
#
# 20141208 AJL - Uncomment the following line to enable checking if amplitude attenuation error ratio too high or too low - gain may be incorrect, ignore magnitudes
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_6= -D USE_MAGNITUDE_AMP_ATTEN_CHECK
#
# 20150408 AJL - Uncomment the following line to﻿correct mB to better agree with mb_NEIC.
# ADDED TO PROP FILE, CHECK WORKS  DEV_DEFS_9= -D USE_MB_CORRECTION
#
# =========================================================================================================================



# Options specific for GCC
export CC = gcc
#
# GNU
GNU_SOURCE=-D _GNU_SOURCE
# IMPORTANT: try the following if the compiler you are using does not support GNU extensions!
#GNU_SOURCE=
#
#
CCFLAGS_BASIC =  -std=c99 -Wall $(GNU_SOURCE) `xml2-config --cflags` $(TTIME_DEFS) $(DEV_DEFS)
# 20160513 - -mcmodel=large flag may be needed on Linux per suggestion by mehmety@boun.edu.tr
#CCFLAGS_BASIC = -mcmodel=large -std=c99 -Wall $(GNU_SOURCE) `xml2-config --cflags` $(TTIME_DEFS) $(DEV_DEFS)

#
# all warnings
#export CCFLAGS_BASIC = -Wextra $(CCFLAGS_BASIC)
#
#
# basic
#export CCFLAGS = $(CCFLAGS_BASIC)
#
# optimized
export CCFLAGS = -O3 $(CCFLAGS_BASIC)
#
# profile
#export CCFLAGS= -pg -p $(CCFLAGS_BASIC)
#
# debug - gdb, valgrind, ...
#export CCFLAGS = -g $(CCFLAGS_BASIC)
#export CCFLAGS = -g -O3 $(CCFLAGS_BASIC) -g -O3
# gdb:
# gdb --args exe_name
# gdb-apple --args exe_name
# lldb -- exe_name
# valgrind:
# !!! Note: 20121213 AJL - valgrind only works if remove -c ./plot_warning_report_seedlink_runtime.bash from seedlink_monitor command line!
#           20141126 AJL - Above may not be necessary with Valgrind-3.10.0
# valgrind --leak-check=yes --dsymutil=yes exe_name <args>
# w/o -O3 : valgrind --leak-check=yes --dsymutil=yes --track-origins=yes exe_name <args>
# valgrind --leak-check=full --show-reachable=yes --dsymutil=yes exe_name <args>
# valgrind --tool=callgrind --dsymutil=yes exe_name <args>  # callgrind_annotate [options] --tree=both --inclusive=yes callgrind.out.<pid>



export LIB_OBJS = ../timedomain_processing/timedomain_processing.o ../timedomain_processing/timedomain_processing_data.o  ../timedomain_processing/timedomain_processing_report.o \
../timedomain_processing/timedomain_processes.o ../timedomain_processing/timedomain_filter.o ../timedomain_processing/timedomain_memory.o \
../timedomain_processing/ttimes.o ../timedomain_processing/location.o ../timedomain_processing/loc2xml.o \
../vector/vector.o \
../statistics/statistics.o \
../alomax_matrix/alomax_matrix.o ../alomax_matrix/alomax_matrix_svd.o ../alomax_matrix/polarization.o ../alomax_matrix/eigv.o ../matrix_statistics/matrix_statistics.o \
../picker/PickData.o ../picker/FilterPicker5_Memory.o  ../picker/FilterPicker5.o \
../net/net_lib.o ../response/response_lib.o ../math/complex.o ../octtree/octtree.o ../ran1/ran1.o ../feregion/feregion.o \
../settings/settings.o ../settings/strmap.o ../miniseed_utils/miniseed_utils.o


MODULES = \
	math \
	vector \
	statistics \
	alomax_matrix \
	matrix_statistics \
	picker \
	net \
	response \
	ran1 \
	octtree \
	feregion \
	settings \
	miniseed_utils \
	timedomain_processing \
	fmamp \
	applications

all:
	@for x in $(MODULES); \
	do \
		(echo ------; cd $$x; echo Making $@ in:; pwd; \
		make -f Makefile); \
	done
	@echo
	@echo "Notes:"
	@echo "   Build requires Libxml2 (XML C parser and toolkit) available from http://xmlsoft.org/ and from http://www.macports.org/"
	@echo "   Velocity model / travel-time grid: $(TTIME_DEFS)"
	@echo

clean:
	@for x in $(MODULES); \
	do \
		(cd $$x; echo Cleaning in:; pwd; \
		make -f Makefile clean); \
	done
