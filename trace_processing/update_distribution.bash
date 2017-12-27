#!/bin/bash

BUILD=NO
TEST_MSPROCESS=NO	# !! Requires BUILD!  Must activate to get example output in distribution.
TEST_SEEDLINK=NO	# !! Requires BUILD!

VERSION=1.2.2  # !!!IMPORTANT: must match WARNING_MONITOR_VERSION in timedomain_processing.h
DISTRIBUTION_NAME=early-est-${VERSION}
echo
echo "Update distribution -------------------------------------------"
echo "Version: ${VERSION}     Distribution name: ${DISTRIBUTION_NAME}"
echo "---------------------------------------------------------------"

PLOT_MAP_MECHANISM_TYPE=fmamp_polarity	# sets mechanism type to plot on map, hash or fmamp, see also plot_event_info._GMT5.gmt, plot_warning_map_GMT5.gmt, processEvents.py

echo
echo "Create distribution directory -------------------------------------------"
mkdir ../temp_dir
rm -rf ../temp_dir/*
echo "-------------------------------------------"

echo
echo "Copy files -------------------------------------------"
rm -rf ${DISTRIBUTION_NAME}/
cp -pRL * ../temp_dir/
/bin/mv ../temp_dir/ ${DISTRIBUTION_NAME}/
# clean up a bit
rm ${DISTRIBUTION_NAME}/*.tar.gz
rm ${DISTRIBUTION_NAME}/*.tgz
rm ${DISTRIBUTION_NAME}/Makefile-trace_processing.mk
rm -r ${DISTRIBUTION_NAME}/${DISTRIBUTION_NAME}/
#rm -r ${DISTRIBUTION_NAME}/doc
#rm -r ${DISTRIBUTION_NAME}/logo
rm -r ${DISTRIBUTION_NAME}/nbproject
#rm -r ${DISTRIBUTION_NAME}/archive
rm -r ${DISTRIBUTION_NAME}/seedlink_out
rm -r ${DISTRIBUTION_NAME}/seedlink_plots
#rm -r ${DISTRIBUTION_NAME}/picker/Luca_Elia
#rm -r ${DISTRIBUTION_NAME}/picker/pick_TP4
rm -r ${DISTRIBUTION_NAME}/picker/test_data
rm -r ${DISTRIBUTION_NAME}/timedomain_processing/ttimes/*
cp -p timedomain_processing/ttimes/ttimes_ak135_0-800_10.h ${DISTRIBUTION_NAME}/timedomain_processing/ttimes/
cp -p timedomain_processing/ttimes/ttimes_ak135_0-800_10_times_phases.h ${DISTRIBUTION_NAME}/timedomain_processing/ttimes/
cp -p timedomain_processing/ttimes/ttimes_ak135_0-800_10_toang_phases.h ${DISTRIBUTION_NAME}/timedomain_processing/ttimes/
cp -p timedomain_processing/ttimes/tvel_ak135.h ${DISTRIBUTION_NAME}/timedomain_processing/ttimes/
cp -p timedomain_processing/ttimes/latlon_deep.h ${DISTRIBUTION_NAME}/timedomain_processing/ttimes/
# copy doc files
cp -p /Users/anthony/work/mseed_processing/README ${DISTRIBUTION_NAME}/
cp -p /Users/anthony/work/early-est/doc/${DISTRIBUTION_NAME}_users_guide.pdf ${DISTRIBUTION_NAME}/
# copy program helper files
mkdir ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/tsunami_warning/miniseed_process.prop ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/seedlink_monitor.prop ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/process_events.prop ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/plot_warning_report_seedlink_runtime.bash ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/*_GMT5.gmt ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/PB2002_steps.dat.txt.*.xy ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/*.ras ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/plates.lonlat ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/www/projects/early-est/warning*.html ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/www/projects/early-est/alarm.mp3 ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/www/projects/early-est/*.pdf ${DISTRIBUTION_NAME}/work
# station corrections
mkdir ${DISTRIBUTION_NAME}/work/sta_corr
cp -pr /Users/anthony/work/mseed_processing/sta_corr/* ${DISTRIBUTION_NAME}/work/sta_corr
# copy HASH program
cp -pRL /Users/anthony/opt/HASH_v1.2 ${DISTRIBUTION_NAME}
cp -p /Users/anthony/work/mseed_processing/hash_driver_early_est.bash ${DISTRIBUTION_NAME}/work
# copy fmamp files
cp -p /Users/anthony/work/mseed_processing/fmamp_driver_early_est*.bash ${DISTRIBUTION_NAME}/work
# copy python files
mkdir ${DISTRIBUTION_NAME}/work/python
cp -pr /Users/anthony/work/mseed_processing/python/* ${DISTRIBUTION_NAME}/work/python
rm ${DISTRIBUTION_NAME}/work/python/evaluate/data/underwater/percent_prob_underwater.grd		# huge file!
mkdir ${DISTRIBUTION_NAME}/work/tsunamiLearn_regressors
cp -pr /Users/anthony/work/early-est/tsunamiLearn/out/regressors/* ${DISTRIBUTION_NAME}/work/tsunamiLearn_regressors
# copy example files
mkdir ${DISTRIBUTION_NAME}/work/seedlink_out
mkdir ${DISTRIBUTION_NAME}/work/msprocess_out
mkdir ${DISTRIBUTION_NAME}/work/msprocess_plots
cp -p /Users/anthony/work/mseed_processing/msprocess_data/Honshu_2011_0_90deg.mseed ${DISTRIBUTION_NAME}/work/msprocess_out
cp -p /Users/anthony/work/mseed_processing/*station_coordinates*.csv ${DISTRIBUTION_NAME}/work
cp -p /Users/anthony/work/mseed_processing/*gainfile*.csv ${DISTRIBUTION_NAME}/work

echo "-------------------------------------------"

cd ${DISTRIBUTION_NAME}/

export MYBIN=$(pwd)/work
export PATH=${MYBIN}:${PATH}

# build =======================================================================

if [ ${BUILD} =  YES ]; then
echo
echo "Clean and make ${DISTRIBUTION_NAME} -------------------------------------------"
make clean
make
echo "-------------------------------------------"

echo
echo "Clean and make HASH_v1.2 -------------------------------------------"
cd HASH_v1.2/
make clean
make
cd ..
echo "-------------------------------------------"
fi


# test =======================================================================

if [ ${TEST_MSPROCESS} =  YES ]; then
echo
echo "Testing miniseed_process -------------------------------------------"
unset FTP_USER_PW_HOST
unset FTP_DIR_WARNING
cd work
COMMAND="${MYBIN}/miniseed_process msprocess_out/Honshu_2011_0_90deg.mseed -g miniseed_process_station_coordinates.csv -a -pz gainfile.csv -pz-query-host service.iris.edu -pz-query irisws/resp/1/query -pz-query-type IRIS_WS_RESP -sta-query-host service.iris.edu -sta-query fdsnws/station/1/query -sta-query-type FDSN_WS_STATION -timeseries-query-hosturl http://alomax.net/webtools/sgweb.html?http://service.iris.edu -timeseries-query irisws/timeseries/1/query -timeseries-query-type IRIS_WS_TIMESERIES -timeseries-query-sladdr rtserve.iris.washington.edu:18000 -timeseries-query-hosturl http://alomax.net/webtools/sgweb.html?http://service.iris.edu -timeseries-query irisws/timeseries/1/query -timeseries-query-type IRIS_WS_TIMESERIES -timeseries-query-sladdr 137.227.224.97:18000 -mwp -mb -mwpd -report-delay 1200 -v -alarm -i AU_XMIS"
echo ${COMMAND}
${COMMAND}
# HASH focal mechanism calculation =======================
COMMAND="python python/processEvents.py MECHANISM msprocess_plots/ msprocess_out/Honshu_2011_0_90deg.mseed.out/ msprocess_out/ msprocess_out/Honshu_2011_0_90deg.mseed.out/ Honshu_2011_0_90deg ${PLOT_MAP_MECHANISM_TYPE}"
echo ${COMMAND}
${COMMAND}
# event processing =======================
COMMAND="python python/processEvents.py EVENTS msprocess_plots/ msprocess_out/Honshu_2011_0_90deg.mseed.out/ msprocess_out/ msprocess_out/Honshu_2011_0_90deg.mseed.out/ Honshu_2011_0_90deg ${PLOT_MAP_MECHANISM_TYPE}"
echo ${COMMAND}
${COMMAND}
# main report processing =======================
COMMAND="./plot_warning_report_GMT5.gmt msprocess_plots/ msprocess_out/Honshu_2011_0_90deg.mseed.out/ msprocess_out/Honshu_2011_0_90deg.mseed.out/ Honshu_2011_0_90deg 0.2 ${PLOT_MAP_MECHANISM_TYPE}"
echo ${COMMAND}
${COMMAND}
COMMAND="$PS_VIEWER msprocess_plots/Honshu_2011_0_90deg_Monitor.pdf"
echo ${COMMAND}
${COMMAND}
read -rsp $'Press any key to continue...\n' -n1 key
rm -r msprocess_out/Honshu_2011_0_90deg.mseed.out
cp -p msprocess_plots/Honshu_2011_0_90deg_Monitor.pdf msprocess_plots/Honshu_2011_0_90deg_Monitor_ORIGINAL.pdf
cp -p msprocess_plots/Honshu_2011_0_90deg_Monitor.ps msprocess_plots/Honshu_2011_0_90deg_Monitor_ORIGINAL.ps
#rm -r ./msprocess_plots
rm ./test.picks
cd ..
echo "-------------------------------------------"
fi

if [ ${TEST_SEEDLINK} =  YES ]; then
echo
echo "Testing seedlink_monitor -------------------------------------------"
unset FTP_USER_PW_HOST
unset FTP_DIR_WARNING
cd work
mkdir seedlink_out
mkdir seedlink_out/iris
mkdir seedlink_plots
COMMAND="${MYBIN}/seedlink_monitor \
\
discovery.rm.ingv.it:39962 -S MN_AQU:BHZ,MN_BLY:BHZ,MN_BNI:BHZ,MN_CEL:BHZ,MN_CLTB:BHZ,MN_CUC:BHZ,MN_DIVS:BHZ,MN_IDI:BHZ,MN_PDG:BHZ,MN_RTC:BHZ,MN_TIP:BHZ,MN_TIR:BHZ,MN_TRI:BHZ,MN_TUE:BHZ,MN_VLC:BHZ,MN_VSL:BHZ,MN_VTS:BHZ,MN_WDD:BHZ,IV_PTCC:BHZ,IV_SALO:BHZ,IV_BOB:BHZ,IV_MONC:BHZ,IV_MGAB:BHZ,IV_MRVN:BHZ,IV_RAFF:BHZ \
\
-t -3600 -nt 120 -nd 60 -locs --,00,10,01 \
\
137.227.224.97:18000 -S IU_*:00BHZ,CU_*:BHZ,GT_*:BHZ,IC_*:00BHZ,IM_*:BHZ,US_*:00BHZ,IW_SMCO:BHZ,GS_KAN11:HHZ,GS_OK025:HHZ  \
\
geosbud.ipgp.fr:18000 -S G_ATD:BHZ,G_CAN:BHZ,G_CLF:BHZ,G_CRZF:BHZ,G_DZM:BHZ,G_FDF:BHZ,G_FOMA:BHZ,G_HDC:BHZ,G_INU:BHZ,G_IVI:BHZ,G_KIP:BHZ,G_PAF:BHZ,G_PEL:BHZ,G_PPTF:BHZ,G_RER:BHZ,G_SPB:BHZ,G_SSB:BHZ,G_TAM:BHZ,G_TAOE:BHZ,G_TRIS:BHZ,G_UNM:BHZ \
\
-c ./plot_warning_report_seedlink_runtime.bash -g realtime_station_coordinates.csv -o seedlink_out/iris -a -v -pz gainfile.csv -pz-query-host service.iris.edu -pz-query irisws/resp/1/query -pz-query-type IRIS_WS_RESP -sta-query-host service.iris.edu -sta-query fdsnws/station/1/query -sta-query-type FDSN_WS_STATION -timeseries-query-hosturl http://alomax.net/webtools/sgweb.html?http://service.iris.edu -timeseries-query irisws/timeseries/1/query -timeseries-query-type IRIS_WS_TIMESERIES -timeseries-query-sladdr rtserve.iris.washington.edu:18000 -timeseries-query-hosturl http://alomax.net/webtools/sgweb.html?http://service.iris.edu -timeseries-query irisws/timeseries/1/query -timeseries-query-type IRIS_WS_TIMESERIES -timeseries-query-sladdr 137.227.224.97:18000 -mwp -mb -mwpd -si -alarm  -agency-id my.domain.name"
echo ${COMMAND}
echo "-------------------------------------------"
echo "NOTE!!  IMPORTANT!!  After several 'New report generated' messages, press Ctrl-C twice rapidly to stop program."
echo "-------------------------------------------"
${COMMAND}
COMMAND="$PS_VIEWER seedlink_plots/t50.pdf"
echo ${COMMAND}
${COMMAND}
read -rsp $'Press any key to continue...\n' -n1 key
rm -r ./seedlink_out
rm -r ./seedlink_plots
cd ..
fi



#if [ ${BUILD} =  YES ]; then
echo
echo "Clean -------------------------------------------"
make clean
cd HASH_v1.2/
make clean
cd ..
echo "-------------------------------------------"
#fi


# construct distribution archive =========================================

echo
echo "Make gzipped archive -------------------------------------------"
cd ..
ARCHIVE_NAME=${DISTRIBUTION_NAME}.tgz
tar czf ${ARCHIVE_NAME} ${DISTRIBUTION_NAME}
echo "cp -p ${ARCHIVE_NAME} /Users/anthony/www/projects/early-est/software"
cp -p ${ARCHIVE_NAME} /Users/anthony/www/projects/early-est/software
cp -p /Users/anthony/work/early-est/doc/${DISTRIBUTION_NAME}_users_guide.pdf /Users/anthony/www/projects/early-est/early-est_users_guide.pdf
echo "mv -f ${ARCHIVE_NAME} /Users/anthony/work/early-est"
mv -f ${ARCHIVE_NAME} /Users/anthony/work/early-est
echo "-------------------------------------------"

echo
echo "Remove temporary distribution directories -------------------------------------------"
rm -r ${DISTRIBUTION_NAME}
echo "-------------------------------------------"


echo
echo "Finished updating distribution -------------------------------------------"
echo "Version: ${VERSION}     Distribution name: ${DISTRIBUTION_NAME}"
echo "--------------------------------------------------------------------------"
