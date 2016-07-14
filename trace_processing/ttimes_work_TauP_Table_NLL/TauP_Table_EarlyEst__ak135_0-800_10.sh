echo "============"
echo "Script to create spherical, layered-model, travel-time grids for EarlyEst - Global mode"
echo "see http://early-est.alomax.net"
echo
echo "IMPORTANT:  Requires:"
echo "   1. Java - the command \"java\" must be on your path"
echo "   2. TauP Toolkit classes must be on your java classpath- see: http://www.seis.sc.edu/software/TauP"
echo "   3. The additional class edu.sc.seis.TauP.TauP_Table_NLL must be on your java classpath"
echo "     - this class is available at http://alomax.net/nlloc/java"

echo ============ global

# EDIT HERE !! ##############################################
MODEL=ak135

echo MODEL: ${MODEL}
OUTPATH=./temp
echo OUTPATH: ${OUTPATH}
rm -r ${OUTPATH}/*
mkdir ${OUTPATH}

# EDIT HERE !! ##############################################
OUTROOT=${OUTPATH}/ttimes_${MODEL}_0-800_10_TauP

java edu.sc.seis.TauP.TauP_Create -nd ${MODEL}.nd

# NLL_GRID=<num cells in distance>,<distance min in deg>,<distance max in deg>,<num cells in depth>,<depth min in km>,<depth max in km>
#  Set the parameters with regards to the area (cylindrical slice of the spherical earth) that you want to study.  The default is for the whole earth.  The grid is 2D, so it is not critical if you make it larger than needed.
#  IMPORTANT: NUM_DEPTH_INDEX is NUM_DEPTH + 2 (distance, values available flag).
#NLL_GRID=1801,0.0,180.0,281,0.0,700.0
# EDIT HERE !! ##############################################
NUM_DIST=1801
NUM_DEPTH=81
NUM_DEPTH_INDEX=83
NLL_GRID=${NUM_DIST},0.0,180.0,${NUM_DEPTH},0.0,800.0

echo $NLL_GRID


# TauP User's Guide: Currently there is no support for PKPab, PKPbc, or PKPdf phase names. They lead to increased algorith- mic complexity that at this point seems unwarranted. Currently, in regions where triplications develop, the triplicated phase will have multiple arrivals at a given distance. So, PKPab and PKPbc are both labeled just PKP while PKPdf is called PKIKP.
# EDIT HERE !! ##############################################
NUM_PHASES=21
P_PHASES="P,p,Pg,Pn"
S_PHASES="S,s,Sg,Sn"
# IMPORTANT! - order here must match ttimes.c
#PHASE_LIST="${P_PHASES} Pdiff Pg  PKiKP   PKPab PKPbc PKPdf   PKKPab PKKPbc PKKPdf   PKPPKPab PKPPKPbc PKPPKPdf  PP ${S_PHASES} Sdiff Sg   SKSac SKSdf SKPdf   pP sP PcP ScP ScS   PKKSab PKKSbc PKKSdf   SKKSac SKKSdf"
PHASE_LIST="${P_PHASES} Pdiff Pg  PKiKP   PKP PKIKP   PKKP   PKPPKP  PP ${S_PHASES} Sdiff Sg   SKS SKP   pP sP PcP ScP ScS   PKKS   SKKS"
# !! NOTE: PKKP gives no times! ???
# !! NOTE: P'P' gives no times, use PKPPKP
echo ${PHASE_LIST}


cat  << END > ${OUTROOT}_times_phases.h
static double dist_time[${NUM_PHASES}][${NUM_DIST}][${NUM_DEPTH_INDEX}] = {
END
cat  << END > ${OUTROOT}_toang_phases.h
static double dist_toang[${NUM_PHASES}][${NUM_DIST}][${NUM_DEPTH_INDEX}] = {
END

for PHASE in ${PHASE_LIST}; do

	echo "Processing ${PHASE}"

	java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph ${PHASE} -o ${OUTROOT}
	java edu.sc.seis.TauP.TauP_Table_NLL -earlyest_angles ${NLL_GRID}  -mod ${MODEL} -ph ${PHASE} -o ${OUTROOT}
	cat  "${OUTROOT}_times_phases.h_part_${PHASE%%,*}" >> ${OUTROOT}_times_phases.h
	cat  "${OUTROOT}_toang_phases.h_part_${PHASE%%,*}" >> ${OUTROOT}_toang_phases.h

done

cat  << END >> ${OUTROOT}_times_phases.h
};
END
cat  << END >> ${OUTROOT}_toang_phases.h
};
END


echo "Finished!"
echo "cp -p ${OUTROOT}*.h ../timedomain_processing/ttimes/"

