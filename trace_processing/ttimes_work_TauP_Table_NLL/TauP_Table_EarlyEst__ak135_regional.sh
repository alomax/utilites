echo "============"
echo "Script to create spherical, layered-model, travel-time grids for EarlyEst - Global mode"
echo "see http://early-est.alomax.net"
echo
echo "IMPORTANT:  Requires:"
echo "   1. Java - the command \"java\" must be on your path"
echo "   2. TauP Toolkit classes must be on your java classpath- see: http://www.seis.sc.edu/software/TauP"
echo "   3. The additional class edu.sc.seis.TauP.TauP_Table_NLL must be on your java classpath"
echo "     - this class is available at http://alomax.net/nlloc/java"

echo ============ regional

# EDIT HERE !! ##############################################
MODEL=ak135

echo MODEL: ${MODEL}
OUTPATH=./temp
echo OUTPATH: ${OUTPATH}

mkdir ${OUTPATH}

# EDIT HERE !! ##############################################
OUTROOT=${OUTPATH}/ttimes_${MODEL}_regional

java edu.sc.seis.TauP.TauP_Create -tvel ${MODEL}.nd

# NLL_GRID=<num cells in distance>,<distance min in deg>,<distance max in deg>,<num cells in depth>,<depth min in km>,<depth max in km>
#  Set the parameters with regards to the area (cylindrical slice of the spherical earth) that you want to study.  The default is for the whole earth.  The grid is 2D, so it is not critical if you make it larger than needed.
#  IMPORTANT: NUM_DEPTH_INDEX is NUM_DEPTH + 2 (distance, values available flag).
#NLL_GRID=1801,0.0,180.0,281,0.0,700.0
# EDIT HERE !! ##############################################
NUM_DIST=2001
NUM_DEPTH=81
NUM_DEPTH_INDEX=83
NLL_GRID=${NUM_DIST},0.0,20.0,${NUM_DEPTH},0.0,80.0

echo $NLL_GRID


# EDIT HERE !! ##############################################
NUM_PHASES=2

cat  << END > ${OUTROOT}_times_phases.h
static double dist_time[${NUM_PHASES}][${NUM_DIST}][${NUM_DEPTH_INDEX}] = {
END
cat  << END > ${OUTROOT}_toang_phases.h
static double dist_toang[${NUM_PHASES}][${NUM_DIST}][${NUM_DEPTH_INDEX}] = {
END

echo First arriving P and S

java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph P,p,Pn,Pdiff,PKP,PKiKP,PKIKP -o ${OUTROOT}
java edu.sc.seis.TauP.TauP_Table_NLL -earlyest_angles ${NLL_GRID}  -mod ${MODEL} -ph P,p,Pn,Pdiff,PKP,PKiKP,PKIKP -o ${OUTROOT}
cat  "${OUTROOT}_times_phases.h_part_P" >> ${OUTROOT}_times_phases.h
cat  "${OUTROOT}_toang_phases.h_part_P" >> ${OUTROOT}_toang_phases.h
java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph S,s,Sn,Sdiff,SKS,SKiKS,SKIKS -o ${OUTROOT}
java edu.sc.seis.TauP.TauP_Table_NLL -earlyest_angles ${NLL_GRID}  -mod ${MODEL} -ph S,s,Sn,Sdiff,SKS,SKiKS,SKIKS -o ${OUTROOT}
cat  "${OUTROOT}_times_phases.h_part_S" >> ${OUTROOT}_times_phases.h
cat  "${OUTROOT}_toang_phases.h_part_S" >> ${OUTROOT}_toang_phases.h

cat  << END >> ${OUTROOT}_times_phases.h
};
END
cat  << END >> ${OUTROOT}_toang_phases.h
};
END




echo Pg and Sg

#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph Pg -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph Sg -o ${OUTROOT}


echo pP and sP

#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sP -o ${OUTROOT}


echo PcP and ScP and PcS and ScS, etc

#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph PcP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pPcP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sPcP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph ScP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pScP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sScP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph PcS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pPcS  -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sPcS  -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph ScS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pScS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sScS -o ${OUTROOT}


echo PSP, SKS, SKP and PKS, etc

#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph PKP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph SKS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pSKP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sSKP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph PKS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pPKS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sPKS -o ${OUTROOT}


echo PP and PS and SP and SS

#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph PP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pPP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sPP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph PS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pPS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sPS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph SP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pSP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sSP -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph SS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph pSS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph sSS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph SSS -o ${OUTROOT}
#java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph SSSS -o ${OUTROOT}

echo 410, 660, etc

### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph p410S -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph P410S -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph Pv410S -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph s410P -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph S410P -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph Sv410P -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph p660S -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph P660S -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph Pv660S -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph s660P -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph S660P -o ${OUTROOT}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -earlyest ${NLL_GRID}  -mod ${MODEL} -ph Sv660P -o ${OUTROOT}


echo "Finished!"
echo "cp -p ${OUTROOT}*.h ../timedomain_processing/ttimes/"

