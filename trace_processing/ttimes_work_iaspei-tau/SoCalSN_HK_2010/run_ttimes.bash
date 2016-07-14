IASPEI_TAU_PATH=/Users/anthony/opt/iaspei-tau

MODEL_GRID_NAME=SoCalSN_HK_2010

DEPTH="0.0 80.0 1.0"
DISTANCES="0.0 20.0 0.01"
# IMPORTANT! - array sizes below (dist_time[], dist_toang[]) must agree with above

rm ttimes.out
rm times.tmp
rm toang.tmp

declare -i N_PHASE=0

#   The files in this anonymous FTP implement the travel time computation
#algorithm described by Buland and Chapman (1983) in "The Computation of
#Seismic Travel Times", BSSA, v. 73, pp. 1271-1302, for the IASPEI phase
#set derived from the IASPEI travel-time model developed by Brian Kennett.
#This phase set currently includes:  P, Pdiff, PKP, PKiKP, pP, pPdiff,
#pPKP, pPKiKP, sP, sPdiff, sPKP, sPKiKP, PP, P'P', S, Sdiff, SKS, pS,
#pSdiff, pSKS, sS, sSdiff, sSKS, SS, S'S', PS, PKS, SP, SKP, SKiKP, PcP,
#PcS, ScP, ScS, PKKP, PKKS, SKKP, and SKKS.

# IMPORTANT! - order here must match ttimes.c
PHASE_LIST="P S"
# !! NOTE: P includes P, Pg, Pb, Pn, and similar for S
# !! NOTE: PKKP gives no times! ???
# !! NOTE: P'P' gives no times, use PKPPKP
# test here
#PHASE_LIST="PKKPab PKKPbc PKKPdf"
#PHASE_LIST="PKKPdf"
#
echo "Phase list: ${PHASE_LIST}"
for PHASE in ${PHASE_LIST} ; do

	echo "Processing: ${PHASE}"
	
	rm ttimes.lst
	rm toang.lst
	
	KEYWORD=basic
	if [ ${PHASE:0:5} == "PKiKP" ] ; then
		KEYWORD="PKiKP"
	elif [ ${PHASE:0:5} == "PKKP" ] ; then
		KEYWORD="PKKP"
	elif [ ${PHASE:0:3} == "SKS" ] ; then
		KEYWORD="SKS"
	elif [ ${PHASE:0:3} == "SKP" ] ; then
		KEYWORD="SKP"
	elif [ ${PHASE} == "PcP" ] ; then
		KEYWORD="PcP"
	elif [ ${PHASE} == "ScP" ] ; then
		KEYWORD="ScP"
	elif [ ${PHASE} == "ScS" ] ; then
		KEYWORD="ScS"
	elif [ ${PHASE:0:4} == "PKKS" ] ; then
		KEYWORD="PKKS"
	elif [ ${PHASE:0:4} == "SKKS" ] ; then
		KEYWORD="SKKS"
	fi

nice ttimes_ajl << EOT1 >> ttimes.out
y
${IASPEI_TAU_PATH}/build-tables/scsn_hk
${PHASE}
${KEYWORD}

n
${DEPTH}
${DISTANCES}
-10. -10. -10.
-10. -10. -10.
EOT1

	cat ttimes.lst >> times.tmp
	cat toang.lst >> toang.tmp
	
	N_PHASE=N_PHASE+1

done

sed 's/-1.00/-1/g' times.tmp > times1.tmp
sed 's/ //g' times1.tmp > times2.tmp
sed 's/-1.00/-1/g' toang.tmp > toang1.tmp
sed 's/ //g' toang1.tmp > toang2.tmp
sed 's/NaN/90/g' toang2.tmp > toang3.tmp

cat << END > times.tmp
static double dist_time[${N_PHASE}][2001][83] = {
END
cat times2.tmp >> times.tmp
cat << END >> times.tmp
};
END

cat << END > toang.tmp
static double dist_toang[${N_PHASE}][2001][83] = {
END
cat toang3.tmp >> toang.tmp
cat << END >> toang.tmp
};
END

mv times.tmp ttimes_${MODEL_GRID_NAME}_times_phases.h_part
mv toang.tmp ttimes_${MODEL_GRID_NAME}_toang_phases.h_part

rm *.tmp

echo "To install phases header files:"
echo "cp ttimes_${MODEL_GRID_NAME}_times_phases.h_part ~/soft/trace_processing/timedomain_processing/ttimes/ttimes_${MODEL_GRID_NAME}_times_phases.h"
echo "cp ttimes_${MODEL_GRID_NAME}_toang_phases.h_part ~/soft/trace_processing/timedomain_processing/ttimes/ttimes_${MODEL_GRID_NAME}_toang_phases.h"
echo
echo "Remember to create and edit ~/soft/trace_processing/timedomain_processing/ttimes/ttimes_${MODEL_GRID_NAME}.h"

