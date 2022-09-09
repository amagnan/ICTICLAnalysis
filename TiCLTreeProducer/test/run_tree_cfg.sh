#! /bin/sh

# 1. Prepare the input file with all the step3 files to process.
#EOSDIRIN=/eos/cms/store/user/amagnan/HGCAL/TiCL
#GEOMDIR=D49_EM_trkFirst/
#PTETA=_

EOSDIRIN=/eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan

#210420
#210325
#EMonlyNew
#EOSDIRIN=/eos/cms/store/user/amagnan/HGCAL/TiCL/EMonlyNewDef
#GEOMDIR=D49_Dummy_new/
GEOMDIR=D49_FineCalo/
#GEOMDIR=D49_DefSC/

#DataType=ChargedPionsFromVtx #CloseByPhotons #CloseByPhotonsFromVtxWithPU
DataType=SinglePhotons/FineCalo/
#ElectronsFromVtx #CloseByPhotonsFromVtx #CloseByPhotonsWithPU #CloseByPhotonsFromVtxWithPU

for PT in 40to400
#for PT in 3 5 10 15 20 30 40 50 75 100 150 200;
#for PT in 3 10 20 50 100 150 200;
#for PT in 40 50 75 100 150 200;
do
    for eta in 21; #17 19 23 25 27;
    do

	if (( ${eta} == "17" )); then
	    if (( ${PT} == "5" )); then
		fillTriplets=1
	    else 
		fillTriplets=0
	    fi
	else
	    fillTriplets=0
	fi

	echo "Filling triplets: $fillTriplets"

	PTETA=_En${PT}_eta${eta}

	for d in $DataType 
	do
	    echo "Processing "$GEOMDIR" "$d
	    # Establish output directory based on sample name.
	    # Check if it exists, create if not.
	    OUT_DIR=$(basename $d)
	    if [ ! -e ${GEOMDIR}/$d ]; then
		echo "Creating directory ${GEOMDIR}/$d"
		mkdir -p ${GEOMDIR}/${d}
	    fi

	    echo step3${PTETA}_\*.root
	    # Gather all files to process and save them in local ASCII file.
	    find $EOSDIRIN/${d}/ -maxdepth 1 -type f -iname step3ticl${PTETA}_\*.root -print | tail -250 > ${GEOMDIR}/files_${OUT_DIR}${PTETA}.txt
           #find $EOSDIRIN/${d}/ -maxdepth 1 -type f -iname step3${PTETA}_\*.root -print > ${GEOMDIR}/files_${OUT_DIR}.txt
	    #find $EOSDIRIN/${d}/ -maxdepth 1 -type f -iname step3_\*.root -print | grep -v "_pt" > ${GEOMDIR}/files_${OUT_DIR}.txt
	    
	    cat ${GEOMDIR}/files_${OUT_DIR}${PTETA}.txt
	    #sleep 60
	    
	    # 2. Loop over the files and lauch N_PARALLEL cmsRun jobs
	    
	    
	    OUT_SUFFIX="_FlatTracksters"
	    N_PARALLEL=10
	    for f in $(cat ${GEOMDIR}/files_${OUT_DIR}${PTETA}.txt); do
		while [ "$(pgrep -c cmsRun)" -ge "$N_PARALLEL" ]; do
		    sleep 10
		    echo -n "."
		done
		#IN_FILES=`awk 'BEGIN{ORS=","}{print "file:"$1}' files_${OUT_DIR}.txt`
		OUT_FILE=${GEOMDIR}/${d}/$(basename $f .root)
		#${x%?} strips out last character... 
		#echo "cmsRun run_pid_cfg.py inputFiles=${IN_FILES%?} outputFile=${OUT_FILE}${OUT_SUFFIX}.root | tee ${OUT_FILE}.log"
		echo "cmsRun run_tree_cfg.py fillTriplets=${fillTriplets} inputFiles=file:$f outputFile=${OUT_FILE}${OUT_SUFFIX}.root | tee ${OUT_FILE}.log"
		cmsRun run_tree_cfg.py fillTriplets=${fillTriplets} inputFiles=file:$f outputFile=${OUT_FILE}${OUT_SUFFIX}.root &> ${OUT_FILE}.log &
		sleep 2
		#echo "rm ${GEOMDIR}/${d}/step3ticl${PTETA}${OUT_SUFFIX}.root; hadd ${GEOMDIR}/${d}/step3ticl${PTETA}${OUT_SUFFIX}.root ${GEOMDIR}/${d}/step3ticl${PTETA}*${OUT_SUFFIX}.root"
		#break # uncommented if want to only process first file?
	    done
	done
    done
done

for d in $DataType
do
    #echo "for PT in 3 5 10 15 20 30 40 50 75 100 150 200; do for eta in 17 19 21 23 25 27; do PTETA=_pt\${PT}_eta\${eta}; rm $GEOMDIR/$d/step3ticl\${PTETA}_FlatTracksters.root; hadd $GEOMDIR/$d/step3ticl\${PTETA}_FlatTracksters.root $GEOMDIR/$d/step3ticl\${PTETA}_run*_FlatTracksters.root; done; done"
    echo "for PT in 50; do for eta in 21; do PTETA=_pt\${PT}_eta\${eta}; rm $GEOMDIR/$d/step3ticl\${PTETA}_FlatTracksters.root; hadd $GEOMDIR/$d/step3ticl\${PTETA}_FlatTracksters.root $GEOMDIR/$d/step3ticl\${PTETA}_run*_FlatTracksters.root; ll $GEOMDIR/$d/*; done; done"
done
