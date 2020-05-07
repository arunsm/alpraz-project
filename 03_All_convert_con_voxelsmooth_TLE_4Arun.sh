source /data_/mica1/03_projects/lorenzo/micasoft/bash/mica_do_cmd
IN=/data_/mica1/04_data/18_londonGradientsBB/
OUT=${IN}/coreg3_max_TLE/
OUT2=${IN}/coreg3_con_max5mm_TLE/
mkdir -v ${OUT2}

TM=${IN}/con-native-spm-TLE/
SUBJECTS_DIR=/data_/mica1/04_data/18_londonGradientsBB/TLE_HS/
for i in `ls -d ${IN}/meanfunc-native-TLE/mean??????_*.nii`
do 
	echo ${i}
	 
	b=`basename ${i}`
	pre=${b:4:6}
	echo ${pre}
	# prefix should now look like 'ADWI01'
	a=${b%_0000.nii}
	a=${a:4}
	echo $a
	
	FSFOLDER=`ls -d ${SUBJECTS_DIR}/co${pre}*`
	FSFOLDER=`basename ${FSFOLDER}`
	echo ${FSFOLDER}

	echo "${OUT}/${b%_0000.nii}_bbr.reg "
	
	bbregister --mov ${i} \
	--s ${FSFOLDER} --bold \
	--reg ${OUT}/${b%_0000.nii}_bbr.reg 


	echo "estimate registration and map mean func time series to cortex"
	 
	for tmap in `ls ${TM}*${a}*con*.nii`
	do
	
		map=`basename ${tmap}`
		mri_vol2surf --mov ${TM}/${map} \
		--reg ${OUT}/${b%_0000.nii}_bbr.reg \
		--projfrac-max .2 .8 .1 \
		--cortex --hemi lh --interp trilinear \
		--o ${OUT2}/${map%.nii}_bbr_lh.mgh
	
		mri_vol2surf --mov ${TM}/${map} \
		--reg ${OUT}/${b%_0000.nii}_bbr.reg \
		--projfrac-max .2 .8 .1 \
		--cortex --hemi rh --interp trilinear \
		--o ${OUT2}/${map%.nii}_bbr_rh.mgh
		
		# this step multiplies by 10000 in left hemisphere! 
		#if you do -- remember to then divide by 10000!
		fscalc ${OUT2}/${map%.nii}_bbr_lh.mgh \
		mul 10000 --o ${OUT2}/${map%.nii}_bbr_lh_10000.mgh

		# this step multiplies by 10000 in right hemisphere! 
		#if you do -- remember to then divide by 10000!
		fscalc ${OUT2}/${map%.nii}_bbr_rh.mgh \
		mul 10000 --o ${OUT2}/${map%.nii}_bbr_rh_10000.mgh
		
		# the following surf2surf applies smoothing - left hemisphere
		#don't do if using already smoothed data
		mri_surf2surf --s ${FSFOLDER} \
		--sval ${OUT2}/${map%.nii}_bbr_lh_10000.mgh \
		--trgsubject fsaverage5 \
		--cortex \
		--noreshape \
		--hemi lh \
		--fwhm-src 5 \
		--tval ${OUT2}/${map%.nii}_bbr_lh_fsa5_10000.mgh
		
		# the following surf2surf applies smoothing - right hemisphere
		#don't do if using already smoothed data
		mri_surf2surf --s ${FSFOLDER} \
		--sval ${OUT2}/${map%.nii}_bbr_rh_10000.mgh \
		--trgsubject fsaverage5 \
		--cortex \
		--noreshape \
		--hemi rh \
		--fwhm-src 5 \
		--tval ${OUT2}/${map%.nii}_bbr_rh_fsa5_10000.mgh

		# the following is a version of surf2surf without smoothing - left hemisphere
		mri_surf2surf --s ${FSFOLDER} \
		--sval ${OUT2}/${map%.nii}_bbr_lh_10000.mgh \
		--trgsubject fsaverage5 \
		--cortex \
		--noreshape \
		--hemi lh \
		--tval ${OUT2}/${map%.nii}_bbr_lh_fsa5_10000_unsm.mgh
	
		# the following is a version of surf2surf without smoothing - right hemisphere
		mri_surf2surf --s ${FSFOLDER} \
		--sval ${OUT2}/${map%.nii}_bbr_rh_10000.mgh \
		--trgsubject fsaverage5 \
		--cortex \
		--noreshape \
		--hemi rh \
		--tval ${OUT2}/${map%.nii}_bbr_rh_fsa5_10000_unsm.mgh
	done 
done
