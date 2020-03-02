# script to compute parcel averages of contrasts for all subjects and sessions

MASK_PATH=/data/joy/BBL/studies/alpraz/parcellations/LausanneNifti/ROIv_scale125_dilated.nii.gz
BASE_DIR=/data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks


SUB_ARRAY=( $(tail -n +2 /data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks/group/n191_quality.csv | cut -d , -f1))
SES_ARRAY=( $(tail -n +2 /data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks/group/n191_quality.csv | cut -d , -f2))

LENGTH=${#SUB_ARRAY[@]}

for (( i=0; i<=$LENGTH-1; i++ ))
	
	do
	#echo $i
	SUB=${SUB_ARRAY[$i]}
	SES=${SES_ARRAY[$i]}
        
        echo $SUB
	echo $SES
	
	# EMOTION-ID

	# contrast1: threatcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast1_threatcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast1_threatcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# contrast3: nonthreatcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast3_nonthreatcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast3_nonthreatcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# contrast5: neutralcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast5_neutralcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast5_neutralcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# contrast7: anystimcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast7_anystimcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_contrast7_anystimcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# individual subject masks
        OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_maskStd_lausanne_ROIv_scale125_dilated.txt"
        MASK_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionid/norm/"${SUB}_${SES}_task-emotionid_maskStd.nii.gz"
        # call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
        3dROIstats -mask $MASK_PATH $MASK_IMAGE_PATH > $OUTPUT_PATH

	# EMOTION-REC

	# contrast1: threatcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast1_threatcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast1_threatcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# contrast3: nonthreatcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast3_nonthreatcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast3_nonthreatcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# contrast5: neutralcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast5_neutralcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast5_neutralcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# contrast7: anystimcorrect	
	OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast7_anystimcorrectStd_lausanne_ROIv_scale125_dilated.txt"
	CONTRAST_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_contrast7_anystimcorrectStd.nii.gz"
	# call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
	3dROIstats -mask $MASK_PATH $CONTRAST_IMAGE_PATH > $OUTPUT_PATH

	# individual subject masks
        OUTPUT_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_maskStd_lausanne_ROIv_scale125_dilated.txt"
        MASK_IMAGE_PATH=$BASE_DIR/$SUB/$SES/task-emotionrec/norm/"${SUB}_${SES}_task-emotionrec_maskStd.nii.gz"
        # call AFNI's 3dROIstats program to compute average contrast values over all ROIs in the 234-node Lausanne parecllation
        3dROIstats -mask $MASK_PATH $MASK_IMAGE_PATH > $OUTPUT_PATH

done

