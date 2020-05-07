# resample Lausanne nifti to resolution of PET maps (3mm MNI space)
MASK=/data/joy/BBL/studies/alpraz/parcellations/LausanneNifti/ROIv_scale125_dilated.nii.gz
RESIZED_MASK=/data/joy/BBL/studies/alpraz/parcellations/LausanneNifti/ROIv_scale125_dilated_resized_3mm.nii.gz
MASTER=/data/joy/BBL/studies/alpraz/PETatlas/5HT1a_WAY_HC36.nii

3dresample -master $MASTER -prefix $RESIZED_MASK -input $MASK

# extract Lausanne parcel values for PET maps
FNAME_ARRAY=(/data/joy/BBL/studies/alpraz/PETatlas/*)
LENGTH=${#FNAME_ARRAY[@]}

for (( i=0; i<=$LENGTH-1; i++ ))

	do
	CURRENTFNAME=${FNAME_ARRAY[$i]}
	echo $CURRENTFNAME
	OUTPUT_PATH="${CURRENTFNAME}_lausanne_ROIv_scale125_dilated_resized_3mm.txt"
	
	3dROIstats -mask $RESIZED_MASK $CURRENTFNAME > $OUTPUT_PATH
done
