#!/bin/bash

# A script for obtaining dicom info from a directory of dicoms, structured: {subjID}/{scanID}/{scantype}/dicoms/{dicom}.nii.gz
# Modified to get dicom info for recently downloaded XNAT dicoms 2018/11/15; modifications shown by commented out lines

# obtain scan and session labels

# Modified_TT
# scans=/data/jux/BBL/studies/reward/rawData/*/*/
scans=/data/joy/BBL/studies/alpraz/rawData/*/*/

for sc in $scans; 

	# Modified_TT
	# do ses=$(echo $sc|cut -d'/' -f9); 
	# subID=$(echo $sc|cut -d'/' -f8);
	do
	echo $sc
	ses=$(echo $sc|cut -d'/' -f9); 
	subID=$(echo $sc|cut -d'/' -f8);
	
	echo $ses
	echo $subID
# USE SINGULARITY HERE TO RUN HEUDICONV FOR DICOM INFO
# note to replace axu with your chead name instead

	# Modified_TT
	# singularity run -B /data/jux/BBL/studies/reward/rawData:/home/axu/base /data/joy/BBL/applications/heudiconv/heudiconv-latest.simg -d /home/axu/base/{subject}/{session}/*/*/*.dcm -o /home/axu/base/output -f convertall -s ${subID} -ss ${ses}  -c none --overwrite;
	
	singularity run -B /data/joy/BBL/studies/alpraz/rawData:/home/amahad/data /data/joy/BBL/applications/heudiconv/heudiconv-latest.simg -d /home/amahad/data/{subject}/{session}/*.dcm -o /home/amahad/data/output -f convertall -s ${subID} -ss ${ses}  -c none --overwrite;

done 
