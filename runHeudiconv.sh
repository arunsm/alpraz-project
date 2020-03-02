#!/bin/bash

#$ -j y
#$ -cwd
#$ -l h_vmem=10.5G,s_vmem=10G

# obtain scan and session labels
scans=/data/joy/BBL/studies/alpraz/rawData/*/*/

for sc in $scans; 
	do 
	ses=$(echo $sc|cut -d'/' -f9); 
	subID=$(echo $sc|cut -d'/' -f8);
	
	echo $ses
	echo $subID

# USE SINGULARITY HERE TO RUN HEUDICONV FOR BIDS FORMATTING

	singularity run \
		-B /data/joy/BBL/studies/alpraz/rawData:/mnt \
		/data/joy/BBL/applications/heudiconv/heudiconv-latest.simg \
		-d /mnt/{subject}/{session}/*.dcm \
		-o /mnt/NIFTI \
		-f /mnt/alpraz_heuristic.py \
		-s ${subID} -ss ${ses} \
		-c dcm2niix -b --overwrite

done 
