#!/bin/bash

HOME_DIR=/home/amahad
BASE_DIR=/data/joy/BBL/studies/alpraz/rawData
TOOLS_DIR=/data/joy/BBL/applications/bids_apps

subject=${1}
#subject=013541
echo "job running"

singularity run --cleanenv -B ${BASE_DIR}:/mnt ${TOOLS_DIR}/fmriprep.simg \
/mnt/NIFTI/ /mnt/derivatives \
participant \
 -w /tmp \
 --participant-label ${subject} \
 --fs-license-file ${HOME_DIR}/license.txt \
 --use-syn-sdc \
 --fs-no-reconall
