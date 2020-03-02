#!/bin/bash

# -l h_vmem=16G,s_vmem=15.5G

SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=/data/joy/BBL/applications/bids_apps/xcpEngine.simg
 
${SNGL} run -B /data:/home/amahad/data ${SIMG} \
 -c /home/amahad/data/joy/BBL/studies/alpraz/rawData/derivatives/cohortFile.csv \
 -d /home/amahad/data/joy/BBL/studies/alpraz/rawData/derivatives/xcpEngine-master/designs/task_alpraz.dsn \
 -o /home/amahad/data/joy/BBL/studies/alpraz/rawData/derivatives/xcp_output_allTasks  \
 -r /home/amahad 
  
 
