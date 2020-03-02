#!/bin/bash
#$ -l h_vmem=19.6G,s_vmem=19.5G
#$ -cwd

singularity build /data/joy/BBL/studies/alpraz/xcpEngine-master/xcpEngine.simg docker://pennbbl/xcpengine:latest
