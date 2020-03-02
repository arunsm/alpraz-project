#!/bin/bash
unset PYTHONPATH
export PATH=/data/joy/BBL/applications/miniconda3/bin:$PATH
source activate py2k
python /data/joy/BBL/studies/alpraz/scripts/downloadAlpraz.py
