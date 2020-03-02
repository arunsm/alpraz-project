#!/bin/bash/

#$ -j y
#$ -cwd

matlab -nodisplay -r "parameterAnalysis; exit"
