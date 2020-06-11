#!/bin/bash/

#$ -j y
#$ -cwd

matlab -nodisplay -r "calculateControlEnergy; exit"
