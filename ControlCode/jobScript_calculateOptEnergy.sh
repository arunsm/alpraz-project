#!/bin/bash/

#$ -j y
#$ -cwd

matlab -nodisplay -r "calculateOptEnergy; exit"
