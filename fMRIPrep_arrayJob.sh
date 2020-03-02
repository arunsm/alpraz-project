#!/bin/bash

#$ -j y
#$ -cwd
#$ -l h_vmem=25.5G,s_vmem=25G
#$ -M arun.mu@gmail.com

SUBJECT_LIST=/data/joy/BBL/studies/alpraz/rawData/subjectList2.txt

mapfile -t ARRAY < ${SUBJECT_LIST}
LENGTH=${#ARRAY[@]}
echo Number of Subjects: $LENGTH 
echo SGE_TASK_ID: $SGE_TASK_ID

INDX=`expr $SGE_TASK_ID - 1`;

if [[ $INDX -ge $LENGTH ]]; then
 echo Array index greater than number of elements
else
 SUB=${ARRAY[$INDX]}
 echo Subject $SUB
 ./fMRIPrep_singleSubject.sh ${SUB}
fi
