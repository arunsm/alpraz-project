#assumes have already run created working directory wd with design files, as well as mask and merged images in merged_images directory.
#designs are expected to be named "design".mat, .con, .grp, etc etc
#this version operates on grid by calling flameo_grid_effort.sh, probably unnecessary for neff as relatively small groups;

#set number of subjects, working directory, and copes to operate on.
nsubj=n45  #eg n45

#wd=/import/monstrum/neff//group_results_n1445/frac2back/voxelwise_analyses/n479_6EVbehav_td_ps4_zeroinc
wd=/import/monstrum/neff/analyses/3rdlev/n45_tasksubvalz_all45mean/

lev2_cope=cope1
#only specify 1 level 2 cope
copes="cope1 cope2" #1stlevel copes to get "cope1 cope2 cope3" etc

nickname="tasksubvalz_20150427"

echo "working directory is $wd"
echo ""
echo "nickname is $nickname"
echo ""
if [ ! -d "$wd" ]; then
	echo "working directory not found!! exiting!"
	exit 1
fi

echo "*********"
echo "for all analyses, lev2 cope is $lev2_cope"
echo "*******"

for c in $copes; do
	echo ""
	echo "working on $c"

cfile=$(ls -d /import/monstrum/neff/analyses/merged_images/${nsubj}_${nickname}_lev2${lev2_cope}_${c}.nii.gz)
vfile=$(ls -d /import/monstrum/neff/analyses/merged_images/${nsubj}_${nickname}_lev2${lev2_cope}_var${c}.nii.gz)

mask=/import/monstrum/neff/analyses/merged_images/${nickname}${nsubj}_final_maskoverlap.nii.gz

outdir=${nickname}_${nsubj}_lev2${lev2_cope}_${c}  #makes output dir within working directory defined above


outpath=$(ls -d $wd/$outdir 2> /dev/null)


#flameo --cope=$cfile --varcope=$vfile --mask=$mask --dm=${wd}/design.mat --tc=${wd}/design.con  --cs=${wd}/design.grp --runmode=flame1 --logdir=${wd}/${outdir} --npo

qsub -V -q long.q -S /bin/bash -o ~/sge_out -e ~/sge_out /import/monstrum/neff/progs/bash/imaging/flameo_grid_effort.sh $cfile $vfile $mask $wd $outdir


# if have fstat designs, add the above after the .grp: --fc=${wd}/design.fts

done
