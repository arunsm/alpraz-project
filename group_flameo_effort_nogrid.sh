#assumes have already run created working directory wd with design files, as well as mask and merged images in merged_images directory.
#designs are expected to be named "design".mat, .con, .grp, etc etc
#specify nsubj wd nickname(1st/2ndlev) and first and second level copes below
#this version runs without grid (actually, now it runs with grid or without, depending on what lines are commented)

#set number of subjects, working directory, nickname, and copes to operate on.
nsubj=n52 #eg n45

#wd=/import/monstrum/neff/analyses/3rdlev/n30_taskchosenvalz_all30mean_cainsmot/
#wd=/import/monstrum/neff/analyses/3rdlev/n45_tasksubvalz_ctsz_within_pcthard/
#wd=/import/monstrum/neff/analyses/3rdlev/n8_hardeasyhardsvzeasysvz_all8mean/
#wd=/import/monstrum/neff/analyses/3rdlev/n30_tasksubvalz_all30mean/
#wd=/import/monstrum/neff/analyses/3rdlev/n43_zerotwo2C_taskpressmoneyz_2002then0220all43mean
#wd=/import/monstrum/neff/analyses/3rdlev/n28_zerotwo_taskchosenvalz_2002then0220all28mean
#wd=/import/monstrum/neff/analyses/3rdlev/n44v2r1r2_taskpresszmoneyz_all44v2r1r2mean
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n44v2r1r2_taskpresszmoneyz_all44v2r1r2_ct_sz_ctmot_szmot
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n44v2r1r2_taskpresszmoneyz_all44v2r1r2_ct_sz_ctcnbiq_szcnbiq
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n44v2r1r2_taskpresszmoneyz_all44v2r1r2_ct_sz_ctcnbiq_szcnbiq
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n35v2r1r2_taskpresszmoneyz_n35v2r1r2_ct_sz
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n35v2r1r2_taskpresszmoneyz_n35v2r1r2_ct_sz_ctmot_szmot
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n35v2r1r2_taskpresszmoneyz_n35v2r1r2_ct_sz_ctcnbiq_szcnbiq
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n35v2r1r2_taskpresszmoneyz_n35v2r1r2_ct_sz_ctsaps_szsaps
#wd=/import/monstrum/neff/analyses/3rdlev/v2/n35v2r1r2_taskpresszmoneyz_n35v2r1r2_ct_sz_ctmbrkl_szmbrkl
wd=/import/monstrum/neff/analyses/3rdlev/v2/n52v2r1r2_taskpresszmoneyz_all52v2r1r2mean

lev2_cope="cope1 cope2 cope3 cope4" #[does this multi cope loop work ok? i think so]
copes="cope1 cope2 cope3 cope4" #1stlevel copes to get "cope1 cope2 cope3" etc

#nickname="tasksubvalz_20150427" #note nickname is critical - gets right 1st/2ndlevel analysis
#nickname="taskdconfz_20150706"
#nickname="taskpresszmoneyz_20150706"
#nickname="taskdconfrtz_20160224"  
#nickname="tasksubvalrtz_20160224"
#nickname="taskpressmoneyrtz_20160224"
#nickname="taskchosenvalz_20160303"
#nickname="hardeasyhardsvzeasysvz_20160301"
#nickname="tasksubvalzhyper_20160307"
#nickname="zerotwo_????_taskpresszmoneyz_20150706"
#nickname="zerotwo2C_????_taskpresszmoneyz_20150706"
#nickname="zerotwo_????_taskchosenvalz_20160303"
nickname="v2_2E4C_taskpresszmoneyz_20160308"

echo "working directory is $wd"
echo ""
echo "nickname is $nickname"
echo ""
if [ ! -d "$wd" ]; then
	echo "working directory not found!! exiting!"
	exit 1
fi

for cc in $lev2_cope; do

echo "*********"
echo "lev2 cope(s) is/are $lev2_cope"
echo "*******"

for c in $copes; do
	echo ""
	echo "working on level1 $c in level2 $cc"

cfile=$(ls -d /import/monstrum/neff/analyses/merged_images/${nsubj}_${nickname}_lev2${cc}_${c}.nii.gz)
vfile=$(ls -d /import/monstrum/neff/analyses/merged_images/${nsubj}_${nickname}_lev2${cc}_var${c}.nii.gz)
#cfile=$(ls -d /import/monstrum/neff/analyses/merged_images/${nsubj}_${nickname}_lev2${lev2_cope}_${c}.nii.gz)
#vfile=$(ls -d /import/monstrum/neff/analyses/merged_images/${nsubj}_${nickname}_lev2${lev2_cope}_var${c}.nii.gz)

mask=/import/monstrum/neff/analyses/merged_images/${nickname}${nsubj}_final_maskoverlap.nii.gz

outdir=${nickname}_${nsubj}_lev2${cc}_${c}  #makes output dir within working directory defined above
#outdir=${nickname}_${nsubj}_lev2${lev2_cope}_${c}

outpath=$(ls -d $wd/$outdir 2> /dev/null)

#for grid leave next 3 lines operative, to not use grid comment out 3 lines and open up the flameo command below that;
#commandvar="flameo --cope=$cfile --varcope=$vfile --mask=$mask --dm=${wd}/design.mat --tc=${wd}/design.con  --cs=${wd}/design.grp --runmode=flame1 --logdir=${wd}/${outdir} --npo"
#echo $commandvar > tmpflameo.sh
#qsub -V -q long.q -S /bin/bash -o ~/sge_out -e ~/sge_out /import/monstrum/neff/progs/bash/imaging/tmpflameo.sh

nohup flameo --cope=$cfile --varcope=$vfile --mask=$mask --dm=${wd}/design.mat --tc=${wd}/design.con  --cs=${wd}/design.grp --runmode=flame1 --logdir=${wd}/${outdir} --npo &  #without grid

done #end of c level1 loop

done #end of cc level2 loop

#if have fstat designs, add the above after the .grp: --fc=${wd}/design.fts
