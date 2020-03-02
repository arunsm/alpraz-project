
#goal of this is to merge cope & varcope images and make mask for group analysis with flameo
#designed to run on 2ndlevel neff analyses
#set bblid list, nickname, copes

#######
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_effortv1_n45.txt)  #this is list of bbblids to merge for flameo use later
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_effortv1_n29nosvconfcorr.txt) #29 subjects that on average lack sv-dconf correl
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_effortv1_n30okforchosenval.txt) #30 subjects with reasonable hard/easy props
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_effortv1_n8okforhardeasy.txt) #8 subj
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n29_effortv1_2002.txt) #29 subjects with 2002 order
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n14_effortv1_0220.txt) #14 subjects with 0220 order
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n29_effortv1_2002.txt /import/monstrum/neff/analyses/bblidlists/bblid_n14_effortv1_0220.txt) #43 subjects 29 with 2002 then 14 with 0220 order
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n16_effortv1chosenval_2002.txt /import/monstrum/neff/analyses/bblidlists/bblid_n12_effortv1chosenval_0220.txt) #28 subjects for runspecific chosenval (intersection of n30 with n29 plus n14);
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n44_v2r2.txt) #all subjects with run2 have run1, but converse is not true;
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n49_v2r1.txt)
#bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n35_v2r2.txt) #subjects with run2, excluding those with high motion and weirdness
bblids=$(cat /import/monstrum/neff/analyses/bblidlists/bblid_n52_v2r2.txt) #subjects with run2, including the additional subjects sheila ran after aylin left


#for effort_v1 is n45
outdir=/import/monstrum/neff/analyses/merged_images/
#feat=level2_fe_tasksubvalz_20150427 #2nd level flameo folder
#nickname="tasksubvalz_20150427"
#nickname="taskdconfz_20150706"
#nickname="taskpresszmoneyz_20150706"
#nickname="taskdconfrtz_20160224"  
#nickname="tasksubvalrtz_20160224"
#nickname="taskpressmoneyrtz_20160224"
#nickname="tasksubvalzhyper_20160307"
#nickname="taskchosenvalz_20160303"
#nickname="hardeasyhardsvzeasysvz_20160301"
#nickname="zerotwo2C_????_taskpresszmoneyz_20150706"
#nickname="zerotwo_????_taskchosenvalz_20160303"
nickname="v2_2E4C_taskpresszmoneyz_20160308"


feat=level2_fe_$nickname #2nd level flameo folder
#feat name can use wild cards (maybe not in current coding where nickname is used to define featname and also used separately).

lev2_cope="cope1 cope2 cope3 cope4" #eg, group mean. note flameo directory structure is flameodirectory/lev1cope/lev2cope.nii.gz
#only specify 1 level 2 cope at a time here? I tried to let it loop here;
copes="cope1 cope2 cope3 cope4" #level1 copes to merge-- can be multiple copes e.g., "cope1 cope2 cope3 cope7"
mask="mask"
varcopeflag=1   #if=1, also will merge varcopes 
maskflag=1 # if =1, also will merge masks   
#could automate mask decision using code below....
#if [ ! -e "$outdir/n${nsubj}_mask.nii.gz" ]; then
#echo ""
#echo "making mask"
#fi


logdir=/import/monstrum/neff/progs/bash/imaging/logs/
#######

for cc in $lev2_cope; do

rm -f $logdir/effort_missing_copes_for_merge.txt
echo "feat is $feat"

nsubj=$(echo $bblids | wc |  awk '{print $2}')
echo "$nsubj subjects"

echo "*********"
echo "lev2 copes are: $lev2_cope"
echo "*******"

for c in $copes; do
	echo ""
	echo ""
	echo "********************"
	echo "now working on $c"
	echo "********************"
	rm -f copelist_tmp.txt
	rm -f varcopelist_tmp.txt

	for b in $bblids; do
		echo ""
		echo $b
	
		#image=$(ls -d /import/monstrum/neff/subjects/${b}*/level2/*${feat}*/${c}/${lev2_cope}.nii.gz)
		image=$(ls -d /import/monstrum/neff/subjects/${b}*/level2/*${feat}*/${c}/${cc}.nii.gz)
		#don't use {c}*  bc breaks if >9 copes and want cope1

		if [ -z "$image" ]; then
			echo "COPE MISSING"
			echo -e"$b \t $c" >> $logdir/effort_missing_copes_for_merge.txt
			exit 1
		fi
		echo $image
		echo $image >> copelist_tmp.txt
		
		if [ "$varcopeflag" == 1 ]; then
			echo "merging varcopes also"	
			#varimage=$(ls -d /import/monstrum/neff/subjects/${b}*/level2/*${feat}*/${c}/var${lev2_cope}.nii.gz)	
			varimage=$(ls -d /import/monstrum/neff/subjects/${b}*/level2/*${feat}*/${c}/var${cc}.nii.gz)
		#note be careful some analyses use ${c}_2mm. nii.gz and others use ${c}.nii.gz 

			if [ -z "$varimage" ]; then
				echo "VARCOPE MISSING"
				echo -e"$b \t var${c}" >> $logdir/effort_missing_copes_for_merge.txt
				exit 1
			fi
		echo $varimage
		echo $varimage >> varcopelist_tmp.txt
		fi

		if [ "$maskflag" == 1 ]; then
			echo "merging masks also"	
			maskimage=$(ls -d /import/monstrum/neff/subjects/${b}*/level2/*${feat}*/${c}/mask.nii.gz)	 

			if [ -z "$maskimage" ]; then
				echo "MASK MISSING"
				echo -e"$b \t mask" >> $logdir/effort_missing_copes_for_merge.txt
				exit 1
			fi
		echo $maskimage
		echo $maskimage >> masklist_tmp.txt
		fi


	done	


	echo ""
	echo "********************"
	echo "now merging $nsubj images for $c"
	#fslmerge -t ${outdir}/n${nsubj}_${nickname}_lev2${lev2_cope}_${c} $(cat copelist_tmp.txt) 
	fslmerge -t ${outdir}/n${nsubj}_${nickname}_lev2${cc}_${c} $(cat copelist_tmp.txt) 
	
	if [ "$varcopeflag" == 1 ]; then
	
		echo ""
		echo "now merging $nsubj images for var${c}"
		#fslmerge -t ${outdir}/n${nsubj}_${nickname}_lev2${lev2_cope}_var${c} $(cat varcopelist_tmp.txt) 
		fslmerge -t ${outdir}/n${nsubj}_${nickname}_lev2${cc}_var${c} $(cat varcopelist_tmp.txt) 
		echo "********************"
	fi

	if [ "$maskflag" == 1 ]; then
	
		echo ""
		echo "now merging $nsubj images for mask"
		fslmerge -t ${outdir}/n${nsubj}_mask $(cat masklist_tmp.txt) 
		fslmaths ${outdir}/n${nsubj}_mask -Tmin -bin -mas /import/monstrum/Applications/fsl5/data/standard/MNI152_T1_2mm_brain_mask ${outdir}/${nickname}n${nsubj}_final_maskoverlap
		echo "********************"
	fi

done

rm -f copelist_tmp.txt
rm -f varcopelist_tmp.txt
rm -f masklist_tmp.txt

#echo "output at ${outdir}/n${nsubj}_${nickname}_lev2${lev2_cope}_${c}"
echo "output at ${outdir}/n${nsubj}_${nickname}_lev2${cc}_${c}"

done #end of 2nd level cope loop?
