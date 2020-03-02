%make a 2nd or 3rd level EV design matrix (for design.mat) for neff flameo
%analysis

addpath /import/monstrum/Users/danwolf/matlab/all2file

%bblidimg=load('/import/monstrum/neff/analyses/bblidlists/bblid_effortv1_n45.txt');%get bblidlist from /import/monstrum/neff/progs/bash/imaging/effort_merge_images_mask.sh;
%bblidimg=load('/import/monstrum/neff/analyses/bblidlists/bblid_n44_v2r2.txt'); %NOTE this list is in unix number order not normal numerical order!
bblidimg=load('/import/monstrum/neff/analyses/bblidlists/bblid_n35_v2r2.txt'); %NOTE this list is in unix number order not normal numerical order!


%datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n45demobehavlong_7-14-15.txt');%this is the matrix with all the subject variables;
%datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n44_fslorder_demodxcains_5-20-16.txt');%this is the matrix with all the subject variables;
%data here for v2 5-20-16 has: 1bblid	2age	3sex	4smoke	5education	6pared	7dx0c2s	race	ethnicity	handedness	fr_soc_mot	fr_recvoc_mot	13fr_global_mot	fr_flat_affect	fr_psymotor_slow	fr_soc_anh_past	fr_soc_anh_int	fr_soc_anh	fr_soc_anh_amot	fr_soc_all	fr_recvoc_anh_past	fr_recvoc_anh_int	fr_recvoc_anh	fr_recvoc_anh_amot	fr_recvoc_all	fr_phys_anh_past	fr_phys_anh_int	fr_phys_anh	fr_phys_all	fr_anh_amot_all	fr_anh_int_past_all	fr_anh_int_fut_all	fr_anh_freq_past_all	fr_anh_past_all	fr_anh_global	fr_mot_anticip_anh	fr_total_scale_avg	fr_global_avg_psych	fr_anh_avol_global_exp
%datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n35_v2_cnbiqnopmat.txt'); % this is the matrix with zscores for each cnb task and a cnbiq score
%datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n35_v2_cainsmot.txt');
%data here for v2 5-24-16 has: bblid	er40_Zpctacc	cpw_Zpctacc lnb_Zpctacc	pcetdw_Zacc2	pcpt_Zpctacc	spvrt_Zpctacc	pmat_Zpctacc	cnbiq	cnbiqnopmat dx0c2s
%datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n35_v2_kavrank.txt');
datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n35_v2_sapstotal.txt'); % saps_global	
%datastruct=importdata('/import/monstrum/neff/analyses/matlab_analysis/n35_v2_mbrkl.txt'); % mbrkl

data=datastruct.data; 
%cnbiq=data(:,11);
%cainsmot=data(:,13);
%kavrank=data(:,3);
%mbrkl=data(:,3);
sapstotal=data(:,3);

%load in additional clinical data and motion data
%datastruct2=importdata('/import/monstrum/neff/analyses/matlab_analysis/n45clindata79_1-21-16.txt');%45 rows=subjects, 79 var columns;


%rtcorrmat=importdata('n45rtcorrmat_2-29-16.txt');rtcorrmat=rtcorrmat.data;
%this is correl across all trials of reaction time and model metrics:
%bblid then:svrelative4,svconfidence4,svchoserelnot4,svchosen4,choice4,releff4,relmag4
%note the svconfidence-rt corr captures a QA feature, expect neg corr at -.1189 or stronger for one-tail p.05, 6 people miss this, including some with 50/50 choice pattern and some with nearly all-hard or all-easy choice pattern

%svcorrdconf=load('n45svcorrdconf_2-29-16.txt');
%bblid, corr between svrel (hard-easy) and dconf;
%poscorind=find(svcorrdconf(:,2)>=0);
%negcorind=find(svcorrdconf(:,2)<0);
%for n45, mean svcorrdconf is  -0.3362; 
%NOTE rtcorrmat and svcorrdconf are sorted by bblid

%merge data and data2;
%[mat1,mat2]=bblidmatch(data2,data);
%data=[mat2, mat1]; %NOTE ORIGINAL DATA NOW EXTENDED TO INCLUDE MANY MORE VARIABLES, with the newer 79 variables tacked onto the end.  

%use bblidmatch to make sure fmri order and matlab order are the same; (so far they are the same);
%[mat1,mat2]=bblidmatch(bblidimg,data);
%data=mat2; %data now reordered if necessary

%define other variables;

%SET COVARIATES
covarmat=[sapstotal];%also do kadj and log10kadj and mbrkl?
szcovarmat=size(covarmat);
numnan=length(find(isnan(covarmat)));
nanind=find(isnan(covarmat(:))==1);
if numnan>0;
disp(['There were this many NaNs in the covariate matrix: ' num2str(numnan)]);
disp(['Warning! Replacing NaNs in nanzscored covariate matrix with zeros']);
end;%end numnan test;
covarmat=nanzscore(covarmat); covarmat(nanind)=0;

%SET FLAGS
numgp=2; %1 or 2? %
covariate=1; %0 if only running group, 1 if including covariate(s)
withingroup=1;
%0 if covariate demeaned across both groups, 1 if demeaned separately within each group;

%exclude data based on missing values (set to group missing right now)
%excludeind=find(isnan(data(:,4)));
%data(excludeind,:)=''; 

%exclude based on high motion (>.3)
%motionind=find(data(:,11)>.3);
%data(motionind,:)=''; 

%exclude based on fmri_exclude column
%exclude=find(data(:,3)==1);
%data=data(setdiff(1:size(data,1),exclude),:);

%dx EVs
%group=data(:,7); %NC=0,SZ=2 for cainsmot?
group=data(:,2); %NC=0,SZ=2 for cnb, kavrank, saps total
%group=zeros(45,1); group(negcorind)=2;%for now make 2 groups be poscor,negcor
ncind=find(group==0);
szind=find(group==2);
ncinclude=zeros(length(group),1); 
ncinclude(ncind)=1;
szinclude=zeros(length(group),1); 
szinclude(szind)=1;

clear design; %clear prior design;

%if numgp==1; %need to distinguish between treating both groups as 1 and just using 1 of the 2 gruops;
%data(ncind,:)='';
%end;

if numgp==2 & covariate==0;
design=[ncinclude szinclude];% run most basic group only no covariates model
end;

if numgp==1 & covariate==0;
design=ones(size(data,1),1);
end;

if numgp==1 & covariate==1;
design=[ones(size(data,1),1) covarmat];
end;
    
if covariate==1 & numgp==2 & withingroup==0;
design=[ncinclude szinclude covarmat];%covariates across all subjects
end;

if covariate==1 & numgp==2 & withingroup==1;%covariate separately by group
covarmat(nanind)=NaN;%put back nans to properly demean within group
design=zeros(size(data,1),2+2*szcovarmat(2));
design(:,[1 2])=[ncinclude szinclude];
for i=1:szcovarmat(2);
design(ncind,3+2*(i-1))=nanzscore(covarmat(ncind,i));
design(szind,4+2*(i-1))=nanzscore(covarmat(szind,i));
design(find(isnan(design)==1))=0;
end; %end i loop
end;

designincl=design(:,:);
 
%demean all but the group EVs;
szdes=size(designincl);
designincl2=designincl;

%if covariate==1;
%for i=1+numgp:szdes(2);
%designincl2(:,i)=designincl(:,i)-nanmean(designincl(:,i)); %THIS DOESN'T WORK!!
%end;
%end;

numnan=length(find(isnan(designincl2)));
disp(['There were this many NaNs in the design matrix: ' num2str(numnan)]);
designincl2(isnan(designincl2))=0;
disp('Your final design matrix is called designincl2');

%create design.grp
designgrp=ones(length(designincl2),1);

%print out design.mat template
designincl2


ind=find(group==2);
for i=1:size(data,2);dwcorr1(i,1)=nancorr2(mbrkl(ind),data(ind,i));end;
tmpdata=data; tmpstd=nanstd(data);badind=find(tmpstd==0 | isnan(tmpstd)==1);
for i=setdiff([1:1:size(data,2)],badind);
[b,~,stats]=glmfit(nanzscore([group log10kadj]), nanzscore(tmpdata(:,i)),'normal','constant','off');
dwstats3(i,:)=stats.beta;
end;

