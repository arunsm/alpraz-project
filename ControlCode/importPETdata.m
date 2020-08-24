PETlabels = {'5HT1a', '5HT1b', '5HT2a', 'D1', 'D2', 'DAT', 'FDOPA', 'GABAa', 'NAT', 'SERT'};

parameters

% load PET data
PETdir = '../../data/PETatlas/';

X = importdata(strcat(PETdir, '5HT1a_WAY_HC36.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT1a_WAY_HC36 = X.data'; PET_5HT1a_WAY_HC36(idx_brainstem) = [];
X = importdata(strcat(PETdir, '5HT1b_P943_HC22.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT1b_P943_HC22 = X.data'; PET_5HT1b_P943_HC22(idx_brainstem) = [];
X = importdata(strcat(PETdir, '5HT2a_ALT_HC19.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT2a_ALT_HC19 = X.data'; PET_5HT2a_ALT_HC19(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'D1_SCH23390_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_D1_SCH23390_c11 = X.data'; PET_D1_SCH23390_c11(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'D2_RACLOPRIDE_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_D2_RACLOPRIDE_c11 = X.data'; PET_D2_RACLOPRIDE_c11(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'DAT_DATSPECT.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_DAT_DATSPECT = X.data'; PET_DAT_DATSPECT(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'FDOPA_f18.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_FDOPA_f18 = X.data'; PET_FDOPA_f18(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'GABAa_FLUMAZENIL_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_GABAa_FLUMAZENIL_c11 = X.data'; PET_GABAa_FLUMAZENIL_c11(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'NAT_MRB_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_NAT_MRB_c11 = X.data'; PET_NAT_MRB_c11(idx_brainstem) = [];
X = importdata(strcat(PETdir, 'SERT_DASB_HC30.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_SERT_DASB_HC30 = X.data'; PET_SERT_DASB_HC30(idx_brainstem) = [];