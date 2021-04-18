% directory to store results
resultsDir = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_QA_betas/';

idx_brainstem = 234; % index of brain stem parcel    

% load Lausanne parcel info
load LindenYeoPurity/yeo7netlabelsLaus125EJC.mat;
subSystemLabels = {'Visual', 'Somatomator', 'DorsalAttention', 'VentralAttention', 'Limbic', 'FrontoparietalControl', 'DefaultMode', 'Subcortical'};
subcorticalIndices = find(finalLabels == 8);
subcorticalIndices(subcorticalIndices==idx_brainstem) = [];
nifti = load_nii('ROIv_scale125_dilated.nii.gz'); % load Lausanne nifti
X = importdata('../../data/lausanne2008/LausanneParcelNames.xlsx');
LausanneParcelNames = X.textdata;
nSubSystems = numel(subSystemLabels);

contrastLabels = {'contrast1_threatcorrectStd', 'contrast3_nonthreatcorrectStd', ...
    'contrast5_neutralcorrectStd'};
nContrasts = numel(contrastLabels);

parcelCoverageThresh = 0.5;
nNodes = 233;
nTimeSteps = 1001;

pathToDemographics = '../../data/Alpraz_subjectDemographics.xlsx';
pathToMotionFile = '../../data/subjectMotion.csv';
subjectDemographics = readtable(pathToDemographics);
subjectMotion = readtable(pathToMotionFile);