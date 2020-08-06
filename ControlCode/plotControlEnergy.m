%% loading relevant files

clear all;

addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/LausanneSurfaceFigures'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/LausanneCoordinates'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/mult_comp_perm_t1'));

resultsDir = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_QA_betas/';

allControlEnergies_emotionid = readtable(strcat(resultsDir, 'allControlEnergies_emotionid.csv'));
allControlEnergies_emotionrec = readtable(strcat(resultsDir, 'allControlEnergies_emotionrec.csv'));
load(strcat(resultsDir, 'allControlTrajectories_emotionid.mat'));
load(strcat(resultsDir, 'allControlTrajectories_emotionrec.mat'));
load(strcat(resultsDir, 'structuralAdjacencyMatrix.mat'));

load LindenYeoPurity/yeo7netlabelsLaus125EJC.mat;
subSystemLabels = {'Visual', 'Somatomator', 'DorsalAttention', 'VentralAttention', 'Limbic', 'FrontoparietalControl', 'DefaultMode', 'Subcortical'};
subcorticalIndices = find(finalLabels == 8);
nifti = load_nii('ROIv_scale125_dilated.nii.gz');

X = importdata('../../data/lausanne2008/LausanneParcelNames.xlsx');
LausanneParcelNames = X.textdata;

nSubSystems = numel(subSystemLabels);
contrastLabels = {'contrast1_threatcorrectStd', 'contrast3_nonthreatcorrectStd', ...
    'contrast5_neutralcorrectStd'};
nContrasts = numel(contrastLabels);

parcelCoverageThresh = 0.5;
nNodes = 234;

%% Table 1 - clinical and demographic variables of cohort

pathToDemographics = '../../data/Alpraz_subjectDemographics.xlsx';
subjectDemographics = readtable(pathToDemographics);

group = subjectDemographics.Group;

fprintf('percent female controls = %.1f\n', 100*sum(subjectDemographics.sex_M0F1(group==0)/sum(group==0)))
fprintf('proportion female controls: %iF/%iM\n', sum(subjectDemographics.sex_M0F1(group==0)), sum(group==0)-sum(subjectDemographics.sex_M0F1(group==0)))
fprintf('percent female relatives = %.1f\n', 100*sum(subjectDemographics.sex_M0F1(group==1)/sum(group==1)))
fprintf('proportion female controls: %iF/%iM\n', sum(subjectDemographics.sex_M0F1(group==1)), sum(group==1)-sum(subjectDemographics.sex_M0F1(group==1)))
x = table([sum(subjectDemographics.sex_M0F1(group==0)) sum(subjectDemographics.sex_M0F1(group==1)); ...
     sum(group==0)-sum(subjectDemographics.sex_M0F1(group==0)) sum(group==1)-sum(subjectDemographics.sex_M0F1(group==1))]);
[~, pValue] = fishertest(x)

fprintf('percent right-handed = %.1f\n', 100-(100*sum(subjectDemographics.hand_R0L1(group==0)/sum(group==0))))
fprintf('proportion right-handed: %iR/%iL\n', sum(group==0)-sum(subjectDemographics.hand_R0L1(group==0)), sum(subjectDemographics.hand_R0L1(group==0)))
fprintf('percent right-handed = %.1f\n', 100-(100*sum(subjectDemographics.hand_R0L1(group==1)/sum(group==1))))
fprintf('proportion right-handed: %iR/%iL\n', sum(group==1)-sum(subjectDemographics.hand_R0L1(group==1)), sum(subjectDemographics.hand_R0L1(group==1)))
x = table([sum(subjectDemographics.hand_R0L1(group==0)) sum(subjectDemographics.hand_R0L1(group==1)); ...
     sum(group==0)-sum(subjectDemographics.hand_R0L1(group==0)) sum(group==1)-sum(subjectDemographics.hand_R0L1(group==1))]);
[~, pValue] = fishertest(x)

fprintf('percent non-smokers = %.1f\n', 100*sum(subjectDemographics.smoke_Y0N1(group==0)/sum(group==0)))
fprintf('proportion smoke: %iN/%iY\n', sum(subjectDemographics.smoke_Y0N1(group==0)), sum(group==0)-sum(subjectDemographics.smoke_Y0N1(group==0)))
fprintf('percent non-smokers = %.1f\n', 100*sum(subjectDemographics.smoke_Y0N1(group==1)/sum(group==1)))
fprintf('proportion smoke: %iN/%iY\n', sum(subjectDemographics.smoke_Y0N1(group==1)), sum(group==1)-sum(subjectDemographics.smoke_Y0N1(group==1)))
x = table([sum(subjectDemographics.smoke_Y0N1(group==0)) sum(subjectDemographics.smoke_Y0N1(group==1)); ...
     sum(group==0)-sum(subjectDemographics.smoke_Y0N1(group==0)) sum(group==1)-sum(subjectDemographics.smoke_Y0N1(group==1))]);
[~, pValue] = fishertest(x)

fprintf('mean age (std) = %.1f(%.1f)\n', mean(subjectDemographics.AgeAtFMRI(group==0)), std(subjectDemographics.AgeAtFMRI(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.AgeAtFMRI(group==0)), max(subjectDemographics.AgeAtFMRI(group==0)))
fprintf('mean age (std) = %.1f(%.1f)\n', mean(subjectDemographics.AgeAtFMRI(group==1)), std(subjectDemographics.AgeAtFMRI(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.AgeAtFMRI(group==1)), max(subjectDemographics.AgeAtFMRI(group==1)))
pValue = ranksum(subjectDemographics.AgeAtFMRI(group==0), subjectDemographics.AgeAtFMRI(group==1))

fprintf('mean education (std) = %.1f(%.1f)\n', mean(subjectDemographics.educ(group==0)), std(subjectDemographics.educ(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.educ(group==0)), max(subjectDemographics.educ(group==0)))
fprintf('mean education (std) = %.1f(%.1f)\n', mean(subjectDemographics.educ(group==1)), std(subjectDemographics.educ(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.educ(group==1)), max(subjectDemographics.educ(group==1)))
[~, pValue] = ttest2(subjectDemographics.educ(group==0), subjectDemographics.educ(group==1))

fprintf('mean parental education (std) = %.1f(%.1f)\n', mean(subjectDemographics.par_ed(group==0), 'omitnan'), std(subjectDemographics.par_ed(group==0), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.par_ed(group==0)), max(subjectDemographics.par_ed(group==0)))
fprintf('mean parental education (std) = %.1f(%.1f)\n', mean(subjectDemographics.par_ed(group==1), 'omitnan'), std(subjectDemographics.par_ed(group==1), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.par_ed(group==1)), max(subjectDemographics.par_ed(group==1)))
[~, pValue] = ttest2(subjectDemographics.par_ed(group==0), subjectDemographics.par_ed(group==1))

fprintf('mean height = %.1f(%.1f)\n', mean(subjectDemographics.height(group==0)), std(subjectDemographics.height(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.height(group==0)), max(subjectDemographics.height(group==0)))
fprintf('mean height (std) = %.1f(%.1f)\n', mean(subjectDemographics.height(group==1)), std(subjectDemographics.height(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.height(group==1)), max(subjectDemographics.height(group==1)))
[~, pValue] = ttest2(subjectDemographics.height(group==0), subjectDemographics.height(group==1))

fprintf('mean weight = %.1f(%.1f)\n', mean(subjectDemographics.weight(group==0)), std(subjectDemographics.weight(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.weight(group==0)), max(subjectDemographics.weight(group==0)))
fprintf('mean weight (std) = %.1f(%.1f)\n', mean(subjectDemographics.weight(group==1)), std(subjectDemographics.weight(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.weight(group==1)), max(subjectDemographics.weight(group==1)))
[~, pValue] = ttest2(subjectDemographics.weight(group==0), subjectDemographics.weight(group==1))

fprintf('mean BMI = %.1f(%.1f)\n', mean(subjectDemographics.BMI(group==0)), std(subjectDemographics.BMI(group==0)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.BMI(group==0)), max(subjectDemographics.BMI(group==0)))
fprintf('mean BMI (std) = %.1f(%.1f)\n', mean(subjectDemographics.BMI(group==1)), std(subjectDemographics.BMI(group==1)))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.BMI(group==1)), max(subjectDemographics.BMI(group==1)))
[~, pValue] = ttest2(subjectDemographics.BMI(group==0), subjectDemographics.BMI(group==1))

fprintf('mean SISTOTAL = %.1f(%.1f)\n', mean(subjectDemographics.SISTOTAL(group==0), 'omitnan'), std(subjectDemographics.SISTOTAL(group==0), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.SISTOTAL(group==0)), max(subjectDemographics.SISTOTAL(group==0)))
fprintf('mean SISTOTAL (std) = %.1f(%.1f)\n', mean(subjectDemographics.SISTOTAL(group==1), 'omitnan'), std(subjectDemographics.SISTOTAL(group==1), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(subjectDemographics.SISTOTAL(group==1)), max(subjectDemographics.SISTOTAL(group==1)))
[~, pValue] = ttest2(subjectDemographics.SISTOTAL(group==0), subjectDemographics.SISTOTAL(group==1))

% using different table for STAI_TRAIT and alpraz_levels
allSubjectInfo_emotionID = readtable('../../data/Alpraz_emotionid.xlsx');
subjectIDs = subjectDemographics.bblid;

STAI_TRAIT = zeros(numel(subjectIDs), 1);
for i = 1:numel(subjectIDs)
    currentSubjectID = subjectIDs(i)
    STAI_TRAIT(i) = mean(allSubjectInfo_emotionID.STAI_TRAIT(allSubjectInfo_emotionID.bblid == currentSubjectID)); % averaging anxiety values over sessions
end
fprintf('mean STAI_TRAIT = %.1f(%.1f)\n', mean(STAI_TRAIT(group==0)), std(STAI_TRAIT(group==0)))
fprintf('range: %.1f-%.1f\n', min(STAI_TRAIT(group==0)), max(STAI_TRAIT(group==0)))
fprintf('mean STAI_TRAIT (std) = %.1f(%.1f)\n', mean(STAI_TRAIT(group==1)), std(STAI_TRAIT(group==1)))
fprintf('range: %.1f-%.1f\n', min(STAI_TRAIT(group==1)), max(STAI_TRAIT(group==1)))
pValue = ranksum(STAI_TRAIT(group==0), STAI_TRAIT(group==1))

drug = allSubjectInfo_emotionID.drug;
alpraz_levels = allSubjectInfo_emotionID.alpraz_levels(drug == 0); % using alpraz_levels in session where drug was administered
fprintf('mean alpraz_levels = %.1f(%.1f)\n', mean(alpraz_levels(group==0), 'omitnan'), std(alpraz_levels(group==0), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(alpraz_levels(group==0)), max(alpraz_levels(group==0)))
fprintf('mean alpraz_levels (std) = %.1f(%.1f)\n', mean(alpraz_levels(group==1), 'omitnan'), std(alpraz_levels(group==1), 'omitnan'))
fprintf('range: %.1f-%.1f\n', min(alpraz_levels(group==1)), max(alpraz_levels(group==1)))
[~, pValue] = ttest2(alpraz_levels(group==0), alpraz_levels(group==1))

%% Figure 1

%% plot structural matrix

figure; set(gcf, 'color', 'w'); hold on;
imagesc(A);
colormap('parula');
colorbar;
axis off;
ax = gca;
ax.FontSize = 12;

% h = refline(0, 234-109); h.Color = 'r';
% h = refline(0, 234-115); h.Color = 'r';
% h = refline(0, 234-227); h.Color = 'r';
% h = refline(0, 234-234); h.Color = 'r';

%% Figure depicting Lausanne parcellation

parcelValues = zeros(234, 1);
[~, ax, ph, ~] = fcn_lausannesurf(parcelValues, white, [-20 100]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');

%% plot brain activation maps w/ and w/o alpraz for individual subject

brainStates_alpraz = zeros(234, 1);
xf = allControlTrajectories_emotionid.xf{2};
parcelsToInclude_idx = allControlTrajectories_emotionid.parcelsToInclude_idx{1};
brainStates_alpraz(parcelsToInclude_idx) = xf;

[~, ax, ph, ~] = fcn_lausannesurf(brainStates_alpraz, redbluecmap);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');

%% Supplementary Figure 1 - plot average brain mask coverage
M = importdata('../../data/slabCoverage_combinedMaskStd_emotionid.txt');
slabCoverage = M.data';
slabCoverage = double(slabCoverage > parcelCoverageThresh);
[~, ax, ph, ~] = fcn_lausannesurf(slabCoverage, parula);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');

%% Figure 2

%% creating mixed models to examine effects of clinical and demographic variables on persistence_allNodes during emotionid

contrast = allControlEnergies_emotionid.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i}
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % fitting model after checking for normality
    persistence_allNodes = allControlEnergies_emotionid_currentContrast.persistence_allNodes;
    %h = kstest(zscore(persistence_allNodes));
    %figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'persistence_allNodes ~ SISTOTAL + STAI_TRAIT + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    save(strcat(resultsDir, 'mixedModel_emotionid_', currentContrast, '.mat'), 'mixedModel_emotionid_currentContrast');
    printFileName = strcat(resultsDir, 'mixedModel_emotionid_', currentContrast, '.xlsx');
    printLinearMixedModel(mixedModel_emotionid_currentContrast, printFileName);
end

%% creating mixed models to examine effects of clinical and demographic variables on persistence_allNodes during emotionrec

contrast = allControlEnergies_emotionrec.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    % fitting model
    persistence_allNodes = allControlEnergies_emotionrec_currentContrast.persistence_allNodes;
    %h = kstest(zscore(persistence_allNodes));
    %figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'persistence_allNodes ~ SISTOTAL + STAI_TRAIT + gender + age + (1|subjectID)', 'FitMethod', 'ML')
    save(strcat(resultsDir, 'mixedModel_emotionrec_', currentContrast, '.mat'), 'mixedModel_emotionrec_currentContrast');
    printFileName = strcat(resultsDir, 'mixedModel_emotionrec_', currentContrast, '.xlsx');
    printLinearMixedModel(mixedModel_emotionrec_currentContrast, printFileName);
end

%% comparing test statistics from mixed model against null models - emotionid

nIterations = 500;
contrast = allControlEnergies_emotionid.contrast;

% loading results of null model
resultsDir_nullModel = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_QA_spinTest/';
allControlEnergies_emotionid_nullModel = readtable(strcat(resultsDir_nullModel, 'allControlEnergies_emotionid.csv'));

for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % fitting model after checking for normality
    %persistence_allNodes = allControlEnergies_emotionid_currentContrast.persistence_allNodes;
    %h = kstest(zscore(persistence_allNodes));
    %figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + STAI_TRAIT*drug + SISTOTAL + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    [betas_emotionid_currentContrast, ~, ~] = fixedEffects(mixedModel_emotionid_currentContrast);
    
    betas_emotionid_nullModel_currentContrast = zeros(numel(betas_emotionid_currentContrast), nIterations);
    pValues_emotionid_nullModel_currentContrast = zeros(numel(betas_emotionid_currentContrast), 1);
    for j = 1:nIterations
        fprintf('iteration %d\n', j);
        allControlEnergies_emotionid_nullModel_currentContrast = allControlEnergies_emotionid_nullModel(strcmp(contrast, currentContrast), :);
        
        % converting categorical variables
        allControlEnergies_emotionid_nullModel_currentContrast.group = categorical(allControlEnergies_emotionid_nullModel_currentContrast.group);
        allControlEnergies_emotionid_nullModel_currentContrast.drug = categorical(allControlEnergies_emotionid_nullModel_currentContrast.drug);
        allControlEnergies_emotionid_nullModel_currentContrast.gender = categorical(allControlEnergies_emotionid_nullModel_currentContrast.gender);
        
        currentModelFormula = strcat('persistence_allNodes_', num2str(j), ' ~ drug + group + drug*group + STAI_TRAIT*drug + SISTOTAL + gender + age + (1|subjectID)');
        mixedModel_emotionid_nullModel_currentContrast = fitlme(allControlEnergies_emotionid_nullModel_currentContrast, currentModelFormula, 'FitMethod', 'ML');
        [betas, ~, ~] = fixedEffects(mixedModel_emotionid_nullModel_currentContrast);
        betas_emotionid_nullModel_currentContrast(:, j) = betas;
    end
    
    for j = 1:numel(betas_emotionid_currentContrast)
        pValues_emotionid_nullModel_currentContrast(i) = sum(abs(betas_emotionid_currentContrast(j)) > abs(betas_emotionid_nullModel_currentContrast(j, :)))/nIterations;
    end
    
    save(strcat(resultsDir, 'betas_emotionid_nullModel_', currentContrast, '.mat'), 'tStats_emotionid_nullModel_currentContrast');
    save(strcat(resultsDir, 'pValues_emotionid_nullModel_', currentContrast, '.mat'), 'pValues_emotionid_nullModel_currentContrast');
end

%% comparing test statistics from mixed model against null models - emotionrec

contrast = allControlEnergies_emotionrec.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i}
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(contrast, currentContrast), :);
    
    % converting categorical variables
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    % fitting model
    %persistence_allNodes = allControlEnergies_emotionrec_currentContrast.persistence_allNodes;
    %h = kstest(zscore(persistence_allNodes));
    %figure, histogram(persistence_allNodes); title(strcat('KS test: ', num2str(h)));
    mixedModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'persistence_allNodes ~ drug + group + drug*group + STAI_TRAIT*drug + SISTOTAL + gender + age + (1|subjectID)', 'FitMethod', 'REML')
    save(strcat(resultsDir, 'mixedModel_emotionrec_', currentContrast, '.mat'), 'mixedModel_emotionrec_currentContrast');
end

%% plot persistence during emotion id +/- Alpraz for controls and relatives; note that drug(0,1)=(alpraz,placebo)

drug = allControlEnergies_emotionid.drug;
group = allControlEnergies_emotionid.group;
contrast = allControlEnergies_emotionid.contrast;

f = figure; set(gcf, 'color', 'w');
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 6 5];

x = [allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))];

groups = [ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    2*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    3*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    4*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    5*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    6*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    7*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    8*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    9*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    10*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    11*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    12*ones(size(allControlEnergies_emotionid.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))))];

positions = [1, 1.2, 1.4, 1.6, 2, 2.2, 2.4, 2.6, 3, 3.2, 3.4, 3.6];

boxplot(x, groups, 'positions', positions, 'OutlierSize', 8, 'Symbol', 'k.');
set(gca, 'XTickLabel', {'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A'});

%set(gca,'xtick',[])
color = repmat({[0.4660, 0.6740, 0.1880], [0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560], [0.4940, 0.1840, 0.5560]}, 1, 3); % right to left
faceAlpha = repmat([0.75, 0.25], 1, 6); % right to left
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    currentFaceAlpha = faceAlpha(j);
    patch(get(h(j), 'XData'), get(h(j), 'YData'), color{j}, 'FaceAlpha', currentFaceAlpha);
end

ylabel('persistence energy')
%ylim([0.2, 0.65]);
set(gca, 'FontSize', 20);

%% plot persistence during emotion rec +/- Alpraz for controls and relatives; note that drug(0,1)=(alpraz,placebo)

drug = allControlEnergies_emotionrec.drug;
group = allControlEnergies_emotionrec.group;
contrast = allControlEnergies_emotionrec.contrast;

f = figure; set(gcf, 'color', 'w');
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 6 5];

x = [allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')); ...
    allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))];

groups = [ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    2*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    3*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    4*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast1_threatcorrectStd')))); ...
    5*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    6*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    7*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    8*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast3_nonthreatcorrectStd')))); ...
    9*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    10*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==0 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    11*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==1 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd')))); ...
    12*ones(size(allControlEnergies_emotionrec.persistence_allNodes(drug==0 & group==1 & strcmp(contrast, 'contrast5_neutralcorrectStd'))))];

positions = [1, 1.2, 1.4, 1.6, 2, 2.2, 2.4, 2.6, 3, 3.2, 3.4, 3.6];

boxplot(x, groups, 'positions', positions, 'OutlierSize', 8, 'Symbol', 'k.');
set(gca, 'XTickLabel', {'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A', 'P', 'A'});

%set(gca,'xtick',[])
color = repmat({[0.4660, 0.6740, 0.1880], [0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560], [0.4940, 0.1840, 0.5560]}, 1, 3); % right to left
faceAlpha = repmat([0.75, 0.25], 1, 6); % right to left
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    currentFaceAlpha = faceAlpha(j);
    patch(get(h(j), 'XData'), get(h(j), 'YData'), color{j}, 'FaceAlpha', currentFaceAlpha);
end

ylabel('persistence energy')
%ylim([0.2, 0.65]);
set(gca, 'FontSize', 20);

%% calculating and plotting average control impact across all subjects

%% emotionid
avgeControlImpact_emotionid_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    group = allControlTrajectories_emotionid_currentContrast.group;
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionid = controlInput_emotionid{i};
        controlInput_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    avgeControlImpact_emotionid_currentContrast_allSubjects = mean(controlInput_emotionid_allSubjects); % averaging over all subjects
    avgeControlImpact_emotionid_currentContrast_allSubjects(isnan(avgeControlImpact_emotionid_currentContrast_allSubjects)) = 0; % setting NaN values (outside slab) to 0
    avgeControlImpact_emotionid_allSubjects(:, c) = avgeControlImpact_emotionid_currentContrast_allSubjects;
    
    emotionid_controlImpact_nodeNames = cell(234, 2);
    for i = 1:234
        emotionid_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionid_allSubjects(i, c);
        emotionid_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
    end
    
    emotionid_controlImpact_nodeNames = sortrows(emotionid_controlImpact_nodeNames, 'descend'); % sorting nodes by descending control impact value
    emotionid_controlImpact_nodeNames = cell2table(emotionid_controlImpact_nodeNames, 'VariableNames', {'controlImpact_sorted', 'nodeName_Lausanne'});
    writetable(emotionid_controlImpact_nodeNames, strcat(resultsDir, 'emotionid_controlImpact_nodeNames_', currentContrast, '.csv'));
end

avgeControlImpact_emotionid_allSubjects = mean(avgeControlImpact_emotionid_allSubjects, 2); % averaging over all contrasts

% emotionid_controlImpact_nodeNames = cell(234, 2);
% for i = 1:234
%     emotionid_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionid_allSubjects(i);
%     emotionid_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
% end
% 
% emotionid_controlImpact_nodeNames = sortrows(emotionid_controlImpact_nodeNames, 'descend'); % sorting nodes by descending control impact value
% emotionid_controlImpact_nodeNames = cell2table(emotionid_controlImpact_nodeNames, 'VariableNames', {'controlImpact_sorted', 'nodeName_Lausanne'});
% writetable(emotionid_controlImpact_nodeNames, strcat(resultsDir, 'emotionid_controlImpact_nodeNames.csv'));

% [f, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionid_allSubjects, redbluecmap, [0 1.5]);
% view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
% axis(ax(1), 'off'); axis(ax(2), 'off');
% savePath = strcat(resultsDir, 'emotionid_controlImpactPlots_cortex1.svg'); saveas(f(1), savePath); close(f(1));
% savePath = strcat(resultsDir, 'emotionid_controlImpactPlots_cortex2.svg'); saveas(f(2), savePath); close(f(2));
% 
% [f, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionid_allSubjects, redbluecmap, [0 1.5]);
% view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
% view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
% axis(ax(1), 'off'); axis(ax(2), 'off');
% savePath = strcat(resultsDir, 'emotionid_controlImpactPlots_cortex3.svg'); saveas(f(1), savePath); close(f(1));
% savePath = strcat(resultsDir, 'emotionid_controlImpactPlots_cortex4.svg'); saveas(f(2), savePath); close(f(2));
% 
% figure; [f1, f2] = plot_subcortvol(avgeControlImpact_emotionid_allSubjects(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, redbluecmap, 0, 1.5);
% savePath1 = strcat(resultsDir, 'emotionid_controlImpactPlots_subcortex1.svg'); saveas(f1, savePath1); 
% savePath2 = strcat(resultsDir, 'emotionid_controlImpactPlots_subcortex2.svg'); saveas(f2, savePath2); 
% close(f1); close(f2);

%% emotionrec

avgeControlImpact_emotionrec_allSubjects = zeros(nNodes, nContrasts);

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlImpact_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = controlInput_emotionrec{i};
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    avgeControlImpact_emotionrec_currentContrast_allSubjects = mean(controlInput_emotionrec_allSubjects);
    avgeControlImpact_emotionrec_currentContrast_allSubjects(isnan(avgeControlImpact_emotionrec_currentContrast_allSubjects)) = 0;
    
    avgeControlImpact_emotionrec_allSubjects(:, c) = avgeControlImpact_emotionrec_currentContrast_allSubjects;
end

avgeControlImpact_emotionrec_allSubjects = mean(avgeControlImpact_emotionrec_allSubjects, 2);

emotionrec_controlImpact_nodeNames = cell(234, 2);
for i = 1:234
    emotionrec_controlImpact_nodeNames{i, 1} = avgeControlImpact_emotionrec_allSubjects(i);
    emotionrec_controlImpact_nodeNames{i, 2} = LausanneParcelNames{i};
end

emotionrec_controlImpact_nodeNames = sortrows(emotionrec_controlImpact_nodeNames, 'descend');
emotionrec_controlImpact_nodeNames = cell2table(emotionrec_controlImpact_nodeNames, 'VariableNames', {'controlImpact_sorted', 'nodeName_Lausanne'});
writetable(emotionrec_controlImpact_nodeNames, strcat(resultsDir, 'emotionrec_controlImpact_nodeNames.csv'));

[f, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionrec_allSubjects, redbluecmap, [0 1.5]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, 'emotionrec_controlImpactPlots_cortex1.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, 'emotionrec_controlImpactPlots_cortex2.svg'); saveas(f(2), savePath); close(f(2));

[f, ax, ph, ~] = fcn_lausannesurf(avgeControlImpact_emotionrec_allSubjects, redbluecmap, [0 1.5]);
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, 'emotionrec_controlImpactPlots_cortex3.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, 'emotionrec_controlImpactPlots_cortex4.svg'); saveas(f(2), savePath); close(f(2));

figure; [f1, f2] = plot_subcortvol(avgeControlImpact_emotionrec_allSubjects(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, redbluecmap, 0, 1.5);
savePath1 = strcat(resultsDir, 'emotionrec_controlImpactPlots_subcortex1.svg'); saveas(f1, savePath1); 
savePath2 = strcat(resultsDir, 'emotionrec_controlImpactPlots_subcortex2.svg'); saveas(f2, savePath2); 
close(f1); close(f2);

%% Figure 3

%% plotting SISTOTAL against persistence_allNodes for nonthreat emotionid

contrast = allControlEnergies_emotionid.contrast;
allControlEnergies_emotionid_nonthreat = allControlEnergies_emotionid(strcmp(contrast, 'contrast3_nonthreatcorrectStd'), :);
age = allControlEnergies_emotionid_nonthreat.age;
gender = allControlEnergies_emotionid_nonthreat.gender;
persistence_allNodes_nonthreat = allControlEnergies_emotionid_nonthreat.persistence_allNodes;
SISTOTAL = allControlEnergies_emotionid_nonthreat.SISTOTAL;

figure; set(gcf, 'color', 'white');
plot(SISTOTAL, persistence_allNodes_nonthreat, 'k.', 'MarkerSize', 20);
h = lsline; h.LineWidth = 2;
xlabel('SISTOTAL');
ylabel('persistence energy - nonthreat');
set(gca, 'FontSize', 20);
% [rho, pValue] = partialcorr(SISTOTAL, persistence_allNodes_nonthreat, [age gender], 'Rows', 'Complete');
% text(30, 0.25, strcat('\rho=', num2str(rho, 2)), 'FontSize', 14);
% text(30, 0.225, strcat('p=', num2str(pValue, 2)), 'FontSize', 14);
ylim([0.2, 0.65]);
box off;

%% plotting STAI_TRAIT against persistence_allNodes for neutral emotionid

allControlEnergies_emotionid_neutral = allControlEnergies_emotionid(strcmp(contrast, 'contrast5_neutralcorrectStd'), :);
persistence_allNodes_neutral = allControlEnergies_emotionid_neutral.persistence_allNodes;
age = allControlEnergies_emotionid_neutral.age;
gender = allControlEnergies_emotionid_neutral.gender;
drug = allControlEnergies_emotionid_neutral.drug;
STAI_TRAIT = allControlEnergies_emotionid_neutral.STAI_TRAIT;

figure; set(gcf, 'color', 'w'); set(gca, 'FontSize', 20); hold on;
plot(STAI_TRAIT(drug==1), persistence_allNodes_neutral(drug==1), 'b.', 'MarkerSize', 20); % placebo
plot(STAI_TRAIT(drug==0), persistence_allNodes_neutral(drug==0), 'r.', 'MarkerSize', 20); % alprazolam
lsline;
h = lsline; h(1).LineWidth = 2; h(2).LineWidth = 2;
xlabel('trait anxiety');
ylabel('persistence energy - neutral');
legend('placebo', 'alprazolam', 'location', 'northeast'); legend boxoff;
% [rho, pValue] = partialcorr(STAI_TRAIT(drug==0), persistence_allNodes_neutral(drug==0), [age(drug==0), gender(drug==0)], 'Rows', 'Complete');
% text(60, 0.25, strcat('\rho=', num2str(rho, 2)), 'Color', 'b', 'FontSize', 14);
% text(60, 0.225, strcat('p=', num2str(pValue, 2)), 'Color', 'b', 'FontSize', 14);
% [rho, pValue] = partialcorr(STAI_TRAIT(drug==1), persistence_allNodes_neutral(drug==1), [age(drug==1), gender(drug==1)], 'Rows', 'Complete');
% text(60, 0.30, strcat('\rho=', num2str(rho, 2)), 'Color', 'r', 'FontSize', 14);
% text(60, 0.275, strcat('p=', num2str(pValue, 2)), 'Color', 'r', 'FontSize', 14);
ylim([0.2, 0.65]);

%% mixed models for accuracy and reaction time versus persistence

% emotionid
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(allControlEnergies_emotionid.contrast, currentContrast), :);
    
    allControlEnergies_emotionid_currentContrast.drug = categorical(allControlEnergies_emotionid_currentContrast.drug);
    allControlEnergies_emotionid_currentContrast.group = categorical(allControlEnergies_emotionid_currentContrast.group);
    allControlEnergies_emotionid_currentContrast.gender = categorical(allControlEnergies_emotionid_currentContrast.gender);
    
    % calculating efficiency = accuracy/reactionTime
    efficiency_corrthreat = allControlEnergies_emotionid_currentContrast.pctcorr_threat./allControlEnergies_emotionid_currentContrast.rtmdn_threatcorr;
    efficiency_corrnonthreat = allControlEnergies_emotionid_currentContrast.pctcorr_nonthreat./allControlEnergies_emotionid_currentContrast.rtmdn_nonthreatcorr;
    efficiency_corrneutral = allControlEnergies_emotionid_currentContrast.pctcorr_neutral./allControlEnergies_emotionid_currentContrast.rtmdn_neutralcorr;
    
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrthreat);
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrnonthreat);
    allControlEnergies_emotionid_currentContrast = addvars(allControlEnergies_emotionid_currentContrast, efficiency_corrneutral);
    
    switch currentContrast
        case 'contrast1_threatcorrectStd'
            %figure, histogram(allControlEnergies_emotionid_currentContrast.pctcorr_threat); title('emotionid - pctcorr_threat');
            %figure, histogram(allControlEnergies_emotionid_currentContrast.rtmdn_threatcorr); title('emotionid - rtmdn_threatcorr');
            figure, histogram(allControlEnergies_emotionid_currentContrast.efficiency_corrthreat); title('emotionid - efficiency_corrthreat');
            accuracyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'pctcorr_threat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            reactionTimeModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'rtmdn_threatcorr ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast3_nonthreatcorrectStd'
            %figure, histogram(allControlEnergies_emotionid_currentContrast.pctcorr_nonthreat); title('emotionid - pctcorr_nonthreat');
            %figure, histogram(allControlEnergies_emotionid_currentContrast.rtmdn_nonthreatcorr); title('emotionid - rtmdn_nonthreatcorr');
            figure, histogram(allControlEnergies_emotionid_currentContrast.efficiency_corrnonthreat); title('emotionid - efficiency_corrnonthreat');
            accuracyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'pctcorr_nonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            reactionTimeModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'rtmdn_nonthreatcorr ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrnonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast5_neutralcorrectStd'
            %figure, histogram(allControlEnergies_emotionid_currentContrast.pctcorr_neutral); title('emotionid - pctcorr_neutral');
            %figure, histogram(allControlEnergies_emotionid_currentContrast.rtmdn_neutralcorr); title('emotionid - rtmdn_neutralcorr');
            figure, histogram(allControlEnergies_emotionid_currentContrast.efficiency_corrneutral); title('emotionid - efficiency_corrneutral');
            accuracyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'pctcorr_neutral ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            reactionTimeModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'rtmdn_neutralcorr ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            efficiencyModel_emotionid_currentContrast = fitlme(allControlEnergies_emotionid_currentContrast, 'efficiency_corrneutral ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    end
    
    save(strcat(resultsDir, 'accuracyModel_emotionid_', currentContrast, '.mat'), 'accuracyModel_emotionid_currentContrast');
    save(strcat(resultsDir, 'reactionTimeModel_emotionid_', currentContrast, '.mat'), 'reactionTimeModel_emotionid_currentContrast');
    save(strcat(resultsDir, 'efficiencyModel_emotionid_', currentContrast, '.mat'), 'efficiencyModel_emotionid_currentContrast');
end

% emotionrec
for i = 1:nContrasts
    currentContrast = contrastLabels{i};
    allControlEnergies_emotionrec_currentContrast = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, currentContrast), :);
    
    allControlEnergies_emotionrec_currentContrast.drug = categorical(allControlEnergies_emotionrec_currentContrast.drug);
    allControlEnergies_emotionrec_currentContrast.group = categorical(allControlEnergies_emotionrec_currentContrast.group);
    allControlEnergies_emotionrec_currentContrast.gender = categorical(allControlEnergies_emotionrec_currentContrast.gender);
    
    % calculating efficiency = accuracy/reactionTime
    efficiency_corrthreat = allControlEnergies_emotionrec_currentContrast.pctcorr_threat./allControlEnergies_emotionrec_currentContrast.rtmdn_threatcorr;
    efficiency_corrnonthreat = allControlEnergies_emotionrec_currentContrast.pctcorr_nonthreat./allControlEnergies_emotionrec_currentContrast.rtmdn_nonthreatcorr;
    efficiency_corrneutral = allControlEnergies_emotionrec_currentContrast.pctcorr_neutral./allControlEnergies_emotionrec_currentContrast.rtmdn_neutralcorr;
    
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrthreat);
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrnonthreat);
    allControlEnergies_emotionrec_currentContrast = addvars(allControlEnergies_emotionrec_currentContrast, efficiency_corrneutral);

    switch currentContrast
        case 'contrast1_threatcorrectStd'
            %figure, histogram(allControlEnergies_emotionrec_currentContrast.pctcorr_threat); title('emotionrec - pctcorr_threat');
            %figure, histogram(allControlEnergies_emotionrec_currentContrast.rtmdn_threatcorr); title('emotionrec - rtmdn_threatcorr');
            figure, histogram(allControlEnergies_emotionrec_currentContrast.efficiency_corrthreat); title('emotionrec - efficiency_corrthreat');
            accuracyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'pctcorr_threat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            reactionTimeModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'rtmdn_threatcorr ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            efficiencyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'efficiency_corrthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast3_nonthreatcorrectStd'
            %figure, histogram(allControlEnergies_emotionrec_currentContrast.pctcorr_nonthreat); title('emotionrec - pctcorr_nonthreat');
            %figure, histogram(allControlEnergies_emotionrec_currentContrast.rtmdn_nonthreatcorr); title('emotionrec - rtmdn_nonthreatcorr');            
            figure, histogram(allControlEnergies_emotionrec_currentContrast.efficiency_corrnonthreat); title('emotionrec - efficiency_corrnonthreat');
            accuracyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'pctcorr_nonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            reactionTimeModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'rtmdn_nonthreatcorr ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            efficiencyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'efficiency_corrnonthreat ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
        case 'contrast5_neutralcorrectStd'
            %figure, histogram(allControlEnergies_emotionrec_currentContrast.pctcorr_neutral); title('emotionrec - pctcorr_neutral');
            %figure, histogram(allControlEnergies_emotionrec_currentContrast.rtmdn_neutralcorr); title('emotionrec - rtmdn_neutralcorr');
            figure, histogram(allControlEnergies_emotionrec_currentContrast.efficiency_corrneutral); title('emotionrec - efficiency_corrneutral');
            accuracyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'pctcorr_neutral ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            reactionTimeModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'rtmdn_neutralcorr ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
            efficiencyModel_emotionrec_currentContrast = fitlme(allControlEnergies_emotionrec_currentContrast, 'efficiency_corrneutral ~ persistence_allNodes + group + drug + drug*group + gender + age + (1|subjectID)', 'FitMethod', 'ML');
    end
    
    save(strcat(resultsDir, 'accuracyModel_emotionrec_', currentContrast, '.mat'), 'accuracyModel_emotionrec_currentContrast');
    save(strcat(resultsDir, 'reactionTimeModel_emotionrec_', currentContrast, '.mat'), 'reactionTimeModel_emotionrec_currentContrast');
    save(strcat(resultsDir, 'efficiencyModel_emotionrec_', currentContrast, '.mat'), 'efficiencyModel_emotionrec_currentContrast');
end

%% plot task efficiency versus persistence energy for all contrasts in emotion ID

for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionrec(strcmp(allControlEnergies_emotionrec.contrast, currentContrast), :);
    
    contrastString = currentContrast(11:end); contrastString = strrep(contrastString, 'correctStd', ''); 
    accuracy_contrastString = strcat('pctcorr_', contrastString); rtmdn_contrastString = strcat('rtmdn_', contrastString, 'corr');
    accuracy = allControlEnergies_emotionid_currentContrast.(accuracy_contrastString);
    rtmdn = allControlEnergies_emotionid_currentContrast.(rtmdn_contrastString);
    efficiency = accuracy./rtmdn;
    persistence = allControlEnergies_emotionid_currentContrast.persistence_allNodes;
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'w'); hold on;
    set(gca, 'FontSize', 20);
    plot(persistence, efficiency, 'k.', 'MarkerSize', 20);
    h = lsline; h.LineWidth = 2; 
    xlabel('persistence energy');
    ylabel('efficiency');
    xlim([0.2, 0.65]); %ylim([0.2, 1]);
    
    savePath = strcat(resultsDir, 'emotionrec_efficiency_persistenceEnergy_', currentContrast, '.eps'); saveas(f, savePath); close(f);
end

%% Figure 4

%% import AHBA gene expression, PET neurotransmitter profiles

% calculate GABA gene expression from Allen human brain atlas data
load('/Users/ArunMahadevan/Documents/BBL/studies/alpraz/GeneExpression/ParcellatedGeneExpressionLausanne125.mat');
GABRA1 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA1'));
GABRA2 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA2'));
GABRA3 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA3'));
GABRA5 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRA5'));
GABRB1 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRB1'));
GABRB2 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRB2'));
GABRB3 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRB3'));
GABRG2 = lausanneParcelExpression(:, strcmp(gene_names, 'GABRG2'));

PETlabels = {'5HT1a', '5HT1b', '5HT2a', 'D1', 'D2', 'DAT', 'FDOPA', 'GABAa', 'NAT', 'SERT'};

% load PET data
PETdir = '/Users/ArunMahadevan/Documents/BBL/studies/alpraz/PETatlas/';
X = importdata(strcat(PETdir, '5HT1a_WAY_HC36.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT1a_WAY_HC36 = X.data';
X = importdata(strcat(PETdir, '5HT1b_P943_HC22.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT1b_P943_HC22 = X.data';
X = importdata(strcat(PETdir, '5HT2a_ALT_HC19.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_5HT2a_ALT_HC19 = X.data';
X = importdata(strcat(PETdir, 'D1_SCH23390_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_D1_SCH23390_c11 = X.data';
X = importdata(strcat(PETdir, 'D2_RACLOPRIDE_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_D2_RACLOPRIDE_c11 = X.data';
X = importdata(strcat(PETdir, 'DAT_DATSPECT.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_DAT_DATSPECT = X.data';
X = importdata(strcat(PETdir, 'FDOPA_f18.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_FDOPA_f18 = X.data';
X = importdata(strcat(PETdir, 'GABAa_FLUMAZENIL_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_GABAa_FLUMAZENIL_c11 = X.data';
X = importdata(strcat(PETdir, 'NAT_MRB_c11.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_NAT_MRB_c11 = X.data';
X = importdata(strcat(PETdir, 'SERT_DASB_HC30.nii_lausanne_ROIv_scale125_dilated_resized_3mm.txt'));
PET_SERT_DASB_HC30 = X.data';

%% plot matrix of PET maps arranged by Yeo7 systems

PET_maps = [PET_5HT1a_WAY_HC36 PET_5HT1b_P943_HC22 PET_5HT2a_ALT_HC19 PET_D1_SCH23390_c11 PET_D2_RACLOPRIDE_c11 PET_DAT_DATSPECT PET_FDOPA_f18 PET_GABAa_FLUMAZENIL_c11 PET_NAT_MRB_c11 PET_SERT_DASB_HC30];
PET_maps_Yeo7 = rearrangeMatrix_Yeo7(PET_maps);
nMaps = size(PET_maps_Yeo7, 2);
f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
imagesc(PET_maps_Yeo7); cb = colorbar; cb.Location = 'westoutside';
set(gca, 'FontSize', 20);
xlim([0.5, nMaps+0.5]); ylim([0, nNodes]);
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(7), refLines(8)]), 'Subcortical');

saveas(f, strcat(resultsDir, 'PET_allMaps_Yeo7.svg'));

%% plot GABAa spatial maps

[f, ax, ph, ~] = fcn_lausannesurf(PET_GABAa_FLUMAZENIL_c11, redbluecmap, [0 100]);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, 'PET_GABAa_FLUMAZENIL_c11_cortex1.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, 'PET_GABAa_FLUMAZENIL_c11_cortex2.svg'); saveas(f(2), savePath); close(f(2));

[f, ax, ph, ~] = fcn_lausannesurf(PET_GABAa_FLUMAZENIL_c11, redbluecmap, [0 100]);
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, 'PET_GABAa_FLUMAZENIL_c11_cortex3.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, 'PET_GABAa_FLUMAZENIL_c11_cortex4.svg'); saveas(f(2), savePath); close(f(2));

figure; [f1, f2] = plot_subcortvol(PET_GABAa_FLUMAZENIL_c11(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, redbluecmap, 0, 100);
savePath1 = strcat(resultsDir, 'PET_GABAa_FLUMAZENIL_c11_subcortex1.svg'); saveas(f1, savePath1); 
savePath2 = strcat(resultsDir, 'PET_GABAa_FLUMAZENIL_c11_subcortex2.svg'); saveas(f2, savePath2); 
close(f1); close(f2);

%% average control input difference w/ and w/o drug vs PET neurotransmitter maps - emotionid

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    drug = allControlTrajectories_emotionid_currentContrast.drug;
    group = allControlTrajectories_emotionid_currentContrast.group;
    
    % populate matrix of control input for [nSubjects*2 x nNodes]
    controlInputs_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInputs_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    controlInputDiff_emotionid_currentContrast_drug = abs(controlInputs_emotionid_allSubjects(drug==0, :) - controlInputs_emotionid_allSubjects(drug==1, :))'; % calculating absolute value of difference in control input w/ and w/o drug
    controlInputDiff_emotionid_currentContrast_drug_Yeo7 = rearrangeMatrix_Yeo7(controlInputDiff_emotionid_currentContrast_drug);
    
    nSubjects = size(controlInputDiff_emotionid_currentContrast_drug_Yeo7, 2);
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    imagesc(controlInputDiff_emotionid_currentContrast_drug_Yeo7); cb = colorbar; cb.Location = 'westoutside';
    set(gca, 'FontSize', 20);
    xlim([1, nSubjects]); ylim([1, nNodes]);
    refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
    h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([0, refLines(1)]), 'Visual');
    h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(1), refLines(2)]), 'Somatomator');
    h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
    h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(3), refLines(4)]), 'Ventral Attention');
    h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(4), refLines(5)]), 'Limbic');
    h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
    h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(6), refLines(7)]), 'Default Mode');
    h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(7), refLines(8)]), 'Subcortical');
    
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDiff_heatmap_drug_vs_noDrug.svg'));
    
    % compute correlations between control input differences w/ and w/o drug against PET maps
    controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs = [corr(controlInputDiff_emotionid_currentContrast_drug, PET_5HT1a_WAY_HC36, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_5HT1b_P943_HC22, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_5HT2a_ALT_HC19, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_D1_SCH23390_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_D2_RACLOPRIDE_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_DAT_DATSPECT, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_FDOPA_f18, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_NAT_MRB_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionid_currentContrast_drug, PET_SERT_DASB_HC30, 'Type', 'Spearman', 'Rows', 'Complete')];
    
    controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs = atanh(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % Fisher z-transform
    %[h, pValues] = ttest(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % using one-sample t-test to check for significance of subject-wise correlations
    [pValues, ~, ~, ~, ~] = mult_comp_perm_t1(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs);
    h = double((pValues < 0.05));
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    ylim([-0.4, 0.4]);
    
    boxplot(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs, 'Labels', PETlabels, 'LabelOrientation', 'inline', 'OutlierSize', 8, 'Symbol', 'k.');
    refline(0, 0);
    for i = 1:numel(h)
        if h(i) == 1
            text(i, 0, '*', 'Color', 'r', 'FontSize', 20);
        end
    end
    
    ylabel('Fisher z (Spearman \rho)');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDiff_drug_vs_noDrug_PETatlases.eps'));
    
end

%% average control input w/ and w/o drug vs PET neurotransmitter maps - emotionrec

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    drug = allControlTrajectories_emotionrec_currentContrast.drug;
    group = allControlTrajectories_emotionrec_currentContrast.group;
    
    % populate matrix of control input for [nSubjects*2 x nNodes]
    controlInputs_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps; % calculating total control input for current node as the sum of integral of time-varying control input
        controlInputs_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    controlInputDiff_emotionrec_currentContrast_drug = abs(controlInputs_emotionrec_allSubjects(drug==0, :) - controlInputs_emotionrec_allSubjects(drug==1, :))'; % calculating absolute value of difference in control input w/ and w/o drug
    controlInputDiff_emotionrec_currentContrast_drug_Yeo7 = rearrangeMatrix_Yeo7(controlInputDiff_emotionrec_currentContrast_drug);
    
    nSubjects = size(controlInputDiff_emotionrec_currentContrast_drug_Yeo7, 2);
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    imagesc(controlInputDiff_emotionrec_currentContrast_drug_Yeo7); cb = colorbar; cb.Location = 'westoutside';
    set(gca, 'FontSize', 20);
    xlim([1, nSubjects]); ylim([1, nNodes]);
    refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
    h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([0, refLines(1)]), 'Visual');
    h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(1), refLines(2)]), 'Somatomator');
    h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
    h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(3), refLines(4)]), 'Ventral Attention');
    h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(4), refLines(5)]), 'Limbic');
    h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
    h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(6), refLines(7)]), 'Default Mode');
    h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nSubjects, mean([refLines(7), refLines(8)]), 'Subcortical');
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDiff_heatmap_drug_vs_noDrug.svg'));
    
    % compute correlations between control input differences w/ and w/o drug against PET maps
    controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs = [corr(controlInputDiff_emotionrec_currentContrast_drug, PET_5HT1a_WAY_HC36, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_5HT1b_P943_HC22, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_5HT2a_ALT_HC19, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_D1_SCH23390_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_D2_RACLOPRIDE_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_DAT_DATSPECT, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_FDOPA_f18, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_NAT_MRB_c11, 'Type', 'Spearman', 'Rows', 'Complete'), ...
        corr(controlInputDiff_emotionrec_currentContrast_drug, PET_SERT_DASB_HC30, 'Type', 'Spearman', 'Rows', 'Complete')];
    
    controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs = atanh(controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs); % Fisher z-transform
    %[h, pValues] = ttest(controlInputDiff_emotionid_currentContrast_drug_PETatlasCorrs); % using one-sample t-test to check for significance of subject-wise correlations
    [pValues, ~, ~, ~, ~] = mult_comp_perm_t1(controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs);
    h = double((pValues < 0.05));
    
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    ylim([-0.4, 0.4]);
    
    boxplot(controlInputDiff_emotionrec_currentContrast_drug_PETatlasCorrs, 'Labels', PETlabels, 'LabelOrientation', 'inline', 'OutlierSize', 8, 'Symbol', 'k.');
    refline(0, 0);
    for i = 1:numel(h)
        if h(i) == 1
            text(i, 0, '*', 'Color', 'r', 'FontSize', 20);
        end
    end
    
    ylabel('Fisher z (Spearman \rho)');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDiff_drug_vs_noDrug_PETatlases.eps'));
    
end

%% plot control input beta coefficients vs PET neurotransmitter maps - emotionid

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionid_currentContrast = allControlTrajectories_emotionid(strcmp(allControlTrajectories_emotionid.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionid = allControlTrajectories_emotionid_currentContrast.controlInputs_persistence; % extracting control input
    parcelsToInclude_emotionid = allControlTrajectories_emotionid_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionid);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionid_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionid{i};
        current_controlInput_emotionid = trapz(controlInput_emotionid{i}.^2)/nTimeSteps;
        controlInput_emotionid_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionid;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    % looping over nodes, building a model for each
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = abs(controlInput_emotionid_allSubjects(:, i)); % taking absolute value of total control input for current node
        controlInput_currentNode = zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % skipping node if one or more subjects does not have parcel coverage in slab
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionid_currentContrast, controlInput_currentNode); % creating table with additional variable
            allControlTrajectories_currentNode.drug = 1 - allControlTrajectories_currentNode.drug; % NOTE: inverting drug indicator to drug(0, 1) = (placebo, alpraz) for ease of interpretation
            % converting to categorical variables
            allControlTrajectories_currentNode.drug = categorical(allControlTrajectories_currentNode.drug);
            allControlTrajectories_currentNode.group = categorical(allControlTrajectories_currentNode.group);
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'ML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group_1')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug_1')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group_1:drug_1'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group_1')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug_1')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group_1:drug_1'));
        end
    end
    
    %betas_groupMain = abs(betas_groupMain); betas_drugMain = abs(betas_drugMain); betas_groupDrugInteraction = abs(betas_groupDrugInteraction); % taking absolute values of coefficients
    
    %%%%% plotting coefficients of drug main effect versus PET neurotransmitter profiles %%%%%
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDrugEffect_PET_GABAa_FLUMAZENIL_c11.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction (z-score)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputGroupDrugInteraction_PET_GABAa_FLUMAZENIL_c11.eps'));
    
    %%%% plotting coefficients of drug main effect versus GABA receptor expression from Allen atlas %%%%
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_drugMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputDrugEffect_AHBA_GABA.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_groupDrugInteraction, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of group x drug interaction (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionid_', currentContrast, '_controlInputGroupDrugInteraction_AHBA_GABA.eps'));
end

%% plot control input beta coefficients vs PET neurotransmitter maps - emotionrec

nTimeSteps = 1001;
for c = 1:nContrasts
    currentContrast = contrastLabels{c};
    allControlTrajectories_emotionrec_currentContrast = allControlTrajectories_emotionrec(strcmp(allControlTrajectories_emotionrec.contrast, currentContrast), :); % extracting table for current contrast
    
    controlInput_emotionrec = allControlTrajectories_emotionrec_currentContrast.controlInputs_persistence; % extracting control impact
    parcelsToInclude_emotionrec = allControlTrajectories_emotionrec_currentContrast.parcelsToInclude_idx; % extracting parcel indices in imaging slab
    nIterations = numel(controlInput_emotionrec);
    
    % populate matrix of control impact for [nSubjects*2 x nNodes]
    controlInput_emotionrec_allSubjects = NaN(nIterations, nNodes);
    for i = 1:nIterations
        current_parcelsToInclude_idx = parcelsToInclude_emotionrec{i};
        current_controlInput_emotionrec = trapz(controlInput_emotionrec{i}.^2)/nTimeSteps;
        controlInput_emotionrec_allSubjects(i, current_parcelsToInclude_idx) = current_controlInput_emotionrec;
    end
    
    % evaluate effects of drug and group on control impact of each node
    betas_groupMain = NaN(nNodes, 1); pValues_groupMain = NaN(nNodes, 1);
    betas_drugMain = NaN(nNodes, 1); pValues_drugMain = NaN(nNodes, 1);
    betas_groupDrugInteraction = NaN(nNodes, 1); pValues_groupDrugInteraction = NaN(nNodes, 1);
    
    for i = 1:nNodes
        %fprintf('node%d\n', i);
        controlInput_currentNode = abs(controlInput_emotionrec_allSubjects(:, i));
        controlInput_currentNode =  zscore(controlInput_currentNode);
        if sum(isnan(controlInput_currentNode)) > 0 % skipping node if one or more subjects does not have parcel coverage in slab
            continue;
        else
            allControlTrajectories_currentNode = addvars(allControlTrajectories_emotionrec_currentContrast, controlInput_currentNode); % creating table with additional variable
            allControlTrajectories_currentNode.drug = 1 - allControlTrajectories_currentNode.drug; % NOTE: inverting drug indicator to drug(0, 1) = (placebo, alpraz) for ease of interpretation
            % converting to categorical variables
            allControlTrajectories_currentNode.drug = categorical(allControlTrajectories_currentNode.drug);
            allControlTrajectories_currentNode.group = categorical(allControlTrajectories_currentNode.group);
            controlInputModel = fitlme(allControlTrajectories_currentNode, 'controlInput_currentNode ~ drug + group + drug*group + (1|subjectID)', 'FitMethod', 'ML'); % fitting mixed model
            coefficientNames = controlInputModel.CoefficientNames; coefficients = controlInputModel.Coefficients.Estimate; pValues = controlInputModel.Coefficients.pValue; % extracting coefficients and p-values
            betas_groupMain(i) = coefficients(strcmp(coefficientNames, 'group_1')); betas_drugMain(i) = coefficients(strcmp(coefficientNames, 'drug_1')); betas_groupDrugInteraction(i) = coefficients(strcmp(coefficientNames, 'group_1:drug_1'));
            pValues_groupMain(i) = pValues(strcmp(coefficientNames, 'group_1')); pValues_drugMain(i) = pValues(strcmp(coefficientNames, 'drug_1')); pValues_groupDrugInteraction(i) = pValues(strcmp(coefficientNames, 'group_1:drug_1'));
        end
    end
    
    %betas_groupMain = abs(betas_groupMain); betas_drugMain = abs(betas_drugMain); betas_groupDrugInteraction = abs(betas_groupDrugInteraction); % taking absolute values of coefficients
    
%     pValues_groupMain_FDR = mafdr(pValues_groupMain, 'BHFDR', true);
%     pValues_drugMain_FDR = mafdr(pValues_drugMain, 'BHFDR', true);
%     pValues_groupDrugInteraction_FDR = mafdr(pValues_groupDrugInteraction, 'BHFDR', true);
%     
%     sum(pValues_groupMain_FDR < 0.05/234)
%     sum(pValues_drugMain_FDR < 0.05/234)
%     sum(pValues_groupDrugInteraction_FDR < 0.05/234)
    
    % plotting coefficients of group main effect versus PET
    % neurotransmitter profiles
%     f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
%     set(gca, 'FontSize', 20);
%     
%     xlim([0, 1]);
%     ylim([10, 80]);
%     
%     plot(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
%     [rho, pVal_corr] = corr(betas_groupMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
%     text(0.1, 78, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
%     
%     lsline;
%     xlabel('coefficients of group main effect');
%     ylabel('GABA receptor density');
%     saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupEffect_GABA.eps'));
    
    % plotting coefficients of drug main effect versus PET
    % neurotransmitter profiles
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDrugEffect_PET_GABAa_FLUMAZENIL_c11.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([10, 80]);
    
    plot(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'k.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, PET_GABAa_FLUMAZENIL_c11, 'Type', 'Spearman', 'Rows', 'complete');
    text(0.8, 20, strcat('\rho=', num2str(rho)), 'Color', 'k', 'FontSize', 12);
    text(0.8, 15, strcat('p=', num2str(pVal_corr)), 'Color', 'k', 'FontSize', 12);
    
    lsline;
    xlabel('coefficients of group x drug interaction (z-scores)');
    ylabel('GABA receptor density');
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupDrugInteraction_PET_GABAa_FLUMAZENIL_c11.eps'));

    %%%% plotting coefficients of drug main effect versus GABA receptor expression from Allen atlas %%%%
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_drugMain, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_drugMain, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_drugMain, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of drug main effect (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputDrugEffect_AHBA_GABA.eps'));
    
    % plotting coefficients of drug x group interaction versus GABA receptor expression
    f = figure('Visible', 'off'); set(gcf, 'color', 'white'); hold on;
    set(gca, 'FontSize', 20);
    xlim([-1, 1]); ylim([0, 1]);
    
    plot(betas_groupDrugInteraction, GABRA1, 'r.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA1, 'Rows', 'complete');
    text(-0.9, 0.98, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'r', 'FontSize', 12);
    
    plot(betas_groupDrugInteraction, GABRA2, 'b.', 'MarkerSize', 20);
    [rho, pVal_corr] = corr(betas_groupDrugInteraction, GABRA2, 'Rows', 'complete');
    text(-0.9, 0.95, strcat('\rho=', num2str(rho), ', p=', num2str(pVal_corr)), 'Color', 'b', 'FontSize', 12);
   
    lsline;
    xlabel('coefficients of group x drug interaction (z-score)');
    ylabel('GABA receptor expression');
    legend('\alpha1', '\alpha2', 'location', 'northeastoutside'); legend boxoff;
    
    saveas(f, strcat(resultsDir, 'emotionrec_', currentContrast, '_controlInputGroupDrugInteraction_AHBA_GABA.eps'));
end

%% plot coefficient maps for neutral emotionrec

betas_drugMain(isnan(betas_drugMain)) = 0;

[f, ax, ph, ~] = fcn_lausannesurf(betas_drugMain, redbluecmap);
view(ax(1), [-90, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [90, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, 'betas_drugMain_cortex1.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, 'betas_drugMain_cortex2.svg'); saveas(f(2), savePath); close(f(2));

[f, ax, ph, ~] = fcn_lausannesurf(betas_drugMain, redbluecmap);
view(ax(1), [-270, 0]); lighting(ax(1), 'gouraud'); camlight(ax(1), 'headlight'); material(ph(1), 'dull');
view(ax(2), [270, 0]); lighting(ax(2), 'gouraud'); camlight(ax(2), 'headlight'); material(ph(2), 'dull');
axis(ax(1), 'off'); axis(ax(2), 'off');
savePath = strcat(resultsDir, 'betas_drugMain_cortex3.svg'); saveas(f(1), savePath); close(f(1));
savePath = strcat(resultsDir, 'betas_drugMain_cortex4.svg'); saveas(f(2), savePath); close(f(2));

figure; [f1, f2] = plot_subcortvol(betas_drugMain(subcorticalIndices), subcorticalIndices, subcorticalIndices, nifti, redbluecmap);
savePath1 = strcat(resultsDir, 'betas_drugMain_subcortex1.svg'); saveas(f1, savePath1); 
savePath2 = strcat(resultsDir, 'betas_drugMain_subcortex2.svg'); saveas(f2, savePath2); 
close(f1); close(f2);


%% structural null model analysis

nullResultsDir = 'Results/minimalControl/avge_FD_thresh_0.5_parcelCoverageThresh_0.5_EuclideanNorm_allNodes_nullModel/';

allControlEnergies_emotionid_null = readtable(strcat(nullResultsDir, 'allControlEnergies_emotionid.csv'));
allControlEnergies_emotionrec_null = readtable(strcat(nullResultsDir, 'allControlEnergies_emotionrec.csv'));

contrast = allControlEnergies_emotionid.contrast;

for i = 1:nContrasts
    currentContrast = contrastLabels{i}
    allControlEnergies_emotionid_currentContrast = allControlEnergies_emotionid(strcmp(contrast, currentContrast), :);
    allControlEnergies_emotionid_currentContrast_null = allControlEnergies_emotionid_null(strcmp(contrast, currentContrast), :);
    
    drug = allControlEnergies_emotionid_currentContrast.drug;
    group = allControlEnergies_emotionid_currentContrast.group;
    
    [~, p, ci, stats] = ttest(allControlEnergies_emotionid_currentContrast.persistence_allNodes(drug==1), allControlEnergies_emotionid_currentContrast_null.persistence_allNodes(drug==1))
    figure, boxplot([allControlEnergies_emotionid_currentContrast.persistence_allNodes(drug==1), allControlEnergies_emotionid_currentContrast_null.persistence_allNodes(drug==1)]);
end

%% 

% plotting control trajectories with drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionid_drug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionid_drug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories without drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionid_noDrug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(1001, mean([0, refLines(1)]), 'Visual');
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionid_noDrug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories with drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionrec_drug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionrec_drug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);

% plotting control trajectories without drug
f = figure('visible', 'off'); set(gcf, 'color', 'w'); hold on;
imagesc(mean(U_opt_stability_emotionrec_noDrug, 3)'); cb = colorbar; cb.Location = 'westoutside';
xlabel('Time (a.u.)');
ylabel('Control input by region');
refLines = cumsum([sum(finalLabels==1), sum(finalLabels==2), sum(finalLabels==3), sum(finalLabels==4), sum(finalLabels==5), sum(finalLabels==6), sum(finalLabels==7), sum(finalLabels==8)]);
xlim([0, nTimeSteps]); ylim([0, nNodes]);
h = refline(0, refLines(1)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([0, refLines(1)]), 'Visual');
h = refline(0, refLines(2)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(1), refLines(2)]), 'Somatomator');
h = refline(0, refLines(3)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(2), refLines(3)]), 'Dorsal Attention');
h = refline(0, refLines(4)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(3), refLines(4)]), 'Ventral Attention');
h = refline(0, refLines(5)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(4), refLines(5)]), 'Limbic');
h = refline(0, refLines(6)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(5), refLines(6)]), 'Frontoparietal Control');
h = refline(0, refLines(7)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(6), refLines(7)]), 'Default Mode');
h = refline(0, refLines(8)); h.Color = 'r'; h.LineWidth = 2; text(nTimeSteps, mean([refLines(7), refLines(8)]), 'Subcortical');
ax = gca;
ax.FontSize = 20;
currentSaveFileName = strcat(resultsDir, 'stateTrajectories_emotionrec_noDrug_', currentContrastLabel, '.svg');
saveas(f, currentSaveFileName);
close(f);